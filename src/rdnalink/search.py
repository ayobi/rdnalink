"""One MMseqs2 nucleotide search per marker window, as in the rdnalink test (protocol note 2).

- Exhaustive: --max-seqs and --max-accept equal the number of distinct reference windows, so every candidate with a
  k-mer hit is aligned (a capped prefilter lost hits that belong in the top ranks in the test's pilots).
- A hit counts when the alignment covers >= 80% of the shorter of the two windows. Identity = identical
  positions / alignment columns. Hits are ranked by identity, then identical positions, then target ID; bit scores
  and E-values are kept but never used for ranking (MMseqs2's bit scores can differ by a few bits between runs of
  the same search). Queries are written in sorted order, so a rerun gives the same rows in the same order.
- Identical windows are searched once and their hits copied to every read that carries them.
- The best --keep hits per query and marker are stored (300 by default, as the test's calls read), plus any hits
  tied with the last one. A query that is also in the reference keeps its self hit, marked in the 'self' column.

The reference is a window prefix, or a folder built by `rdnalink ref`, whose distinct windows and indexes are then
used as they are (the same ones this module would build). Writes <out>/<marker>.hits.tsv.gz and <out>/search.log;
a finished chunk is kept, so a stopped run resumes.
"""
import collections
import gzip
import os
import shutil
import subprocess
import sys

from .reference import resolve, target_db
from .seqio import read_fasta, read_id, sha256

MARKERS = ('ssu', 'its', 'lsu')
COLS = 'query,target,nident,alnlen,qstart,qend,qlen,tstart,tend,tlen,bits,evalue'
HEADER = 'query\ttarget\trank\tidentity\tnident\talnlen\tqlen\ttlen\tqcov\ttcov\tbits\tevalue\tself\n'


def _run(cmd, log, hint=''):
    log.write('$ ' + ' '.join(cmd) + '\n')
    log.flush()
    p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if p.returncode:
        sys.exit(f'failed ({p.returncode}): {" ".join(cmd)}; see {log.name}' + (f'\n{hint}' if hint else ''))


def _windows(prefix, marker):
    for path in (f'{prefix}.{marker}.fasta', f'{prefix}.{marker}.fasta.gz'):
        if os.path.exists(path):
            return {read_id(h): s.upper() for h, s in read_fasta(path)}
    sys.exit(f'no {prefix}.{marker}.fasta(.gz); run rdnalink windows first')


def distinct_windows(t):
    """{read: window} -> {representative: [reads carrying that window]}; the representative is the first read ID
    in sort order."""
    by_seq = collections.defaultdict(list)
    for r, s in t.items():
        by_seq[s].append(r)
    members = {}
    for rs in by_seq.values():
        rs.sort()
        members[rs[0]] = rs
    return members


def build_target(d, t, members, mmseqs, threads, kmer, log):
    """MMseqs2 database and index of the distinct windows in folder d; skipped when already built."""
    os.makedirs(d, exist_ok=True)
    tdb = os.path.join(d, 'tdb')
    if not os.path.exists(tdb + '.done'):
        with open(os.path.join(d, 'target.fasta'), 'w') as fh:
            for rep in sorted(members):
                fh.write(f'>{rep}\n{t[rep]}\n')
        _run([mmseqs, 'createdb', os.path.join(d, 'target.fasta'), tdb, '--dbtype', '2'], log)
        _run([mmseqs, 'createindex', tdb, os.path.join(d, 'idxtmp'), '--search-type', '3', '--strand', '1',
              '--threads', str(threads)] + (['-k', str(kmer)] if kmer else []), log,
             hint='If the log says the database does not fit into memory: the default k-mer size (15) needs about '
                  '10 GB of RAM; use --kmer 13.')
        shutil.rmtree(os.path.join(d, 'idxtmp'), ignore_errors=True)
        open(tdb + '.done', 'w').close()
    return tdb


def _rank_and_write(raw, fo, members_t, members_q, keep, min_cov):
    # Queries are written in sorted order: MMseqs2 returns them in an order that varies with thread scheduling.
    lines = {}

    def flush(uq, rows):
        rows.sort(key=lambda x: (-x[0], -x[3], x[2]))
        take, n = [], 0
        for r in rows:
            tie = take and (r[0], r[3]) == (take[-1][0], take[-1][3])
            if n >= keep and not tie:
                break
            take.append(r)
            n += len(members_t[r[2]])
        exp = sorted(((r[0], r[3], t, r) for r in take for t in members_t[r[2]]), key=lambda x: (-x[0], -x[1], x[2]))
        last = exp[keep - 1][:2] if len(exp) >= keep else None
        exp = [x for i, x in enumerate(exp) if i < keep or (last and x[:2] == last)]
        lines[uq] = [f'{qread}\t{t}\t{i}\t{idn:.4f}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]:.4f}\t{r[8]:.4f}\t'
                     f'{r[1]}\t{r[9]}\t{int(t == qread)}\n'
                     for qread in members_q[uq] for i, (idn, _, t, r) in enumerate(exp, 1)]
    cur, rows = None, []
    with open(raw) as fh:
        for line in fh:
            v = line.rstrip('\n').split('\t')
            if v[0] != cur:
                if cur is not None:
                    flush(cur, rows)
                cur, rows = v[0], []
            ql, tl = int(v[6]), int(v[9])
            qc = (abs(int(v[5]) - int(v[4])) + 1) / ql
            tc = (abs(int(v[8]) - int(v[7])) + 1) / tl
            if (qc if ql <= tl else tc) < min_cov:
                continue
            nid, aln = int(v[2]), int(v[3])
            rows.append((100 * nid / aln, float(v[10]), v[1], nid, aln, ql, tl, qc, tc, v[11]))
        if cur is not None:
            flush(cur, rows)
    for uq in sorted(lines):
        fo.writelines(lines[uq])


def search(query_prefix, ref_prefix, out, mmseqs='mmseqs', threads=4, keep=300, min_cov=0.8, sens=5.7,
           chunk=200, markers=MARKERS, kmer=0):
    if shutil.which(mmseqs) is None and not os.path.exists(mmseqs):
        sys.exit(f'{mmseqs} not found on PATH')
    ref_windows, man = resolve(ref_prefix)
    if man:
        if kmer and kmer != man['mmseqs']['kmer']:
            sys.exit(f'{ref_prefix} was indexed with --kmer {man["mmseqs"]["kmer"] or "default"}; '
                     f'search it with that, or rebuild it')
        kmer = man['mmseqs']['kmer']
    os.makedirs(out, exist_ok=True)
    log = open(os.path.join(out, 'search.log'), 'a')
    summary = {}
    for m in markers:
        q, t = _windows(query_prefix, m), _windows(ref_windows, m)
        members_t = distinct_windows(t)
        rep_of = {t[rep]: rep for rep in members_t}
        qseq = {}
        members_q = collections.defaultdict(list)
        for r, s in q.items():
            u = rep_of.get(s, 'Q' + r)                      # a query identical to a reference window shares its ID
            qseq[u] = s
            members_q[u].append(r)
        uq = sorted(members_q)
        d = os.path.join(out, 'work', m)
        os.makedirs(d, exist_ok=True)
        if man:
            if sha256(f'{ref_windows}.{m}.fasta') != man['markers'][m]['windows_sha256']:
                sys.exit(f'{ref_windows}.{m}.fasta has changed since {ref_prefix} was built; rebuild it')
            tdb = target_db(ref_prefix, m)
            if not os.path.exists(tdb + '.done'):
                sys.exit(f'{tdb} is missing; rebuild {ref_prefix}')
        else:
            tdb = build_target(d, t, members_t, mmseqs, threads, kmer, log)
        cap = str(len(members_t) + 1)
        params = ['--search-type', '3', '-s', str(sens), '--max-seqs', cap, '--max-accept', cap, '-e', '1e6',
                  '--strand', '1', '-a', '--threads', str(threads)] + (['-k', str(kmer)] if kmer else [])
        parts = []
        for ci in range(0, len(uq), chunk):
            ch = uq[ci:ci + chunk]
            part = os.path.join(d, f'chunk_{ci // chunk:05d}.tsv.gz')
            parts.append(part)
            if os.path.exists(part):
                continue
            cd = os.path.join(d, f'c{ci // chunk:05d}')
            shutil.rmtree(cd, ignore_errors=True)
            os.makedirs(cd)
            with open(os.path.join(cd, 'q.fasta'), 'w') as fh:
                for u in ch:
                    fh.write(f'>{u}\n{qseq[u]}\n')
            _run([mmseqs, 'createdb', os.path.join(cd, 'q.fasta'), os.path.join(cd, 'qdb'), '--dbtype', '2'], log)
            _run([mmseqs, 'search', os.path.join(cd, 'qdb'), tdb, os.path.join(cd, 'aln'), os.path.join(cd, 'tmp')]
                 + params, log)
            raw = os.path.join(cd, 'aln.tsv')
            _run([mmseqs, 'convertalis', os.path.join(cd, 'qdb'), tdb, os.path.join(cd, 'aln'), raw,
                  '--format-output', COLS, '--threads', str(threads)], log)
            with gzip.open(part + '.tmp', 'wt', compresslevel=4) as fo:
                _rank_and_write(raw, fo, members_t, {u: members_q[u] for u in ch}, keep, min_cov)
            os.replace(part + '.tmp', part)
            shutil.rmtree(cd)
        res = os.path.join(out, f'{m}.hits.tsv.gz')
        with open(res, 'wb') as fo:
            fo.write(gzip.compress(HEADER.encode()))
            for part in parts:
                with open(part, 'rb') as fi:
                    shutil.copyfileobj(fi, fo)
        summary[m] = (len(q), len(uq), len(t), len(members_t))
    log.close()
    return summary
