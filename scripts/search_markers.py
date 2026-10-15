#!/usr/bin/env python3
"""Phase 0, protocol section 5: one MMseqs2 search per marker (SSU, ITS, LSU).

Queries are the windows of the reads of the four test groups (protocol section 3) that have all three windows of
>= 300 nt (section 4, Amendment 1). The reference is every window of that marker in the whole long-read set. For
each query the script keeps the hits whose alignment covers >= 80% of the shorter of the two windows, with
identity = identical positions / alignment columns (gaps included, as BLAST's pident), ranked by identity, then
the number of identical positions, then target ID. It stores the best --keep hits (default 1,000) with their rank. Scoring reads ranks
<= 300 (section 5); the deeper ranks are there for the depth check, so a group that fails it is re-called from
ranks <= 1,000 of the same search instead of searching again.

Exhaustive by default. A capped search (--max-seqs N) passes only the N candidates with the best k-mer scores to
alignment, and the pilots showed that no cap up to 30,000 gives a stable identity-ranked list (13 October 2026).
So every candidate with a k-mer hit is aligned. MMseqs2's own no-prefilter mode (--prefilter-mode 2) is not used:
it aligns every pair on the zero diagonal, which misaligns windows that start at different offsets.

To save work, identical sequences are searched once: targets with the same window sequence are one database entry
and their hits are copied to every read that carries them, and so are queries. The result is the same as searching
every read, because identical sequences give identical alignments; ties are still broken by read ID. Queries run in
chunks; a finished chunk is kept, so a stopped run resumes where it left off.

Self hits (query == target) stay in the table and are marked; every fold drops them with the held-out reads, and
K1 drops them explicitly.

  python3 scripts/search_markers.py --windows windows --fasta ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta \\
      --out search --threads 32

Check the inputs first (read IDs, window orientation, query counts against Amendment 1 note 1), without searching:
  python3 scripts/search_markers.py ... --check
Pilot (the first N/4 queries of each group; writes to search/pilot/):
  python3 scripts/search_markers.py ... --pilot 400

Writes, per marker m: search/{m}.hits.tsv.gz, and search/queries.tsv (query reads, group, window lengths),
search/counts.tsv and search/manifest.txt (MMseqs2 version, parameters, input and output SHA-256, counts, times).
"""
import argparse
import collections
import concurrent.futures
import datetime
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time

MARKERS = ['ssu', 'its', 'lsu']
GROUPS = [('Glomeromycota', 'p'), ('Agaricales', 'o'), ('Hypocreales', 'o'), ('Rhabditida', 'o')]
# reads with all three windows >= 300 nt per group, Amendment 1 note 1 (section 3); the script checks them
EXPECTED = {'Glomeromycota': 5546, 'Agaricales': 9436, 'Hypocreales': 3318, 'Rhabditida': 2901}
TAX = re.compile(r'(?:^|[;|,])\s*([kpcofgs])__([^;|,]*)')        # as in eukaryome_inventory.py
MIN_WINDOW = 300
COLS = 'query,target,nident,alnlen,qstart,qend,qlen,tstart,tend,tlen,bits,evalue'
HEADER = 'query\ttarget\trank\tidentity\tnident\talnlen\tqlen\ttlen\tqcov\ttcov\tbits\tevalue\tself\n'
WINDOW_FILES = {'ssu': 'win.ssu.fasta', 'its': 'win.full.fasta', 'lsu': 'win.lsu.fasta'}   # cut_windows.py names


def read_id(header):
    """The read ID shared by the long-read FASTA and the window FASTAs: the header up to the first space or ';'.
    The script stops if any window ID is not a read ID, so a wrong guess here cannot pass silently."""
    return re.split(r'[\s;]', header, maxsplit=1)[0]


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def fasta(path):
    head, seq = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if head is not None:
                    yield head, ''.join(seq)
                head, seq = line[1:].rstrip('\n'), []
            else:
                seq.append(line.strip())
    if head is not None:
        yield head, ''.join(seq)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def find_window(d, m):
    for f in (WINDOW_FILES[m], WINDOW_FILES[m] + '.gz'):
        if os.path.exists(os.path.join(d, f)):
            return os.path.join(d, f)
    sys.exit(f'no {WINDOW_FILES[m]}(.gz) in {d}; pass --{m} explicitly')


def group_of(header):
    tax = {}
    for r, v in TAX.findall(header):
        tax.setdefault(r, v.strip())
    for g, r in GROUPS:
        if tax.get(r) == g:
            return g
    return ''


def ssu_orientation(path, n=2000):
    """Share of SSU windows ending in the conserved 18S 3' end on each strand (a check that windows are oriented)."""
    motif, rc = re.compile('GTAGGTGAACC'), re.compile('GGTTCACCTAC')
    fw = rv = k = 0
    for _, s in fasta(path):
        fw += bool(motif.search(s[-200:].upper()))
        rv += bool(rc.search(s[:200].upper()))
        k += 1
        if k >= n:
            break
    return fw, rv, k


def run(cmd, log):
    log.write('$ ' + ' '.join(cmd) + '\n')
    log.flush()
    p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if p.returncode:
        sys.exit(f'failed ({p.returncode}): {" ".join(cmd)}; see {log.name}')


def filter_chunk(raw, part, members_t, members_q, keep, min_cov, keep_within=0.0, max_keep=0):
    """Hits of one chunk of unique queries: coverage filter, ranking, expansion to reads; writes part (gzip TSV)
    and returns {query read: hits kept} and the number of alignments dropped by the coverage filter."""
    nhits, dropped = {}, 0
    tmp = part + '.tmp'
    with open(raw) as fh, gzip.open(tmp, 'wt', compresslevel=4) as fo:
        def flush(uq, rows):
            # rows: (identity, bits, rep, nident, alnlen, qlen, tlen, qcov, tcov, evalue) for unique targets
            # rank by identity, then identical positions (longer alignments first), then read ID; bit scores and
            # e-values are reported but not used, as they shift slightly with the composition of the database
            rows.sort(key=lambda x: (-x[0], -x[3], x[2]))
            best = rows[0][0] if rows else 0.0
            cap = max(keep, max_keep)
            take, n = [], 0
            for r in rows:                                  # enough unique targets to fill `keep` reads (and, with
                near = keep_within > 0 and r[0] >= best - keep_within   # keep_within, every read that close to the
                tie = take and (r[0], r[3]) == (take[-1][0], take[-1][3])  # best, up to max_keep), plus every tie
                if (n >= keep and not near and not tie) or (n >= cap and not tie):  # at the boundary
                    break
                take.append(r)
                n += len(members_t[r[2]])
            exp = sorted(((r[0], r[3], t, r) for r in take for t in members_t[r[2]]),
                         key=lambda x: (-x[0], -x[1], x[2]))
            lim = keep
            if keep_within > 0:
                lim = max(keep, sum(1 for x in exp if x[0] >= best - keep_within))
            exp = exp[:min(lim, cap)]
            for qread in members_q[uq]:
                for i, (idn, _, t, r) in enumerate(exp, 1):
                    fo.write(f'{qread}\t{t}\t{i}\t{idn:.4f}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]:.4f}\t{r[8]:.4f}\t'
                             f'{r[1]}\t{r[9]}\t{int(t == qread)}\n')
                nhits[qread] = len(exp)
        cur, rows = None, []
        for line in fh:
            v = line.rstrip('\n').split('\t')
            if v[0] != cur:
                if cur is not None:
                    flush(cur, rows)
                cur, rows = v[0], []
            ql, tl = int(v[6]), int(v[9])
            qc = (int(v[5]) - int(v[4]) + 1) / ql
            tc = (abs(int(v[8]) - int(v[7])) + 1) / tl
            if (qc if ql <= tl else tc) < min_cov:
                dropped += 1
                continue
            nid, aln = int(v[2]), int(v[3])
            rows.append((100 * nid / aln, float(v[10]), v[1], nid, aln, ql, tl, qc, tc, v[11]))
        if cur is not None:
            flush(cur, rows)
    for uq in members_q:                                    # queries with no hit at all
        for qread in members_q[uq]:
            nhits.setdefault(qread, 0)
    os.replace(tmp, part)
    with open(part + '.counts', 'w') as fh:
        fh.write(f'#dropped\t{dropped}\n')
        for qread, n in nhits.items():
            fh.write(f'{qread}\t{n}\n')
    os.remove(raw)
    return nhits, dropped


def read_counts(part):
    nhits, dropped = {}, 0
    with open(part + '.counts') as fh:
        for line in fh:
            k, v = line.rstrip('\n').split('\t')
            if k == '#dropped':
                dropped = int(v)
            else:
                nhits[k] = int(v)
    return nhits, dropped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--windows', required=True, help='folder of the window FASTAs from cut_windows.py')
    for m in MARKERS:
        ap.add_argument(f'--{m}', help=f'{m.upper()} window FASTA (default: found in --windows)')
    ap.add_argument('--fasta', required=True, help='the long-read FASTA, for the taxonomy of each read')
    ap.add_argument('--out', default='search')
    ap.add_argument('--markers', default='ssu,its,lsu', help='markers to search (default all three)')
    ap.add_argument('--threads', type=int, default=os.cpu_count())
    ap.add_argument('--mmseqs', default='mmseqs')
    ap.add_argument('--keep', type=int, default=1000, help='hits stored per query and marker (scoring uses 300)')
    ap.add_argument('--max-seqs', type=int, default=0, help='cap the k-mer candidates per query (0 = exhaustive)')
    ap.add_argument('--keep-within', type=float, default=0.0, help='also store every hit within this many points '
                    'of the best identity (e.g. 0.2, the tie of section 6), up to --max-keep hits')
    ap.add_argument('--max-keep', type=int, default=20000)
    ap.add_argument('--no-dedupe', action='store_true', help='search every read separately (for checking only)')
    ap.add_argument('--evalue', default='1e6', help='MMseqs2 -e (default: no e-value filter; the protocol has none)')
    ap.add_argument('--chunk', type=int, default=200, help='unique queries per MMseqs2 run')
    ap.add_argument('--sens', type=float, default=5.7, help='MMseqs2 -s')
    ap.add_argument('--min-cov', type=float, default=0.8)
    ap.add_argument('--strand', type=int, default=1, help='MMseqs2 --strand (1 = forward; windows are oriented)')
    ap.add_argument('--kmer', type=int, default=0, help='MMseqs2 -k (default: MMseqs2 picks 15 for nucleotides, '
                    'whose index needs > 8 GB of RAM; 13 for small test machines only)')
    ap.add_argument('--check', action='store_true', help='run the input checks and write queries.tsv, then stop')
    ap.add_argument('--pilot', type=int, default=0, help='search only the first N/4 queries of each group')
    ap.add_argument('--no-expected', action='store_true', help='skip the check against the Amendment 1 counts')
    a = ap.parse_args()
    markers = [m for m in a.markers.split(',') if m]

    out = os.path.join(a.out, 'pilot') if a.pilot else a.out
    os.makedirs(out, exist_ok=True)
    log = open(os.path.join(out, 'mmseqs.log'), 'a')
    t0 = time.time()
    paths = {m: getattr(a, m) or find_window(a.windows, m) for m in MARKERS}

    grp, nfa = {}, set()
    for h, _ in fasta(a.fasta):
        nfa.add(read_id(h))
        g = group_of(h)
        if g:
            grp[read_id(h)] = g
    print(f'{len(nfa)} reads in {a.fasta}; {len(grp)} in the four test groups', flush=True)

    wseq = {m: {} for m in MARKERS}
    for m in MARKERS:
        n = 0
        for h, s in fasta(paths[m]):
            wseq[m][read_id(h)] = s.upper()
            n += 1
        if n != len(wseq[m]):
            sys.exit(f'{m}: {n} windows but {len(wseq[m])} distinct read IDs; adapt read_id()')
        print(f'{m}: {len(wseq[m])} windows, {len(set(wseq[m].values()))} distinct sequences ({paths[m]})',
              flush=True)
    unknown = set().union(*[set(wseq[m]) for m in MARKERS]) - nfa
    if unknown:
        sys.exit(f'{len(unknown)} window IDs are not read IDs of the FASTA, e.g. {sorted(unknown)[:3]}: '
                 'adapt read_id()')
    fw, rv, k = ssu_orientation(paths['ssu'])
    print(f'SSU orientation check: 18S 3\' end at the window end in {fw} of {k}, reverse complement at the start '
          f'in {rv}', flush=True)
    if rv > 0.01 * k and a.strand == 1:
        sys.exit('some SSU windows look reverse-complemented; check cut_windows.py output or use --strand 2')

    q = sorted(r for r, g in grp.items() if all(len(wseq[m].get(r, '')) >= MIN_WINDOW for m in MARKERS))
    per_group = collections.Counter(grp[r] for r in q)
    print('queries per group:', dict(per_group), flush=True)
    if not a.no_expected and dict(per_group) != EXPECTED:
        sys.exit(f'query counts differ from Amendment 1 note 1 ({EXPECTED}); stop and check the windows')
    with open(os.path.join(out, 'queries.tsv'), 'w') as fh:
        fh.write('read\tgroup\t' + '\t'.join(f'{m}_len' for m in MARKERS) + '\n')
        for r in q:
            fh.write(f'{r}\t{grp[r]}\t' + '\t'.join(str(len(wseq[m][r])) for m in MARKERS) + '\n')
    if a.check:
        print(f'checks passed; {len(q)} queries written to {out}/queries.tsv')
        return
    if a.pilot:
        by_g = collections.defaultdict(list)
        for r in q:
            by_g[grp[r]].append(r)
        q = sorted(x for g in by_g for x in by_g[g][:max(1, a.pilot // len(by_g))])

    ver = subprocess.run([a.mmseqs, 'version'], capture_output=True, text=True).stdout.strip()
    counts, outputs, times, uniq = [], {}, {}, {}
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=1)        # filters chunk i while MMseqs2 runs i+1
    for m in markers:
        tm = time.time()
        d = os.path.join(out, 'work', m)
        os.makedirs(d, exist_ok=True)
        # unique target sequences (representative = smallest read ID) and unique query sequences
        by_seq = collections.defaultdict(list)
        for r, s in wseq[m].items():
            by_seq[s].append(r)
        members_t, rep_of_seq = {}, {}
        for s, rs in by_seq.items():
            rs.sort()
            if a.no_dedupe:
                for r in rs:
                    members_t[r] = [r]
            else:
                members_t[rs[0]] = rs
            rep_of_seq[s] = rs[0]
        qby = collections.defaultdict(list)
        for r in q:
            qby[r if a.no_dedupe else rep_of_seq[wseq[m][r]]].append(r)
        uq = sorted(qby)
        uniq[m] = (len(members_t), len(uq))
        tdb = os.path.join(d, 'tdb')
        if not os.path.exists(tdb + '.done'):
            with open(os.path.join(d, 'target.fasta'), 'w') as fh:
                for rep in sorted(members_t):
                    fh.write(f'>{rep}\n{wseq[m][rep]}\n')
            run([a.mmseqs, 'createdb', os.path.join(d, 'target.fasta'), tdb, '--dbtype', '2'], log)
            idx = [a.mmseqs, 'createindex', tdb, os.path.join(d, 'idxtmp'), '--search-type', '3', '--strand',
                   str(a.strand), '--threads', str(a.threads)] + (['-k', str(a.kmer)] if a.kmer else [])
            run(idx, log)
            open(tdb + '.done', 'w').close()
        cap = a.max_seqs or len(members_t) + 1
        params = ['--search-type', '3', '-s', str(a.sens), '--max-seqs', str(cap), '--max-accept', str(cap),
                  '-e', a.evalue, '--strand', str(a.strand), '-a', '--threads', str(a.threads)]
        if a.kmer:
            params += ['-k', str(a.kmer)]
        chunks = [uq[i:i + a.chunk] for i in range(0, len(uq), a.chunk)]
        print(f'{m}: {len(q)} queries ({len(uq)} distinct) against {len(wseq[m])} windows ({len(members_t)} '
              f'distinct); {len(chunks)} chunks; {"exhaustive" if not a.max_seqs else f"cap {a.max_seqs}"}',
              flush=True)
        futures, parts = [], []
        for ci, ch in enumerate(chunks):
            part = os.path.join(d, f'chunk_{ci:05d}.tsv.gz')
            parts.append(part)
            if os.path.exists(part) and os.path.exists(part + '.counts'):
                continue
            cd = os.path.join(d, f'c{ci:05d}')
            shutil.rmtree(cd, ignore_errors=True)
            os.makedirs(cd)
            with open(os.path.join(cd, 'q.fasta'), 'w') as fh:
                for rep in ch:
                    fh.write(f'>{rep}\n{wseq[m][rep]}\n')
            run([a.mmseqs, 'createdb', os.path.join(cd, 'q.fasta'), os.path.join(cd, 'qdb'), '--dbtype', '2'], log)
            run([a.mmseqs, 'search', os.path.join(cd, 'qdb'), tdb, os.path.join(cd, 'aln'), os.path.join(cd, 'tmp')]
                + params, log)
            raw = os.path.join(d, f'c{ci:05d}.aln.tsv')
            run([a.mmseqs, 'convertalis', os.path.join(cd, 'qdb'), tdb, os.path.join(cd, 'aln'), raw,
                 '--format-output', COLS, '--threads', str(a.threads)], log)
            shutil.rmtree(cd)
            for f in futures:                                    # surface a filtering error early
                if f.done():
                    f.result()
            futures.append(pool.submit(filter_chunk, raw, part, members_t, {u: qby[u] for u in ch}, a.keep,
                                       a.min_cov, a.keep_within, a.max_keep))
            el = time.time() - tm
            print(f'  {m} chunk {ci + 1}/{len(chunks)} aligned ({el / 60:.1f} min so far)', flush=True)
        for f in futures:
            f.result()
        # assemble: header, then the chunks in order (concatenated gzip members are one valid gzip file)
        res = os.path.join(out, f'{m}.hits.tsv.gz')
        nhits, dropped = {}, 0
        with open(res, 'wb') as fo:
            fo.write(gzip.compress(HEADER.encode()))
            for part in parts:
                with open(part, 'rb') as fi:
                    shutil.copyfileobj(fi, fo)
                nh, dr = read_counts(part)
                nhits.update(nh)
                dropped += dr
        by_group = collections.defaultdict(list)
        for r in q:
            by_group[grp[r]].append(nhits.get(r, 0))
        for g, xs in sorted(by_group.items()):
            xs = sorted(xs)
            counts.append({'marker': m, 'group': g, 'queries': len(xs), 'no_hit': sum(x == 0 for x in xs),
                           'median_hits': xs[len(xs) // 2], 'full_300': sum(x >= 300 for x in xs),
                           f'full_{a.keep}': sum(x >= a.keep for x in xs)})
        times[m] = round(time.time() - tm)
        outputs[m] = (res, ' '.join(params))
        print(f'{m}: {times[m]} s; {sum(nhits.values())} hits kept; {dropped} alignments under {a.min_cov:.0%} '
              f'coverage dropped', flush=True)

    with open(os.path.join(out, 'counts.tsv'), 'w') as fh:
        keys = list(counts[0])
        fh.write('\t'.join(keys) + '\n')
        for c in counts:
            fh.write('\t'.join(str(c[k]) for k in keys) + '\n')
    with open(os.path.join(out, 'manifest.txt'), 'w') as fh:
        fh.write(f'date\t{datetime.datetime.now().isoformat(timespec="seconds")}\n')
        fh.write(f'script\tsearch_markers.py\t{sha256(os.path.abspath(__file__))}\n')
        fh.write(f'mmseqs2\t{ver}\n')
        for m in markers:
            fh.write(f'mmseqs2 search parameters\t{m}\t{outputs[m][1]}\n')
        fh.write(f'search\t{"exhaustive: every candidate with a k-mer hit is aligned" if not a.max_seqs else "capped"}'
                 f'; identical sequences searched once ({", ".join(f"{m} {uniq[m][0]} distinct targets, {uniq[m][1]} distinct queries" for m in markers)})\n')
        fh.write('ids\tqueries and targets written as read_id() (header up to the first space or ;)\n')
        fh.write(f'filter\tcoverage of the shorter window >= {a.min_cov}; identity = nident / alnlen\n')
        fh.write(f'ranking\tidentity, then identical positions, then target ID; {a.keep} kept per query (scoring '
                 f'uses 300)' + (f', plus every hit within {a.keep_within} points of the best, up to {a.max_keep}'
                                 if a.keep_within else '') + '\n')
        if a.pilot:
            fh.write(f'pilot\tfirst {a.pilot // 4} queries of each group\n')
        fh.write(f'input\t{a.fasta}\t{sha256(a.fasta)}\n')
        for m in MARKERS:
            fh.write(f'input\t{paths[m]}\t{sha256(paths[m])}\n')
        for m in markers:
            fh.write(f'output\t{outputs[m][0]}\t{sha256(outputs[m][0])}\n')
        fh.write(f'queries\t{len(q)}\t' + ', '.join(f'{g} {per_group[g]}' for g, _ in GROUPS) + '\n')
        for m in markers:
            fh.write(f'time\t{m}\t{times[m]} s\n')
        fh.write(f'time\ttotal\t{round(time.time() - t0)} s\n')
    print(f'done in {round(time.time() - t0)} s; see {out}/counts.tsv and {out}/manifest.txt')


if __name__ == '__main__':
    main()
