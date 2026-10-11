#!/usr/bin/env python3
"""Phase 0 search (protocol section 5): one MMseqs2 nucleotide search per marker.

Queries are the test-group reads with all three windows >= 300 nt (section 4). The reference is every window
cut for that marker, from the whole long-read set. A hit counts only if the alignment covers >= 80% of the shorter
of the two windows. Identity = identical positions / alignment columns (MMseqs2 pident, seq-id-mode 0; its
nident column is 0 in nucleotide searches without backtraces). Self-hits are dropped. Hits are ranked by identity,
then bit score, then target id, and up to 1,000 are stored per query and marker: the protocol uses the first 300,
and the stored 1,000 serve the depth-check rerun (section 5) without a new search.

Identical windows are searched once. Each marker's windows are collapsed to unique sequences, unique query
sequences are searched against unique reference sequences, and every hit is expanded back to all reads that carry
that exact sequence. Identical sequences give identical alignments, so this changes no hit's identity, coverage or
bit score, only the time. It also makes the 1,000 prefilter slots go further. (MMseqs2 E-values shrink with the
smaller database; only hits with E <= 1e-3 are reported, a bound no alignment covering 80% of a 300 nt window at
the identities used here comes near.)

  python3 search_markers.py --windows ../windows/win --out full [--threads 8] [--kmer 13] [--max-queries 300]

Reads <windows>.ssu.fasta, <windows>.full.fasta (ITS) and <windows>.lsu.fasta. Writes, in --out:
  queries.tsv               query id, group, SSU / ITS / LSU window lengths
  hits.<marker>.tsv.gz      query, target, rank, identity, alnlen, qlen, tlen, cov_short, bits, evalue
  summary.tsv               per group and marker: queries, with hits, full lists, saturated lists, ...
  manifest.json             MMseqs2 version, commands, parameters, unique counts, hashes, counts, timings
--max-queries N searches a random sample of N queries (seed 20261007), for timing.
"""
import argparse, gzip, hashlib, json, os, random, shutil, statistics, subprocess, sys, time
from collections import Counter, defaultdict

MARKERS = [('ssu', 'ssu'), ('its', 'full'), ('lsu', 'lsu')]     # marker name, window file suffix
TEST = [('Glomeromycota', 'p__Glomeromycota;'), ('Agaricales', 'o__Agaricales;'),
        ('Hypocreales', 'o__Hypocreales;'), ('Rhabditida', 'o__Rhabditida;')]
MIN_LEN = 300            # section 4
MIN_COV = 0.80           # section 5: of the shorter window
TOP = 300                # section 5: hits kept for the protocol
STORE = 1000             # stored for the depth-check rerun
MAX_SEQS = 1000          # MMseqs2 prefilter candidates per query (unique sequences)
TIE = 0.002              # 0.2 percentage points (unique-best-name rule, section 6)
SEED = 20261007
FORMAT = 'query,target,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen,bits,evalue'


def key(h):
    return h.split()[0].split(';')[0].split('|')[0]


def group_of(h):
    for name, pat in TEST:
        if pat in h:
            return name
    return None


def read_fasta(path):
    head, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if head is not None:
                    yield head, ''.join(chunks)
                head, chunks = line[1:].rstrip('\n'), []
            else:
                chunks.append(line.strip())
    if head is not None:
        yield head, ''.join(chunks)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def run(cmd, log, timings, label):
    t = time.time()
    log.write('$ ' + ' '.join(cmd) + '\n'); log.flush()
    r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    timings[label] = round(time.time() - t, 1)
    if r.returncode != 0:
        sys.exit(f'failed ({r.returncode}): {" ".join(cmd)}  - see {log.name}')


def expand(rows_rep, members, q):
    """rows_rep: hits of q's sequence to unique targets, sorted by (-identity, -bits). Expand each unique target to
    all reads carrying it (minus q itself), order ties by target id, and keep the first STORE."""
    out, i = [], 0
    while i < len(rows_rep) and len(out) < STORE:
        j, grp = i, []
        while j < len(rows_rep) and rows_rep[j][0] == rows_rep[i][0] and rows_rep[j][1] == rows_rep[i][1]:
            ident, bits, t_rep, aln, ql, tl, cov, ev = rows_rep[j]
            for mbr in members[t_rep]:
                if mbr != q:
                    grp.append((mbr, ident, aln, ql, tl, cov, bits, ev))
            j += 1
        grp.sort(key=lambda r: r[0])
        out.extend(grp[:STORE - len(out)])
        i = j
    return out


def filter_hits(m8, rep_queries, members, out_gz, stats_m, short_ref):
    """Stream MMseqs2 output (one block per unique query sequence), filter by coverage, expand, rank, store."""
    per_q, seen = {}, set()
    raw = kept_lines = 0

    def flush(r, rows_rep):
        nonlocal kept_lines
        rows_rep.sort(key=lambda x: (-x[0], -x[1]))
        for q in rep_queries[r]:
            rows = expand(rows_rep, members, q)
            for i, (t, ident, aln, ql, tl, cov, bits, ev) in enumerate(rows, 1):
                out.write(f'{q}\t{t}\t{i}\t{ident:.6f}\t{aln}\t{ql}\t{tl}\t{cov:.4f}\t{bits:g}\t{ev:g}\n')
            n = len(rows)
            best = rows[0][1] if rows else None
            sat = n >= TOP and rows[TOP - 1][1] >= best - TIE
            per_q[q] = (n, best, sat, bool(rows) and rows[0][0] in short_ref)
            kept_lines += n

    with open(m8) as fh, gzip.open(out_gz, 'wt') as out:
        out.write('query\ttarget\trank\tidentity\talnlen\tqlen\ttlen\tcov_short\tbits\tevalue\n')
        cur, rows = None, []
        for line in fh:
            c = line.rstrip('\n').split('\t')
            r, t = c[0], c[1]
            raw += 1
            if r != cur:
                if cur is not None:
                    flush(cur, rows)
                    seen.add(cur)
                if r in seen:
                    sys.exit(f'{m8}: hits for {r} are not contiguous; sort the file by query first')
                cur, rows = r, []
            ident, aln = float(c[2]) / 100, int(c[3])     # pident: matches / alignment columns, 3 decimals
            qs, qe, ql, ts, te, tl = (int(x) for x in c[4:10])
            qcov = (abs(qe - qs) + 1) / ql
            tcov = (abs(te - ts) + 1) / tl
            cov = qcov if ql <= tl else tcov
            if cov < MIN_COV or aln == 0:
                continue
            rows.append((ident, float(c[10]), t, aln, ql, tl, cov, float(c[11])))
        if cur is not None:
            flush(cur, rows)
    for qs_ in rep_queries.values():
        for q in qs_:
            per_q.setdefault(q, (0, None, False, False))
    stats_m['raw_alignments'] = raw
    stats_m['stored_hits'] = kept_lines
    return per_q


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--windows', required=True, help='prefix of the window FASTAs, e.g. ../windows/win')
    ap.add_argument('--out', required=True)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--split-memory-limit', default='0', help='MMseqs2 --split-memory-limit (e.g. 6G); 0 = all')
    ap.add_argument('--kmer', type=int, default=13,
                    help='MMseqs2 -k [13]; its default for nucleotides (15) needs a k-mer table of over 8 GB')
    ap.add_argument('--max-queries', type=int, default=0)
    ap.add_argument('--markers', default='ssu,its,lsu', help='comma-separated subset, for testing')
    ap.add_argument('--mmseqs', default='mmseqs')
    ap.add_argument('--keep-tmp', action='store_true')
    a = ap.parse_args()
    markers = [(m, s) for m, s in MARKERS if m in a.markers.split(',')]

    os.makedirs(a.out, exist_ok=True)
    for d in ('fasta', 'db', 'res', 'tmp'):
        os.makedirs(os.path.join(a.out, d), exist_ok=True)
    P = lambda *x: os.path.join(a.out, *x)
    timings, manifest = {}, {'started': time.strftime('%Y-%m-%d %H:%M:%S'), 'parameters': {
        'min_window_len': MIN_LEN, 'min_cov_shorter': MIN_COV, 'top': TOP, 'stored': STORE,
        'max_seqs_unique': MAX_SEQS, 'tie_pp': TIE * 100, 'seed': SEED, 'max_queries': a.max_queries,
        'threads': a.threads, 'split_memory_limit': a.split_memory_limit, 'kmer': a.kmer,
        'dereplication': 'exact sequence, per marker; hits expanded to every read carrying the sequence'}}

    # 1. windows -> unique sequences per marker (rep = first read with that exact sequence)
    lens, groups, inputs, rep_of, members = {}, {}, {}, {}, {}
    t = time.time()
    for m, suf in MARKERS:
        path = f'{a.windows}.{suf}.fasta'
        inputs[path] = sha256(path)
        lens[m], rep_of[m], members[m] = {}, {}, defaultdict(list)
        by_hash = {}
        with open(P('fasta', f'ref.{m}.fasta'), 'w') as out:
            for h, s in read_fasta(path):
                k = key(h)
                lens[m][k] = len(s)
                if k not in groups:
                    groups[k] = group_of(h)
                d = hashlib.blake2b(s.upper().encode(), digest_size=16).digest()
                r = by_hash.get(d)
                if r is None:
                    by_hash[d] = r = k
                    out.write(f'>{k}\n{s}\n')
                rep_of[m][k] = r
                members[m][r].append(k)
    timings['dereplicate'] = round(time.time() - t, 1)
    manifest['inputs'] = inputs
    queries = sorted(k for k, g in groups.items()
                     if g and all(lens[m].get(k, 0) >= MIN_LEN for m, _ in MARKERS))
    all_unique = {m: len({rep_of[m][q] for q in queries}) for m, _ in MARKERS}
    print(f'all {len(queries)} queries: unique sequences ' + ', '.join(f'{m} {all_unique[m]}' for m, _ in MARKERS))
    manifest['all_queries_unique'] = all_unique
    if a.max_queries and a.max_queries < len(queries):
        queries = sorted(random.Random(SEED).sample(queries, a.max_queries))
    with open(P('queries.tsv'), 'w') as out:
        out.write('query\tgroup\tssu_len\tits_len\tlsu_len\n')
        for k in queries:
            out.write(f'{k}\t{groups[k]}\t' + '\t'.join(str(lens[m][k]) for m, _ in MARKERS) + '\n')
    gcount = Counter(groups[k] for k in queries)
    print('windows: ' + ', '.join(f'{m} {len(lens[m])} ({len(members[m])} unique)' for m, _ in MARKERS))
    print(f'queries: {len(queries)} (' + ', '.join(f'{g} {gcount[g]}' for g, _ in TEST) + ')')

    # 2. MMseqs2, one search per marker, unique query sequences against unique reference sequences
    ver = subprocess.run([a.mmseqs, 'version'], capture_output=True, text=True).stdout.strip()
    manifest['mmseqs_version'] = ver
    manifest['commands'], manifest['unique'] = {}, {}
    summary_rows, per_marker = [], {}
    for m, _ in markers:
        rep_queries = defaultdict(list)
        for q in queries:
            rep_queries[rep_of[m][q]].append(q)
        with open(P('fasta', f'query.{m}.fasta'), 'w') as out:
            for h, s in read_fasta(P('fasta', f'ref.{m}.fasta')):
                if h in rep_queries:
                    out.write(f'>{h}\n{s}\n')
        manifest['unique'][m] = {'reference_windows': len(lens[m]), 'reference_unique': len(members[m]),
                                 'queries': len(queries), 'query_unique': len(rep_queries)}
        log = open(P(f'mmseqs.{m}.log'), 'w')
        cmds = [
            [a.mmseqs, 'createdb', P('fasta', f'ref.{m}.fasta'), P('db', f'ref_{m}'), '--dbtype', '2'],
            [a.mmseqs, 'createdb', P('fasta', f'query.{m}.fasta'), P('db', f'q_{m}'), '--dbtype', '2'],
            [a.mmseqs, 'search', P('db', f'q_{m}'), P('db', f'ref_{m}'), P('res', m), P('tmp', m),
             '--search-type', '3', '--strand', '1', '--max-seqs', str(MAX_SEQS), '--alignment-mode', '3',
             '--seq-id-mode', '0', '-e', '1e-3', '-k', str(a.kmer), '--threads', str(a.threads),
             '--split-memory-limit', a.split_memory_limit],
            [a.mmseqs, 'convertalis', P('db', f'q_{m}'), P('db', f'ref_{m}'), P('res', m), P('res', f'{m}.m8'),
             '--format-output', FORMAT, '--threads', str(a.threads)],
        ]
        manifest['commands'][m] = [' '.join(c) for c in cmds]
        for i, c in enumerate(cmds):
            run(c, log, timings, f'{m}.{["createdb_ref", "createdb_query", "search", "convertalis"][i]}')
        log.close()
        short_ref = {k for k, L in lens[m].items() if L < MIN_LEN}
        st = {}
        t = time.time()
        per_q = filter_hits(P('res', f'{m}.m8'), rep_queries, members[m], P(f'hits.{m}.tsv.gz'), st, short_ref)
        timings[f'{m}.filter'] = round(time.time() - t, 1)
        per_marker[m] = st
        print(f'{m}: {len(rep_queries)} unique queries against {len(members[m])} unique windows; '
              f'{st["raw_alignments"]} alignments, {st["stored_hits"]} hits stored '
              f'(search {timings[m + ".search"]:.0f} s)')
        for g in [x for x, _ in TEST] + ['ALL']:
            ks = [k for k in queries if g == 'ALL' or groups[k] == g]
            if not ks:
                continue
            ns = [per_q[k][0] for k in ks]
            bests = [per_q[k][1] for k in ks if per_q[k][1] is not None]
            summary_rows.append([g, m, len(ks), sum(n > 0 for n in ns), sum(n >= TOP for n in ns),
                                 sum(per_q[k][2] for k in ks), sum(per_q[k][3] for k in ks),
                                 int(statistics.median(ns)),
                                 f'{100 * statistics.median(bests):.2f}' if bests else ''])
        if not a.keep_tmp:
            os.remove(P('res', f'{m}.m8'))
            shutil.rmtree(P('tmp', m), ignore_errors=True)

    with open(P('summary.tsv'), 'w') as out:
        hdr = ['group', 'marker', 'queries', 'with_hits', f'full_{TOP}', 'saturated', 'best_hit_short_window',
               'median_hits_stored', 'median_best_identity_pct']
        out.write('\t'.join(hdr) + '\n')
        print('\t'.join(hdr))
        for r in summary_rows:
            out.write('\t'.join(map(str, r)) + '\n')
            print('\t'.join(map(str, r)))
    manifest['counts'] = {'queries': len(queries), 'queries_by_group': dict(gcount), 'per_marker': per_marker}
    manifest['timings_s'] = timings
    manifest['outputs'] = {f: sha256(P(f)) for f in ['queries.tsv', 'summary.tsv'] +
                           [f'hits.{m}.tsv.gz' for m, _ in markers]}
    manifest['finished'] = time.strftime('%Y-%m-%d %H:%M:%S')
    with open(P('manifest.json'), 'w') as out:
        json.dump(manifest, out, indent=1)
    print(f'saturated = list full and the {TOP}th hit within {TIE * 100:g} pp of the best; '
          f'best_hit_short_window = best hit is a window under {MIN_LEN} nt')


if __name__ == '__main__':
    main()
