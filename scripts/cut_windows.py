#!/usr/bin/env python3
"""Cut SSU / ITS / LSU windows from full-operon reads, from either of two boundary finders.

SSU = everything upstream of the ITS, ITS = ITS1 + 5.8S + ITS2, LSU = everything downstream of the ITS
(protocol section 4). Only the way the two outer ITS bounds are found differs:

  itsxfull  bounds from ITSxRust's full-ITS output (out.full.fasta, --plain-ids). Works whenever ITSxRust
            found the full ITS, from the 4-anchor chain or from the SSU + LSU partial chain; ITSxRust's own
            --region ssu/lsu export needs the 4-anchor chain. --check-ssu / --check-lsu compare the cut
            flanks with ITSxRust's own export where both exist (they should be identical).
  barrnap   bounds from barrnap >= 1.0 (--kingdom fun: Rfam SSU_rRNA_eukarya RF01960 and LSU_rRNA_eukarya
            RF02543, both eukaryote-wide): ITS = 18S end to 28S start. Column 6 of its GFF is the E-value;
            the lowest E-value hit of each model and strand is the anchor.
            --join also uses the other 18S / 28S hits on the anchor strand (>= 50 nt each): the ITS starts
            after the last 18S hit that ends before the anchor 28S, and ends before the first 28S hit after
            that. This joins a gene split into several hits (an intron, or a local alignment that stops
            early). A read with an 18S hit after its 28S anchor, or a 28S hit before its 18S anchor, holds
            more than one operon and is set aside as multi_operon.
            --fix-ssu-end: if the conserved 18S 3' end (GTAGGTGAACCTGC, <= 2 mismatches) lies inside barrnap's
            ITS rather than just before it, the ITS starts 27 nt after the motif instead (an insert, usually an
            intron, near the end of the 18S that --fast leaves unaligned). Such reads get source barrnap+motif.
  hybrid    itsxfull bounds where ITSxRust found the full ITS, barrnap bounds for the rest.
  compare   boundary agreement between two runs' log files, on reads both cut.
  explain   the reads on which two runs disagree by more than --min-diff nt, with barrnap's hits for each.

  python3 cut_windows.py itsxfull --reads pilot.fasta --full out.full.fasta --check-ssu out.ssu.fasta \\
      --check-lsu out.lsu.fasta --out cut_itsx
  python3 cut_windows.py barrnap  --reads pilot.fasta --gff pilot.barrnap.gff [--join] [--fix-ssu-end] --out cut_barrnap
  python3 cut_windows.py hybrid   --reads pilot.fasta --full out.full.fasta --gff pilot.barrnap.gff [--join] \\
      [--fix-ssu-end] --out cut_hybrid
  python3 cut_windows.py compare  --logs cut_itsx.log.tsv cut_barrnap.log.tsv
  python3 cut_windows.py explain  --logs cut_itsx.log.tsv cut_barrnap.log.tsv --gff pilot.barrnap.gff \\
      --reads pilot.fasta [--min-diff 20] [--show 15] --out explain

Writes <out>.ssu.fasta, <out>.full.fasta and <out>.lsu.fasta, with the reads' original headers, and
<out>.log.tsv with one line per read: id, length, status, source, strand, its_start, its_end (0-based,
half-open, on the oriented read).
"""
import argparse, math, re, statistics, sys
from collections import Counter, defaultdict

COMP = str.maketrans('ACGTUNacgtunRYKMSWBDHVrykmswbdhv', 'TGCAANtgcaanYRMKSWVHDByrmkswvhdb')
SSU_NAMES = ('18S_rRNA',)
LSU_NAMES = ('28S_rRNA',)
MIN_JOIN = 50                                       # nt; shorter extra hits are not joined
SSU_END_MOTIF = 'GTAGGTGAACCTGC'                    # conserved 18S 3' end (ITS1 primer site)
MOTIF_TO_ITS = 27                                   # nt from motif start to ITS start (ITSxRust's convention, pilot median)
MOTIF_TOL, MOTIF_MAXMM = 10, 2
TEST = [('Glomeromycota', 'p__Glomeromycota;'), ('Agaricales', 'o__Agaricales;'),
        ('Hypocreales', 'o__Hypocreales;'), ('Rhabditida', 'o__Rhabditida;')]


def rc(s):
    return s.translate(COMP)[::-1]


def key(header):
    return header.split()[0].split(';')[0].split('|')[0]


def group_of(header):
    for name, pat in TEST:
        if pat in header:
            return name
    m = re.search(r'p__([^;]*)', header)
    return 'other: ' + (m.group(1) if m else '(none)')


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


def wrap(s, n=80):
    return '\n'.join(s[i:i + n] for i in range(0, len(s), n))


def motif_hits(seq, lo, hi):
    lo, hi = max(0, lo), min(len(seq) - len(SSU_END_MOTIF), hi)
    out = []
    for p in range(lo, hi + 1):
        mm = 0
        for x, y in zip(seq[p:p + len(SSU_END_MOTIF)], SSU_END_MOTIF):
            if x != y:
                mm += 1
                if mm > MOTIF_MAXMM:
                    break
        if mm <= MOTIF_MAXMM:
            out.append((p, mm))
    return out


def fix_ssu_end(oriented, i, j):
    """Move the ITS start to the 18S 3' end motif when it lies inside barrnap's ITS (an insert, usually an
    intron, near the end of the 18S that --fast does not align). Returns the new ITS start, or None."""
    h = motif_hits(oriented.upper(), i - MOTIF_TO_ITS - MOTIF_TOL, j - len(SSU_END_MOTIF))
    if any(abs(p + MOTIF_TO_ITS - i) <= MOTIF_TOL for p, _ in h):
        return None                                  # barrnap's call already sits at the 18S end
    later = [(mm, p) for p, mm in h if i + MOTIF_TOL < p + MOTIF_TO_ITS < j]
    if not later:
        return None
    return min(later)[1] + MOTIF_TO_ITS


def bounds_itsxfull(reads, full_path):
    full = {key(h): s.upper() for h, s in read_fasta(full_path)}
    out = {}
    for k, seq in reads.items():
        f = full.get(k)
        if not f:
            out[k] = ('no_its', 'itsxrust', '', None, None)
            continue
        s = seq.upper()
        i = s.find(f)
        if i >= 0:
            out[k] = ('ok', 'itsxrust', '+', i, i + len(f))
            continue
        i = rc(s).find(f)
        out[k] = ('ok', 'itsxrust', '-', i, i + len(f)) if i >= 0 else \
                 ('its_not_in_read', 'itsxrust', '', None, None)
    return out


def evalue(x):
    try:
        return float(x)
    except ValueError:
        return 1.0                                  # '.' or blank: treat as the weakest possible hit


def load_gff(gff_path, reads=None, quiet=False):
    feats = defaultdict(list)                       # key -> [(name, start, end, strand, evalue)], forward strand
    n_lines, names = 0, Counter()
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                continue
            n_lines += 1
            m = re.search(r'Name=([^;]+)', c[8])
            if not m:
                continue
            names[m.group(1)] += 1
            feats[key(c[0])].append((m.group(1), int(c[3]), int(c[4]), c[6], evalue(c[5])))
    if not quiet:
        matched = sum(k in feats for k in reads) if reads is not None else len(feats)
        total = len(reads) if reads is not None else len(feats)
        sys.stderr.write(f'barrnap GFF: {n_lines} feature lines; names: '
                         + (', '.join(f'{n} {v}' for n, v in names.most_common()) or 'none')
                         + f'; {matched} of {total} reads have a feature\n')
        if not n_lines:
            sys.stderr.write('  the GFF has no features: check the end of the barrnap log for an error\n')
    return feats


def oriented_hits(feat_list, L):
    """{strand: {'SSU': [(s, e, ev)], 'LSU': [...]}} with 1-based coordinates on the oriented read."""
    hits = {'+': {'SSU': [], 'LSU': []}, '-': {'SSU': [], 'LSU': []}}
    for name, s, e, strand, ev in feat_list:
        marker = 'SSU' if name in SSU_NAMES else 'LSU' if name in LSU_NAMES else None
        if marker is None or strand not in '+-':
            continue
        if strand == '-':                            # GFF is 1-based on the forward strand; move to the RC
            s, e = L - e + 1, L - s + 1
        hits[strand][marker].append((s, e, ev))
    return hits


def bounds_barrnap(reads, gff_path, join=False, fix=False):
    feats = load_gff(gff_path, reads)
    out = {}
    for k, seq in reads.items():
        hits = oriented_hits(feats.get(k, []), len(seq))
        cands = []
        for strand in '+-':
            ssu, lsu = hits[strand]['SSU'], hits[strand]['LSU']
            if not ssu or not lsu:
                continue
            a = min(ssu, key=lambda h: h[2])         # anchors: lowest E-value hit of each model
            b = min(lsu, key=lambda h: h[2])
            if not a[1] < b[0]:
                continue
            i, j, st = a[1], b[0] - 1, 'ok'          # ITS = after 18S end .. before 28S start (0-based, half-open)
            if join:
                long_ssu = [h for h in ssu if h[1] - h[0] + 1 >= MIN_JOIN]
                long_lsu = [h for h in lsu if h[1] - h[0] + 1 >= MIN_JOIN]
                if any(h[0] > b[1] for h in long_ssu) or any(h[1] < a[0] for h in long_lsu):
                    st = 'multi_operon'
                i = max([a[1]] + [h[1] for h in long_ssu if h[1] < b[0]])
                j = min([b[0]] + [h[0] for h in long_lsu if h[0] > i]) - 1
            worst = max(a[2], b[2], 1e-300)
            cands.append((-math.log10(worst), strand, i, j, st))
        if cands:
            _, strand, i, j, st = max(cands)
            src = 'barrnap'
            if st == 'ok' and fix:
                n = fix_ssu_end(seq if strand == '+' else rc(seq), i, j)
                if n is not None:
                    i, src = n, 'barrnap+motif'
            out[k] = (st, src, strand, i, j) if st == 'ok' else (st, src, '', None, None)
            continue
        have = {m for strand in '+-' for m in ('SSU', 'LSU') if hits[strand][m]}
        if have >= {'SSU', 'LSU'}:
            st = 'strand_or_order'
        else:
            st = 'missing_' + '_'.join(m for m in ('SSU', 'LSU') if m not in have)
        out[k] = (st, 'barrnap', '', None, None)
    return out


def write_windows(reads, heads, b, out):
    win = {'ssu': {}, 'full': {}, 'lsu': {}}
    status = Counter()
    with open(f'{out}.log.tsv', 'w') as log:
        log.write('id\tlength\tstatus\tsource\tstrand\tits_start\tits_end\n')
        for k, (st, src, strand, i, j) in b.items():
            status[(st, src)] += 1
            log.write(f'{k}\t{len(reads[k])}\t{st}\t{src}\t{strand}\t'
                      f'{"" if i is None else i}\t{"" if j is None else j}\n')
            if st != 'ok':
                continue
            s = reads[k] if strand == '+' else rc(reads[k])
            win['ssu'][k], win['full'][k], win['lsu'][k] = s[:i], s[i:j], s[j:]
    for r, d in win.items():
        with open(f'{out}.{r}.fasta', 'w') as fh:
            for k, s in d.items():
                if s:
                    fh.write(f'>{heads[k]}\n{wrap(s)}\n')
    return win, status


def read_log(path):
    out = {}
    with open(path) as fh:
        cols = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            r = dict(zip(cols, line.rstrip('\n').split('\t')))
            if r.get('status') != 'ok':
                continue
            L, i, j = int(r['length']), int(r['its_start']), int(r['its_end'])
            if r['strand'] == '-':                   # back to forward-strand coordinates
                i, j = L - j, L - i
            out[r['id']] = (r['strand'], i, j, L)
    return out


def compare(paths):
    a, b = read_log(paths[0]), read_log(paths[1])
    both = sorted(set(a) & set(b))
    print(f'compare: {len(a)} cut in {paths[0]}, {len(b)} in {paths[1]}, {len(both)} in both, '
          f'{len(set(b) - set(a))} only in the second, {len(set(a) - set(b))} only in the first')
    if not both:
        return
    same = [k for k in both if a[k][0] == b[k][0]]
    print(f'  same strand: {len(same)} of {len(both)}')
    for label, idx in (('ITS start (SSU/ITS boundary)', 1), ('ITS end (ITS/LSU boundary)', 2)):
        d = sorted(abs(a[k][idx] - b[k][idx]) for k in both)
        p95 = d[min(len(d) - 1, int(0.95 * len(d)))]
        w5 = sum(x <= 5 for x in d); w20 = sum(x <= 20 for x in d)
        print(f'  {label}: median |diff| {statistics.median(d):g} nt, 95th percentile {p95} nt, '
              f'within 5 nt {100 * w5 / len(d):.1f}%, within 20 nt {100 * w20 / len(d):.1f}%, max {d[-1]} nt')


def explain(paths, gff_path, reads_path, min_diff, show, out):
    a, b = read_log(paths[0]), read_log(paths[1])
    groups = {}
    if reads_path:
        for h, _s in read_fasta(reads_path):
            groups[key(h)] = group_of(h)
    feats = load_gff(gff_path, quiet=True)
    rows, summ = [], Counter()
    for k in sorted(set(a) & set(b)):
        ds, de = b[k][1] - a[k][1], b[k][2] - a[k][2]
        if abs(ds) <= min_diff and abs(de) <= min_diff:
            continue
        L = a[k][3]
        f = sorted((s, e, n.split('_')[0], st, ev) for n, s, e, st, ev in feats.get(k, [])
                   if n in SSU_NAMES + LSU_NAMES)
        n18 = sum(x[2] == '18S' for x in f); n28 = sum(x[2] == '28S' for x in f)
        strands = {x[3] for x in f}
        summ['reads'] += 1
        summ['start off'] += abs(ds) > min_diff
        summ['end off'] += abs(de) > min_diff
        summ['>1 18S hit'] += n18 > 1
        summ['>1 28S hit'] += n28 > 1
        summ['hits on both strands'] += len(strands) > 1
        summ['ITS by B longer'] += (b[k][2] - b[k][1]) > (a[k][2] - a[k][1])
        summ['read > 6 kb'] += L > 6000
        summ['group: ' + groups.get(k, '?')] += 1
        hit_txt = ' | '.join(f'{n}{st}{s}-{e} E={ev:.0e}' for s, e, n, st, ev in f)
        rows.append((k, groups.get(k, ''), L, a[k][1], a[k][2], b[k][1], b[k][2], ds, de, n18, n28, hit_txt))
    with open(f'{out}.tsv', 'w') as fh:
        fh.write('id\tgroup\tlength\tA_its_start\tA_its_end\tB_its_start\tB_its_end\tdiff_start\tdiff_end\t'
                 'n18S\tn28S\tbarrnap_hits (forward strand, 1-based)\n')
        for r in rows:
            fh.write('\t'.join(map(str, r)) + '\n')
    print(f'explain: {summ["reads"]} reads differ by > {min_diff} nt at either boundary '
          f'(A = {paths[0]}, B = {paths[1]}; coordinates on the forward strand, 0-based)')
    for kk in ['start off', 'end off', 'ITS by B longer', '>1 18S hit', '>1 28S hit', 'hits on both strands',
               'read > 6 kb']:
        print(f'  {kk}: {summ[kk]}')
    for kk, v in sorted(((kk, v) for kk, v in summ.items() if kk.startswith('group: ')), key=lambda x: -x[1]):
        print(f'  {kk}: {v}')
    print(f'first {min(show, len(rows))} (all in {out}.tsv):')
    for r in rows[:show]:
        print(f'  {r[0]} {r[1]} len {r[2]}  A {r[3]}-{r[4]}  B {r[5]}-{r[6]}  diff {r[7]:+d}/{r[8]:+d}  {r[11]}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['itsxfull', 'barrnap', 'hybrid', 'compare', 'explain'])
    ap.add_argument('--reads')
    ap.add_argument('--full', help='ITSxRust full-ITS FASTA (itsxfull and hybrid modes)')
    ap.add_argument('--gff', help='barrnap GFF (barrnap, hybrid and explain modes)')
    ap.add_argument('--join', action='store_true', help='join split 18S / 28S hits; set multi-operon reads aside')
    ap.add_argument('--fix-ssu-end', action='store_true',
                    help="barrnap reads: move the ITS start to the 18S 3' end motif when it lies inside barrnap's ITS")
    ap.add_argument('--check-ssu'); ap.add_argument('--check-lsu')
    ap.add_argument('--logs', nargs=2, help='two <out>.log.tsv files (compare and explain modes)')
    ap.add_argument('--min-diff', type=int, default=20)
    ap.add_argument('--show', type=int, default=15)
    ap.add_argument('--out')
    a = ap.parse_args()

    if a.mode in ('compare', 'explain'):
        if not a.logs:
            sys.exit(f'--logs is required in {a.mode} mode')
        if a.mode == 'compare':
            compare(a.logs)
        else:
            if not (a.gff and a.out):
                sys.exit('--gff and --out are required in explain mode')
            explain(a.logs, a.gff, a.reads, a.min_diff, a.show, a.out)
        return
    if not (a.reads and a.out):
        sys.exit('--reads and --out are required')
    if a.mode in ('itsxfull', 'hybrid') and not a.full:
        sys.exit(f'--full is required in {a.mode} mode')
    if a.mode in ('barrnap', 'hybrid') and not a.gff:
        sys.exit(f'--gff is required in {a.mode} mode')

    heads, reads = {}, {}
    for h, s in read_fasta(a.reads):
        heads[key(h)] = h
        reads[key(h)] = s
    if a.mode == 'itsxfull':
        b = bounds_itsxfull(reads, a.full)
    elif a.mode == 'barrnap':
        b = bounds_barrnap(reads, a.gff, a.join, a.fix_ssu_end)
    else:
        bi, bb = bounds_itsxfull(reads, a.full), bounds_barrnap(reads, a.gff, a.join, a.fix_ssu_end)
        b = {k: bi[k] if bi[k][0] == 'ok' else bb[k] for k in reads}

    win, status = write_windows(reads, heads, b, a.out)
    flags = (' --join' if a.join else '') + (' --fix-ssu-end' if a.fix_ssu_end else '')
    print(f'{a.mode}{flags}: {len(reads)} reads; ' + ', '.join(
        f'{st}' + (f' ({src})' if a.mode == 'hybrid' or src.endswith('motif') else '') + f' {v}'
        for (st, src), v in status.most_common()))
    for label, path, region in (('SSU', a.check_ssu, 'ssu'), ('LSU', a.check_lsu, 'lsu')):
        if not path:
            continue
        ref = {key(h): s.upper() for h, s in read_fasta(path)}
        both = [k for k in ref if k in win[region]]
        same = sum(ref[k] == win[region][k].upper() for k in both)
        print(f'check {label}: {len(ref)} in ITSxRust export, {len(both)} also cut here, {same} identical')


if __name__ == '__main__':
    main()
