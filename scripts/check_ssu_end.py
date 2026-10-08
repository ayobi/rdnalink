#!/usr/bin/env python3
"""Referee for SSU/ITS boundary disagreements: which call sits at the real 3' end of the 18S?

The last ~26 nt of the 18S gene, GTAGGTGAACCTGC + GGAAGGATCATT (the ITS1 primer site), are conserved across
eukaryotes. On each read (oriented as the logs give it) this finds GTAGGTGAACCTGC with <= 2 mismatches near
the two calls, and asks whether each call's ITS start lies 26 nt (+/- 10) after a motif hit.

  python3 check_ssu_end.py pilot.fasta cut_itsx.log.tsv cut_barrnapfastj.log.tsv disagree.ids
  python3 check_ssu_end.py pilot.fasta --single cut_hybridfastj.log.tsv

--single checks one run's barrnap-cut reads (source barrnap or barrnap+motif): is the motif just before the
ITS start, inside the ITS (an unaligned insert near the 18S end), or not found?

Control: every read whose two calls agree within 5 nt at the SSU/ITS edge (shows how often the motif is
found at all, per group). Test: the reads listed in disagree.ids.
"""
import re, statistics, sys
from collections import Counter, defaultdict

MOTIF = 'GTAGGTGAACCTGC'
OFFSET, TOL, MAXMM = 26, 10, 2
COMP = str.maketrans('ACGTUNacgtun', 'TGCAANtgcaan')
TEST = [('Glomeromycota', 'p__Glomeromycota;'), ('Agaricales', 'o__Agaricales;'),
        ('Hypocreales', 'o__Hypocreales;'), ('Rhabditida', 'o__Rhabditida;')]


def rc(s):
    return s.translate(COMP)[::-1]


def key(h):
    return h.split()[0].split(';')[0].split('|')[0]


def group_of(h):
    for name, pat in TEST:
        if pat in h:
            return name
    m = re.search(r'p__([^;]*)', h)
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


def read_log(path):
    out = {}
    with open(path) as fh:
        cols = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            r = dict(zip(cols, line.rstrip('\n').split('\t')))
            if r.get('status') == 'ok':
                out[r['id']] = (r['strand'], int(r['its_start']))
    return out


def hits(seq, lo, hi):
    lo, hi = max(0, lo), min(len(seq) - len(MOTIF), hi)
    out = []
    for p in range(lo, hi + 1):
        mm = 0
        for a, b in zip(seq[p:p + len(MOTIF)], MOTIF):
            if a != b:
                mm += 1
                if mm > MAXMM:
                    break
        if mm <= MAXMM:
            out.append((p, mm))
    return out


def supported(h, its_start):
    return any(abs(its_start - (p + OFFSET)) <= TOL for p, _ in h)


def single(reads_path, log_path):
    calls = {}
    with open(log_path) as fh:
        cols = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            r = dict(zip(cols, line.rstrip('\n').split('\t')))
            if r.get('status') == 'ok' and r.get('source', '').startswith('barrnap'):
                calls[r['id']] = (r['strand'], int(r['its_start']), int(r['its_end']), r['source'])
    res = defaultdict(Counter); inside = []
    for h, s in read_fasta(reads_path):
        k = key(h)
        if k not in calls:
            continue
        strand, i, j, src = calls[k]
        o = s.upper() if strand == '+' else rc(s.upper())
        hh = hits(o, i - OFFSET - TOL, j - len(MOTIF))
        if any(abs(i - (p + OFFSET)) <= TOL for p, _ in hh):
            v = 'at ITS start'
        elif any(i + TOL < p + OFFSET < j for p, _ in hh):
            v = 'inside ITS'
            inside.append((k, group_of(h), i, j, ' '.join(f'{p}(mm{mm})' for p, mm in hh)))
        else:
            v = 'not found'
        for g in (group_of(h), 'ALL'):
            res[g][v] += 1; res[g][src] += 1
    print(f'{log_path}: barrnap-cut reads, where the 18S 3\' end motif lies')
    for g in ['ALL'] + sorted((g for g in res if g != 'ALL'), key=lambda g: -res[g]['at ITS start'] - res[g]['inside ITS'] - res[g]['not found']):
        c = res[g]
        print(f'  {g}: at ITS start {c["at ITS start"]}, inside ITS {c["inside ITS"]}, not found {c["not found"]}'
              f'  (moved by motif: {c["barrnap+motif"]})')
    if inside:
        print('motif inside the ITS: id, group, ITS start, ITS end, motif hits')
        for r in inside[:20]:
            print('  ' + '\t'.join(map(str, r)))


def main():
    if len(sys.argv) > 2 and sys.argv[2] == '--single':
        single(sys.argv[1], sys.argv[3])
        return
    reads_path, log_a, log_b, ids_path = sys.argv[1:5]
    a, b = read_log(log_a), read_log(log_b)
    test = {l.strip() for l in open(ids_path) if l.strip()}
    reads, groups = {}, {}
    for h, s in read_fasta(reads_path):
        k = key(h)
        if k in a and k in b:
            reads[k], groups[k] = s.upper(), group_of(h)

    ctrl = defaultdict(Counter); offs = []
    res = defaultdict(Counter); rows = []
    for k, s in reads.items():
        (sa, ia), (sb, ib) = a[k], b[k]
        if sa != sb:
            continue
        o = s if sa == '+' else rc(s)
        g = groups[k]
        if k in test:
            h = hits(o, min(ia, ib) - OFFSET - 200, max(ia, ib) + 50)
            va, vb = supported(h, ia), supported(h, ib)
            verdict = 'both' if va and vb else 'A only' if va else 'B only' if vb else 'neither'
            res[g][verdict] += 1; res['ALL'][verdict] += 1
            rows.append((k, g, ia, ib, verdict, ' '.join(f'{p}(mm{mm}, ends {p + OFFSET})' for p, mm in h)))
        elif abs(ia - ib) <= 5:
            h = hits(o, ia - OFFSET - 60, ia + 20)
            near = [ia - p for p, _ in h if abs(ia - (p + OFFSET)) <= TOL]
            ctrl[g]['reads'] += 1; ctrl['ALL']['reads'] += 1
            if near:
                ctrl[g]['found'] += 1; ctrl['ALL']['found'] += 1
                offs.append(near[0])

    print(f'control (calls agree within 5 nt): motif found {OFFSET} +/- {TOL} nt before the ITS start')
    order = ['ALL'] + sorted((g for g in ctrl if g != 'ALL'), key=lambda g: -ctrl[g]['reads'])
    for g in order[:14]:
        c = ctrl[g]
        print(f'  {g}: {c["found"]} of {c["reads"]} ({100 * c["found"] / c["reads"]:.1f}%)')
    if offs:
        print(f'  offset from motif start to ITS start: median {statistics.median(offs):g} nt')
    print(f'\ndisagreeing reads (A = {log_a}, B = {log_b}): which ITS start sits at the 18S 3\' end')
    for g in ['ALL'] + sorted((g for g in res if g != 'ALL'), key=lambda g: -sum(res[g].values())):
        c = res[g]
        print(f'  {g}: A only {c["A only"]}, B only {c["B only"]}, both {c["both"]}, neither {c["neither"]}')
    print('\nper read: id, group, A ITS start, B ITS start, verdict, motif hits (oriented read, 0-based)')
    for r in rows:
        print('  ' + '\t'.join(map(str, r)))


if __name__ == '__main__':
    main()
