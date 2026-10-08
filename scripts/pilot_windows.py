#!/usr/bin/env python3
"""Phase 0 tooling pilot: does ITSxRust cut SSU / ITS / LSU windows across the test groups?

Reads the pilot FASTA (EUKARYOME headers) and ITSxRust's outputs, and reports per group how many reads
gave all three windows at the protocol's minimum length (300 nt each). This checks the window tool only;
it touches no fold, threshold or outcome.

  python3 pilot_windows.py pilot/pilot.fasta pilot/out [--top 12]  > pilot/windows_report.tsv

expects pilot/out.full.fasta, pilot/out.ssu.fasta and pilot/out.lsu.fasta (written with --plain-ids).
"""
import re, statistics, sys
from collections import defaultdict

TEST = [('Glomeromycota', 'p__Glomeromycota;'), ('Agaricales', 'o__Agaricales;'),
        ('Hypocreales', 'o__Hypocreales;'), ('Rhabditida', 'o__Rhabditida;')]
MIN_LEN = 300


def lengths(path):
    out, rid, n = {}, None, 0
    try:
        fh = open(path)
    except FileNotFoundError:
        sys.stderr.write(f'missing {path}\n')
        return out
    with fh:
        for line in fh:
            if line.startswith('>'):
                if rid is not None:
                    out[rid] = n
                rid, n = line[1:].split()[0].split(';')[0].split('|')[0], 0
            else:
                n += len(line.strip())
    if rid is not None:
        out[rid] = n
    return out


def group_of(header):
    for name, pat in TEST:
        if pat in header:
            return name
    m = re.search(r'p__([^;]*)', header)
    return 'other: ' + (m.group(1) if m else '(none)')


def main():
    args = sys.argv[1:]
    top = None
    if '--top' in args:
        i = args.index('--top'); top = int(args[i + 1]); del args[i:i + 2]
    pilot, prefix = args[0], args[1]
    reads = {}
    with open(pilot) as fh:
        for line in fh:
            if line.startswith('>'):
                h = line[1:].strip()
                reads[h.split()[0].split(';')[0]] = group_of(h)
    full, ssu, lsu = (lengths(f'{prefix}.{r}.fasta') for r in ('full', 'ssu', 'lsu'))
    stats = defaultdict(lambda: {'n': 0, 'its': 0, 'ssu': 0, 'lsu': 0, 'all3': 0, 'L': defaultdict(list)})
    for rid, g in reads.items():
        for key in (g, 'ALL'):
            s = stats[key]
            s['n'] += 1
            s['its'] += rid in full
            s['ssu'] += rid in ssu
            s['lsu'] += rid in lsu
            ok = all(rid in d and d[rid] >= MIN_LEN for d in (full, ssu, lsu))
            s['all3'] += ok
            for name, d in (('ssu', ssu), ('its', full), ('lsu', lsu)):
                if rid in d:
                    s['L'][name].append(d[rid])
    others = sorted((k for k in stats if k.startswith('other')), key=lambda k: -stats[k]['n'])
    order = [t[0] for t in TEST] + (others[:top] if top else others) + ['ALL']
    print('group\treads\tITS_found\tSSU_found\tLSU_found\tall_three_ge300\tpct_all_three\t'
          'median_SSU\tmedian_ITS\tmedian_LSU')
    for k in order:
        s = stats[k]
        if not s['n']:
            continue
        med = [str(int(statistics.median(s['L'][x]))) if s['L'][x] else '' for x in ('ssu', 'its', 'lsu')]
        print(f"{k}\t{s['n']}\t{s['its']}\t{s['ssu']}\t{s['lsu']}\t{s['all3']}\t"
              f"{100 * s['all3'] / s['n']:.1f}\t" + '\t'.join(med))


if __name__ == '__main__':
    main()
