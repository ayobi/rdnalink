#!/usr/bin/env python3
"""Phase 0, protocol section 9 (K1): within-read marker conflicts, on every query read and on the spike-ins.

Each marker reads ranks <= 300 of the search table with the read itself dropped (note 2); a chimera spike-in drops
both source reads, a label swap its own read. A marker is confident for genus G when its best remaining hit is in
the named genus G, that identity is >= the genus threshold t_m(5), and it leads the best remaining hit of any other
genus (named or not) by >= 1 percentage point. If no hit of another genus is left, the lead is the best identity
minus the identity at the last rank when the list was full (a lower bound), and unlimited when it was not.

t_m(5) is the 5th lower percentile (numpy linear) of each read's identity to its nearest other read of the same
named genus, pooled over the group's usable genera (make_folds.py), with the full reference (only the read itself
dropped; same censoring rule as score_k2.py).

  F1 internal conflict  two markers confident for different named genera
  F2 label conflict     all three markers confident for one named genus other than the read's label, for reads
                        whose label is a named genus (the swapped label for a swap spike-in)

Recall: chimeras (both kinds) caught by F1, label swaps caught by F2, per group and pooled. Precision sample: 50 F1
and 50 F2 drawn at random (seed 20261007) from the real flags of all four groups (all of a type if fewer than 50);
precision_sample.tsv lists each read with its 20 best remaining hits per marker for the tree check.

  python3 scripts/score_k1.py --search search --fasta ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta \\
      --folds folds --spikeins spikeins/spikeins.tsv --out scores/k1
"""
import argparse
import collections
import datetime
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_calls import taxonomy  # noqa: E402
from score_k2 import EPS, GROUPS, MARKERS, SEED, TOP, load_hits, sha256  # noqa: E402

LEAD = 1.0
N_PREC = 50
N_HITS_TREE = 20


def marker_call(idn, tg, genus, held, full, t5):
    """(genus or None, best identity, lead, confident)"""
    best, g0, other = None, None, None
    for i, t in zip(idn, tg):
        if t in held:
            continue
        if best is None:
            best, g0 = float(i), genus.get(t)
            continue
        if genus.get(t) != g0:
            other = float(i)
            break
    if best is None:
        return None, None, None, False
    if other is not None:
        lead = best - other
    elif full and len(idn):
        lead = best - float(idn[-1])
    else:
        lead = float('inf')
    conf = g0 is not None and best >= t5 - EPS and lead >= LEAD - EPS
    return g0, best, lead, conf


def flags(calls, label):
    conf = {m: c[0] for m, c in calls.items() if c[3]}
    f1 = len(set(conf.values())) >= 2
    gs = {calls[m][0] for m in MARKERS}
    f2 = (label is not None and len(conf) == 3 and len(gs) == 1 and next(iter(gs)) != label)
    return f1, f2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--search', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--folds', required=True)
    ap.add_argument('--spikeins', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--groups', nargs='*', default=GROUPS)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sp, ge = taxonomy(a.fasta)
    genus = {r: g for r, g in ge.items() if g}
    group = {}
    with open(os.path.join(a.search, 'queries.tsv')) as fh:
        next(fh)
        for line in fh:
            r, g = line.rstrip('\n').split('\t')[:2]
            group[r] = g
    keep = set(group)
    hits = {}
    for m in MARKERS:
        hits[m], _ = load_hits(os.path.join(a.search, f'{m}.hits.tsv.gz'), keep, TOP, set())
    with_hits = {r for r in keep if all(r in hits[m] for m in MARKERS)}
    use_ge = collections.defaultdict(set)
    with open(os.path.join(a.folds, 'units.tsv')) as fh:
        next(fh)
        for line in fh:
            g_, rk, n = line.rstrip('\n').split('\t')[:3]
            if rk == 'genus':
                use_ge[g_].add(n)
    full = {m: {r: len(hits[m][r][1]) >= TOP for r in with_hits} for m in MARKERS}

    # genus thresholds t_m(5), full reference minus the read itself
    t5 = {}
    for grp in a.groups:
        pool = [r for r in with_hits if group[r] == grp and ge[r] in use_ge[grp]]
        t5[grp] = {}
        for m in MARKERS:
            vals = []
            for r in pool:
                idn, tg = hits[m][r]
                v = next((float(i) for i, t in zip(idn, tg) if t != r and genus.get(t) == ge[r]), None)
                if v is None:
                    v = float(idn[-1]) if full[m][r] and len(idn) else 0.0
                vals.append(v)
            t5[grp][m] = float(np.percentile(vals, 5)) if vals else float('inf')

    def score(windows, held, grp, label):
        calls = {}
        for m in MARKERS:
            src = windows[m]
            idn, tg = hits[m][src]
            calls[m] = marker_call(idn, tg, genus, held, full[m][src], t5[grp][m])
        return calls, flags(calls, label)

    # real flags
    real, real_rows = collections.defaultdict(list), []
    for r in sorted(with_hits):
        grp = group[r]
        if grp not in a.groups:
            continue
        calls, (f1, f2) = score({m: r for m in MARKERS}, {r}, grp, ge[r])
        row = [grp, r, ge[r] or '']
        for m in MARKERS:
            g0, b, ld, c = calls[m]
            row += [g0 or '', '' if b is None else f'{b:.4f}', '' if ld is None else f'{ld:.4f}', int(c)]
        row += [int(f1), int(f2)]
        real_rows.append(row)
        if f1:
            real['F1'].append(r)
        if f2:
            real['F2'].append(r)

    # spike-ins
    spk_rows, rec = [], collections.defaultdict(lambda: [0, 0])
    with open(a.spikeins) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            x = dict(zip(head, line.rstrip('\n').split('\t')))
            grp = x['group']
            if grp not in a.groups:
                continue
            w = {'ssu': x['ssu_read'], 'its': x['its_read'], 'lsu': x['lsu_read']}
            if not all(w[m] in with_hits for m in MARKERS):
                continue
            held = set(w.values())
            label = x['label_genus'] if x['kind'] == 'swap' else x['genus_x']
            calls, (f1, f2) = score(w, held, grp, label)
            caught = f1 if x['kind'] != 'swap' else f2
            typ = 'chimera' if x['kind'] != 'swap' else 'swap'
            rec[(grp, typ)][0] += caught
            rec[(grp, typ)][1] += 1
            spk_rows.append([x['spike_id'], grp, x['kind'], label] +
                            [calls[m][0] or '' for m in MARKERS] + [int(calls[m][3]) for m in MARKERS] +
                            [int(f1), int(f2), int(caught)])

    # precision sample for the tree check
    rng = random.Random(SEED)
    sample = []
    for typ in ('F1', 'F2'):
        pool = sorted(real[typ])
        pick = pool if len(pool) <= N_PREC else rng.sample(pool, N_PREC)
        for r in sorted(pick):
            for m in MARKERS:
                idn, tg = hits[m][r]
                top = [(t, float(i)) for i, t in zip(idn, tg) if t != r][:N_HITS_TREE]
                sample.append([typ, group[r], r, ge[r] or '', m] + [f'{t}:{i:.2f}' for t, i in top])

    def write(fn, head, rows):
        with open(os.path.join(a.out, fn), 'w') as fh:
            fh.write('\t'.join(head) + '\n')
            for r in rows:
                fh.write('\t'.join(map(str, r)) + '\n')
    mh = []
    for m in MARKERS:
        mh += [f'{m}_genus', f'{m}_best', f'{m}_lead', f'{m}_confident']
    write('flags.tsv', ['group', 'read', 'label_genus'] + mh + ['F1', 'F2'], real_rows)
    write('spikeins_scored.tsv', ['spike_id', 'group', 'kind', 'label_genus'] + [f'{m}_genus' for m in MARKERS] +
          [f'{m}_confident' for m in MARKERS] + ['F1', 'F2', 'caught'], spk_rows)
    write('thresholds.tsv', ['group'] + [f't5_{m}' for m in MARKERS],
          [[g] + [f'{t5[g][m]:.4f}' for m in MARKERS] for g in a.groups])
    write('precision_sample.tsv', ['flag', 'group', 'read', 'label_genus', 'marker', 'best_hits...'], sample)
    rows = []
    for typ in ('chimera', 'swap'):
        tot = [0, 0]
        for g in a.groups:
            c, n = rec[(g, typ)]
            tot[0] += c
            tot[1] += n
            rows.append([g, typ, c, n, f'{100 * c / n:.1f}' if n else ''])
        rows.append(['all', typ, tot[0], tot[1], f'{100 * tot[0] / tot[1]:.1f}' if tot[1] else ''])
    write('recall.tsv', ['group', 'spike_type', 'caught', 'spikeins', 'recall_pct'], rows)
    fl = [[g, sum(1 for x in real_rows if x[0] == g), sum(1 for r in real['F1'] if group[r] == g),
           sum(1 for r in real['F2'] if group[r] == g)] for g in a.groups]
    write('flag_counts.tsv', ['group', 'reads', 'F1', 'F2'], fl)
    with open(os.path.join(a.out, 'manifest.txt'), 'w') as fh:
        fh.write(f'date\t{datetime.datetime.now().isoformat(timespec="seconds")}\n')
        fh.write(f'script\tscore_k1.py\t{sha256(os.path.abspath(__file__))}\n')
        for p in [os.path.join(a.search, f'{m}.hits.tsv.gz') for m in MARKERS] + [a.spikeins,
                                                                                 os.path.join(a.folds, 'units.tsv')]:
            fh.write(f'input\t{p}\t{sha256(p)}\n')
        fh.write(f'seed\t{SEED}\nlead\t{LEAD}\nranks\t{TOP}\n')
    print('recall (pooled): ' + ', '.join(f'{r[1]} {r[4]}%' for r in rows if r[0] == 'all'))
    print(f'real flags: F1 {len(real["F1"])}, F2 {len(real["F2"])}; precision sample written; see {a.out}')


if __name__ == '__main__':
    main()
