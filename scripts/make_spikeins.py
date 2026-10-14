#!/usr/bin/env python3
"""Phase 0, protocol section 9: build the K1 spike-ins before anything is scored (note 2, section 5).

Per test group, with one random generator seeded 20261007 and the groups taken in the order of section 3:

  chim_a  200 chimeras: SSU + ITS of read X, LSU of read Y
  chim_b  200 chimeras: SSU of read X, ITS + LSU of read Y
  swap    200 label swaps: read X with its genus label replaced by another named genus of the same order

X and Y are query reads of the group (search/queries.tsv) with a named genus (section 2) in a named order; Y is in
the same order as X and in another named genus. In each set the X reads are drawn without replacement, and each Y
is drawn after its X. A swap's new genus is drawn from the named genera of X's order in the whole long-read set.

No search is needed: every window of a spike-in is the window of a query read, whose hits are already in the
table. When a spike-in is scored, its source reads are dropped from its hits (note 2, section 3).

  python3 scripts/make_spikeins.py --fasta ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta \\
      --queries search/queries.tsv --out spikeins

The FASTA can also be its headers alone (grep '>' ... | gzip), which gives the same output.
Writes spikeins/spikeins.tsv and spikeins/manifest.txt (inputs, seed, counts).
"""
import argparse
import collections
import datetime
import gzip
import hashlib
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_calls import PLACEHOLDER, read_id, taxonomy  # noqa: E402

SEED = 20261007
GROUPS = ['Glomeromycota', 'Agaricales', 'Hypocreales', 'Rhabditida']
N = 200
TAX = re.compile(r'(?:^|[;|,])\s*([kpcofgs])__([^;|,]*)')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def orders(path):
    o = {}
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as fh:
        for line in fh:
            if line.startswith('>'):
                h = line[1:].rstrip('\n')
                t = {}
                for r, v in TAX.findall(h):
                    t.setdefault(r, v.strip())
                v = t.get('o', '')
                o[read_id(h)] = v if v and v[0].isupper() and not PLACEHOLDER.search(v) else None
    return o


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True, help='long-read FASTA, or its headers (gzip allowed)')
    ap.add_argument('--queries', required=True, help='search/queries.tsv')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    _, ge = taxonomy(a.fasta)
    order = orders(a.fasta)
    group = {}
    with open(a.queries) as fh:
        next(fh)
        for line in fh:
            r, g = line.rstrip('\n').split('\t')[:2]
            group[r] = g
    # named genera of each named order, over the whole long-read set (for label swaps)
    set_genera = collections.defaultdict(set)
    for r, g in ge.items():
        if g and order.get(r):
            set_genera[order[r]].add(g)
    rng = random.Random(SEED)
    rows, counts = [], []
    for grp in GROUPS:
        reads = sorted(r for r, g in group.items() if g == grp and ge.get(r) and order.get(r))
        by_order = collections.defaultdict(list)
        for r in reads:
            by_order[order[r]].append(r)
        # X needs a Y: another named genus of its order among the group's query reads
        gen_in = {o: {ge[r] for r in rs} for o, rs in by_order.items()}
        x_chim = [r for r in reads if len(gen_in[order[r]]) >= 2]
        x_swap = [r for r in reads if len(set_genera[order[r]]) >= 2]
        for kind in ('chim_a', 'chim_b'):
            xs = rng.sample(x_chim, N)
            for i, x in enumerate(xs, 1):
                ys = [r for r in by_order[order[x]] if ge[r] != ge[x]]
                y = rng.choice(ys)
                ssu, its, lsu = (x, x, y) if kind == 'chim_a' else (x, y, y)
                rows.append([f'{grp}_{kind}_{i:03d}', grp, kind, ssu, its, lsu, order[x], ge[x], ge[y], ge[x]])
        xs = rng.sample(x_swap, N)
        for i, x in enumerate(xs, 1):
            new = rng.choice(sorted(set_genera[order[x]] - {ge[x]}))
            rows.append([f'{grp}_swap_{i:03d}', grp, 'swap', x, x, x, order[x], ge[x], '', new])
        counts.append((grp, len(reads), len(x_chim), len(x_swap)))
    os.makedirs(a.out, exist_ok=True)
    out = os.path.join(a.out, 'spikeins.tsv')
    with open(out, 'w') as fh:
        fh.write('spike_id\tgroup\tkind\tssu_read\tits_read\tlsu_read\torder\tgenus_x\tgenus_y\tlabel_genus\n')
        for r in rows:
            fh.write('\t'.join(r) + '\n')
    with open(os.path.join(a.out, 'manifest.txt'), 'w') as fh:
        fh.write(f'date\t{datetime.datetime.now().isoformat(timespec="seconds")}\n')
        fh.write(f'script\tmake_spikeins.py\t{sha256(os.path.abspath(__file__))}\n')
        fh.write(f'seed\t{SEED}\n')
        fh.write(f'input\t{a.queries}\t{sha256(a.queries)}\n')
        fh.write(f'output\t{out}\t{sha256(out)}\n')
        fh.write('group\tquery_reads_named_genus_and_order\teligible_chimera_x\teligible_swap_x\n')
        for c in counts:
            fh.write('\t'.join(map(str, c)) + '\n')
    print(f'{len(rows)} spike-ins written to {out}; sha256 {sha256(out)}')


if __name__ == '__main__':
    main()
