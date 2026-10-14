#!/usr/bin/env python3
"""Phase 0, protocol section 7: fix the test units and the known-kept split before anything is scored (note 2).

Per test group, from the query reads (search/queries.tsv) and the read labels (section 2):
  usable species  named species with >= 3 query reads (species leave-out folds)
  usable genera   named genera with >= 3 query reads (genus leave-out folds)
  known-kept      the query reads of usable species split into 5 folds within each species. Reads identical in
                  all three windows are one unit and go to the same fold. Two splits are written:
                    original   identical reads together (section 7 as written)
                    assembly   also the copies of one genome assembly (ID and ID_dupN) together (note 2, s. 4)
                  Units are shuffled with random.Random('20261007-<group>-<species>') and each goes to the fold
                  with the fewest reads so far (lowest fold number on a tie).

  python3 scripts/make_folds.py --windows windows --fasta ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta \\
      --queries search/queries.tsv --out folds

Writes folds/units.tsv (group, rank, name, query reads) and folds/known_kept.tsv (group, read, species, unit,
fold_original, fold_assembly), and folds/manifest.txt.
"""
import argparse
import collections
import datetime
import hashlib
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_calls import read_id, taxonomy  # noqa: E402

SEED = 20261007
GROUPS = ['Glomeromycota', 'Agaricales', 'Hypocreales', 'Rhabditida']
WINDOW_FILES = {'ssu': 'win.ssu.fasta', 'its': 'win.full.fasta', 'lsu': 'win.lsu.fasta'}
DUP = re.compile(r'_dup\d+$')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def windows(path, keep):
    seqs, cur, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if cur in keep:
                    seqs[cur] = ''.join(buf).upper()
                cur, buf = read_id(line[1:].rstrip('\n')), []
            else:
                buf.append(line.strip())
        if cur in keep:
            seqs[cur] = ''.join(buf).upper()
    return seqs


def split(units, rng_key):
    """units: list of read lists for one species -> {read: fold}"""
    us = sorted(sorted(u) for u in units)
    random.Random(rng_key).shuffle(us)
    size, fold = [0] * 5, {}
    for u in us:
        k = min(range(5), key=lambda j: (size[j], j))
        for r in u:
            fold[r] = k
        size[k] += len(u)
    return fold


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--windows', required=True)
    ap.add_argument('--fasta', required=True, help='long-read FASTA or its headers (gzip allowed)')
    ap.add_argument('--queries', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    sp, ge = taxonomy(a.fasta)
    group = {}
    with open(a.queries) as fh:
        next(fh)
        for line in fh:
            r, g = line.rstrip('\n').split('\t')[:2]
            group[r] = g
    seq = {m: windows(os.path.join(a.windows, f), set(group)) for m, f in WINDOW_FILES.items()}
    missing = [r for r in group if not all(r in seq[m] for m in seq)]
    if missing:
        sys.exit(f'{len(missing)} query reads lack a window, e.g. {missing[:3]}')
    os.makedirs(a.out, exist_ok=True)
    unit_rows, kk_rows, counts = [], [], []
    for grp in GROUPS:
        q = sorted(r for r, g in group.items() if g == grp)
        n_sp = collections.Counter(sp[r] for r in q if sp[r])
        n_ge = collections.Counter(ge[r] for r in q if ge[r])
        use_sp = sorted(s for s, n in n_sp.items() if n >= 3)
        use_ge = sorted(g for g, n in n_ge.items() if n >= 3)
        for s in use_sp:
            unit_rows.append((grp, 'species', s, n_sp[s]))
        for g in use_ge:
            unit_rows.append((grp, 'genus', g, n_ge[g]))
        uset = set(use_sp)
        n_units = {}
        for s in use_sp:
            rs = [r for r in q if sp[r] == s]
            by_seq = collections.defaultdict(list)
            for r in rs:
                by_seq[(seq['ssu'][r], seq['its'][r], seq['lsu'][r])].append(r)
            orig = list(by_seq.values())
            # assembly: merge identical-read units that share an assembly
            parent = {r: r for r in rs}

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x
            for u in orig:
                for r in u[1:]:
                    parent[find(r)] = find(u[0])
            by_asm = collections.defaultdict(list)
            for r in rs:
                by_asm[DUP.sub('', r)].append(r)
            for u in by_asm.values():
                for r in u[1:]:
                    parent[find(r)] = find(u[0])
            asm = collections.defaultdict(list)
            for r in rs:
                asm[find(r)].append(r)
            asm = list(asm.values())
            f_o = split(orig, f'{SEED}-{grp}-{s}')
            f_a = split(asm, f'{SEED}-{grp}-{s}')
            uid = {r: i for i, u in enumerate(sorted(sorted(u) for u in orig)) for r in u}
            for r in rs:
                kk_rows.append((grp, r, s, f'{s}#{uid[r]}', f_o[r], f_a[r]))
            n_units[s] = (len(orig), len(asm))
        counts.append((grp, len(q), len(use_sp), len(use_ge), sum(1 for r in q if sp[r] in uset),
                       sum(x for x, _ in n_units.values()), sum(y for _, y in n_units.values())))
    with open(os.path.join(a.out, 'units.tsv'), 'w') as fh:
        fh.write('group\trank\tname\tquery_reads\n')
        for r in unit_rows:
            fh.write('\t'.join(map(str, r)) + '\n')
    kk = os.path.join(a.out, 'known_kept.tsv')
    with open(kk, 'w') as fh:
        fh.write('group\tread\tspecies\tidentical_unit\tfold_original\tfold_assembly\n')
        for r in kk_rows:
            fh.write('\t'.join(map(str, r)) + '\n')
    with open(os.path.join(a.out, 'manifest.txt'), 'w') as fh:
        fh.write(f'date\t{datetime.datetime.now().isoformat(timespec="seconds")}\n')
        fh.write(f'script\tmake_folds.py\t{sha256(os.path.abspath(__file__))}\n')
        fh.write(f'seed\t{SEED}\ninput\t{a.queries}\t{sha256(a.queries)}\n')
        for m, f in WINDOW_FILES.items():
            p = os.path.join(a.windows, f)
            fh.write(f'input\t{p}\t{sha256(p)}\n')
        for f in ('units.tsv', 'known_kept.tsv'):
            p = os.path.join(a.out, f)
            fh.write(f'output\t{p}\t{sha256(p)}\n')
        fh.write('group\tquery_reads\tusable_species\tusable_genera\tknown_kept_reads\tunits_original\t'
                 'units_assembly\n')
        for c in counts:
            fh.write('\t'.join(map(str, c)) + '\n')
    for c in counts:
        print('\t'.join(map(str, c)))
    print(f'written to {a.out}/')


if __name__ == '__main__':
    main()
