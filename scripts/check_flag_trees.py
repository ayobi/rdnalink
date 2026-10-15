#!/usr/bin/env python3
"""Phase 0, protocol section 9, K1 precision: is a sampled flag real? (definitions: protocol note 3)

For each flag drawn by score_k1.py (precision_sample.tsv) and each marker: the read's window and the windows of its
20 best remaining hits are aligned with MAFFT (--auto) and a tree is built with IQ-TREE (ModelFinder, 1,000
UFBoot replicates, --seed 20261007, one thread). The tree is read unrooted, as bipartitions:

  placement   the smallest side of a bipartition with UFBoot >= 70 that holds the read and >= 1 other sequence;
              if the other sequences on that side with a named genus all share one genus G (and there is >= 1),
              the marker places the read in G; otherwise it places it nowhere. With fewer than 4 sequences
              (the read and fewer than 3 remaining hits) no tree is built: an unrooted tree that small has no
              internal branch, so the marker places the read nowhere (note 5)
  F1 real     two markers place the read in different named genera
  F2 real     all three markers place the read in named genera other than its label

  python3 scripts/check_flag_trees.py --sample scores/k1/precision_sample.tsv --windows windows \\
      --fasta ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta --out scores/k1_trees --jobs 16

Writes scores/k1_trees/<flag>_<read>_<marker>.{fasta,aln,treefile,...}, placements.tsv, precision.tsv, manifest.txt.
"""
import argparse
import collections
import concurrent.futures
import datetime
import hashlib
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_calls import read_id, taxonomy  # noqa: E402

SEED = 20261007
SUPPORT = 70.0
MIN_TREE = 4        # an unrooted tree of 3 or fewer sequences has no internal branch (note 5)
MARKERS = ['ssu', 'its', 'lsu']
WINDOW_FILES = {'ssu': 'win.ssu.fasta', 'its': 'win.full.fasta', 'lsu': 'win.lsu.fasta'}


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
                    seqs[cur] = ''.join(buf)
                cur, buf = read_id(line[1:].rstrip('\n')), []
            else:
                buf.append(line.strip())
        if cur in keep:
            seqs[cur] = ''.join(buf)
    return seqs


def parse_newick(s):
    """-> list of (leaf set below an internal node, support or None), and the set of all leaves"""
    s = s.strip().rstrip(';')
    stack, splits, leaves = [[]], [], set()
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == '(':
            stack.append([])
            i += 1
        elif c == ',':
            i += 1
        elif c == ')':
            i += 1
            j = i
            while j < n and s[j] not in ',();:':
                j += 1
            label = s[i:j]
            i = j
            if i < n and s[i] == ':':
                i += 1
                while i < n and s[i] not in ',()':
                    i += 1
            kids = stack.pop()
            below = set().union(*kids) if kids else set()
            sup = None
            if label:
                try:
                    sup = float(label.split('/')[-1])
                except ValueError:
                    sup = None
            stack[-1].append(below)
            splits.append((below, sup))
        else:
            j = i
            while j < n and s[j] not in ',():':
                j += 1
            name = s[i:j]
            i = j
            if i < n and s[i] == ':':
                i += 1
                while i < n and s[i] not in ',()':
                    i += 1
            leaves.add(name)
            stack[-1].append({name})
    return splits, leaves


def placement(treefile, read, genus):
    if treefile is None:                    # too few sequences for a tree: no bipartition, places nowhere (note 5)
        return None, 0
    splits, leaves = parse_newick(open(treefile).read())
    best = None
    for below, sup in splits:
        if sup is None or sup < SUPPORT - 1e-9 or not below or below == leaves:
            continue
        side = below if read in below else leaves - below
        if len(side) < 2 or len(side) == len(leaves):
            continue
        if best is None or len(side) < len(best):
            best = side
    if best is None:
        return None, 0
    gs = {genus.get(x) for x in best if x != read and genus.get(x)}
    return (next(iter(gs)) if len(gs) == 1 else None), len(best)


def run_one(job):
    mafft, iqtree, prefix, seqs = job
    tree = prefix + '.treefile'
    if os.path.exists(tree):
        return tree
    with open(prefix + '.fasta', 'w') as fh:
        for r, s in seqs:
            fh.write(f'>{r}\n{s}\n')
    if len(seqs) < MIN_TREE:
        return None
    with open(prefix + '.aln', 'w') as fh:
        subprocess.run([mafft, '--auto', '--quiet', prefix + '.fasta'], stdout=fh, check=True)
    subprocess.run([iqtree, '-s', prefix + '.aln', '-B', '1000', '-T', '1', '--seed', str(SEED), '-pre', prefix,
                    '-quiet', '-redo'], check=True)
    return tree


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sample', required=True)
    ap.add_argument('--windows', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--mafft', default='mafft')
    ap.add_argument('--iqtree', default='iqtree2')
    ap.add_argument('--jobs', type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    _, ge = taxonomy(a.fasta)
    genus = {r: g for r, g in ge.items() if g}
    rows = []
    with open(a.sample) as fh:
        next(fh)
        for line in fh:
            v = line.rstrip('\n').split('\t')
            flag, grp, r, label, m = v[:5]
            hits = [x.rsplit(':', 1)[0] for x in v[5:] if x]
            rows.append((flag, grp, r, label, m, hits))
    need = collections.defaultdict(set)
    for flag, grp, r, label, m, hits in rows:
        need[m].add(r)
        need[m].update(hits)
    seqs = {m: windows(os.path.join(a.windows, WINDOW_FILES[m]), need[m]) for m in MARKERS}
    jobs, keys = [], []
    for flag, grp, r, label, m, hits in rows:
        prefix = os.path.join(a.out, f'{flag}_{r}_{m}')
        jobs.append((a.mafft, a.iqtree, prefix, [(x, seqs[m][x]) for x in [r] + hits]))
        keys.append((flag, grp, r, label, m))
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.jobs) as ex:
        trees = list(ex.map(run_one, jobs))
    place = {}
    with open(os.path.join(a.out, 'placements.tsv'), 'w') as fh:
        fh.write('flag\tgroup\tread\tlabel_genus\tmarker\tplaced_genus\tclade_size\n')
        for (flag, grp, r, label, m), tree in zip(keys, trees):
            g, size = placement(tree, r, genus)
            place[(flag, r, m)] = g
            fh.write(f'{flag}\t{grp}\t{r}\t{label}\t{m}\t{g or ""}\t{size}\n')
    real = collections.defaultdict(lambda: [0, 0])
    reads = sorted({(flag, grp, r, label) for flag, grp, r, label, m in keys})
    with open(os.path.join(a.out, 'flags_checked.tsv'), 'w') as fh:
        fh.write('flag\tgroup\tread\tlabel_genus\t' + '\t'.join(f'{m}_placed' for m in MARKERS) + '\treal\n')
        for flag, grp, r, label in reads:
            pl = [place.get((flag, r, m)) for m in MARKERS]
            if flag == 'F1':
                ok = len({g for g in pl if g}) >= 2
            else:
                ok = all(g and g != label for g in pl)
            real[flag][0] += ok
            real[flag][1] += 1
            fh.write(f'{flag}\t{grp}\t{r}\t{label}\t' + '\t'.join(g or '' for g in pl) + f'\t{int(ok)}\n')
    tot = [sum(x[0] for x in real.values()), sum(x[1] for x in real.values())]
    with open(os.path.join(a.out, 'precision.tsv'), 'w') as fh:
        fh.write('flag\treal\tchecked\tprecision_pct\n')
        for flag in ('F1', 'F2', 'all'):
            c, n = tot if flag == 'all' else real[flag]
            fh.write(f'{flag}\t{c}\t{n}\t{100 * c / n:.1f}\n' if n else f'{flag}\t0\t0\t\n')
    ver = subprocess.run([a.iqtree, '--version'], capture_output=True, text=True).stdout.splitlines()[:1]
    mv = subprocess.run([a.mafft, '--version'], capture_output=True, text=True)
    with open(os.path.join(a.out, 'manifest.txt'), 'w') as fh:
        fh.write(f'date\t{datetime.datetime.now().isoformat(timespec="seconds")}\n')
        fh.write(f'script\tcheck_flag_trees.py\t{sha256(os.path.abspath(__file__))}\n')
        fh.write(f'iqtree\t{ver[0] if ver else ""}\nmafft\t{(mv.stderr or mv.stdout).strip()}\n')
        fh.write(f'input\t{a.sample}\t{sha256(a.sample)}\nseed\t{SEED}\nsupport\t{SUPPORT}\n')
    print(open(os.path.join(a.out, 'precision.tsv')).read())


if __name__ == '__main__':
    main()
