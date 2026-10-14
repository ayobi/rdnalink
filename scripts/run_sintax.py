#!/usr/bin/env python3
"""Phase 0, protocol section 6, B2: vsearch SINTAX on the ITS window, one run per fold against that fold's ITS
reference (every ITS window of the long-read set except the fold's held-out reads).

Folds are those of make_folds.py and score_k2.py: species leave-out (all reads of the species leave the
reference), genus leave-out (all reads of the genus), and the 5 known-kept folds of both splits (original;
assembly, where every copy of a fold read's assembly also leaves). Reference labels are EUKARYOME's own
(tax=d:,p:,c:,o:,f:,g:,s:<genus>_<epithet>); placeholder names stay as labels, and score_k2.py counts a call to
one as no name (section 2). vsearch runs with --threads 1 and --randseed 20261007, so each run is reproducible;
--jobs runs folds in parallel. A finished fold is kept, so a stopped run resumes.

  python3 scripts/run_sintax.py --windows windows --fasta ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta \\
      --folds folds --queries search/queries.tsv --out sintax --jobs 24 [--vsearch vsearch] [--drop-reads FILE]

Writes sintax/<group>/<kind>__<unit>.tsv (vsearch --tabbedout) and sintax/manifest.txt.
"""
import argparse
import concurrent.futures
import datetime
import hashlib
import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_calls import read_id  # noqa: E402

SEED = 20261007
GROUPS = ['Glomeromycota', 'Agaricales', 'Hypocreales', 'Rhabditida']
TAX = re.compile(r'(?:^|[;|,])\s*([kpcofgs])__([^;|,]*)')
DUP = re.compile(r'_dup\d+$')
RANKS = [('k', 'd'), ('p', 'p'), ('c', 'c'), ('o', 'o'), ('f', 'f'), ('g', 'g'), ('s', 's')]


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def clean(x):
    return re.sub(r'[,:;\s]+', '_', x.strip()) or 'NA'


def labels(path):
    """read -> SINTAX tax string, from the long-read FASTA headers"""
    import gzip
    out = {}
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            h = line[1:].rstrip('\n')
            t = {}
            for r, v in TAX.findall(h):
                t.setdefault(r, v.strip())
            parts = []
            for r, s in RANKS:
                v = t.get(r, '')
                if not v:
                    continue
                if r == 's':
                    v = f"{t.get('g', 'NA')}_{v}"
                parts.append(f'{s}:{clean(v)}')
            out[read_id(h)] = ','.join(parts)
    return out


def read_windows(path):
    seqs, cur, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if cur is not None:
                    seqs[cur] = ''.join(buf)
                cur, buf = read_id(line[1:].rstrip('\n')), []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = ''.join(buf)
    return seqs


SEQS, LAB = {}, {}            # set before the worker pool forks, so workers share them


def run_fold(job):
    vsearch, out, queries, held, tmpdir = job
    seqs, lab = SEQS, LAB
    if os.path.exists(out):
        return out, 'kept'
    with tempfile.TemporaryDirectory(dir=tmpdir) as d:
        q, ref = os.path.join(d, 'q.fa'), os.path.join(d, 'ref.fa')
        with open(q, 'w') as fh:
            for r in queries:
                fh.write(f'>{r}\n{seqs[r]}\n')
        with open(ref, 'w') as fh:
            for r, s in seqs.items():
                if r not in held:
                    fh.write(f'>{r};tax={lab[r]}\n{s}\n')
        part = out + '.part'
        subprocess.run([vsearch, '--sintax', q, '--db', ref, '--tabbedout', part, '--strand', 'plus',
                        '--threads', '1', '--randseed', str(SEED), '--quiet'], check=True)
        os.replace(part, out)
    return out, 'done'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--windows', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--folds', required=True)
    ap.add_argument('--queries', required=True, help='search/queries.tsv')
    ap.add_argument('--out', required=True)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--vsearch', default='vsearch')
    ap.add_argument('--drop-reads')
    ap.add_argument('--groups', nargs='*', default=GROUPS)
    ap.add_argument('--tmp', default=None, help='scratch directory for the per-fold references (default: --out)')
    a = ap.parse_args()
    drop = set()
    if a.drop_reads:
        drop = {line.split()[0] for line in open(a.drop_reads) if line.strip()}
    lab = labels(a.fasta)
    seqs = {r: s for r, s in read_windows(os.path.join(a.windows, 'win.full.fasta')).items() if r not in drop}
    SEQS.update(seqs)
    LAB.update(lab)
    print(f'{len(seqs)} ITS windows in the reference', flush=True)
    sp_reads, ge_reads, copies = {}, {}, {}
    for r, t in lab.items():
        g = re.search(r'(?:^|,)g:([^,]*)', t)
        s = re.search(r'(?:^|,)s:([^,]*)', t)
        if g:
            ge_reads.setdefault(g.group(1), set()).add(r)
        if s:
            sp_reads.setdefault(s.group(1), set()).add(r)
        copies.setdefault(DUP.sub('', r), set()).add(r)
    units, kk = [], []
    with open(os.path.join(a.folds, 'units.tsv')) as fh:
        next(fh)
        units = [line.rstrip('\n').split('\t')[:3] for line in fh]
    with open(os.path.join(a.folds, 'known_kept.tsv')) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        kk = [dict(zip(head, line.rstrip('\n').split('\t'))) for line in fh]
    # the queries of a leave-out fold: the group's query reads of that species / genus
    group = {}
    with open(a.queries) as fh:
        next(fh)
        for line in fh:
            r, g = line.rstrip('\n').split('\t')[:2]
            if r not in drop:
                group[r] = g
    tag = {r: lab.get(r, '') for r in group}
    jobs = []
    tmp = a.tmp or a.out
    for grp in a.groups:
        os.makedirs(os.path.join(a.out, grp), exist_ok=True)
        q = [r for r in group if group[r] == grp]
        for g_, rank, nm in units:
            if g_ != grp:
                continue
            # species units are 'Genus_epithet' (compare_calls names); SINTAX labels use the same form
            if rank == 'species':
                qs = [r for r in q if re.search(rf'(?:^|,)s:{re.escape(nm)}(?:,|$)', tag[r])]
                held = sp_reads.get(nm, set())
            else:
                qs = [r for r in q if re.search(rf'(?:^|,)g:{re.escape(nm)}(?:,|$)', tag[r])]
                held = ge_reads.get(nm, set())
            if qs:
                out = os.path.join(a.out, grp, f'{rank}__{clean(nm)}.tsv')
                jobs.append((a.vsearch, out, sorted(qs), held | drop, tmp))
        for split in ('original', 'assembly'):
            for k in range(5):
                qs = sorted(x['read'] for x in kk if x['group'] == grp and int(x[f'fold_{split}']) == k)
                held = set(qs)
                if split == 'assembly':
                    for r in qs:
                        held |= copies.get(DUP.sub('', r), set())
                if qs:
                    out = os.path.join(a.out, grp, f'kept_{split}__{k}.tsv')
                    jobs.append((a.vsearch, out, qs, held | drop, tmp))
    print(f'{len(jobs)} folds to run', flush=True)
    done = 0
    ctx = __import__('multiprocessing').get_context('fork')
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.jobs, mp_context=ctx) as ex:
        for out, state in ex.map(run_fold, jobs):
            done += 1
            if done % 25 == 0 or done == len(jobs):
                print(f'  {done}/{len(jobs)} folds ({state}: {out})', flush=True)
    ver = subprocess.run([a.vsearch, '--version'], capture_output=True, text=True)
    with open(os.path.join(a.out, 'manifest.txt'), 'w') as fh:
        fh.write(f'date\t{datetime.datetime.now().isoformat(timespec="seconds")}\n')
        fh.write(f'script\trun_sintax.py\t{sha256(os.path.abspath(__file__))}\n')
        fh.write(f'vsearch\t{(ver.stderr or ver.stdout).splitlines()[0]}\n')
        fh.write(f'parameters\t--sintax --strand plus --threads 1 --randseed {SEED}\n')
        fh.write(f'input\t{os.path.join(a.windows, "win.full.fasta")}\t'
                 f'{sha256(os.path.join(a.windows, "win.full.fasta"))}\n')
        for f in ('units.tsv', 'known_kept.tsv'):
            fh.write(f'input\t{os.path.join(a.folds, f)}\t{sha256(os.path.join(a.folds, f))}\n')
        fh.write(f'folds\t{len(jobs)}\ndrop_reads\t{a.drop_reads}\t{len(drop)}\n')


if __name__ == '__main__':
    main()
