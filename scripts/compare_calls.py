#!/usr/bin/env python3
"""Does the MMseqs2 prefilter cap change what the calls see? Compare two search_markers.py runs of the same
queries that differ only in --max-seqs, on the quantities the protocol's calls are built from, each read off the
top 300 hits (section 5) after dropping self hits:

  best       best identity, and the unique best species name (section 6: the species of the best hit, if every
             hit within 0.2 points of it has that species; otherwise none)
  conspec    identity to the nearest other read of the query's own species (what the thresholds t_m(q) are
             computed from)
  sp_out     after dropping the query's species (species leave-out): best identity and unique best species name
  g_out      after dropping the query's genus (genus leave-out): best identity and unique best genus name

A query agrees when run A gives the same values as run B (identities to 4 decimals). Raw hit lists are not
compared: conserved SSU and LSU have thousands of targets within a few tenths of a point, so the targets that fill
the top 300 change with the cap even when no call input does.

Names follow protocol section 2 (as eukaryome_inventory.py): a named genus has no placeholder marks; a named
species is a named genus plus a lower-case epithet.

  python3 scripts/compare_calls.py --fasta ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta \\
      search/pilot search_ms10k/pilot
"""
import argparse
import collections
import gzip
import re

MARKERS = ['ssu', 'its', 'lsu']
TOP = 300
TIE = 0.2
TAX = re.compile(r'(?:^|[;|,])\s*([kpcofgs])__([^;|,]*)')
EPITHET = re.compile(r'^[a-z][a-z-]+$')
BINOMIAL = re.compile(r'^[A-Z][a-z]+[_ ][a-z][a-z-]+$')
PLACEHOLDER = re.compile(r'(?i)(^|[_ .])(sp|spp|cf|aff|nov|unidentified|uncultured|unclassified|incertae|'
                         r'sedis|environmental|metagenome)([_ .]|$)|\d|\.')


def read_id(header):
    return re.split(r'[\s;]', header, maxsplit=1)[0]


def taxonomy(path):
    sp, ge = {}, {}
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            h = line[1:].rstrip('\n')
            t = {}
            for r, v in TAX.findall(h):
                t.setdefault(r, v.strip())
            g, s = t.get('g', ''), t.get('s', '')
            named_g = bool(g) and g[0].isupper() and not PLACEHOLDER.search(g)
            if s and EPITHET.match(s) and named_g:
                s = f'{g}_{s}'
            r = read_id(h)
            ge[r] = g if named_g else None
            sp[r] = s if (s and BINOMIAL.match(s) and not PLACEHOLDER.search(s)) else None
    return sp, ge


def load(path):
    hits = collections.defaultdict(list)                       # query -> [(identity, target)] in rank order
    with gzip.open(path, 'rt') as fh:
        head = fh.readline().rstrip('\n').split('\t')
        iq, it, ii = head.index('query'), head.index('target'), head.index('identity')
        for line in fh:
            v = line.rstrip('\n').split('\t')
            q, t = read_id(v[iq]), read_id(v[it])
            if q == t or len(hits[q]) >= TOP:               # self hit (the 'self' column of early runs missed it)
                continue
            hits[q].append((round(float(v[ii]), 4), t))
    return hits


def unique_best(hl, name):
    if not hl:
        return None, None
    best = hl[0][0]
    n = name.get(hl[0][1])
    if n is None or any(name.get(t) != n for idn, t in hl if idn >= best - TIE):
        return best, None
    return best, n


def summary(q, hl, sp, ge):
    out = {'best': unique_best(hl, sp)}
    s, g = sp.get(q), ge.get(q)
    if s:
        con = [idn for idn, t in hl if sp.get(t) == s]
        out['conspec'] = con[0] if con else None
        out['sp_out'] = unique_best([(i, t) for i, t in hl if sp.get(t) != s], sp)
    if g:
        out['g_out'] = unique_best([(i, t) for i, t in hl if ge.get(t) != g], ge)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('run_a')
    ap.add_argument('run_b')
    a = ap.parse_args()
    sp, ge = taxonomy(a.fasta)
    group = {}
    with open(f'{a.run_b}/queries.tsv') as fh:
        next(fh)
        for line in fh:
            r, g = line.split('\t')[:2]
            group[r] = g
    keys = ['best', 'conspec', 'sp_out', 'g_out']
    print('marker\tgroup\tqueries\t' + '\t'.join(f'{k}_agree' for k in keys) + '\tall_agree')
    worst = 100.0
    for m in MARKERS:
        A, B = load(f'{a.run_a}/{m}.hits.tsv.gz'), load(f'{a.run_b}/{m}.hits.tsv.gz')
        st = collections.defaultdict(collections.Counter)
        for q in B:
            sa, sb = summary(q, A.get(q, []), sp, ge), summary(q, B[q], sp, ge)
            c = st[group.get(q, '?')]
            c['queries'] += 1
            ok_all = True
            for k in keys:
                if k in sb:
                    c[k + '_n'] += 1
                    same = sa.get(k) == sb[k]
                    c[k] += same
                    ok_all &= same
            c['all'] += ok_all
        for g in sorted(st):
            c = st[g]
            cells = [f"{c[k]}/{c[k + '_n']}" if c[k + '_n'] else '-' for k in keys]
            pct = 100 * c['all'] / c['queries']
            worst = min(worst, pct)
            print(f"{m}\t{g}\t{c['queries']}\t" + '\t'.join(cells) + f"\t{c['all']} ({pct:.1f}%)")
    print(f'lowest share of queries with every call input the same: {worst:.1f}%')


if __name__ == '__main__':
    main()
