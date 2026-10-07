#!/usr/bin/env python3
"""Phase 0, step 1 (v3): inventory a EUKARYOME FASTA (e.g. the v2.1 long-read SSU-ITS-LSU set).

Counts sequences, genera and species per taxonomic group so the Phase 0 test groups can be chosen
before any test runs. v2 reads taxon values up to the next ';' (spaces allowed), separates named
species and genera from placeholders (EUKARYOME's 'unclassified' and dotted codes such as
'Densosporales.fam01.gen01'), joins EUKARYOME's epithet-only species field to its genus, and ranks groups at genus level as well as species level, because
the long-read set is mostly environmental. Standard library only.

  python3 eukaryome_inventory.py General_EUK_longread_v2.1.fasta --out inventory

Writes to the --out folder:
  headers_sample.txt   the first 20 headers
  lengths.tsv          length summary
  species_values.tsv   the 60 commonest species-field values other than 'unclassified', with counts
  groups.tsv           per phylum, class and order: sequences, genera, usable genera, named genera,
                       named species, usable named species, median length
  candidates.tsv       orders ranked as Phase 0 test groups (Glomeromycota listed separately)
  report.txt           a short summary to paste back
"""
import argparse, os, re, statistics, sys
from collections import Counter, defaultdict

TAX = re.compile(r'(?:^|[;|,])\s*([kpcofgs])__([^;|,]*)')
BINOMIAL = re.compile(r'^[A-Z][a-z]+[_ ][a-z][a-z-]+$')
EPITHET = re.compile(r'^[a-z][a-z-]+$')
PLACEHOLDER = re.compile(r'(?i)(^|[_ .])(sp|spp|cf|aff|nov|unidentified|uncultured|unclassified|incertae|'
                         r'sedis|environmental|metagenome)([_ .]|$)|\d|\.')


def read_fasta(path):
    head, n = None, 0
    with open(path, 'r', errors='replace') as fh:
        for line in fh:
            if line.startswith('>'):
                if head is not None:
                    yield head, n
                head, n = line[1:].rstrip('\n'), 0
            else:
                n += len(line.strip())
    if head is not None:
        yield head, n


def parse_tax(header):
    tax = {}
    for r, v in TAX.findall(header):
        v = v.strip()
        if v and r not in tax:
            tax[r] = v
    return tax


def named_species(s):
    return bool(s) and bool(BINOMIAL.match(s)) and not PLACEHOLDER.search(s)


def named_genus(g):
    return bool(g) and g[0].isupper() and not PLACEHOLDER.search(g)


def pct(xs, q):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(round(q / 100 * (len(xs) - 1))))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fasta')
    ap.add_argument('--out', default='inventory')
    ap.add_argument('--min-seqs', type=int, default=3, help='sequences a taxon needs to count as usable')
    a = ap.parse_args()
    k = a.min_seqs
    os.makedirs(a.out, exist_ok=True)

    sample, lengths, parsed = [], [], 0
    sp_values = Counter()
    sp_kind = Counter()
    rank_of = {'p': {}, 'c': {}, 'o': {}}                         # taxon -> parent lineage
    seqs = {r: Counter() for r in 'pco'}
    lens = {r: defaultdict(list) for r in 'pco'}
    genus_n = {r: defaultdict(Counter) for r in 'pco'}           # taxon -> genus -> n
    species_n = {r: defaultdict(Counter) for r in 'pco'}         # taxon -> named species -> n

    for head, n in read_fasta(a.fasta):
        if len(sample) < 20:
            sample.append(head)
        lengths.append(n)
        t = parse_tax(head)
        if t:
            parsed += 1
        L = {r: t.get(r, '') for r in 'kpcofgs'}
        s = L['s']
        # EUKARYOME stores the epithet alone ('mosseae'): join it to a named genus
        if s and s.lower() != 'unclassified' and EPITHET.match(s) and named_genus(L['g']):
            s = f"{L['g']}_{s}"
        if not s or s.lower() == 'unclassified':
            sp_kind['unclassified'] += 1
        else:
            sp_values[s] += 1
            sp_kind['named' if named_species(s) else 'other'] += 1
        for r in 'pco':
            key = L[r] or '(none)'
            seqs[r][key] += 1
            lens[r][key].append(n)
            rank_of[r].setdefault(key, (L['p'], L['c']))
            if L['g']:
                genus_n[r][key][L['g']] += 1
            if named_species(s):
                species_n[r][key][s] += 1

    with open(os.path.join(a.out, 'headers_sample.txt'), 'w') as fh:
        fh.write('\n'.join(sample) + '\n')
    if not lengths:
        sys.exit('no sequences read')
    if parsed == 0:
        sys.exit('no k__/p__/.../s__ taxonomy found in the headers; send headers_sample.txt')

    with open(os.path.join(a.out, 'lengths.tsv'), 'w') as fh:
        fh.write('n\tmin\tp5\tmedian\tp95\tmax\n')
        fh.write(f'{len(lengths)}\t{min(lengths)}\t{pct(lengths, 5)}\t{int(statistics.median(lengths))}\t'
                 f'{pct(lengths, 95)}\t{max(lengths)}\n')
    with open(os.path.join(a.out, 'species_values.tsv'), 'w') as fh:
        fh.write('species_value\tsequences\tnamed\n')
        for v, c in sp_values.most_common(60):
            fh.write(f'{v}\t{c}\t{named_species(v)}\n')

    def stats(r, key):
        g = genus_n[r][key]; sp = species_n[r][key]
        return dict(seqs=seqs[r][key], genera=len(g), genera_ok=sum(c >= k for c in g.values()),
                    named_genera=sum(named_genus(x) for x in g),
                    named_genera_ok=sum(named_genus(x) and c >= k for x, c in g.items()),
                    named_species=len(sp), named_species_ok=sum(c >= k for c in sp.values()),
                    median_len=int(statistics.median(lens[r][key])))

    cols = ['seqs', 'genera', 'genera_ok', 'named_genera', 'named_genera_ok', 'named_species',
            'named_species_ok', 'median_len']
    head_cols = ['sequences', 'genera', f'genera_ge{k}', 'named_genera', f'named_genera_ge{k}',
                 'named_species', f'named_species_ge{k}', 'median_length']
    orders = []
    with open(os.path.join(a.out, 'groups.tsv'), 'w') as fh:
        fh.write('rank\ttaxon\tphylum\tclass\t' + '\t'.join(head_cols) + '\n')
        for r, name in (('p', 'phylum'), ('c', 'class'), ('o', 'order')):
            for key, _ in seqs[r].most_common():
                st = stats(r, key); p, c = rank_of[r][key]
                fh.write(f'{name}\t{key}\t{p}\t{c}\t' + '\t'.join(str(st[x]) for x in cols) + '\n')
                if r == 'o':
                    orders.append((key, p, c, st))

    skip = {'(none)', 'unclassified'}
    cand = [o for o in orders if o[0] not in skip and not PLACEHOLDER.search(o[0])
            and o[1] != 'Glomeromycota' and o[3]['genera_ok'] >= 5]
    cand.sort(key=lambda o: (-o[3]['named_species_ok'], -o[3]['named_genera_ok'], -o[3]['genera_ok']))
    with open(os.path.join(a.out, 'candidates.tsv'), 'w') as fh:
        fh.write('order\tphylum\tclass\t' + '\t'.join(head_cols) + '\n')
        for o, p, c, st in cand[:60]:
            fh.write(f'{o}\t{p}\t{c}\t' + '\t'.join(str(st[x]) for x in cols) + '\n')

    gl = stats('p', 'Glomeromycota') if 'Glomeromycota' in seqs['p'] else None
    report = [
        f'sequences: {len(lengths)}; with parsed taxonomy: {parsed}',
        f'length median {int(statistics.median(lengths))} (5-95%: {pct(lengths, 5)}-{pct(lengths, 95)})',
        f'phyla: {len(seqs["p"])}; classes: {len(seqs["c"])}; orders: {len(seqs["o"])}',
        f'species field: {sp_kind["unclassified"]} unclassified, {sp_kind["named"]} binomial names, '
        f'{sp_kind["other"]} other values (see species_values.tsv)',
    ]
    if gl:
        report.append(f'Glomeromycota: {gl["seqs"]} sequences, {gl["genera"]} genera ({gl["named_genera"]} named, '
                      f'{gl["named_genera_ok"]} named with >= {k} seqs), {gl["named_species"]} named species '
                      f'({gl["named_species_ok"]} with >= {k} seqs)')
    report.append('top candidate orders (named species >= k / named genera >= k / all genera >= k, sequences):')
    report += [f'  {o} ({p}): {st["named_species_ok"]} / {st["named_genera_ok"]} / {st["genera_ok"]}, {st["seqs"]}'
               for o, p, c, st in cand[:15]]
    open(os.path.join(a.out, 'report.txt'), 'w').write('\n'.join(report) + '\n')
    print('\n'.join(report))


if __name__ == '__main__':
    main()
