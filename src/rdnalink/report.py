"""One row per read: what each marker says, whether the markers agree, and flags.

For each marker the read's hits (best --top ranks, the read itself left out) give:
  best        identity of the best hit
  species     unique best species: the species of the best hit if every hit within 0.2 points of it carries that
              same named species, else none (a hit to a read without a named species carries no name)
  genus       the same for genera
  lead        best identity minus the best hit of any other genus (named or not); if no other genus is among the
              hits, best minus the last rank when the list is full, unlimited when it is not
  confident   genus named, best >= the marker's --min-identity and lead >= --lead

Agreement and flags:
  linked_species / linked_genus   the name all three markers give, when they agree and each best >= --min-identity
  its_species                     ITS alone: its unique best species when best >= the ITS --min-identity
  conflict                        two confident markers name different genera (a possible chimera)
  label_conflict                  all three markers confidently name one genus other than the read's own label
  paralog                         short5.8S when the ITS carries the short-5.8S hairpin and the read is placed in
                                  Rhizophagus; short5.8S_motif when it carries it but is placed elsewhere or nowhere

These are the definitions of the rdnalink test (protocol sections 6 and 9, note 2). The test calibrates its
thresholds per group; until its results are in, --min-identity defaults to 0 (no threshold) and --lead to 1 point.
"""
import collections
import gzip
import math
import os

from .motifs import paralog_flag
from .seqio import read_fasta, read_id
from .taxonomy import names

MARKERS = ('ssu', 'its', 'lsu')
TIE = 0.2
EPS = 1e-9


def load_hits(path, top, drop_self=True):
    """query -> [(identity, target)] in rank order (ranks <= top), and query -> whether its list was full."""
    hits, full = collections.defaultdict(list), {}
    with gzip.open(path, 'rt') as fh:
        head = fh.readline().rstrip('\n').split('\t')
        iq, it, ir, ii = (head.index(c) for c in ('query', 'target', 'rank', 'identity'))
        for line in fh:
            v = line.rstrip('\n').split('\t')
            if int(v[ir]) > top:
                continue
            q, t = v[iq], v[it]
            full[q] = int(v[ir]) >= top
            if drop_self and q == t:
                continue
            hits[q].append((float(v[ii]), t))
    return hits, full


def unique_best(hl, name):
    if not hl:
        return None, None
    best, nm = hl[0][0], name.get(hl[0][1])
    if nm is None:
        return None, best
    for i, t in hl[1:]:
        if i < best - TIE - EPS:
            break
        if name.get(t) != nm:
            return None, best
    return nm, best


def genus_lead(hl, genus, full):
    if not hl:
        return None, None, None
    best, g0 = hl[0][0], genus.get(hl[0][1])
    for i, t in hl[1:]:
        if genus.get(t) != g0:
            return g0, best, best - i
    return g0, best, (best - hl[-1][0]) if full else math.inf


def report(hits_dir, query_prefix, ref_prefix, out_tsv, top=300, lead=1.0, min_identity=None, drop_self=True):
    min_identity = min_identity or {}
    sp, ge = {}, {}
    for m in MARKERS:
        for h, _ in read_fasta(f'{ref_prefix}.{m}.fasta'):
            k = read_id(h)
            if k not in sp:
                s, g, _ = names(h)
                sp[k], ge[k] = s, g
    sp = {k: v for k, v in sp.items() if v}
    ge = {k: v for k, v in ge.items() if v}
    labels, its_seq, order = {}, {}, []
    for m in MARKERS:
        for h, s in read_fasta(f'{query_prefix}.{m}.fasta'):
            k = read_id(h)
            if k not in labels:
                labels[k] = names(h)
                order.append(k)
            if m == 'its':
                its_seq[k] = s
    H, F = {}, {}
    for m in MARKERS:
        if not os.path.exists(f'{hits_dir}/{m}.hits.tsv.gz'):
            raise SystemExit(f'no {hits_dir}/{m}.hits.tsv.gz; run rdnalink search first')
        H[m], F[m] = load_hits(f'{hits_dir}/{m}.hits.tsv.gz', top, drop_self)
    cols = ['read', 'label_species', 'label_genus']
    for m in MARKERS:
        cols += [f'{m}_best', f'{m}_species', f'{m}_genus', f'{m}_lead', f'{m}_confident']
    cols += ['linked_species', 'linked_genus', 'its_species', 'conflict', 'label_conflict', 'paralog']
    counts = collections.Counter()
    with open(out_tsv, 'w') as fo:
        fo.write('\t'.join(cols) + '\n')
        for r in order:
            lab_sp, lab_g, _ = labels[r]
            row, call = [r, lab_sp or '', lab_g or ''], {}
            for m in MARKERS:
                hl = H[m].get(r, [])
                s_nm, best = unique_best(hl, sp)
                g_nm, _ = unique_best(hl, ge)
                g0, _, ld = genus_lead(hl, ge, F[m].get(r, False))
                ok_id = best is not None and best >= min_identity.get(m, 0.0) - EPS
                conf = g0 is not None and ok_id and ld is not None and ld >= lead - EPS
                call[m] = (s_nm, g_nm, best, ok_id, g0 if conf else None)
                row += ['' if best is None else f'{best:.2f}', s_nm or '', g_nm or '',
                        '' if ld is None else ('inf' if math.isinf(ld) else f'{ld:.2f}'), int(conf)]
            lsp = {call[m][0] for m in MARKERS}
            lg = {call[m][1] for m in MARKERS}
            all_ok = all(call[m][3] for m in MARKERS)
            linked_sp = next(iter(lsp)) if len(lsp) == 1 and None not in lsp and all_ok else ''
            linked_g = next(iter(lg)) if len(lg) == 1 and None not in lg and all_ok else ''
            its_sp = call['its'][0] if call['its'][0] and call['its'][3] else ''
            conf = {m: call[m][4] for m in MARKERS if call[m][4]}
            conflict = ''
            if len(set(conf.values())) >= 2:
                conflict = ';'.join(f'{m}:{g}' for m, g in conf.items())
            label_conflict = ''
            if lab_g and len(conf) == 3 and len(set(conf.values())) == 1 and next(iter(conf.values())) != lab_g:
                label_conflict = next(iter(conf.values()))
            placed = linked_g or call['its'][1] or ''
            par = paralog_flag(its_seq.get(r, ''), placed)
            row += [linked_sp, linked_g, its_sp, conflict, label_conflict, par]
            fo.write('\t'.join(map(str, row)) + '\n')
            counts['reads'] += 1
            counts['linked_species'] += bool(linked_sp)
            counts['its_species'] += bool(its_sp)
            counts['conflict'] += bool(conflict)
            counts['label_conflict'] += bool(label_conflict)
            if par:
                counts[f'paralog:{par}'] += 1
    return counts
