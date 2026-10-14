#!/usr/bin/env python3
"""Phase 0, protocol sections 6-8 (K2), for the hit-table methods B1 (ITS best hit) and L (linked). B2 (SINTAX)
needs its own run per fold (run_sintax.py) and is added with --sintax; without it the decision is provisional.

B2 (--sintax): species S if SINTAX's species is a named species with bootstrap >= c (c from 0.50 to 1.00 in steps
of 0.01); a placeholder label is no name. K2 (section 8) is decided on per-read rates at species level, at known-kept
90% or, where a method cannot reach it, at the highest known-kept that L, B1 and B2 all reach (k2_decision.txt).

Reads the search table (search_markers.py) under protocol note 2:
  - ranks <= 300 (or <= 1,000 for a group that fails the depth check), self hit included in the ranks;
  - a fold drops its held-out reads from those ranks (the query is always one of them); lower ranks never move up.

Folds (section 7), per group, as fixed by make_folds.py:
  species leave-out  every usable species (named, >= 3 query reads): all reads of that species in the whole set
                     leave the reference; its query reads are queried. A call naming any named species is a false
                     known (FK).
  genus leave-out    the same for every usable genus; FK if the genus call or the species call names anything.
  known kept         the 5-fold split of make_folds.py (--split original: identical reads together, as section 7;
                     --split assembly: also the copies of one assembly, every copy of it leaving the reference with
                     its fold, note 2 section 4). Kept (KK) if the call is the read's own species (genus rank: its
                     own genus).

Thresholds t_m(q) (section 6), per fold from training reads only: the q-th lower percentile (numpy linear) of each
training read's identity to its nearest other read of the same named species (genus rank: genus), pooled over the
group's usable species (genera). Training reads are the group's query reads of usable species (genera) that the fold
does not hold out. A read with no such read among its ranks scores the identity of its last rank if its list was
full (an upper bound, note 2), and 0 if it was not (no such read aligns). q runs from 0.5 to 50 in steps of 0.5.

Unique best name (section 6): the name of the best remaining hit, if every remaining hit within 0.2 points of it
carries that same name; a hit to a read without a named species (genus) carries no name.

  python3 scripts/score_k2.py --search search --fasta ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta \\
      --folds folds --out scores/k2_original [--split assembly] [--drop-reads FILE]

Writes scores/k2/: sweep.tsv (every group, rank, method, q: FK and KK, per read and averaged over species),
k2_summary.tsv (FK at known-kept 80/90/95%), k2_decision_rows.tsv, k2_decision.txt, depth_check.tsv, near_ties.tsv, manifest.txt.
"""
import argparse
import collections
import datetime
import gzip
import hashlib
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_calls import read_id, taxonomy  # noqa: E402

MARKERS = ['ssu', 'its', 'lsu']
GROUPS = ['Glomeromycota', 'Agaricales', 'Hypocreales', 'Rhabditida']
TOP, DEEP = 300, 1000
TIE = 0.2
EPS = 1e-9
SEED = 20261007
QS = [round(0.5 * i, 1) for i in range(1, 101)]            # 0.5 .. 50
CS = [round(0.5 + 0.01 * i, 2) for i in range(51)]          # SINTAX bootstrap cut-off 0.50 .. 1.00
KK_LEVELS = [80.0, 90.0, 95.0]
DUP = re.compile(r'_dup\d+$')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load_hits(path, keep, depth, drop):
    """query -> (identities, targets) for ranks <= depth, self included, and query -> targets identical over the
    whole window (identity 100, alignment = both window lengths). Reads in drop leave the lists after the rank
    cut, so lower ranks never move up (Amendment 1, original protocol)."""
    out, same = {}, collections.defaultdict(set)
    rows = collections.defaultdict(list)
    with gzip.open(path, 'rt') as fh:
        head = fh.readline().rstrip('\n').split('\t')
        iq, it, ir, ii, ia, ql, tl = (head.index(c) for c in ('query', 'target', 'rank', 'identity', 'alnlen',
                                                              'qlen', 'tlen'))
        for line in fh:
            v = line.split('\t', 9)
            q = read_id(v[iq])
            if q not in keep or int(v[ir]) > depth:
                continue
            t = sys.intern(read_id(v[it]))
            if t in drop:
                continue
            i = float(v[ii])
            rows[q].append((i, t))
            if i >= 100.0 - EPS and v[ia] == v[ql] == v[tl]:
                same[q].add(t)
    for q, x in rows.items():
        out[q] = (np.array([i for i, _ in x], dtype=np.float64), [t for _, t in x])
    return out, same


def clean(x):
    return re.sub(r'[,:;\s]+', '_', x.strip()) or 'NA'


def read_sintax(path):
    """read -> {'g': (label, bootstrap), 's': (label, bootstrap)} from vsearch --tabbedout"""
    out = {}
    with open(path) as fh:
        for line in fh:
            v = line.rstrip('\n').split('\t')
            d = {}
            for part in v[1].split(',') if len(v) > 1 and v[1] else []:
                m = re.match(r'([a-z]):(.*)\(([0-9.]+)\)$', part)
                if m:
                    d[m.group(1)] = (m.group(2), float(m.group(3)))
            out[read_id(v[0])] = d
    return out


def unique_best(idn, tg, name, held):
    """(name or None, best identity or None, n remaining, within-0.2 set reaches the end) after dropping held."""
    best, nm, ok, n, last_in = None, None, True, 0, False
    for i, t in zip(idn, tg):
        if t in held:
            continue
        n += 1
        if best is None:
            best, nm = i, name.get(t)
            ok = nm is not None
            last_in = True
            continue
        if i >= best - TIE - EPS:
            last_in = True
            if name.get(t) != nm:
                ok = False
        else:
            last_in = False
            break
    return (nm if ok else None), best, n, last_in


def nearest_same(idn, tg, name, own, held, me):
    for i, t in zip(idn, tg):
        if t != me and t not in held and name.get(t) == own:
            return i
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--search', required=True)
    ap.add_argument('--fasta', required=True, help='long-read FASTA or its headers (gzip allowed)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--folds', required=True, help='make_folds.py output directory')
    ap.add_argument('--split', choices=['original', 'assembly'], default='original',
                    help='known-kept split: section 7 as written, or note 2 section 4 (assembly copies together)')
    ap.add_argument('--drop-reads', help='file of read IDs removed from queries and reference (e.g. Amendment 1 '
                                         'original protocol: the barrnap-cut reads)')
    ap.add_argument('--groups', nargs='*', default=GROUPS)
    ap.add_argument('--sintax', help='run_sintax.py output directory (adds B2)')
    ap.add_argument('--depth', type=int, default=0, help='force ranks <= N for every group (testing only)')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    log = open(os.path.join(a.out, 'run.log'), 'w')

    def say(*x):
        msg = ' '.join(str(y) for y in x)
        print(msg, flush=True)
        log.write(msg + '\n')

    sp, ge = taxonomy(a.fasta)
    drop = set()
    if a.drop_reads:
        drop = {line.split()[0] for line in open(a.drop_reads) if line.strip()}
    group = {}
    with open(os.path.join(a.search, 'queries.tsv')) as fh:
        next(fh)
        for line in fh:
            r, g = line.rstrip('\n').split('\t')[:2]
            if r not in drop:
                group[r] = g
    keep = set(group)
    hits = {}
    for m in MARKERS:
        hits[m], _ = load_hits(os.path.join(a.search, f'{m}.hits.tsv.gz'), keep, max(DEEP, a.depth or 0),
                                      drop)
        say(f'{m}: hit lists for {len(hits[m])} queries')
    with_hits = {r for r in keep if all(r in hits[m] for m in MARKERS)}
    if len(with_hits) < len(keep):
        say(f'queries without a hit list on every marker (left out; expected only for a pilot): '
            f'{len(keep) - len(with_hits)}')
    # read sets over the whole long-read set, for removing a species/genus/assembly from the reference
    reads_of_sp, reads_of_ge, copies = collections.defaultdict(set), collections.defaultdict(set), \
        collections.defaultdict(set)
    for r in sp:
        if r in drop:
            continue
        if sp[r]:
            reads_of_sp[sp[r]].add(r)
        if ge[r]:
            reads_of_ge[ge[r]].add(r)
        copies[DUP.sub('', r)].add(r)
    name = {'species': {r: s for r, s in sp.items() if s and r not in drop},
            'genus': {r: g for r, g in ge.items() if g and r not in drop}}

    units = []
    with open(os.path.join(a.folds, 'units.tsv')) as fh:
        next(fh)
        for line in fh:
            g, rk, n = line.rstrip('\n').split('\t')[:3]
            units.append((g, rk, n))
    kk_fold = {}
    with open(os.path.join(a.folds, 'known_kept.tsv')) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        ir, if_ = head.index('read'), head.index(f'fold_{a.split}')
        for line in fh:
            v = line.rstrip('\n').split('\t')
            kk_fold[v[ir]] = int(v[if_])
    sweep_rows, summ_rows, depth_rows, tie_rows = [], [], [], []
    methods = ['B1', 'B2', 'L'] if a.sintax else ['B1', 'L']
    named_sp = {clean(x) for x in name['species'].values()}
    named_ge = {clean(x) for x in name['genus'].values()}
    for grp in a.groups:
        q_all = sorted(r for r in with_hits if group[r] == grp)
        use_sp = sorted(n for g, rk, n in units if g == grp and rk == 'species')
        use_ge = sorted(n for g, rk, n in units if g == grp and rk == 'genus')
        q_sp = [r for r in q_all if sp[r] in set(use_sp)]
        q_ge = [r for r in q_all if ge[r] in set(use_ge)]
        fold_of = {r: kk_fold[r] for r in q_sp}                 # from make_folds.py (note 2)

        # ---- the folds as (kind, unit, query reads, held-out set)
        folds = []
        for s in use_sp:
            folds.append(('species', s, [r for r in q_sp if sp[r] == s], reads_of_sp[s]))
        for g in use_ge:
            folds.append(('genus', g, [r for r in q_ge if ge[r] == g], reads_of_ge[g]))
        for k in range(5):
            qs = [r for r in q_sp if fold_of[r] == k]
            held = set(qs)
            if a.split == 'assembly':
                for r in qs:
                    held |= copies[DUP.sub('', r)]
            folds.append(('kept', str(k), qs, held))

        # ---- depth check (section 5): ranks <= 300, list full, fewer than 10 left
        fail, nfq = 0, 0
        for kind, unit, qs, held in folds:
            for r in qs:
                nfq += 1
                bad = False
                for m in MARKERS:
                    idn, tg = hits[m][r]
                    top = tg[:TOP]
                    if len(top) == TOP and sum(1 for t in top if t not in held) < 10:
                        bad = True
                fail += bad
        pct = 100.0 * fail / nfq if nfq else 0.0
        depth = a.depth or (DEEP if pct > 1.0 else TOP)
        depth_rows.append((grp, nfq, fail, f'{pct:.2f}', depth))
        say(f'{grp}: {len(q_all)} queries; usable species {len(use_sp)}, genera {len(use_ge)}; depth check '
            f'{fail}/{nfq} fold queries fail ({pct:.2f}%) -> ranks <= {depth}')
        H = {m: {r: (hits[m][r][0][:depth], hits[m][r][1][:depth]) for r in q_all} for m in MARKERS}
        full = {m: {r: len(hits[m][r][1]) >= depth for r in q_all} for m in MARKERS}

        # ---- threshold values: identity to the nearest other read of the same name, given held-out set
        def thr_value(m, r, rank, held):
            nm = name[rank]
            own = nm.get(r)
            idn, tg = H[m][r]
            v = nearest_same(idn, tg, nm, own, held, r)
            if v is None:
                v = float(idn[-1]) if full[m][r] and len(idn) else 0.0
            return v
        base = {rank: {m: {} for m in MARKERS} for rank in ('species', 'genus')}
        pool = {'species': q_sp, 'genus': q_ge}
        for rank in base:
            for m in MARKERS:
                for r in pool[rank]:
                    base[rank][m][r] = thr_value(m, r, rank, set())

        def thresholds(rank, held, kind):
            out = {}
            for m in MARKERS:
                if kind == 'kept':
                    vals = [thr_value(m, r, rank, held) for r in pool[rank] if r not in held]
                else:
                    vals = [v for r, v in base[rank][m].items() if r not in held]
                out[m] = np.percentile(np.array(vals), QS) if vals else np.full(len(QS), np.inf)
            return out

        # ---- calls
        res = collections.defaultdict(list)       # (rank, method, kind) -> [(unit, read, ok array over q)]
        ties = collections.Counter()
        for kind, unit, qs, held in folds:
            if not qs:
                continue
            T = {rank: thresholds(rank, held, kind) for rank in ('species', 'genus')}
            for r in qs:
                ub = {}
                for rank in ('species', 'genus'):
                    for m in MARKERS:
                        idn, tg = H[m][r]
                        nm, best, n, last_in = unique_best(idn, tg, name[rank], held)
                        ub[(rank, m)] = (nm, best)
                        if last_in and full[m][r] and n > 0:
                            ties[(kind, m)] += 1
                for rank in ('species', 'genus'):
                    nm_its, b_its = ub[(rank, 'its')]
                    t = T[rank]
                    b1 = np.zeros(len(QS), bool) if nm_its is None else (b_its >= t['its'] - EPS)
                    names = {ub[(rank, m)][0] for m in MARKERS}
                    if len(names) == 1 and nm_its is not None:
                        lk = np.ones(len(QS), bool)
                        for m in MARKERS:
                            lk &= ub[(rank, m)][1] >= t[m] - EPS
                    else:
                        lk = np.zeros(len(QS), bool)
                    own = (sp if rank == 'species' else ge)[r]
                    for method, called, nm in (('B1', b1, nm_its), ('L', lk, nm_its)):
                        if kind == 'kept':
                            ok = called & (nm == own)
                        else:
                            ok = called                 # a named call in a leave-out fold is a false known
                        res[(rank, method, kind)].append((unit, r, ok))
            # genus leave-out: FK also when the species call names something
        # genus-fold FK at genus rank includes species-level calls (section 7)
        for method in ('B1', 'L'):
            g_rows = res[('genus', method, 'genus')]
            s_rows = {(u, r): ok for u, r, ok in res[('species', method, 'genus')]}
            res[('genus', method, 'genus')] = [(u, r, ok | s_rows[(u, r)]) for u, r, ok in g_rows]
        for (kind, m), n in sorted(ties.items()):
            tie_rows.append((grp, kind, m, n))

        # ---- B2: SINTAX predictions of each fold (run_sintax.py)
        if a.sintax:
            for kind, unit, qs, held in folds:
                if not qs:
                    continue
                fn = f'kept_{a.split}__{unit}.tsv' if kind == 'kept' else f'{kind}__{clean(unit)}.tsv'
                pred = read_sintax(os.path.join(a.sintax, grp, fn))
                for r in qs:
                    p_ = pred.get(r, {})
                    for rank, key, named in (('species', 's', named_sp), ('genus', 'g', named_ge)):
                        lab, boot = p_.get(key, (None, 0.0))
                        nm = lab if lab in named else None
                        called = np.array([nm is not None and boot >= c - EPS for c in CS])
                        own = clean((sp if rank == 'species' else ge)[r] or '')
                        ok = called & (nm == own) if kind == 'kept' else called
                        res[(rank, 'B2', kind)].append((unit, r, ok))
            g_rows = res[('genus', 'B2', 'genus')]
            s_rows = {(u, r): ok for u, r, ok in res[('species', 'B2', 'genus')]}
            res[('genus', 'B2', 'genus')] = [(u, r, ok | s_rows[(u, r)]) for u, r, ok in g_rows]

        # ---- rates and K2 summary
        for rank in ('species', 'genus'):
            for method in methods:
                grid = CS if method == 'B2' else QS
                leave = 'species' if rank == 'species' else 'genus'
                fk_rows, kk_rows = res[(rank, method, leave)], res[(rank, method, 'kept')]
                if not fk_rows or not kk_rows:
                    continue
                fk = np.array([ok for _, _, ok in fk_rows], float)
                kk = np.array([ok for _, _, ok in kk_rows], float)
                unit_fk = collections.defaultdict(list)
                for u, _, ok in fk_rows:
                    unit_fk[u].append(ok)
                unit_kk = collections.defaultdict(list)
                for _, r, ok in kk_rows:
                    unit_kk[sp[r]].append(ok)
                fk_macro = np.mean([np.mean(v, axis=0) for v in unit_fk.values()], axis=0)
                kk_macro = np.mean([np.mean(v, axis=0) for v in unit_kk.values()], axis=0)
                fk_read, kk_read = fk.mean(axis=0), kk.mean(axis=0)
                for i, q in enumerate(grid):
                    sweep_rows.append((grp, rank, method, q, len(fk_rows), len(unit_fk), len(kk_rows),
                                       len(unit_kk), 100 * fk_read[i], 100 * fk_macro[i], 100 * kk_read[i],
                                       100 * kk_macro[i]))
                for avg, f, k in (('read', fk_read, kk_read), ('macro', fk_macro, kk_macro)):
                    for lvl in KK_LEVELS:
                        summ_rows.append((grp, rank, method, avg, lvl, *at_kept(100 * k, 100 * f, lvl)))

    # ---- write
    def write(fn, head, rows):
        with open(os.path.join(a.out, fn), 'w') as fh:
            fh.write('\t'.join(head) + '\n')
            for r in rows:
                fh.write('\t'.join(f'{x:.3f}' if isinstance(x, float) else str(x) for x in r) + '\n')
    write('sweep.tsv', ['group', 'rank', 'method', 'q_or_c', 'fk_reads', 'fk_units', 'kk_reads', 'kk_species',
                        'fk_read_pct', 'fk_macro_pct', 'kk_read_pct', 'kk_macro_pct'], sweep_rows)
    write('k2_summary.tsv', ['group', 'rank', 'method', 'average', 'known_kept_target', 'known_kept_reached',
                             'false_known_pct', 'note'], summ_rows)
    write('depth_check.tsv', ['group', 'fold_queries', 'failing', 'failing_pct', 'ranks_used'], depth_rows)
    write('near_ties.tsv', ['group', 'fold_kind', 'marker', 'fold_calls_with_0.2_set_at_last_rank'], tie_rows)
    dec_rows = []
    for grp in a.groups:
        for avg in ('read', 'macro'):
            curve = {}
            for method in methods:
                rows = [x for x in sweep_rows if x[0] == grp and x[1] == 'species' and x[2] == method]
                if rows:
                    col = (8, 10) if avg == 'read' else (9, 11)
                    curve[method] = (np.array([x[col[1]] for x in rows]), np.array([x[col[0]] for x in rows]))
            if 'L' not in curve or len(curve) < 2:
                continue
            level = min([90.0] + [float(k.max()) for k, _ in curve.values()])
            fk = {m: at_kept(k, f, level)[1] for m, (k, f) in curve.items()}
            base = min((m for m in curve if m != 'L'), key=lambda m: fk[m])
            rr = 1 - fk['L'] / fk[base] if fk[base] > 0 else float('nan')
            dec_rows.append((grp, avg, level, fk.get('B1', ''), fk.get('B2', ''), fk['L'], base, rr))
    write('k2_decision_rows.tsv', ['group', 'average', 'known_kept_level', 'fk_B1_pct', 'fk_B2_pct', 'fk_L_pct',
                                   'best_baseline', 'R'], dec_rows)
    with open(os.path.join(a.out, 'k2_decision.txt'), 'w') as fh:
        rr = {x[0]: x for x in dec_rows if x[1] == 'read'}
        non_amf = [g for g in ('Agaricales', 'Hypocreales', 'Rhabditida') if g in rr]
        n_pass = sum(1 for g in non_amf if rr[g][7] >= 0.25)
        glom = rr.get('Glomeromycota')
        glom_ok = glom is not None and glom[5] <= min(v for v in (glom[3], glom[4]) if v != '') + EPS
        verdict = 'PASS' if n_pass >= 2 and glom_ok else 'FAIL'
        fh.write(f'K2 (section 8, per-read rates, species level): {verdict}\n')
        fh.write(f'  R >= 0.25 in {n_pass} of {len(non_amf)} non-AMF groups; Glomeromycota FK(L) <= best baseline: '
                 f'{glom_ok}\n')
        if not a.sintax:
            fh.write('  PROVISIONAL: B2 (SINTAX) not included; the baseline is B1 only\n')
    with open(os.path.join(a.out, 'manifest.txt'), 'w') as fh:
        fh.write(f'date\t{datetime.datetime.now().isoformat(timespec="seconds")}\n')
        fh.write(f'script\tscore_k2.py\t{sha256(os.path.abspath(__file__))}\n')
        for m in MARKERS:
            p = os.path.join(a.search, f'{m}.hits.tsv.gz')
            fh.write(f'input\t{p}\t{sha256(p)}\n')
        for f in ('units.tsv', 'known_kept.tsv'):
            fh.write(f'input\t{os.path.join(a.folds, f)}\t{sha256(os.path.join(a.folds, f))}\n')
        fh.write(f'split\t{a.split}\ndrop_reads\t{a.drop_reads}\t{len(drop)}\n')
        fh.write(f'seed\t{SEED}\nq\t{QS[0]}..{QS[-1]} step 0.5\n')
    say('done; see', a.out)


def at_kept(kk, fk, lvl):
    """FK at known-kept = lvl, interpolated along the sweep (q order). If lvl is never reached: the highest KK."""
    order = np.argsort(kk, kind='stable')
    kk_s, fk_s = kk[order], fk[order]
    if kk_s[-1] + EPS < lvl:
        i = int(np.argmax(kk))
        return f'{kk[i]:.3f}', float(fk[i]), f'{lvl:g}% not reached'
    if kk_s[0] - EPS > lvl:
        i = int(np.argmin(kk))
        return f'{kk[i]:.3f}', float(fk[i]), f'lowest known-kept {kk[i]:.1f}% is above target'
    return f'{lvl:.3f}', float(np.interp(lvl, kk_s, fk_s)), ''


if __name__ == '__main__':
    main()
