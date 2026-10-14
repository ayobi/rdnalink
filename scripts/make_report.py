#!/usr/bin/env python3
"""Phase 0, protocol section 11: the two-page report for Leho, built from the scoring outputs.

Reads (paths relative to --scores; a run that is missing is shown as not run, and its verdict as pending):
  k2_protocol            score_k2.py --split original --drop-reads <barrnap-cut reads>  (original protocol)
  k2_amended             score_k2.py --split original                                   (Amendment 1 windows)
  k2_amended_assembly    score_k2.py --split assembly                                   (and note 2's split)
  k1_protocol, k1_amended                       score_k1.py, without and with the barrnap-cut reads
  k1_protocol_trees, k1_amended_trees           check_flag_trees.py on each precision sample
plus windows/win.log.tsv.gz (how each read was cut), search/queries.tsv and the protocol's lock files.

Results under the original protocol come first (section 10). Each verdict is computed exactly as the protocol
states it: K2 (section 8) from k2_decision.txt; K1 (section 9) passes if pooled recall >= 80% for chimeras (F1)
and for label swaps (F2) and pooled precision >= 50% (note 3). The report does not choose between the original
and the amended result; it gives both.

Writes report/phase0_report.pdf (two pages), report/phase0_report.md (the same, as text) and
report/numbers.tsv (every number in the report, with its source file).

  python3 scripts/make_report.py --scores scores --out report [--banner "TEST: scrambled labels"]
"""
import argparse
import csv
import datetime
import glob
import gzip
import math
import os
import re
import textwrap
from collections import Counter, OrderedDict

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

GROUPS = ['Glomeromycota', 'Agaricales', 'Hypocreales', 'Rhabditida']
NON_AMF = GROUPS[1:]
METHODS = [('B1', 'ITS best hit', '#2a78d6', '-'), ('B2', 'SINTAX', '#eb6834', '--'), ('L', 'linked', '#1baf7a', '-')]
INK, INK2, GRID, SURFACE = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
K2_RUNS = [('k2_protocol', 'Original protocol (barrnap-cut reads removed, section 7 split)'),
           ('k2_amended', 'Amendment 1 windows, section 7 split'),
           ('k2_amended_assembly', 'Amendment 1 windows, note 2 assembly split')]
K1_RUNS = [('k1_protocol', 'Original protocol'), ('k1_amended', 'Amendment 1 windows')]


# ------------------------------------------------------------------ reading
def tsv(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def f1(x, nd=1):
    v = num(x)
    return '-' if v is None else f'{v:.{nd}f}'


class Numbers:
    """Every number shown, with where it came from."""
    def __init__(self):
        self.rows = []

    def add(self, key, value, source):
        self.rows.append((key, value, source))
        return value

    def write(self, path):
        with open(path, 'w') as fh:
            fh.write('key\tvalue\tsource\n')
            for k, v, s in self.rows:
                fh.write(f'{k}\t{v}\t{s}\n')


def k2_verdict(d):
    path = os.path.join(d, 'k2_decision.txt')
    if not os.path.exists(path):
        return 'not run', ''
    text = open(path).read()
    m = re.search(r':\s*(PASS|FAIL)', text)
    note = 'provisional: B2 (SINTAX) missing' if 'PROVISIONAL' in text else ''
    return (m.group(1) if m else '?'), note


def k1_verdict(d, trees):
    rec = tsv(os.path.join(d, 'recall.tsv'))
    if rec is None:
        return 'not run', {}
    r = {x['spike_type']: num(x['recall_pct']) for x in rec if x['group'] == 'all'}
    prec = tsv(os.path.join(trees, 'precision.tsv'))
    p = {x['flag']: (num(x['precision_pct']), x['real'], x['checked']) for x in prec} if prec else {}
    out = {'recall_chimera': r.get('chimera'), 'recall_swap': r.get('swap'),
           'precision': p.get('all', (None,))[0], 'precision_by_flag': p}
    rec_ok = all(v is not None and v >= 80.0 for v in (r.get('chimera'), r.get('swap')))
    if not prec:
        return ('pending (trees not run)' if rec_ok else 'FAIL'), out
    pr = out['precision']
    return ('PASS' if rec_ok and pr is not None and pr >= 50.0 else 'FAIL'), out


def windows_counts(path):
    c = Counter()
    if not os.path.exists(path):
        return c
    with gzip.open(path, 'rt') as fh:
        next(fh)
        for line in fh:
            v = line.rstrip('\n').split('\t')
            c[(v[2], v[3])] += 1
    return c


def query_counts(path):
    c = Counter()
    if os.path.exists(path):
        with open(path) as fh:
            next(fh)
            for line in fh:
                c[line.split('\t')[1]] += 1
    return c


def locks(protocol):
    out = []
    for f in sorted(glob.glob(os.path.join(protocol, '*.md'))):
        head = open(f).read(600)
        m = re.search(r'Status:\s*LOCKED\s*([0-9-]+)', head)
        title = head.splitlines()[0].lstrip('# ').strip()
        out.append((os.path.basename(f), m.group(1) if m else 'not locked', title))
    return out


# ------------------------------------------------------------------ the figure
def curves(ax_list, sweep, level):
    for ax, grp in zip(ax_list, GROUPS):
        ax.set_facecolor(SURFACE)
        rows = [x for x in sweep if x['group'] == grp and x['rank'] == 'species'] if sweep else []
        ends = []
        for code, label, color, ls in METHODS:
            pts = sorted((num(x['kk_read_pct']), num(x['fk_read_pct'])) for x in rows if x['method'] == code)
            pts = [(a, b) for a, b in pts if a is not None and b is not None]
            if not pts:
                continue
            xs, ys = zip(*pts)
            ax.plot(xs, ys, color=color, lw=1.6, ls=ls, solid_capstyle='round', label=f'{code} {label}')
            ends.append((xs[-1], ys[-1], code, color))
        lv = level.get(grp, 90.0)
        ax.axvline(lv, color=INK2, lw=0.8, ls=':')
        ax.text(lv, 1.0, f' {lv:.0f}%', transform=ax.get_xaxis_transform(), color=INK2, fontsize=6.5, va='top')
        for x, y, code, color in ends:            # direct labels at the high-known-kept end of each line
            ax.annotate(code, (x, y), xytext=(3, 0), textcoords='offset points', fontsize=7, color=INK,
                        va='center')
        ax.set_title(grp, fontsize=8.5, color=INK, loc='left', pad=3)
        xs_all = [num(x['kk_read_pct']) for x in rows if num(x['kk_read_pct']) is not None]
        ax.set_xlim(min([50.0] + [v - 2 for v in xs_all]), 102)
        ax.set_ylim(bottom=0)
        ax.grid(True, color=GRID, lw=0.5)
        ax.tick_params(labelsize=6.5, colors=INK2, length=2)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        for s in ('left', 'bottom'):
            ax.spines[s].set_color(GRID)
        if not rows:
            ax.text(0.5, 0.5, 'not run', transform=ax.transAxes, ha='center', color=INK2, fontsize=8)
    ax_list[0].set_ylabel('false known (% of reads)', fontsize=7, color=INK2)
    for ax in ax_list:
        ax.set_xlabel('known kept (%)', fontsize=7, color=INK2)


# ------------------------------------------------------------------ text helpers
def table(rows, widths):
    out = []
    for r in rows:
        out.append('  '.join(str(c).ljust(w) if i == 0 else str(c).rjust(w) for i, (c, w) in enumerate(zip(r, widths))))
    return '\n'.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scores', default='scores')
    ap.add_argument('--windows-log', default='windows/win.log.tsv.gz')
    ap.add_argument('--queries', default='search/queries.tsv')
    ap.add_argument('--protocol', default='protocol')
    ap.add_argument('--out', default='report')
    ap.add_argument('--banner', default='', help='a line printed in red across both pages (e.g. for a test run)')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    N = Numbers()
    S = lambda *p: os.path.join(a.scores, *p)  # noqa: E731

    # ---- K2
    k2 = OrderedDict()
    for run, label in K2_RUNS:
        verdict, note = k2_verdict(S(run))
        rows = tsv(S(run, 'k2_decision_rows.tsv')) or []
        k2[run] = {'label': label, 'verdict': verdict, 'note': note,
                   'rows': {(x['group'], x['average']): x for x in rows}}
        N.add(f'{run}.verdict', verdict, f'{run}/k2_decision.txt')
        for (g, avg), x in k2[run]['rows'].items():
            for col in ('known_kept_level', 'fk_B1_pct', 'fk_B2_pct', 'fk_L_pct', 'best_baseline', 'R'):
                N.add(f'{run}.{g}.{avg}.{col}', x[col], f'{run}/k2_decision_rows.tsv')
    main_k2 = k2['k2_protocol']

    def k2_table(run, avg='read'):
        rows = [('group', 'at KK %', 'FK B1', 'FK B2', 'FK L', 'R', 'meets rule')]
        for g in GROUPS:
            x = k2[run]['rows'].get((g, avg))
            if not x:
                rows.append((g, '-', '-', '-', '-', '-', 'not run'))
                continue
            base = min(v for v in (num(x['fk_B1_pct']), num(x['fk_B2_pct'])) if v is not None) \
                if any(num(x[c]) is not None for c in ('fk_B1_pct', 'fk_B2_pct')) else None
            fl = num(x['fk_L_pct'])
            if g == 'Glomeromycota':
                ok = 'yes' if base is not None and fl is not None and fl <= base + 1e-9 else 'no'
                rule = f'FK(L) <= best: {ok}'
            else:
                r = num(x['R'])
                rule = 'R undefined' if r is None else ('R >= 0.25: yes' if r >= 0.25 else 'R >= 0.25: no')
            rows.append((g, f1(x['known_kept_level'], 0), f1(x['fk_B1_pct']), f1(x['fk_B2_pct']), f1(x['fk_L_pct']),
                         f1(x['R'], 2), rule))
        return table(rows, (14, 7, 6, 6, 6, 6, 18))

    level = {g: num(main_k2['rows'].get((g, 'read'), {}).get('known_kept_level')) or 90.0 for g in GROUPS}
    sweep = tsv(S('k2_protocol', 'sweep.tsv'))

    # ---- K1
    k1 = OrderedDict()
    for run, label in K1_RUNS:
        verdict, d = k1_verdict(S(run), S(run + '_trees'))
        k1[run] = {'label': label, 'verdict': verdict, **d,
                   'recall': tsv(S(run, 'recall.tsv')) or [], 'counts': tsv(S(run, 'flag_counts.tsv')) or []}
        N.add(f'{run}.verdict', verdict, f'{run}/recall.tsv, {run}_trees/precision.tsv')
        for x in k1[run]['recall']:
            N.add(f'{run}.recall.{x["group"]}.{x["spike_type"]}', f'{x["caught"]}/{x["spikeins"]}', f'{run}/recall.tsv')
        for flag, (p, real, checked) in d.get('precision_by_flag', {}).items():
            N.add(f'{run}.precision.{flag}', f'{real}/{checked}', f'{run}_trees/precision.tsv')
        for x in k1[run]['counts']:
            N.add(f'{run}.flags.{x["group"]}', f'F1 {x["F1"]}, F2 {x["F2"]} of {x["reads"]}', f'{run}/flag_counts.tsv')

    def k1_table(run):
        x = k1[run]
        rec = {(r['group'], r['spike_type']): r for r in x['recall']}
        cnt = {r['group']: r for r in x['counts']}
        rows = [('group', 'chimeras caught', 'swaps caught', 'real F1', 'real F2')]
        for g in GROUPS + ['all']:
            c, s = rec.get((g, 'chimera')), rec.get((g, 'swap'))
            n = cnt.get(g)
            rows.append((g, f'{c["caught"]}/{c["spikeins"]} ({f1(c["recall_pct"])}%)' if c else '-',
                         f'{s["caught"]}/{s["spikeins"]} ({f1(s["recall_pct"])}%)' if s else '-',
                         n['F1'] if n else (sum(int(v['F1']) for v in x['counts']) if g == 'all' and x['counts'] else '-'),
                         n['F2'] if n else (sum(int(v['F2']) for v in x['counts']) if g == 'all' and x['counts'] else '-')))
        out = table(rows, (14, 20, 20, 8, 8))
        p = x.get('precision_by_flag') or {}
        if p:
            parts = [f'{fl} {p[fl][1]}/{p[fl][2]} ({f1(p[fl][0])}%)' for fl in ('F1', 'F2', 'all') if fl in p]
            out += '\nPrecision on trees (real/checked): ' + '; '.join(parts)
        else:
            out += '\nPrecision on trees: not run'
        return out

    # ---- depth check, near ties, excluded reads
    depth = tsv(S('k2_protocol', 'depth_check.tsv')) or tsv(S('k2_amended', 'depth_check.tsv')) or []
    ties = tsv(S('k2_protocol', 'near_ties.tsv')) or []
    wc = windows_counts(a.windows_log)
    qc = query_counts(a.queries)
    cut_ok = sum(n for (st, src), n in wc.items() if st == 'ok')
    barr = sum(n for (st, src), n in wc.items() if st == 'ok' and src.startswith('barrnap'))
    not_cut = sorted(((st, n) for (st, src), n in wc.items() if st != 'ok'), key=lambda x: -x[1])
    N.add('windows.reads', sum(wc.values()), a.windows_log)
    N.add('windows.cut', cut_ok, a.windows_log)
    N.add('windows.barrnap_cut', barr, a.windows_log)
    for st, n in not_cut:
        N.add(f'windows.not_cut.{st}', n, a.windows_log)
    for g in GROUPS:
        N.add(f'queries.{g}', qc.get(g, 0), a.queries)
    for x in depth:
        N.add(f'depth.{x["group"]}', f'{x["failing"]}/{x["fold_queries"]} ({x["failing_pct"]}%), ranks <= {x["ranks_used"]}',
              'k2_protocol/depth_check.tsv')
    tie_total = sum(int(x['fold_calls_with_0.2_set_at_last_rank']) for x in ties)
    N.add('near_ties.fold_calls', tie_total, 'k2_protocol/near_ties.tsv')
    docs = locks(a.protocol)

    # ---- text blocks
    date = datetime.date.today().isoformat()
    head = (f'rdnalink Phase 0: results\nEUKARYOME v2.1 long reads; protocol locked 2026-10-07. '
            f'Report built {date} by scripts/make_report.py.')
    verdicts = ['Verdicts (results under the original protocol first, section 10)']
    for run, x in k2.items():
        verdicts.append(f'  K2  {x["label"]}: {x["verdict"]}' + (f' ({x["note"]})' if x['note'] else ''))
    for run, x in k1.items():
        extra = ''
        if x.get('recall_chimera') is not None:
            pr = x.get('precision')
            extra = (f' (recall: chimeras {f1(x["recall_chimera"])}%, swaps {f1(x["recall_swap"])}%; '
                     + (f'precision {f1(pr)}%)' if pr is not None else 'precision: trees not run)'))
        verdicts.append(f'  K1  {x["label"]}: {x["verdict"]}{extra}')
    rules = ('K2 passes if, at 90% known-kept (or the highest level all methods reach), R = 1 - FK(L)/FK(best of B1, '
             'B2) >= 0.25 in at least two of Agaricales, Hypocreales and Rhabditida, and FK(L) <= FK(best) in '
             'Glomeromycota (species level, per-read rates). K1 passes if spiked chimeras and label swaps are each '
             'caught >= 80% of the time (pooled) and >= 50% of a random sample of real flags hold up on trees.')
    k2_main = ('K2, original protocol: false known (FK, % of reads) at the known-kept level, per read, species '
               'level\n' + k2_table('k2_protocol'))
    k2_other = []
    for run in ('k2_amended', 'k2_amended_assembly'):
        k2_other.append(f'{k2[run]["label"]}: {k2[run]["verdict"]}\n' + k2_table(run))
    k2_macro = 'Original protocol, rates averaged over species (macro; reported, not deciding)\n' + \
        k2_table('k2_protocol', 'macro')
    k1_text = [f'K1, {x["label"]}: {x["verdict"]}\n' + k1_table(run) for run, x in k1.items()]
    depth_text = ['Depth check (note 2: a group with more than 1% of fold queries failing reads ranks <= 1,000)']
    depth_text += [f'  {x["group"]}: {x["failing"]} of {x["fold_queries"]} fold queries fail ({x["failing_pct"]}%); '
                   f'ranks <= {x["ranks_used"]}' for x in depth] or ['  not run']
    depth_text.append(f'Fold calls whose 0.2-point set reaches the last rank read (note 2): {tie_total}')
    excl = [f'Reads and windows: {sum(wc.values()):,} reads; {cut_ok:,} cut into three windows, {barr:,} of them '
            f'by barrnap (removed under the original protocol).']
    if not_cut:
        excl.append('Not cut: ' + '; '.join(f'{st} {n:,}' for st, n in not_cut) + '.')
    excl.append('Query reads per group: ' + '; '.join(f'{g} {qc.get(g, 0):,}' for g in GROUPS) + '.')
    doc_text = ['Protocol, amendments and notes, with the dates they were locked']
    for f, d, t in docs:
        doc_text.append(textwrap.fill(f'{f} (locked {d}): {t}', 112, initial_indent='  ', subsequent_indent='      '))
    note4 = ('Real flags are reads whose markers conflict on per-marker trees. EUKARYOME documents a per-marker '
             'BLASTn check for its (meta)genome and (meta)transcriptome data only, so they are not described as '
             'chimeras its checks missed (note 4). Their accessions are in scores/k1_*/flags.tsv.')

    # ---- markdown
    md = [f'# {head.splitlines()[0]}', '', head.splitlines()[1], '']
    if a.banner:
        md += [f'**{a.banner}**', '']
    md += ['## Verdicts', '', '```', '\n'.join(verdicts), '```', '', rules, '', '## K2', '', '```', k2_main, '```',
           '', '```', '\n\n'.join(k2_other), '```', '', '```', k2_macro, '```', '', '## K1', '', '```',
           '\n\n'.join(k1_text), '```', '', note4, '', '## Checks and exclusions', '', '```', '\n'.join(depth_text),
           '```', '', ' '.join(excl), '', '```', '\n'.join(doc_text), '```', '',
           'Every number above is in numbers.tsv with the file it came from.']
    with open(os.path.join(a.out, 'phase0_report.md'), 'w') as fh:
        fh.write('\n'.join(md) + '\n')
    N.write(os.path.join(a.out, 'numbers.tsv'))

    # ---- PDF: two A4 pages
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})
    mono = {'family': 'DejaVu Sans Mono', 'fontsize': 7.2, 'color': INK, 'va': 'top'}
    with PdfPages(os.path.join(a.out, 'phase0_report.pdf')) as pdf:
        # page 1
        fig = plt.figure(figsize=(8.27, 11.69), facecolor='white')
        fig.text(0.07, 0.965, head.splitlines()[0], fontsize=15, color=INK, weight='bold', va='top')
        fig.text(0.07, 0.94, head.splitlines()[1], fontsize=8, color=INK2, va='top')
        if a.banner:
            fig.text(0.93, 0.965, a.banner, fontsize=9, color='#d03b3b', ha='right', va='top', weight='bold')
        fig.text(0.07, 0.915, '\n'.join(verdicts), **{**mono, 'fontsize': 8})
        fig.text(0.07, 0.835, textwrap.fill(rules, 128), fontsize=7.2, color=INK2, va='top')
        fig.text(0.07, 0.785, k2_main, **mono)
        axs = [fig.add_axes([0.07 + i * 0.225, 0.47, 0.19, 0.17]) for i in range(4)]
        curves(axs, sweep, level)
        handles, labels = axs[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.07, 0.68), ncol=3, frameon=False,
                       fontsize=7.5, labelcolor=INK)
        fig.text(0.07, 0.43, textwrap.fill('Figure. Original protocol, species level, per read: false known against '
                 'known kept as the threshold (q for B1 and L, c for B2) is swept; dotted line, the level K2 is read '
                 'at.', 130), fontsize=7, color=INK2, va='top')
        fig.text(0.07, 0.40, k2_macro, **mono)
        pdf.savefig(fig)
        plt.close(fig)
        # page 2
        fig = plt.figure(figsize=(8.27, 11.69), facecolor='white')
        if a.banner:
            fig.text(0.93, 0.965, a.banner, fontsize=9, color='#d03b3b', ha='right', va='top', weight='bold')
        fig.text(0.07, 0.965, 'K2 with the amendments', fontsize=10, color=INK, weight='bold', va='top')
        fig.text(0.07, 0.945, '\n\n'.join(k2_other), **mono)
        fig.text(0.07, 0.75, 'K1', fontsize=10, color=INK, weight='bold', va='top')
        fig.text(0.07, 0.73, '\n\n'.join(k1_text), **mono)
        fig.text(0.07, 0.445, textwrap.fill(note4, 128), fontsize=7.2, color=INK2, va='top')
        fig.text(0.07, 0.405, 'Checks and exclusions', fontsize=10, color=INK, weight='bold', va='top')
        fig.text(0.07, 0.385, '\n'.join(depth_text), **mono)
        fig.text(0.07, 0.30, '\n'.join(textwrap.fill(x, 128) for x in excl), fontsize=7.2, color=INK, va='top')
        fig.text(0.07, 0.235, '\n'.join(doc_text), **mono)
        fig.text(0.07, 0.04, 'Every number in this report is in numbers.tsv, with the file it came from.',
                 fontsize=7, color=INK2)
        pdf.savefig(fig)
        plt.close(fig)
    print(f'wrote {a.out}/phase0_report.pdf, phase0_report.md and numbers.tsv')
    for line in verdicts:
        print(line)


if __name__ == '__main__':
    main()
