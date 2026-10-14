"""Cut each full-operon read into SSU, ITS and LSU windows.

The method is the one validated for the rdnalink test (Amendment 1 and its note 1):

- Reads are processed in chunks of about 4,000 in file order, because HMMER and Infernal E-values scale with the
  size of the search, and the cut-offs were validated on chunks of that size.
- ITSxRust 0.3.0 (F.hmm; --preset hifi --inc-e 1e-3 --min-anchor-score 15 --max-anchor-evalue 1e-3) finds the full
  ITS. Where it does, SSU = everything upstream and LSU = everything downstream.
- Where it does not, barrnap (--kingdom fun --fast) finds the 18S and 28S and the ITS is what lies between them,
  with two rules. Join: every 18S and 28S hit of >= 50 nt on the anchor strand is used, and a read with an 18S hit
  after its 28S anchor, or a 28S hit before its 18S anchor, holds more than one operon (multi_operon). 18S end: if
  the conserved 18S 3' end (GTAGGTGAACCTGC, <= 2 mismatches) lies inside that ITS, the ITS starts 27 nt after it.
- Windows come out in rRNA orientation. Reads with a window shorter than --min-length are logged, not written.

The boundary functions take the tools' outputs, so they can be tested and rerun without the tools.
"""
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

from .seqio import rc, read_fasta, read_id, write_fasta

ITSXRUST_ARGS = ['--region', 'all', '--preset', 'hifi', '--inc-e', '1e-3', '--min-anchor-score', '15',
                 '--max-anchor-evalue', '1e-3', '--plain-ids']
BARRNAP_ARGS = ['--kingdom', 'fun', '--fast']
SSU_NAMES, LSU_NAMES = ('18S_rRNA',), ('28S_rRNA',)
MIN_JOIN = 50
SSU_END_MOTIF = 'GTAGGTGAACCTGC'
MOTIF_TO_ITS, MOTIF_TOL, MOTIF_MAXMM = 27, 10, 2
MARKERS = ('ssu', 'its', 'lsu')


# ------------------------------------------------------------------ boundaries from the tools' outputs
def bounds_itsxrust(reads, full):
    """reads: {id: seq}; full: {id: full-ITS seq from ITSxRust}. -> {id: (status, source, strand, start, end)}"""
    out = {}
    for k, seq in reads.items():
        f = full.get(k)
        if not f:
            out[k] = ('no_its', 'itsxrust', '', None, None)
            continue
        s = seq.upper()
        i = s.find(f)
        if i >= 0:
            out[k] = ('ok', 'itsxrust', '+', i, i + len(f))
            continue
        i = rc(s).find(f)
        out[k] = ('ok', 'itsxrust', '-', i, i + len(f)) if i >= 0 else ('its_not_in_read', 'itsxrust', '', None, None)
    return out


def _evalue(x):
    try:
        return float(x)
    except ValueError:
        return 1.0


def load_gff(lines):
    feats = defaultdict(list)                       # id -> [(name, start, end, strand, evalue)], forward strand
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9:
            continue
        m = re.search(r'Name=([^;]+)', c[8])
        if m:
            feats[read_id(c[0])].append((m.group(1), int(c[3]), int(c[4]), c[6], _evalue(c[5])))
    return feats


def _motif_hits(seq, lo, hi):
    lo, hi = max(0, lo), min(len(seq) - len(SSU_END_MOTIF), hi)
    out = []
    for p in range(lo, hi + 1):
        mm = 0
        for x, y in zip(seq[p:p + len(SSU_END_MOTIF)], SSU_END_MOTIF):
            if x != y:
                mm += 1
                if mm > MOTIF_MAXMM:
                    break
        if mm <= MOTIF_MAXMM:
            out.append((p, mm))
    return out


def fix_ssu_end(oriented, i, j):
    """New ITS start when the 18S 3' end motif lies inside barrnap's ITS, else None."""
    h = _motif_hits(oriented.upper(), i - MOTIF_TO_ITS - MOTIF_TOL, j - len(SSU_END_MOTIF))
    if any(abs(p + MOTIF_TO_ITS - i) <= MOTIF_TOL for p, _ in h):
        return None
    later = [(mm, p) for p, mm in h if i + MOTIF_TOL < p + MOTIF_TO_ITS < j]
    return min(later)[1] + MOTIF_TO_ITS if later else None


def _oriented_hits(feat_list, n):
    hits = {'+': {'SSU': [], 'LSU': []}, '-': {'SSU': [], 'LSU': []}}
    for name, s, e, strand, ev in feat_list:
        marker = 'SSU' if name in SSU_NAMES else 'LSU' if name in LSU_NAMES else None
        if marker is None or strand not in '+-':
            continue
        if strand == '-':
            s, e = n - e + 1, n - s + 1
        hits[strand][marker].append((s, e, ev))
    return hits


def bounds_barrnap(reads, feats, join=True, fix=True):
    out = {}
    for k, seq in reads.items():
        hits = _oriented_hits(feats.get(k, []), len(seq))
        cands = []
        for strand in '+-':
            ssu, lsu = hits[strand]['SSU'], hits[strand]['LSU']
            if not ssu or not lsu:
                continue
            a = min(ssu, key=lambda h: h[2])
            b = min(lsu, key=lambda h: h[2])
            if not a[1] < b[0]:
                continue
            i, j, st = a[1], b[0] - 1, 'ok'
            if join:
                long_ssu = [h for h in ssu if h[1] - h[0] + 1 >= MIN_JOIN]
                long_lsu = [h for h in lsu if h[1] - h[0] + 1 >= MIN_JOIN]
                if any(h[0] > b[1] for h in long_ssu) or any(h[1] < a[0] for h in long_lsu):
                    st = 'multi_operon'
                i = max([a[1]] + [h[1] for h in long_ssu if h[1] < b[0]])
                j = min([b[0]] + [h[0] for h in long_lsu if h[0] > i]) - 1
            worst = max(a[2], b[2], 1e-300)
            cands.append((-math.log10(worst), strand, i, j, st))
        if cands:
            _, strand, i, j, st = max(cands)
            src = 'barrnap'
            if st == 'ok' and fix:
                n = fix_ssu_end(seq if strand == '+' else rc(seq), i, j)
                if n is not None:
                    i, src = n, 'barrnap+motif'
            out[k] = (st, src, strand, i, j) if st == 'ok' else (st, src, '', None, None)
            continue
        have = {m for s in '+-' for m in ('SSU', 'LSU') if hits[s][m]}
        st = 'strand_or_order' if have >= {'SSU', 'LSU'} else \
            'missing_' + '_'.join(m for m in ('SSU', 'LSU') if m not in have)
        out[k] = (st, 'barrnap', '', None, None)
    return out


def bounds_hybrid(reads, full, feats):
    """ITSxRust where it found the full ITS, barrnap (join, 18S-end fix) for the rest."""
    b = bounds_itsxrust(reads, full)
    rest = {k: s for k, s in reads.items() if b[k][0] != 'ok'}
    b.update(bounds_barrnap(rest, feats, join=True, fix=True))
    return b


def split_read(seq, bound):
    """(ssu, its, lsu) in rRNA orientation, or None when the read was not cut."""
    st, _, strand, i, j = bound
    if st != 'ok':
        return None
    s = seq if strand == '+' else rc(seq)
    return s[:i], s[i:j], s[j:]


# ------------------------------------------------------------------ running the tools
def _run(cmd, log, stdout=None):
    log.write('$ ' + ' '.join(cmd) + '\n')
    log.flush()
    p = subprocess.run(cmd, stdout=stdout or log, stderr=log)
    if p.returncode:
        sys.exit(f'failed ({p.returncode}): {" ".join(cmd)}; see {log.name}')


def tool_versions(itsxrust='itsxrust', barrnap='barrnap'):
    out = {}
    for name, exe in (('itsxrust', itsxrust), ('barrnap', barrnap)):
        try:
            p = subprocess.run([exe, '--version'], capture_output=True, text=True)
            out[name] = (p.stdout or p.stderr).strip().splitlines()[0]
        except (OSError, IndexError):
            out[name] = 'not found'
    return out


def cut(reads_path, out_prefix, hmm, threads=4, chunk=4000, min_length=300, itsxrust='itsxrust',
        barrnap='barrnap', keep_work=False):
    """Run ITSxRust and barrnap chunk by chunk and write <out>.ssu.fasta, <out>.its.fasta, <out>.lsu.fasta and
    <out>.log.tsv. A finished chunk is kept, so a stopped run resumes. Returns status counts."""
    for exe in (itsxrust, barrnap):
        if shutil.which(exe) is None and not os.path.exists(exe):
            sys.exit(f'{exe} not found on PATH')
    work = out_prefix + '.work'
    os.makedirs(work, exist_ok=True)
    log = open(os.path.join(work, 'tools.log'), 'a')
    records = list(read_fasta(reads_path))
    ids = [read_id(h) for h, _ in records]
    dup = [k for k, n in Counter(ids).items() if n > 1]
    if dup:
        sys.exit(f'duplicate read IDs, e.g. {dup[:3]}; IDs must be unique up to the first space, ; or |')
    status = Counter()
    outs = {m: open(f'{out_prefix}.{m}.fasta', 'w') for m in MARKERS}
    with open(f'{out_prefix}.log.tsv', 'w') as lg:
        lg.write('id\tlength\tstatus\tsource\tstrand\tits_start\tits_end\tssu_len\tits_len\tlsu_len\n')
        for ci in range(0, len(records), chunk):
            part = records[ci:ci + chunk]
            cdir = os.path.join(work, f'chunk_{ci // chunk:05d}')
            done = os.path.join(cdir, 'done')
            os.makedirs(cdir, exist_ok=True)
            fa = os.path.join(cdir, 'reads.fasta')
            if not os.path.exists(done):
                with open(fa, 'w') as fh:
                    for h, s in part:
                        write_fasta(fh, read_id(h), s)
                _run([itsxrust, 'extract', '--input', fa, '--hmm', hmm, '--output', os.path.join(cdir, 'itsx'),
                      '--hmmer-cpu', str(threads)] + ITSXRUST_ARGS, log)
                with open(os.path.join(cdir, 'barrnap.gff'), 'w') as gff:
                    _run([barrnap, *BARRNAP_ARGS, '--threads', str(threads), fa], log, stdout=gff)
                open(done, 'w').close()
            reads = {read_id(h): s for h, s in part}
            heads = {read_id(h): h for h, s in part}
            fpath = os.path.join(cdir, 'itsx.full.fasta')
            full = {read_id(h): s.upper() for h, s in read_fasta(fpath)} if os.path.exists(fpath) else {}
            with open(os.path.join(cdir, 'barrnap.gff')) as fh:
                feats = load_gff(fh)
            b = bounds_hybrid(reads, full, feats)
            for k in reads:
                st, src, strand, i, j = b[k]
                w = split_read(reads[k], b[k])
                lens = [len(x) for x in w] if w else ['', '', '']
                if w and min(lens) < min_length:
                    st = 'short_window'
                status[(st, src)] += 1
                lg.write(f'{k}\t{len(reads[k])}\t{st}\t{src}\t{strand}\t{"" if i is None else i}\t'
                         f'{"" if j is None else j}\t' + '\t'.join(map(str, lens)) + '\n')
                if st == 'ok':
                    for m, x in zip(MARKERS, w):
                        write_fasta(outs[m], heads[k], x)
            if not keep_work:
                for f in ('reads.fasta',):
                    p = os.path.join(cdir, f)
                    if os.path.exists(p):
                        os.remove(p)
    for fh in outs.values():
        fh.close()
    log.close()
    return status
