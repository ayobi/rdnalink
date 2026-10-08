# Phase 0 protocol, Amendment 1: how the marker windows are cut

Status: LOCKED 2026-10-08, by the hashes in `LOCK_AMENDMENT1.sha256`.
Author: Aaron O'Brien. Amends section 4 of `PROTOCOL_PHASE0.md` (SHA-256 a81af1fb…, locked 2026-10-07).

No fold, threshold or outcome had been computed when this was written. The only data used are from a pilot of
the window tools, which touches no fold, threshold or outcome. Questions, groups, rules and thresholds are
unchanged.

## 1. Why

Section 4 says ITSxRust finds the ITS boundaries. On the pilot (below), ITSxRust 0.3.0 with its bundled fungal
models found the full ITS in 3,774 of 3,994 reads. Rhabditida got all three windows >= 300 nt in only 92.5% of
reads, short of the 95% per test group a window method must reach (rule 1 in section 4 below).

## 2. Change

Section 4, first bullet, is replaced by:

- **ITSxRust first.** ITSxRust 0.3.0 (`F.hmm`; `--preset hifi --inc-e 1e-3 --min-anchor-score 15
  --max-anchor-evalue 1e-3`) finds the full ITS (ITS1 + 5.8S + ITS2). Where it does, SSU = everything
  upstream of it and LSU = everything downstream (`cut_windows.py itsxfull`). These flanks are identical to
  ITSxRust's own SSU/LSU export on all 2,268 pilot reads where that export exists.
- **barrnap for the rest.** Where ITSxRust finds no full ITS, barrnap 1.10.5 (`--kingdom fun --fast`: Rfam
  SSU_rRNA_eukarya RF01960 and LSU_rRNA_eukarya RF02543, searched by Infernal 1.1.5 in HMM-only mode) finds
  the 18S and 28S, and the ITS is everything between them, with two rules (`cut_windows.py hybrid --join
  --fix-ssu-end`):
  - *Join:* every 18S and 28S hit of >= 50 nt on the strand of the lowest E-value pair is used. The ITS starts
    after the last 18S hit that ends before the anchor 28S, and ends before the first 28S hit after that. A read
    with an 18S hit after its 28S anchor, or a 28S hit before its 18S anchor, holds more than one operon and is
    excluded (`multi_operon`).
  - *18S end:* if the conserved 18S 3' end `GTAGGTGAACCTGC` (<= 2 mismatches) lies inside that ITS, more than
    10 nt after its start, the ITS starts 27 nt after the motif instead.
- Each read's log records which method cut it (`source`: `itsxrust`, `barrnap` or `barrnap+motif`). The
  >= 300 nt rule and the rest of section 4 are unchanged.
- **Original protocol (section 10).** Results under the original section 4 are obtained by removing the
  barrnap-cut reads from queries and reference. The search table is filtered, not rerun, as for folds.

## 3. Pilot

The pilot holds 3,994 reads: 500 drawn at random from each test group plus 2,000 from the whole set (seed
20261007), with 6 duplicates removed. "Kept" means reads with all three windows >= 300 nt.

| Method | Glomeromycota | Agaricales | Hypocreales | Rhabditida | All reads |
|---|---|---|---|---|---|
| ITSxRust alone (section 4 as locked) | 99.8 | 100.0 | 98.2 | 92.5 | 94.0 |
| barrnap `--fast` alone | 99.8 | 100.0 | 98.2 | 98.8 | 99.0 |
| barrnap `--fast --join` alone | 99.8 | 100.0 | 98.2 | 98.6 | 98.6 |
| **Adopted:** ITSxRust, then barrnap `--fast --join --fix-ssu-end` | **99.8** | **100.0** | **98.2** | **99.0** | **98.6** |

Adopted method, all 3,994 reads:
- **Cut:** 3,774 by ITSxRust and 197 by barrnap, of which 0 were moved by the 18S-end rule.
- **Excluded before the length rule (23 reads):** missing 28S 10, missing 18S 9, missing both 2,
  `multi_operon` 1, strand or order conflict 1.
- **Excluded by the length rule:** 33 reads with a window under 300 nt.
- **Kept:** 3,938.

Boundary checks against ITSxRust, on the 3,770 reads both methods cut:

- **barrnap `--fast`:** 248 reads differ by more than 20 nt at either edge. Most were a 28S reported as two hits,
  where the stronger, later hit was taken.
- **With `--join`:** 95 reads differ, 90 of them at the SSU/ITS edge, mostly Cercozoa (68) and Ascomycota outside
  Hypocreales (10). The conserved 18S 3' end sits 27 nt (+/- 10) before ITSxRust's ITS start in 89 of the 95 and
  before barrnap's alone in none; it sits before both in 5, which disagree only at the LSU edge, and before
  neither in 1. So these reads carry an insert near the 3' end of the 18S, most likely an intron, that `--fast`
  leaves unaligned. Control: on the 3,671 reads where the two calls agree within 5 nt, the motif sits before
  the ITS start in 99.0% (every group >= 95.6%), at a median 27 nt.
- **With `--join --fix-ssu-end`:** 90 reads are moved by the motif, and 6 reads differ, all at the ITS/LSU edge
  (5 by 21-27 nt, 1 by 502 nt). The SSU/ITS edge is within 20 nt in 100% of reads (max 11 nt). The ITS/LSU edge
  has a median difference of 4 nt and is within 20 nt in 99.8%.
- **barrnap-cut reads in the adopted method (197):** the 18S end sits at the ITS start in 169, inside the ITS in
  0, and is not found in 28. Of those 28, 24 are Rhabditida, so 24 of that group's 518 reads (4.6%) have SSU/ITS
  boundaries the motif cannot confirm. The report says so.

Rejected:
- **barrnap normal mode** (`cmsearch --rfam -g`): 30 min 53 s for 418 reads, and it found the 28S in only 7 of them.
- **barrnap `--fast` alone:** it misplaces the SSU/ITS boundary where the 18S has an insert, as shown above.

Speed: barrnap `--fast` took 5 min 45 s for the 3,994 pilot reads on 8 threads.

## 4. How the choice was made

The choice rules were written down in sequence, each before the result it applied to:

1. **Before any alternative to ITSxRust was tried:** a method must give all three windows >= 300 nt in >= 95% of each test group.
2. **Before the first barrnap results:** among qualifying methods, the one keeping the most reads wins, and
   `--fast` is used only if normal mode would take more than a day on the full set. Under this rule barrnap
   `--fast` alone (99.0%) would have been chosen.
3. **When the boundary check showed `--fast` misplacing boundaries:** rule 2 was replaced, before any scoring,
   by "ITSxRust where it finds the full ITS (section 4 as locked), barrnap `--fast --join` for the rest". This was
   to be adopted if it cleared rule 1 and the 18S-end check on the 95 disagreeing reads favoured ITSxRust.
   The planned cross-check against normal mode was dropped because normal mode could not find the 28S.
4. **The 18S-end fix** was to be adopted if it brought the reads that differ from ITSxRust by more than 20 nt
   from 95 down to 10 or fewer, with none new, and left every test group >= 95%. It brought them down to 6.

## 5. Files

Pilot commands, logs and reports are in `pilot/` (`report_*.tsv`, `explain_*.tsv`). Scripts (SHA-256):

- `scripts/cut_windows.py` f88f3c5040d2e2d46d29c5b411cfe9d146573c8d45ca5d690a19451e41690065
- `scripts/check_ssu_end.py` d373e4a1c1a25e99ca8d053b814091cc4e45c1a18a5cfa7e848c0a9c8ed1c598
- `scripts/pilot_windows.py` 2380dec4ff563b32e857601bf29c46812b384be5bcd441c1dd0414d3e5976678

Tool versions: ITSxRust 0.3.0, HMMER 3.4, barrnap 1.10.5 (bioconda build 1.10.6), Infernal 1.1.5, seqkit 2.14.0,
Python 3.11.
