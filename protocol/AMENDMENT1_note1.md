# Phase 0 protocol, Amendment 1, note 1: ITSxRust and barrnap run in chunks

Status: LOCKED 2026-10-10, by its hash added to `LOCK_AMENDMENT1.sha256`.
Author: Aaron O'Brien. Adds to `AMENDMENT1_windows.md` (locked 2026-10-08). No question, group, rule or
threshold changes. It is recorded as a dated amendment (protocol section 10) because it changes results. No fold,
threshold or outcome had been computed when it was written.

## 1. Why

HMMER scales E-values with the number of sequences searched, and Infernal (inside barrnap) with the number of
residues. The E-value cut-offs of Amendment 1 (ITSxRust `--inc-e 1e-3 --max-anchor-evalue 1e-3`; barrnap's
`-E 0.001`) were validated on the 3,994-read pilot. One ITSxRust run over all 400,217 reads, 100 times the pilot,
gave:

- the full ITS in 290,268 reads (72.5%), against 94.5% in the pilot and about 91% among its 2,000 random reads;
- the full ITS for only 3,264 of the 3,774 reads that had one in the pilot, although the reads were identical;
- 79 anchor hits per read, against 142 in the pilot.

The short 5.8S anchors fall under the cut-off first. That run was discarded.

## 2. Change

The long-read set is split in file order into 100 chunks of 4,002-4,003 reads (`seqkit split2 -p 100`). Each chunk
is run through ITSxRust and barrnap `--kingdom fun --fast` with the parameters of Amendment 1, as the pilot was,
and the outputs are concatenated before `cut_windows.py`. Chunk order does not matter: E-values depend on the
size of the search, not on what is in it.

## 3. Result (all 400,217 reads)

- All 100 chunks finished cleanly. Wall time was 12 h 0 min on a MacBook Pro (Intel, 8 threads), about 7 min per
  chunk.
- ITSxRust found the full ITS in 365,878 reads (91.4%). Of the 3,774 pilot reads with a full ITS, 3,773 keep it.
- **Cut:** 365,878 by ITSxRust, 31,164 by barrnap and 399 by barrnap with the 18S-end rule (290 of them
  Cercozoa).
- **Excluded before the length rule (2,776 reads):** missing 28S 1,296, missing 18S 684, strand or order
  conflict 439, missing both 344, `multi_operon` 13.
- **Excluded by the length rule:** 5,821 reads with a window under 300 nt.

| Group | Reads | Cut | All three windows >= 300 nt |
|---|---|---|---|
| Glomeromycota | 5,557 | 5,557 | 5,546 (99.8%) |
| Agaricales | 9,454 | 9,441 | 9,436 (99.8%) |
| Hypocreales | 3,388 | 3,323 | 3,318 (97.9%) |
| Rhabditida | 2,919 | 2,906 | 2,901 (99.4%) |
| All reads | 400,217 | 397,441 | 391,620 (97.9%) |

Every test group clears the 95% rule.

Cross-check: barrnap `--fast --join --fix-ssu-end` was also run on every read and compared with ITSxRust on the
365,549 reads both cut.
- Both put the read on the same strand in 365,545 of them.
- SSU/ITS edge: median difference 1 nt, 99.9% within 20 nt (max 2,164 nt).
- ITS/LSU edge: median difference 4 nt, 99.8% within 20 nt (max 2,354 nt).

On the 31,563 barrnap-cut reads of the adopted method, the 18S 3' end sits at the ITS start in 30,642, inside the
ITS in 0, and is not found in 921.

## 4. Files

`windows/manifest.txt` holds the tool versions and the hashes of the reference, the scripts, the concatenated
ITSxRust and barrnap outputs, and the windows. `windows/windows_report.tsv` and `windows/check_fallback.txt` hold
the per-group tables, and `windows/win.log.tsv.gz` holds the per-read cut log.
