# Phase 0 protocol: do linked SSU-ITS-LSU calls beat ITS alone outside AMF?

Status: LOCKED 2026-10-07.
Author: Aaron O'Brien. Drafted 2026-10-07.

## 1. Questions

- **K2 (main).** When a read comes from a taxon the reference does not hold, does a call that uses SSU,
  ITS and LSU from the same read call it a known taxon less often than a call that uses ITS alone, at the
  same rate of keeping known taxa?
- **K1.** When the three markers of one read point to different taxa, or agree on a taxon other than the
  read's label, is that usually a real problem (a chimera or a mislabel)?

If K2 fails, the tool narrows to curation (K1) or stops. If both fail, it stops.

## 2. Data and truth

- Reference and queries: EUKARYOME v2.1 long-read set, `General_EUK_longread_v2.1.fasta` (400,217 reads),
  CC BY 4.0. Its SHA-256 is recorded in `LOCK.sha256` with this file.
- Truth is EUKARYOME's own phylogeny-based classification. The test measures agreement with a curated
  taxonomy, not with cultures; the report says so.
- A **named species** is a named genus plus a lower-case epithet from the species field
  (`g__Funneliformis` + `s__mosseae`). A **named genus** has no placeholder marks: not `unclassified`,
  no dots (`Densosporales.fam01.gen01`), no digits. Placeholders stay in the reference as sequences, but are
  never test units and never count as a correct name.
- The inventory (`inventory/`, script v3) was run before this protocol and is the only look at the data
  allowed before the lock.

## 3. Test groups

| Group | Role | Named species with >= 3 reads | Named genera with >= 3 reads |
|---|---|---|---|
| Glomeromycota (phylum) | Known answer from the AMF benchmark | 28 | 41 |
| Agaricales | Richest fungal order | 298 | 187 |
| Hypocreales | Second fungal phylum, many known complexes | 116 | 112 |
| Rhabditida | Non-fungal (nematodes) | 36 | 65 |

Queries are all reads of these four groups. The reference is always the whole long-read set, so a held-out
read can be called as any taxon in EUKARYOME, not only within its group.

## 4. Marker windows

- ITSxRust (version recorded) finds the ITS1, 5.8S and ITS2 boundaries on every read.
- **SSU** = everything upstream of ITS1; **ITS** = ITS1 + 5.8S + ITS2; **LSU** = everything downstream of ITS2.
- A read enters the tests only if all three windows are found and SSU >= 300 nt, ITS >= 300 nt and
  LSU >= 300 nt. Excluded reads are counted per group and reason, and stay in the reference for markers
  that were found.

## 5. Search

- One search per marker: query windows (the four groups) against reference windows (whole set), MMseqs2 in
  nucleotide mode (version and parameters in the run manifest).
- A hit counts only if the alignment covers >= 80% of the shorter of the two windows. Identity = identical
  positions over aligned positions.
- The top 300 hits per query and marker are kept. Every fold re-calls from this table after dropping the
  held-out reads; no fold re-runs the search.
- **Depth check:** after dropping a fold's held-out reads, a query fails the check if fewer than 10 hits are
  left while its original list was full (300 hits). If more than 1% of a group's fold queries fail, that group's
  search is rerun keeping 1,000 hits before any of its results are scored.

## 6. The calls

All thresholds come from training data only: in every fold the held-out taxon's reads are removed before
thresholds are computed.

- **Marker threshold t_m(q):** for marker m, the q-th lower percentile of each training read's identity to
  its nearest other read of the same named species (pooled over the group's usable species).
- **Unique best name:** in one marker, the species of the best hit, provided every hit within 0.2 percentage
  points of it belongs to that same species; otherwise that marker names no species.
- **B1, ITS best hit:** species S if ITS's unique best name is S and its identity >= t_ITS(q); otherwise novel.
- **B2, SINTAX on ITS:** vsearch SINTAX against the same fold's ITS reference; species S if its bootstrap
  support >= c; otherwise novel.
- **L, linked:** species S if SSU, ITS and LSU each give S as their unique best name, and each best identity
  >= t_m(q), with the same q for all three markers; otherwise novel.
- Genus-level calls use the same rules with genera in place of species.
- q and c are swept (q from 0.5 to 50, c from 0.5 to 1.0); comparisons are made at matched known-kept
  rates, so no single threshold is chosen by hand.
- B2 needs its own SINTAX run per fold (about 900 leave-out folds plus 5 known-kept folds per group); that
  cost is accepted, and it may run in batches overnight.

## 7. Folds and metrics

- **Species leave-out:** for each named species with >= 3 reads, all its reads are removed from the
  reference and queried. Each read called as any named species is a **false known**.
- **Genus leave-out:** for each named genus with >= 3 reads, all its reads are removed and queried. A read
  called as any named genus (or a species in one) is a false known at genus level.
- **Known kept:** the reads of usable named species are split into 5 folds within each species (seed
  20261007), with reads identical in all three windows kept in the same fold. Each fold is queried with its
  reads removed from the reference. A read is kept if it is called as its own species. All methods use the
  same split.
- Rates are per read, and also averaged over species (macro), as in the AMF paper.

## 8. K2 decision

- For each group, method and rank, false-known is read off at known-kept = 90% (interpolated along the sweep),
  and also reported at 80% and 95%.
- **R** = 1 - FK(L) / FK(best of B1 and B2), at species level, known-kept 90%.
- **K2 passes** if R >= 0.25 in at least two of the three non-AMF groups (Agaricales, Hypocreales,
  Rhabditida) and FK(L) <= FK(best baseline) in Glomeromycota. Genus-level results are reported but do not
  decide K2.
- If a method cannot reach 90% known-kept in a group, that group is decided at the highest known-kept both
  methods reach, and the report says so.

## 9. K1 decision

Flags, scored on every query read with the full reference (no leave-out):

- **F1, internal conflict:** two markers' best hits fall in different named genera, each with identity
  >= its genus threshold t_m(5) and leading the best hit of any other genus by >= 1 percentage point.
- **F2, label conflict:** all three markers' best hits agree on a named genus other than the read's own label,
  under the same conditions.

Recall, on spike-ins (built before scoring, kept out of every other test):

- 200 chimeras per group: SSU + ITS from a read of genus X joined to the LSU of a read of genus Y in the same
  order; and 200 more with the SSU of X and the ITS + LSU of Y.
- 200 label swaps per group: a read's genus label replaced by another named genus of the same order.
- Recall = chimeras caught by F1, and swaps caught by F2.

Precision, on real flags: 100 flags drawn at random (50 F1, 50 F2, seed 20261007) from all four groups; if
a type has fewer than 50 flags, all of them are checked.
Each is checked by phylogeny: per marker, the read and its 20 best hits are aligned (MAFFT) and a tree is
built (IQ-TREE, UFBoot 1,000). A flag is **real** if the conflicting placement has UFBoot >= 70 (F1: the
markers place in different genera; F2: the read places outside its labelled genus).

- **K1 passes** if recall >= 80% for both flag types and precision >= 50%.

## 10. Rules

- This file, `LOCK.sha256` (hashes of this file and the FASTA) and the inventory are committed together
  before any fold is scored. Code may be written and fixed after the lock; the questions, groups, rules
  and thresholds may not.
- Any change after the lock is a dated amendment, hashed alongside the original. The report gives results
  under the original protocol first.
- Every run writes a manifest: tool versions, parameters, input hashes, and counts at each step.

## 11. Report

Two pages for Leho: K1 and K2 verdicts, false-known against known-kept for each group and method, the depth
check, excluded reads, and every amendment.
