# Phase 0 protocol, note 3: how the K1 tree check reads a tree

Status: LOCKED 2026-10-14, by the hashes in `LOCK_NOTE3.sha256`.
Author: Aaron O'Brien. Adds to section 9 of `PROTOCOL_PHASE0.md` (locked 2026-10-07). No question, group, rule or
threshold changes. Section 9 says a flag is real "if the conflicting placement has UFBoot >= 70", without saying
how a placement is read from a tree; this note fixes that before any flag is computed. No flag, tree or outcome
had been computed when it was written.

## 1. The trees

For each flag in the precision sample (`score_k1.py`, `precision_sample.tsv`: 50 F1 and 50 F2, or all of a type
if fewer) and each marker: the read's window and the windows of its 20 best remaining hits on that marker (the
read itself dropped), aligned with MAFFT `--auto`; tree by IQ-TREE 2 with ModelFinder, 1,000 UFBoot replicates,
one thread, `--seed 20261007` (`scripts/check_flag_trees.py`).

## 2. Placement

The tree is read unrooted, as bipartitions. The read's placement on a marker is the smallest side of a
bipartition with UFBoot >= 70 that holds the read and at least one other sequence. If the other sequences on that
side that have a named genus (section 2) all share one genus G, and there is at least one, the marker places the
read in G; sequences without a named genus are ignored. Otherwise the marker places the read nowhere.

## 3. Real flags and precision

- **F1** is real if two markers place the read in different named genera.
- **F2** is real if all three markers place the read in named genera other than its label.
- Precision is real flags over flags checked, pooled over F1 and F2 for the decision (section 9: >= 50%), and
  reported per type.

## 4. Limit

Each tree holds only the read's 20 best hits on that marker. When the read's labelled genus is not among them,
"outside its labelled genus" means placed, with support, in a clade of another genus.
