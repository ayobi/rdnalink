# Phase 0 protocol, note 5: code fixes before scoring

Status: LOCKED 2026-10-14, by the hashes in `LOCK_NOTE5.sha256`.
Author: Aaron O'Brien. Records two fixes to the scoring code, made after notes 2 and 3 and before any fold or
flag was scored (section 10: code may be fixed after the lock). No question, group, rule, threshold or definition
changes; each fix makes the code do what a locked note already says. Both were found in tests on the pilot with
labels permuted at random across reads, which say nothing about the outcome.

## 1. Too few sequences for a tree (`check_flag_trees.py`, locked by note 3)

Note 3 builds, for each sampled flag and marker, a tree of the read and its 20 best remaining hits. A marker can
have fewer than three remaining hits (Rhabditida ITS lists are short; note 2, section 7). IQ-TREE stops with an
error on an alignment of fewer than three sequences, and that would stop the whole check. An unrooted tree of
three or fewer sequences has no internal branch, so it has no bipartition with UFBoot support, and note 3 section 2
already says the marker then places the read nowhere. The fixed script builds no tree for fewer than four
sequences and records that marker as placing the read nowhere (`placements.tsv`, clade size 0). Trees of four or
more sequences are built and read as before.

## 2. Reads removed for the original protocol (`score_k2.py`, `score_k1.py`; not hashed by any lock)

Amendment 1 obtains results under the original protocol by removing the barrnap-cut reads from queries and
reference, filtering the search table rather than rerunning it, as for folds; note 2 says ranks below 300 never
move up. `score_k2.py --drop-reads` deleted those reads from the hit lists before cutting them at the ranks read,
so lower ranks moved up, and a full list that lost a dropped read was read as not full. `score_k1.py` had no
`--drop-reads`. Both now keep dropped reads in the lists and skip them wherever the lists are read, exactly as a
fold's held-out reads are skipped, and `score_k1.py` takes `--drop-reads`. With no reads dropped, every output is
byte-identical to the earlier version (24 output files compared on the permuted-label pilot).

## Lock

`LOCK_NOTE5.sha256` holds the hashes of this note and of the three scripts as fixed. `LOCK_NOTE3.sha256` keeps the
hash of the version note 3 locked. `scripts/run_scoring.sh` checks each file against the newest lock that lists
it; only scripts may be listed again by a later lock.
