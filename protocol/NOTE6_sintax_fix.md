# Phase 0 protocol, note 6: a code fix during the SINTAX run

Status: LOCKED 2026-10-15, by the hashes in `LOCK_NOTE6.sha256`.
Author: Aaron O'Brien. Records a fix to `scripts/run_sintax.py`, made after the scoring run had started and
before any K1 or K2 score was computed (section 10: code may be fixed after the lock). No question, group, rule,
threshold or definition changes.

## The fix

The scoring run of 15 October (`scripts/run_scoring.sh all`) stopped in its first SINTAX pass, the one under the
original protocol (`--drop-reads`: the 31,563 barrnap-cut reads), on a dropped read. For the known-kept folds,
`run_sintax.py` took the query reads from `folds/known_kept.tsv` without removing the dropped reads, which have no
ITS window in that pass. Amendment 1 removes them from queries and reference, and the leave-out folds of
`run_sintax.py` and every fold of `score_k2.py` already did so. They now leave the known-kept queries as well.

Folds that had finished before the stop are kept, not rerun: the fix leaves them unchanged. On the pilot with
labels permuted at random across reads, all 352 folds the earlier script finished before stopping are
byte-identical with the fixed script. No SINTAX output had been read and no score computed when the fix was made.

## Lock

`LOCK_NOTE6.sha256` holds the hashes of this note and of the fixed script.
