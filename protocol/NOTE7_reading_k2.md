# Phase 0 protocol, note 7: reading K2 when no known-kept level is shared

Status: LOCKED 2026-10-15, by the hashes in `LOCK_NOTE7.sha256`.
Author: Aaron O'Brien. Written after scoring, on 15 October 2026, once the results were known. It changes no
number and no rule. It records how the report reads K2, because the case below was not foreseen by section 8 or
note 2, and it says so wherever the reading is used.

## What happened

Section 8 reads false known at 90% known kept, or, where a method cannot reach that, at the highest known kept
all methods reach; section 6 says the methods are compared at matched known-kept rates. At species level, per
read, the known-kept range of L lies wholly below those of B1 and B2 in every group and in every run:

| Group | L | B1 | B2 |
|---|---|---|---|
| Glomeromycota | 0.0-0.0% | 41.7-62.5% | 10.0-64.4% |
| Agaricales | 0.1-0.2% | 47.4-73.2% | 46.5-80.2% |
| Hypocreales | 1.9-2.9% | 41.1-50.6% | 21.6-52.9% |
| Rhabditida | 18.5-40.7% | 52.9-83.1% | 56.1-86.8% |

(original protocol, `scores/k2_protocol/sweep.tsv`; the amended runs are alike, `report/numbers.tsv`). No
known-kept level is reached by all three methods, so no matched comparison exists. `score_k2.py` then reads each
baseline at its own lowest known kept (`at_kept`: "lowest known-kept above target") and L at its highest, and
`k2_decision.txt` reports PASS from these unmatched values. L reaches its highest known kept at the loosest
threshold (q = 0.5), so no threshold would close the gap.

## How the report reads K2

The report gives the decision as computed (PASS) and, beside it, the reading at matched rates: in no group can
false known be compared at a known-kept rate all methods reach, so no group meets section 8 and K2 is not
supported. This follows note 2's treatment of an undefined R, under which a group does not count as passing. The
report's third page, added after scoring, shows the ranges and why K1 misses, read from the same score tables.

## Lock

`LOCK_NOTE7.sha256` holds the hashes of this note and of `scripts/make_report.py`, which prints the reading.
