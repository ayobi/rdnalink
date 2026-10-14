# Phase 0 protocol, note 2: the search, the folds, the K1 spike-ins and the scoring definitions

Status: LOCKED 2026-10-14, by the hashes in `LOCK_NOTE2.sha256`.
Author: Aaron O'Brien. Adds to sections 5-9 of `PROTOCOL_PHASE0.md` (locked 2026-10-07) where stated. No
question, group, rule or threshold changes. It is recorded as a dated note (section 10) because it fixes choices
that change results: what the search stores and how scoring reads it, the folds, the spike-ins, and the points
sections 6-9 leave open. No fold of the test had been scored when it was written. What was looked at is listed in
section 8.

## 1. Why the search changed

Section 5 fixes the hit filter (alignment covering >= 80% of the shorter window), identity (identical positions
over alignment columns) and the top 300 hits, and leaves the search parameters to the run manifest. A pilot of
400 queries (the first 100 of each group, 13 October 2026) showed that the parameters change the top 300:

- Capping the candidates MMseqs2 passes from its k-mer prefilter to alignment (`--max-seqs` 2,000, 10,000 and
  30,000) lost hits that belong in the top 300 at every cap, and each larger cap found more.
- Against an exhaustive search, the 30,000 cap left every call input the same (best identity and unique best
  name; identity to the nearest other read of the query's species; best identity and name after species
  leave-out, and after genus leave-out) for only 2-88% of queries on SSU, 37-100% on ITS and 19-83% on LSU,
  depending on the group.
- MMseqs2's own no-prefilter mode (`--prefilter-mode 2`) aligns every pair on the zero diagonal and misaligns
  windows that start at different offsets, so it is not used.

## 2. The search (full run started 13 October 2026)

- `scripts/search_markers.py` (hash in `search/manifest.txt`), MMseqs2 18-8cc5c, nucleotide (`--search-type 3
  -s 5.7 --strand 1`), exhaustive: `--max-seqs` and `--max-accept` equal the number of distinct target
  windows, so every candidate with a k-mer hit is aligned. No E-value filter (`-e 1e6`).
- Identical windows are searched once and their hits copied to every read that carries them. Identical sequences
  give identical alignments, so this changes no hit.
- Hits: coverage >= 80% of the shorter window; identity = identical positions / alignment columns. Hits are
  ranked by identity, then identical positions, then target read ID, so ranks never depend on bit scores or
  E-values, which depend on database composition.
- Stored per query and marker: the top 5,000 hits, every hit within 0.2 points of the best (up to 20,000), and
  every hit tied with the last one stored (`--keep 5000 --keep-within 0.2 --max-keep 20000`).

## 3. What scoring reads from the search

- **Top 300.** The calls of sections 6-9 read ranks <= 300 of this table, as section 5 says; the ranks count the
  self hit. A fold drops its held-out reads from those 300; ranks below 300 never move up.
- **No self hit.** A read is never its own hit. In every fold the query is one of the held-out reads; in K1, "full
  reference" (section 9) means every read except the query itself. A chimera spike-in drops both of its source
  reads, and a label swap drops its own read.
- **Depth check.** Run as section 5 says, over every fold query of a group (leave-out and known-kept), counting a
  query once if any marker fails. Where a group fails, the rerun keeping 1,000 hits is replaced by reading ranks
  <= 1,000 of this table for every result of that group: the search is exhaustive and deterministic, so a rerun
  would return the same ranks. The report says which groups failed and by how much.
- **Near-ties at the last rank.** If the hits within 0.2 points of the best run to the last rank read, section 6
  is applied to the ranks as written. The report gives how many fold calls this affects.

## 4. Folds (section 7)

Fixed before scoring by `scripts/make_folds.py` and hashed with this note (`folds/units.tsv`,
`folds/known_kept.tsv`):

- **Usable** species and genera are named ones (section 2) with >= 3 query reads (reads with all three windows).
- **Leave-out folds** remove every read of the species (genus) in the whole long-read set from the reference, and
  query the group's query reads of it.
- **Known-kept split.** Within each usable species, reads identical in all three windows form one unit. Units are
  shuffled with `random.Random('20261007-<group>-<species>')` and each goes to the fold with the fewest reads so
  far (lowest fold number on a tie).
- **Copies of one assembly (change).** EUKARYOME holds several rDNA copies of some genome assemblies as separate
  reads (`ID`, `ID_dup1`, `ID_dup2`, ...). Copies that differ slightly would otherwise stay in the reference while
  one of them is queried, so a read could be kept because a copy from its own genome is still there. In the
  second split (`fold_assembly`) the copies of one assembly are one unit, and every copy of it, query or not,
  leaves the reference with its fold. Query reads affected: Agaricales 714, Hypocreales 101, Rhabditida 90,
  Glomeromycota 0. As section 10 requires, the report gives known-kept results under the original split
  (`fold_original`) first.

## 5. K1 spike-ins (section 9)

Built by `scripts/make_spikeins.py` (seed 20261007) from `search/queries.tsv` and the read labels, and hashed with
this note: `spikeins/spikeins.tsv`, 2,400 spike-ins. Per group:

- 200 chimeras with SSU + ITS of read X and LSU of read Y, and 200 with SSU of X and ITS + LSU of Y. X and Y are
  query reads of the group with a named genus in a named order; Y is in X's order and another named genus.
- 200 label swaps: read X with its genus label replaced by another named genus of its order, drawn from the named
  genera of that order in the whole long-read set.
- In each set the X reads are drawn without replacement. Every window of a spike-in is a query read's window, so
  no new search is needed.

## 6. Scoring definitions (sections 6-9)

Points the protocol leaves open, fixed here as `scripts/score_k2.py`, `run_sintax.py` and `score_k1.py` implement
them. Code may still be fixed (section 10); these definitions may not.

- **Thresholds.** t_m(q) is the q-th lower percentile (numpy, linear interpolation) of the training reads' values,
  pooled over the group's usable species (genus rank: usable genera). Training reads are the group's query reads
  of usable species (genera) that the fold does not hold out. A read's value is the identity of its nearest other
  read of the same name among its remaining ranks; if there is none, it is the identity at its last rank when
  its list was full (an upper bound) and 0 when it was not (no such read aligns).
- **Sweeps.** q from 0.5 to 50 in steps of 0.5; c from 0.50 to 1.00 in steps of 0.01.
- **Unique best name.** A hit to a read without a named species (genus) carries no name: as best hit it gives no
  name, and within 0.2 points of the best it blocks the name.
- **Genus rank.** Known-kept uses the species split; a read is kept if called as its own genus. A genus-fold read
  is a false known if its genus call or its species call names anything (section 7).
- **Macro averages.** False-known rates are averaged over the held-out species (genera), known-kept rates over
  species.
- **B2.** vsearch SINTAX (`--strand plus --threads 1 --randseed 20261007`, default settings otherwise) against the
  fold's ITS windows labelled with EUKARYOME's own taxonomy (`d:` to `s:<genus>_<epithet>`). Placeholder labels
  stay in the reference; a call to one is no name.
- **K2 decision.** On per-read rates at species level; macro rates are reported beside them. Where L, B1 or B2
  cannot reach 90% known-kept in a group, that group is read at the highest known-kept all three reach. R is
  undefined when the best baseline has no false knowns, and that group does not count as passing.
- **K1 flags.** A marker is confident for genus G when its best remaining hit is in the named genus G, at
  identity >= t_m(5) (genus rank, full reference minus the read), and leads the best remaining hit of any other
  genus, named or not, by >= 1 point. With no hit of another genus left, the lead is the best identity minus the
  identity at the last rank if the list was full, and unlimited if it was not. F1 needs two confident markers with
  different genera; F2 needs all three confident for one genus other than the read's label, and applies only to
  reads labelled with a named genus. Recall is pooled over the four groups for the decision and reported per
  group. The precision sample is drawn with `random.Random(20261007)`.

## 7. Pilot depth check (400 queries, exhaustive, 1,000 hits stored)

`scripts/depth_check.py`, self hits dropped. Counts are of the pilot's 100 queries per group, not of fold queries,
and of queries with no hit left rather than fewer than 10, so they only show the scale.

- Near-ties past rank 300: none in any group or marker.
- After dropping the query's genus, nothing left in the top 300 but hits deeper in the 1,000: Glomeromycota ITS
  2 of 83 named-genus queries, Glomeromycota LSU 2 of 83, no other group or marker.
- Rhabditida ITS lists are short: 17 of 100 queries have no hit but themselves, and the rest a median of 10,
  because few other reads align over 80% of these windows. A deeper list adds nothing there.

## 8. What was looked at before this note was locked

- The search pilots above, and the read counts of folds and spike-ins.
- The scoring scripts were tested on the 400-query pilot. All tests but one used labels permuted at random across
  reads, so their outputs say nothing about the test. One run of `score_k2.py` used the real labels on the pilot
  (14 October 2026); restricted to the 400 pilot queries, it had 1-7 usable species per group, so it reproduced
  no fold of the test. After it, the known-kept split moved into `make_folds.py` (with the per-species seeds
  above), and B2, the K2 decision and the K1 flags were added; the threshold, unique-best-name, genus-rank and
  macro definitions were already in the script and were not changed.

## 9. Files

`search/manifest.txt`; `folds/units.tsv`, `folds/known_kept.tsv` and `folds/manifest.txt`; `spikeins/spikeins.tsv`
and `spikeins/manifest.txt`; `depth_check_pilot.tsv`.
