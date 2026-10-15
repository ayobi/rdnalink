# rdnalink Phase 0: results

EUKARYOME v2.1 long reads; protocol locked 2026-10-07. Report built 2026-10-15 by scripts/make_report.py.

## Verdicts

```
Verdicts (results under the original protocol first, section 10)
  K2  Original protocol (barrnap-cut reads removed, section 7 split): PASS
  K2  Amendment 1 windows, section 7 split: PASS
  K2  Amendment 1 windows, note 2 assembly split: PASS
  K1  Original protocol: FAIL (recall: chimeras 29.8%, swaps 13.5%; precision 16.0%)
  K1  Amendment 1 windows: FAIL (recall: chimeras 31.1%, swaps 14.9%; precision 32.0%)
  K2  at matched known-kept rates (note 7): not supported; no shared level in 4 of 4 groups
```

K2 passes if, at 90% known-kept (or the highest level all methods reach), R = 1 - FK(L)/FK(best of B1, B2) >= 0.25 in at least two of Agaricales, Hypocreales and Rhabditida, and FK(L) <= FK(best) in Glomeromycota (species level, per-read rates). K1 passes if spiked chimeras and label swaps are each caught >= 80% of the time (pooled) and >= 50% of a random sample of real flags hold up on trees.

## K2

```
K2, original protocol: false known (FK, % of reads) at the known-kept level, per read, species level
group           at KK %   FK B1   FK B2    FK L       R          meets rule
Glomeromycota         0     2.3     0.0     0.0       -  FK(L) <= best: yes
Agaricales            0     3.3     8.4     0.0    1.00      R >= 0.25: yes
Hypocreales           3    13.2     8.2     0.0    1.00      R >= 0.25: yes
Rhabditida           41     2.6     7.9     5.8   -1.20       R >= 0.25: no
```

```
Amendment 1 windows, section 7 split: PASS
group           at KK %   FK B1   FK B2    FK L       R          meets rule
Glomeromycota         0     2.3     0.0     0.0       -  FK(L) <= best: yes
Agaricales            0     3.3     8.4     0.0    1.00      R >= 0.25: yes
Hypocreales           2    13.2     8.2     0.0    1.00      R >= 0.25: yes
Rhabditida           40     2.5     6.5     5.5   -1.20       R >= 0.25: no

Amendment 1 windows, note 2 assembly split: PASS
group           at KK %   FK B1   FK B2    FK L       R          meets rule
Glomeromycota         0     2.3     0.0     0.0       -  FK(L) <= best: yes
Agaricales            0     3.3     8.4     0.0    1.00      R >= 0.25: yes
Hypocreales           1    13.2     8.2     0.0    1.00      R >= 0.25: yes
Rhabditida           40     2.5     6.5     5.5   -1.20       R >= 0.25: no
```

```
Original protocol, rates averaged over species (macro; reported, not deciding)
group           at KK %   FK B1   FK B2    FK L       R          meets rule
Glomeromycota         0     2.4     0.0     0.0       -  FK(L) <= best: yes
Agaricales            0     3.5     9.4     0.0    1.00      R >= 0.25: yes
Hypocreales           3    14.8     7.9     0.0    1.00      R >= 0.25: yes
Rhabditida           46     1.9    13.1     8.1   -3.25       R >= 0.25: no
```

## K1

```
K1, Original protocol: FAIL
group                chimeras caught          swaps caught   real F1   real F2
Glomeromycota        157/399 (39.3%)          2/200 (1.0%)         3         0
Agaricales             31/400 (7.8%)          0/200 (0.0%)        96         0
Hypocreales           46/400 (11.5%)          5/200 (2.5%)        40         0
Rhabditida           225/342 (65.8%)        99/183 (54.1%)         2         0
all                 459/1541 (29.8%)       106/783 (13.5%)       141         0
Precision on trees (real/checked): F1 8/50 (16.0%); F2 0/0 (-%); all 8/50 (16.0%)

K1, Amendment 1 windows: FAIL
group                chimeras caught          swaps caught   real F1   real F2
Glomeromycota        156/400 (39.0%)          2/200 (1.0%)         7         0
Agaricales             32/400 (8.0%)          0/200 (0.0%)       115         0
Hypocreales           44/400 (11.0%)          1/200 (0.5%)        43         0
Rhabditida           266/400 (66.5%)       116/200 (58.0%)         2         0
all                 498/1600 (31.1%)       119/800 (14.9%)       167         0
Precision on trees (real/checked): F1 16/50 (32.0%); F2 0/0 (-%); all 16/50 (32.0%)
```

Real flags are reads whose markers conflict on per-marker trees. EUKARYOME documents a per-marker BLASTn check for its (meta)genome and (meta)transcriptome data only, so they are not described as chimeras its checks missed (note 4). Their accessions are in scores/k1_*/flags.tsv.

## After scoring (note 7)

Added after scoring (note 7). These analyses were written once the results were known; they change no number in the pages before, and are read from the same score tables.

Section 8 reads false known at 90% known kept, or at the highest known kept all methods reach, so that methods are compared at matched known-kept rates (section 6). At species level, the known-kept range of L lies wholly below those of B1 and B2 in every group (table): no level is reached by all three, so no matched comparison exists. score_k2.py then reads each baseline at its own lowest known kept and L at its highest, and k2_decision.txt reports PASS from these unmatched values. Read at matched rates, no group meets section 8, and K2 is not supported. L stays this low at every q, so its threshold is not the limit: L names a species only when SSU, ITS and LSU each give the same unique best species, and the SSU and LSU rarely give the same one as the ITS.

```
Known-kept range (%) of each method over its sweep, species level, per read, original protocol
group                      L            B1            B2   shared level
Glomeromycota        0.0-0.0     41.7-62.5     10.0-64.4           none
Agaricales           0.1-0.2     47.4-73.2     46.5-80.2           none
Hypocreales          1.9-2.9     41.1-50.6     21.6-52.9           none
Rhabditida         18.5-40.7     52.9-83.1     56.1-86.8           none
```

K1 misses on recall in every group. A chimera is flagged when two markers name different genera with confidence. The SSU does so for few fungal spike-ins, so chimeras that join the SSU of one read to the ITS and LSU of another (chim_b) are almost never caught in fungi; in Agaricales and Hypocreales the LSU is also rarely confident, which limits chimeras that swap only the LSU (chim_a). Of the 50 sampled real flags checked on trees, 27 had at least two markers placed in a supported clade, the minimum for a flag to be real.

```
Spike-ins caught, by kind, and share of spike-ins on which each marker names a genus with confidence (original protocol)
group              chim_a     chim_b      swaps  confident: SSU / ITS / LSU
Glomeromycota     152/199      5/200      2/200              1% / 92% / 82%
Agaricales         30/200      1/200      0/200              0% / 78% / 20%
Hypocreales        39/200      7/200      5/200              6% / 73% / 22%
Rhabditida        125/172    100/170     99/183             60% / 89% / 76%
```

## Checks and exclusions

```
Depth check (note 2: a group with more than 1% of fold queries failing reads ranks <= 1,000)
  Glomeromycota: 2721 of 5620 fold queries fail (48.42%); ranks <= 1000
  Agaricales: 628 of 10598 fold queries fail (5.93%); ranks <= 1000
  Hypocreales: 214 of 3602 fold queries fail (5.94%); ranks <= 1000
  Rhabditida: 0 of 2003 fold queries fail (0.00%); ranks <= 300
Fold calls whose 0.2-point set reaches the last rank read (note 2): 0
```

Reads and windows: 400,217 reads; 397,441 cut into three windows, 31,563 of them by barrnap (removed under the original protocol). Not cut: missing_LSU 1,296; missing_SSU 684; strand_or_order 439; missing_SSU_LSU 344; multi_operon 13. Query reads per group: Glomeromycota 5,546; Agaricales 9,436; Hypocreales 3,318; Rhabditida 2,901.

```
Protocol, amendments and notes, with the dates they were locked
  AMENDMENT1_note1.md (locked 2026-10-10): Phase 0 protocol, Amendment 1, note 1: ITSxRust and barrnap run in
      chunks
  AMENDMENT1_windows.md (locked 2026-10-08): Phase 0 protocol, Amendment 1: how the marker windows are cut
  NOTE2_search.md (locked 2026-10-14): Phase 0 protocol, note 2: the search, the folds, the K1 spike-ins and the
      scoring definitions
  NOTE3_k1_trees.md (locked 2026-10-14): Phase 0 protocol, note 3: how the K1 tree check reads a tree
  NOTE4_eukaryome_checks.md (locked 2026-10-14): Phase 0 protocol, note 4: what EUKARYOME's own checks cover,
      and how the report describes real K1 flags
  NOTE5_code_fixes.md (locked 2026-10-14): Phase 0 protocol, note 5: code fixes before scoring
  NOTE6_sintax_fix.md (locked 2026-10-15): Phase 0 protocol, note 6: a code fix during the SINTAX run
  NOTE7_reading_k2.md (locked 2026-10-15): Phase 0 protocol, note 7: reading K2 when no known-kept level is
      shared
  PROTOCOL_PHASE0.md (locked 2026-10-07): Phase 0 protocol: do linked SSU-ITS-LSU calls beat ITS alone outside
      AMF?
```

Every number above is in numbers.tsv with the file it came from.
