# rdnalink

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23374303.svg)](https://doi.org/10.5281/zenodo.23374303)

Does linking the SSU, ITS and LSU of one long rDNA read give better species calls than the ITS alone?
This repository holds a pre-registered test of that question (Phase 0) on the EUKARYOME v2.1 long-read set:
the protocol, its dated amendments, the code and the fixed inputs. Phase 0 comes before any tool is built: if
K2 fails, the tool narrows to curation (K1), and if both fail, it is not built.

**Status (October 2026):** protocol locked; full search running; nothing scored yet. A prototype package that
runs the same method on any long reads is in `src/rdnalink`.

## The package (prototype)

`src/rdnalink` turns the test's method into a tool that runs on any set of long reads spanning SSU, ITS and LSU.
It does not touch the test: the protocol, folds and scoring scripts stay as locked, and the package's defaults
follow them.

```
conda env create -f environment.yml && conda activate rdnalink && pip install -e .
rdnalink ref reference.fasta -o ref --hmm F.hmm -t 8            # once per reference
rdnalink windows reads.fasta -o reads --hmm F.hmm -t 8
rdnalink search -q reads -r ref -o hits -t 8
rdnalink report -q reads -r ref --hits hits -o report.tsv
```

`F.hmm` is ITSxRust's fungal profile (`data/hmm/F.hmm` in the [ITSxRust](https://github.com/ayobi/itsxrust)
repository; SHA-256 `a8c38c5c...` at v0.3.0).

- `ref` builds a reference folder once: it cuts the windows, merges identical ones, builds one MMseqs2 index per
  marker and writes `manifest.json` (SHA-256 of the reads and the HMM, tool versions, parameters and counts).
  `search` and `report` take the folder wherever they take a window prefix, and a search against it gives the same
  hits as a search against its windows, without rebuilding the indexes each time. `--windows PREFIX` builds the
  folder from windows already cut. On a machine with less than about 10 GB of RAM, add `--kmer 13`; searches of
  the folder then use it too.
- `windows` cuts the windows as Amendment 1 does: ITSxRust where it finds the full ITS, barrnap with the join and
  18S-end rules for the rest, in chunks of 4,000 reads.
- `search` runs one exhaustive MMseqs2 search per window, as note 2 describes. Hits are ranked by identity, then
  identical positions, then target ID; bit scores and E-values are kept but not used (MMseqs2's bit scores can
  differ by a few bits between runs of the same search).
- `report` writes one row per read. For each marker it gives the best identity, the unique best species and genus,
  and the genus's lead. It then gives the linked species and genus (all three markers agree), the ITS-only
  species, a conflict flag (two markers confidently name different genera), a label conflict (all three
  confidently name a genus other than the read's label), and a paralog flag for the short-5.8S type of
  *Rhizophagus*.

Identity thresholds are off by default (`--min-identity`). The test will say which ones to use.

Checked so far (functional checks, not the test):

- On the 3,994 pilot reads, `windows` gives every read the status, strand and ITS bounds Amendment 1 fixed (3,774
  cut by ITSxRust, 197 by barrnap, 23 rejected), running the real ITSxRust 0.3.0 and barrnap 1.10.5 on macOS.
- On 400 pilot reads outside the four test groups, a search against a `ref` folder and a search against its
  windows give the same hits, ranks and identities, row for row.
- On the 903 culture copies of the AMF benchmark, `windows` cut every copy, and its ITS windows are identical to
  the benchmark's full ITS. `report` flagged exactly the 41 short-5.8S copies of *Rhizophagus*. In 20 synthetic
  chimeras (SSU and ITS of one genus, LSU of another), `report` flagged all 20 as conflicts.

`pytest tests` runs the unit tests, and with `mmseqs` on the PATH also checks a `ref` search against a plain one.
`RDNALINK_HMM=/path/to/F.hmm pytest tests/test_pilot_windows.py` runs the real ITSxRust and barrnap on the 3,994
pilot reads and checks every read's windows against Amendment 1's pilot log (about 8 minutes on a laptop; windows
only, nothing scored).

## Questions

- **K2 (main).** When a read comes from a taxon the reference does not hold, does a call that needs SSU, ITS
  and LSU from the same read to name the same species (L) call it a known taxon less often than ITS alone (best
  hit, B1, or SINTAX, B2), at the same rate of keeping known taxa? K2 passes if, at 90% known-kept, L makes at
  least 25% fewer false-known calls than the better of B1 and B2 in at least two of Agaricales, Hypocreales and
  Rhabditida, and no more in Glomeromycota.
- **K1.** When the three markers of one read point to different genera, or agree on a genus other than the
  read's label, is that usually a real chimera or mislabel? K1 passes if spiked chimeras and label swaps are
  caught at least 80% of the time, and at least half of a random sample of real flags hold up on trees.

## Test groups

| Group | Query reads | Role |
|---|---|---|
| Glomeromycota | 5,546 | Known answer from the AMF culture benchmark |
| Agaricales | 9,436 | Richest fungal order |
| Hypocreales | 3,318 | Second fungal phylum, many known complexes |
| Rhabditida | 2,901 | Non-fungal (nematodes) |

Queries are reads with all three marker windows of at least 300 nt. The reference is always the whole long-read
set.

## Pre-registration

Every change after the lock is a dated document, hashed in a `LOCK*.sha256` file and committed before the step
it governs. The commit history is the record of when each was fixed.

| Date | Document | Fixes |
|---|---|---|
| 2026-10-07 | `protocol/PROTOCOL_PHASE0.md` | Questions, groups, rules and thresholds |
| 2026-10-08 | `protocol/AMENDMENT1_windows.md` | How the marker windows are cut (ITSxRust, barrnap fallback) |
| 2026-10-10 | `protocol/AMENDMENT1_note1.md` | Window tools run in chunks |
| 2026-10-14 | `protocol/NOTE2_search.md` | Exhaustive search, folds, K1 spike-ins, scoring definitions |
| 2026-10-14 | `protocol/NOTE3_k1_trees.md` | How the K1 tree check reads a tree |

To check the locks, put the EUKARYOME file in `ref/eukaryome_v2.1/` and run
`sha256sum -c protocol/LOCK*.sha256` from the repository root (`shasum -a 256 -c` on macOS).

## Layout

| Path | What it holds |
|---|---|
| `protocol/` | Protocol, amendments, notes and their lock files; `eukaryome_inventory.py`, the only look at the data before the lock |
| `inventory/` | Output of the inventory |
| `pilot/` | The 3,994-read pilot of the window tools and its reports |
| `windows/` | Manifest, per-read cut log and reports of the window cut (the window FASTAs are not in git; their hashes are in `windows/manifest.txt`) |
| `search/queries.tsv` | The query reads, their group and window lengths |
| `folds/`, `spikeins/` | Leave-out units, known-kept split and K1 spike-ins, fixed by note 2 |
| `scripts/` | Window cut, search, folds, spike-ins and scoring |
| `src/rdnalink/`, `tests/`, `pyproject.toml`, `environment.yml` | The package (prototype) |

## Pipeline

1. `cut_windows.py` cuts SSU, ITS and LSU windows (Amendment 1).
2. `search_markers.py` runs one MMseqs2 search per marker. The version committed here is the first one; the
   exhaustive version of note 2 is committed with the full run's manifest when the run ends.
3. `make_folds.py` and `make_spikeins.py` fix the folds and the K1 spike-ins.
4. `score_k2.py` scores B1 and L, and adds B2 from `run_sintax.py`; `score_k1.py` scores the flags and draws
   the precision sample; `check_flag_trees.py` checks it on trees. Each takes `--drop-reads` for the original
   protocol (Amendment 1: the barrnap-cut reads leave queries and reference; they are skipped in the hit lists,
   as a fold's held-out reads are, so ranks never move up).
5. `make_report.py` writes the two-page report of section 11 (`report/phase0_report.pdf`, with a text copy and
   `numbers.tsv`, every number with the file it came from).

`bash scripts/run_scoring.sh check` checks the tools, every lock, the search outputs against their manifest and
the queries the folds were built from, and writes the drop list; `bash scripts/run_scoring.sh all [JOBS]` then
runs steps 4 and 5, original protocol first, and can be restarted where it stopped.

Tool versions, parameters and input and output hashes are written to a `manifest.txt` beside every output.

## Data

- EUKARYOME v2.1 long-read set, `General_EUK_longread_v2.1.fasta` (400,217 reads), CC BY 4.0: Tedersoo et al.
  (2024) EUKARYOME: the rRNA gene reference database for identification of all eukaryotes. *Database* 2024,
  baae043, doi:10.1093/database/baae043. Not included; download
  `https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_longread_v2.1.zip` (the FASTA is inside a
  `.7z` archive in it). Its SHA-256 is in `protocol/LOCK.sha256`.
- `pilot/pilot.fasta` holds 3,994 reads of that set, redistributed under the same licence.

## Related

- [amf-novelty-benchmark](https://github.com/ayobi/amf-novelty-benchmark): the AMF culture benchmark behind the
  Glomeromycota group.
- [ITSxRust](https://github.com/ayobi/itsxrust): the ITS boundary finder used to cut the windows.
- EUKARYOME's curation already compares markers: for (meta)genome and (meta)transcriptome data it removes taxa
  whose per-marker BLASTn matches conflict at genus, order or kingdom level (Tedersoo et al. 2024). The conflict
  flag applies the same idea to any long read, with stated thresholds, and K1 measures how often it is right.
  `protocol/NOTE4_eukaryome_checks.md` records what EUKARYOME's documented checks cover.
- Hakimzadeh et al. (2025) Are we throwing away good data? Evaluation of chimera detection algorithms on long-read
  amplicons reveals high false-positive rates across algorithms. *PeerJ* 13:e20456. Of the full-ITS reads each
  de novo tool flagged, 9.8% (uchime_denovo) to 55.0% (chimeras_denovo) matched a EUKARYOME reference at >= 99%
  identity and coverage.

## Licence

To cite: `CITATION.cff`, or the Zenodo record above. Code: MIT (`LICENSE`). Protocols, tables and other outputs: CC BY 4.0 (`DATA_LICENSE.md`). EUKARYOME-derived
reads keep EUKARYOME's CC BY 4.0 terms; cite Tedersoo et al. (2024) when using them.

Aaron O'Brien
