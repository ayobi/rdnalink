# rdnalink

Does linking the SSU, ITS and LSU of one long rDNA read give better species calls than the ITS alone?
This repository holds a pre-registered test of that question (Phase 0) on the EUKARYOME v2.1 long-read set:
the protocol, its dated amendments, the code and the fixed inputs. Phase 0 comes before any tool is built: if
K2 fails, the tool narrows to curation (K1), and if both fail, it is not built.

**Status (October 2026):** protocol locked; full search running; nothing scored yet.

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

## Pipeline

1. `cut_windows.py` cuts SSU, ITS and LSU windows (Amendment 1).
2. `search_markers.py` runs one MMseqs2 search per marker. The version committed here is the first one; the
   exhaustive version of note 2 is committed with the full run's manifest when the run ends.
3. `make_folds.py` and `make_spikeins.py` fix the folds and the K1 spike-ins.
4. `score_k2.py` scores B1 and L, and adds B2 from `run_sintax.py`; `score_k1.py` scores the flags and draws
   the precision sample; `check_flag_trees.py` checks it on trees.

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

## Licence

Code: MIT (`LICENSE`). Protocols, tables and other outputs: CC BY 4.0 (`DATA_LICENSE.md`). EUKARYOME-derived
reads keep EUKARYOME's CC BY 4.0 terms; cite Tedersoo et al. (2024) when using them.

Aaron O'Brien
