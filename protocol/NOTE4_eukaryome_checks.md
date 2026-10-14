# Phase 0 protocol, note 4: what EUKARYOME's own checks cover, and how the report describes real K1 flags

Status: LOCKED 2026-10-14, by the hashes in `LOCK_NOTE4.sha256`.
Author: Aaron O'Brien. Adds to sections 9 and 11 of `PROTOCOL_PHASE0.md` (locked 2026-10-07). No question, group,
rule or threshold changes. Real K1 flags are reads that EUKARYOME's curation kept, so what the report may say
about them depends on what that curation checks. This note records what is documented and fixes the wording
before any flag is computed. No flag, tree or outcome had been computed when it was written.

## 1. What is documented

Sources: Tedersoo et al. (2024), *Database* baae043, which describes EUKARYOME v1.8, with its Supplementary
Text S1 (instructions to annotators and curators) and Table S1. The EUKARYOME website (download and news pages,
read 14 October 2026) describes no curation steps.

- **(Meta)genome and (meta)transcriptome data:** UCHIME on each marker separately, and BLASTn assignments of
  each marker compared against an earlier EUKARYOME version. Taxa with conflicting matches "at genus, order and
  kingdom levels—depending on taxonomic resolution—were removed as potential chimeras". No identity or coverage
  threshold is given, and no rule for which level applies.
- **PacBio long reads** (OTU representatives of published studies): a further round of ITSx and UCHIME v4.2,
  default options, against UNITE, PR2 and SILVA; an indel filter on conserved regions removed 15-25% of them.
- **Nanopore consensus reads:** processed "as described earlier", with 1% excluded.
- **All data:** phylum- or kingdom-level IQ-TREE trees from one alignment of SSU, 5.8S and LSU (Text S1
  recommends removing ITS1 and ITS2), read by curators. Text S1 asks them to record misidentified taxa and,
  optionally ("If you do not feel qualified, do not bother"), suspected chimeras, which "typically have a very
  long unique branch compared with the closest relative".
- **INSDc and JGI sequences:** Table S1 lists 4,256 accessions marked chimeric (508) or low quality (3,748); the
  method is not given.
- **Size:** v1.8 held 32,810 combined SSU-ITS-LSU reads; the v2.1 long-read set holds 400,217. The documents do
  not say how reads added after v1.8 were checked.

So the per-marker comparison closest to F1 is documented only for the (meta)genome and (meta)transcriptome data.
For the long-read amplicons, the documented checks are UCHIME, the indel filter and curators reading
combined-marker trees.

## 2. How the report describes real flags

- A real flag (section 9, note 3) is reported as a read whose markers are placed, on per-marker trees, in
  different genera (F1) or outside its labelled genus (F2). The report does not call real flags chimeras that
  EUKARYOME's checks missed, and does not read few real flags as a sign that those checks already removed them,
  unless EUKARYOME's maintainers state which checks the flagged reads went through.
- If the maintainers say which checks the long reads went through, and at which thresholds, the report gives
  their answer beside the results, with the date it arrived and whether that was before or after scoring.
- Real flags are listed with their accessions in a table for EUKARYOME's maintainers.
- No expectation about the number of real flags is fixed here.
