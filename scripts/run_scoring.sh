#!/usr/bin/env bash
# Phase 0 scoring, start to finish, once the search has ended. Run from the repository root:
#
#   bash scripts/run_scoring.sh check        # checks only: tools, locks, search outputs, folds, drop list
#   bash scripts/run_scoring.sh all [JOBS]   # checks, then SINTAX, K2, K1, flag trees and the report
#
# Results under the original protocol (section 10; Amendment 1: barrnap-cut reads removed from queries and
# reference) are computed first, then with the amendments. A step whose output already exists is skipped, so a
# stopped run can be restarted with the same command; FORCE=1 reruns every step. Stops at the first failure.
set -euo pipefail

MODE=${1:-check}
JOBS=${2:-16}
FASTA=ref/eukaryome_v2.1/General_EUK_longread_v2.1.fasta
SC=scores
DROP=$SC/drop_barrnap.txt
PY=${PYTHON:-python3}

step() { echo; echo "== $(date '+%F %T')  $*"; }
die() { echo "STOP: $*" >&2; exit 1; }
done_or_force() { [ "${FORCE:-0}" != 1 ] && [ -e "$1" ] && { echo "   (done: $1 exists)"; return 0; }; return 1; }

[ -f protocol/PROTOCOL_PHASE0.md ] || die "run from the repository root"
mkdir -p "$SC"

step "1. tools"
for t in "$PY" vsearch mafft iqtree2; do command -v "$t" > /dev/null || die "$t not on PATH"; done
"$PY" -c 'import numpy, matplotlib' || die "python needs numpy and matplotlib"
echo "   $($PY --version 2>&1); $(vsearch --version 2>&1 | head -1); mafft $(mafft --version 2>&1 | head -1); $(iqtree2 --version 2>&1 | grep -m1 -i version)"

step "2. locks (protocol, amendments, notes and the files they hash)"
for f in protocol/LOCK*.sha256; do
    grep -v '^#' "$f" | sha256sum -c --quiet - || die "$f does not match"
    echo "   $f: OK"
done

step "3. search outputs against search/manifest.txt"
[ -f search/manifest.txt ] || die "search/manifest.txt missing: has the search finished?"
awk -F'\t' '$1=="output"{print $3"  "$2}' search/manifest.txt > "$SC/.search_outputs.sha256"
[ -s "$SC/.search_outputs.sha256" ] || die "no output lines in search/manifest.txt"
sha256sum -c --quiet "$SC/.search_outputs.sha256" || die "a search output differs from its manifest hash"
echo "   $(wc -l < "$SC/.search_outputs.sha256") outputs match"
want=$(awk -F'\t' '$1=="input" && $2=="search/queries.tsv"{print $3}' folds/manifest.txt)
have=$(sha256sum search/queries.tsv | cut -d' ' -f1)
[ "$want" = "$have" ] || die "search/queries.tsv is not the one the folds were built from"
echo "   search/queries.tsv is the one folds/ and spikeins/ were built from"

step "4. drop list for the original protocol: reads cut by barrnap (windows/win.log.tsv.gz)"
zcat windows/win.log.tsv.gz | awk -F'\t' 'NR>1 && $3=="ok" && $4 ~ /^barrnap/ {print $1}' > "$DROP"
echo "   $(wc -l < "$DROP") reads in $DROP"

[ "$MODE" = check ] && { step "checks passed; run with 'all' to score"; exit 0; }
[ "$MODE" = all ] || die "unknown mode '$MODE' (check or all)"

step "5. SINTAX (B2), original protocol, then with Amendment 1"
done_or_force sintax_protocol/manifest.txt || \
    "$PY" scripts/run_sintax.py --windows windows --fasta "$FASTA" --folds folds --queries search/queries.tsv \
        --out sintax_protocol --jobs "$JOBS" --drop-reads "$DROP"
done_or_force sintax/manifest.txt || \
    "$PY" scripts/run_sintax.py --windows windows --fasta "$FASTA" --folds folds --queries search/queries.tsv \
        --out sintax --jobs "$JOBS"

step "6. K2: original protocol; Amendment 1; Amendment 1 with the assembly split (note 2)"
done_or_force $SC/k2_protocol/k2_decision.txt || \
    "$PY" scripts/score_k2.py --search search --fasta "$FASTA" --folds folds --split original \
        --drop-reads "$DROP" --sintax sintax_protocol --out $SC/k2_protocol
done_or_force $SC/k2_amended/k2_decision.txt || \
    "$PY" scripts/score_k2.py --search search --fasta "$FASTA" --folds folds --split original \
        --sintax sintax --out $SC/k2_amended
done_or_force $SC/k2_amended_assembly/k2_decision.txt || \
    "$PY" scripts/score_k2.py --search search --fasta "$FASTA" --folds folds --split assembly \
        --sintax sintax --out $SC/k2_amended_assembly

step "7. K1 flags and spike-ins: original protocol; Amendment 1"
done_or_force $SC/k1_protocol/recall.tsv || \
    "$PY" scripts/score_k1.py --search search --fasta "$FASTA" --folds folds --spikeins spikeins/spikeins.tsv \
        --drop-reads "$DROP" --out $SC/k1_protocol
done_or_force $SC/k1_amended/recall.tsv || \
    "$PY" scripts/score_k1.py --search search --fasta "$FASTA" --folds folds --spikeins spikeins/spikeins.tsv \
        --out $SC/k1_amended

step "8. K1 precision on trees (note 3)"
for run in k1_protocol k1_amended; do
    done_or_force $SC/${run}_trees/precision.tsv || \
        "$PY" scripts/check_flag_trees.py --sample $SC/$run/precision_sample.tsv --windows windows \
            --fasta "$FASTA" --out $SC/${run}_trees --jobs "$JOBS"
done

step "9. report (section 11)"
"$PY" scripts/make_report.py --scores $SC --out report
echo; echo "Done. Report: report/phase0_report.pdf (text: report/phase0_report.md; numbers: report/numbers.tsv)"
