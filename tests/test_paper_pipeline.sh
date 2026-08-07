#!/bin/sh

set -eu

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
experiment="$root_dir/Barcode_Experiment"
distance="$root_dir/Barcode_Distance"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/barcode-paper.XXXXXX")
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

count_files() {
    find . ! -name . -prune -type f "$@" -print | awk 'END { print NR + 0 }'
}

cd "$test_dir"

# Section 2.3 used every combination of A/C/E/L/N with ten replicates. Keep
# that complete 32 x 10 design here, but reduce the physical dimensions so the
# scientific integration test remains practical on every CI push. The --paper
# profile itself retains the manuscript's full 68-sample/depth settings.
sh "$experiment" \
    --paper \
    --samples 6 \
    --clusters 12 \
    --variants-per-cluster 3 \
    --fixed-depth 36 \
    --min-depth 4 \
    --max-depth 68 \
    --baseline-length 50 \
    --min-variable-length 35 \
    --max-variable-length 50 \
    --min-errors 1 \
    --max-errors 4 \
    --conserved-sequence ACGTAC \
    --project paper \
    --seed 2020 \
    --no-truth

awk -F '\t' '
    NR == 1 { next }
    {
        expected_a = index($2, "A") ? 1 : 0
        expected_c = index($2, "C") ? 1 : 0
        expected_e = index($2, "E") ? 1 : 0
        expected_l = index($2, "L") ? 1 : 0
        expected_n = index($2, "N") ? 1 : 0
        if ($6 != expected_a || $7 != expected_c || $8 != expected_e ||
            $9 != expected_l || $10 != expected_n)
            exit 1
        if (!($2 in treatment_seen)) {
            treatment_seen[$2] = 1
            treatment_count++
        }
        treatment_rows[$2]++
        if ($3 > maximum_replicate)
            maximum_replicate = $3
    }
    END {
        if (treatment_count != 32 || maximum_replicate != 10)
            exit 1
        for (treatment in treatment_rows)
            if (treatment_rows[treatment] != 60)
                exit 1
        if (!("O" in treatment_seen) || !("ACELN" in treatment_seen))
            exit 1
    }
' paper-manifest.tsv || fail "paper profile did not produce the complete 32-condition, 10-replicate design"

sample_fastas=$(count_files -name 'paper-*-r???-s???.fasta')
[ "$sample_fastas" -eq 1920 ] || fail "paper design produced $sample_fastas sample FASTAs; expected 1920"
reference_fastas=$(count_files -name 'paper-*-r???-reference.fasta')
[ "$reference_fastas" -eq 320 ] || fail "paper design produced $reference_fastas references; expected 320"

sh "$distance" \
    --manifest paper-manifest.tsv \
    --output paper-distances.tsv \
    --summary paper-distance-summary.tsv \
    --matrix-prefix paper \
    --graph paper-distance-summary.svg

result_lines=$(awk 'END { print NR + 0 }' paper-distances.tsv)
[ "$result_lines" -eq 4801 ] || fail "paper analysis has $result_lines rows; expected 4801"
summary_lines=$(awk 'END { print NR + 0 }' paper-distance-summary.tsv)
[ "$summary_lines" -eq 1601 ] || fail "paper summary has $summary_lines rows; expected 1601"
distance_matrices=$(count_files -name 'paper-*-r???-*-matrix.tsv')
[ "$distance_matrices" -eq 1600 ] || fail "direct analyzer produced $distance_matrices matrices; expected 1600"
grep -q '<svg ' paper-distance-summary.svg || fail "paper analysis did not produce an SVG graph"
grep -q 'Generated entirely by POSIX AWK' paper-distance-summary.svg || fail "paper SVG lacks AWK provenance"

awk -F '\t' '
    NR == 1 { next }
    {
        if ($4 != 6 || $5 != 15)
            exit 1
        if ($7 < 0 || $7 > 1 || $8 < 0 || $8 > 1)
            exit 1

        treatment = $2
        metric = $6
        if (metric == "jaccard") {
            if (index(treatment, "E")) {
                error_sum += $7
                error_n++
            } else {
                plain_sum += $7
                plain_n++
            }
        }
        if (metric == "ruzicka") {
            if (index(treatment, "N")) {
                depth_sum += $7
                depth_n++
            } else {
                fixed_sum += $7
                fixed_n++
            }
        }
    }
    END {
        error_mean = error_sum / error_n
        plain_mean = plain_sum / plain_n
        depth_mean = depth_sum / depth_n
        fixed_mean = fixed_sum / fixed_n
        if (!(error_mean > plain_mean))
            exit 1
        if (!(depth_mean > fixed_mean))
            exit 1
        printf "Standalone effect check: Jaccard E %.6f > %.6f; Ruzicka N %.6f > %.6f\n", \
               error_mean, plain_mean, depth_mean, fixed_mean
    }
' paper-distance-summary.tsv || fail "standalone error/depth effects were not recovered"

echo "All standalone AWK pipeline tests passed (32 conditions x 10 replicates; no external analysis software)."
