#!/bin/sh

set -eu

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
experiment="$root_dir/Barcode_Experiment"
post="$root_dir/Barcode_Simulator_Post"
matrix_builder="$root_dir/tests/build_paper_matrices.awk"
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
        replicate_seen[$2 SUBSEP $3] = 1
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

set -- paper-*-r???-s???.fasta
awk -v project=paper -v k=7 -f "$matrix_builder" "$@"

otu_tables=$(count_files -name 'paper-*-r???_otu_table.txt')
[ "$otu_tables" -eq 320 ] || fail "matrix builder produced $otu_tables OTU tables; expected 320"
asv_tables=$(count_files -name 'paper-*-r???_asv_otu_table.txt')
[ "$asv_tables" -eq 320 ] || fail "matrix builder produced $asv_tables ASV tables; expected 320"
kmer_otu_matrices=$(count_files -name 'paper-*-r???_mash_dists.txt')
kmer_asv_matrices=$(count_files -name 'paper-*-r???_asv_mash_dists.txt')
kmer_matrices=$((kmer_otu_matrices + kmer_asv_matrices))
[ "$kmer_matrices" -eq 640 ] || fail "matrix builder produced $kmer_matrices k-mer matrices; expected 640"

sh "$post" distances \
    --meta paper-analysis-meta.tsv \
    --output paper-distances.tsv

result_lines=$(awk 'END { print NR + 0 }' paper-distances.tsv)
[ "$result_lines" -eq 641 ] || fail "paper analysis has $result_lines rows; expected 641"

awk -F '\t' '
    NR == 1 { next }
    {
        if ($3 != 6 || $4 != 15)
            exit 1
        if ($5 == "NA" || $6 == "NA" || $5 < 0 || $5 > 1 || $6 < 0 || $6 > 1)
            exit 1

        treatment = $2
        method = treatment
        sub(/^.*-/, "", method)
        sub(/-(OTU|ASV)$/, "", treatment)

        if (method == "ASV") {
            if (index(treatment, "E")) {
                asv_error_sum += $5
                asv_error_n++
            } else {
                asv_plain_sum += $5
                asv_plain_n++
            }
        }
        if (method == "OTU") {
            if (index(treatment, "N")) {
                otu_depth_sum += $5
                otu_depth_n++
            } else {
                otu_fixed_sum += $5
                otu_fixed_n++
            }
        }
        if ($7 != "NA")
            numeric_mantel++
    }
    END {
        asv_error_mean = asv_error_sum / asv_error_n
        asv_plain_mean = asv_plain_sum / asv_plain_n
        otu_depth_mean = otu_depth_sum / otu_depth_n
        otu_fixed_mean = otu_fixed_sum / otu_fixed_n
        if (!(asv_error_mean > asv_plain_mean))
            exit 1
        if (!(otu_depth_mean > otu_fixed_mean))
            exit 1
        if (numeric_mantel < 320)
            exit 1
        printf "Paper-style effect check: ASV E %.6f > %.6f; OTU N %.6f > %.6f\n", \
               asv_error_mean, asv_plain_mean, otu_depth_mean, otu_fixed_mean
    }
' paper-distances.tsv || fail "paper-style error/depth effects were not recovered"

echo "All paper-style AWK pipeline tests passed (32 conditions x 10 replicates; no R)."
