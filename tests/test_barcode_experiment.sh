#!/bin/sh

set -eu

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
experiment="$root_dir/Barcode_Experiment"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/barcode-experiment.XXXXXX")
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

assert_lines() {
    expected=$1
    file=$2
    actual=$(awk 'END { print NR }' "$file")
    [ "$actual" -eq "$expected" ] || fail "$file has $actual lines; expected $expected"
}

cd "$test_dir"

sh "$experiment" --help | grep -q -- '--factorial' || fail "help output is incomplete"

sh "$experiment" \
    --effects abundance,conserved,errors,lengths,depth \
    --samples 3 \
    --clusters 6 \
    --variants-per-cluster 2 \
    --fixed-depth 5 \
    --min-depth 5 \
    --max-depth 5 \
    --baseline-length 12 \
    --min-variable-length 8 \
    --max-variable-length 10 \
    --min-errors 1 \
    --max-errors 2 \
    --conserved-sequence ACGT \
    --project effects \
    --seed 42

for sample in 001 002 003; do
    assert_lines 10 "effects-ACELN-r001-s${sample}.fasta"
done
assert_lines 12 effects-ACELN-r001-reference.fasta
assert_lines 4 effects-manifest.tsv
assert_lines 16 effects-truth.tsv

awk -F '\t' '
    NR == 1 { next }
    $2 != "ACELN" || $6 != 1 || $7 != 1 || $8 != 1 || $9 != 1 || $10 != 1 { exit 1 }
    $5 != 5 { exit 1 }
' effects-manifest.tsv || fail "effect manifest flags are incorrect"

awk -F '\t' '
    NR == 1 { next }
    $9 < 12 || $9 > 14 { exit 1 }
    $10 < 1 || $10 > 2 { exit 1 }
' effects-truth.tsv || fail "effect truth values are outside configured ranges"

awk 'NR % 2 == 0 && substr($0, 1, 4) != "ACGT" { exit 1 }' \
    effects-ACELN-r001-s001.fasta || fail "conserved prefix is absent"

sh "$experiment" \
    --abundance \
    --samples 1 \
    --clusters 6 \
    --variants-per-cluster 1 \
    --fixed-depth 600 \
    --baseline-length 12 \
    --project abundance \
    --seed 17

awk -F '\t' '
    NR == 1 { next }
    { tier_reads[$8]++; seen[$8 SUBSEP $6] = 1 }
    END {
        for (key in seen) {
            split(key, parts, SUBSEP)
            tier_clusters[parts[1]]++
        }
        high = tier_reads["high"] / tier_clusters["high"]
        middling = tier_reads["middling"] / tier_clusters["middling"]
        low = tier_reads["low"] / tier_clusters["low"]
        if (!(high > middling && middling > low))
            exit 1
    }
' abundance-truth.tsv || fail "abundance tiers do not produce descending per-cluster coverage"

mkdir factorial-a factorial-b
for directory in factorial-a factorial-b; do
    (
        cd "$directory"
        sh "$experiment" \
            --factorial \
            --samples 1 \
            --clusters 6 \
            --variants-per-cluster 1 \
            --fixed-depth 2 \
            --min-depth 2 \
            --max-depth 2 \
            --baseline-length 12 \
            --min-variable-length 10 \
            --max-variable-length 12 \
            --min-errors 1 \
            --max-errors 1 \
            --conserved-sequence AC \
            --project factorial \
            --seed 99
    )
done

fasta_count=$(find factorial-a -type f -name 'factorial-*-s001.fasta' | awk 'END { print NR }')
[ "$fasta_count" -eq 32 ] || fail "factorial produced $fasta_count sample FASTAs; expected 32"
assert_lines 33 factorial-a/factorial-manifest.tsv
assert_lines 65 factorial-a/factorial-truth.tsv

cmp factorial-a/factorial-manifest.tsv factorial-b/factorial-manifest.tsv >/dev/null || \
    fail "same seed did not reproduce the manifest"
cmp factorial-a/factorial-ACELN-r001-s001.fasta \
    factorial-b/factorial-ACELN-r001-s001.fasta >/dev/null || \
    fail "same seed did not reproduce sequence output"

if sh "$experiment" --factorial --errors >/dev/null 2>&1; then
    fail "factorial accepted a conflicting individual effect"
fi

echo "All Barcode Experiment tests passed."
