#!/bin/sh

set -eu

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
distance="$root_dir/Barcode_Distance"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/barcode-distance.XXXXXX")
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

assert_lines() {
    expected=$1
    file=$2
    actual=$(awk 'END { print NR + 0 }' "$file")
    [ "$actual" -eq "$expected" ] || fail "$file has $actual lines; expected $expected"
}

cd "$test_dir"
mkdir inputs

cat > inputs/s1.fasta <<'EOF'
>s1-a1
AA
A
>s1-a2
AAA
>s1-c1
CCC
EOF

cat > inputs/s2.fasta <<'EOF'
>s2-a1
AAA
>s2-g1
GGG
EOF

cat > inputs/s3.fasta <<'EOF'
>s3-c1
CCC
>s3-g1
GGG
>s3-t1
TTT
EOF

cat > inputs/design-manifest.tsv <<'EOF'
file	treatment	replicate	sample
s1.fasta	O	1	S1
s2.fasta	O	1	S2
s3.fasta	O	1	S3
EOF

sh "$distance" --help | grep -q 'Hellinger' || fail "help output does not list the direct metrics"

sh "$distance" \
    --manifest inputs/design-manifest.tsv \
    --output pairs.tsv \
    --summary summary.tsv \
    --matrix-prefix matrix \
    --graph distances.svg

assert_lines 4 pairs.tsv
assert_lines 6 summary.tsv

awk -F '\t' '
    NR == 2 {
        if ($1 != "O-r001" || $2 != "O" || $3 != 1 || $4 != "S1" || $5 != "S2")
            exit 1
        if ($6 != 3 || $7 != 2 || $8 != 1 || $9 != 3)
            exit 1
        if (($10 - 0.666666666667)^2 > 1e-18)
            exit 1
        if (($11 - 0.75)^2 > 1e-18)
            exit 1
        if (($12 - 0.6)^2 > 1e-18)
            exit 1
        if (($13 - 0.367544467966)^2 > 1e-18)
            exit 1
        if (($14 - 0.650115167344)^2 > 1e-18)
            exit 1
    }
' pairs.tsv || fail "direct distance calculations are incorrect"

for metric in jaccard ruzicka bray-curtis cosine hellinger; do
    matrix_file="matrix-O-r001-$metric-matrix.tsv"
    [ -f "$matrix_file" ] || fail "missing $metric matrix"
    assert_lines 4 "$matrix_file"
    awk -F '\t' '
        NR == 1 { next }
        {
            if ($(NR) != 0)
                exit 1
            for (column = 2; column <= NF; column++)
                if ($column < 0 || $column > 1)
                    exit 1
        }
        NR == 2 { forward = $3 }
        NR == 3 && ($2 - forward)^2 > 1e-18 { exit 1 }
    ' "$matrix_file" || fail "$metric matrix is not symmetric with a zero diagonal"
done

grep -q '<svg ' distances.svg || fail "AWK graph is not SVG"
grep -q 'Generated entirely by POSIX AWK' distances.svg || fail "SVG provenance is missing"
grep -q 'Bray-Curtis' distances.svg || fail "SVG does not contain all metrics"

# Outputs must be replaced, not appended, on a deterministic rerun.
sh "$distance" \
    --manifest inputs/design-manifest.tsv \
    --output pairs.tsv \
    --summary summary.tsv \
    --matrix-prefix matrix \
    --graph distances.svg >/dev/null
assert_lines 4 pairs.tsv
assert_lines 6 summary.tsv

# The same analyzer also accepts an ad hoc set of FASTAs without a manifest.
sh "$distance" \
    --dataset direct \
    --no-matrices \
    --no-graph \
    --output direct.tsv \
    --summary direct-summary.tsv \
    inputs/s1.fasta inputs/s2.fasta
assert_lines 2 direct.tsv
assert_lines 6 direct-summary.tsv

# A rejected output/input collision must not truncate the input FASTA.
cp inputs/s1.fasta s1-before.fasta
if sh "$distance" \
    --no-matrices --no-graph \
    --output ./inputs/s1.fasta \
    --summary collision-summary.tsv \
    inputs/s1.fasta inputs/s2.fasta >/dev/null 2>&1; then
    fail "an output path was allowed to overwrite an input FASTA"
fi
cmp inputs/s1.fasta s1-before.fasta >/dev/null || fail "rejected output collision modified its input"

if sh "$distance" --no-matrices --no-graph inputs/s1.fasta >/dev/null 2>&1; then
    fail "a one-sample dataset was accepted"
fi

echo "All direct Barcode Distance tests passed."
