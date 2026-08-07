#!/bin/sh

set -eu

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
post="$root_dir/Barcode_Simulator_Post"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/barcode-post.XXXXXX")
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

sh "$post" --help | grep -q 'all-distances' || fail "help output is incomplete"

cat > case_otu_table.txt <<'EOF'
# Constructed from test data
#OTU ID	S1	S2	S3
otu-1	2	1	0
otu-2	0	1	1
otu-3	0	0	1
EOF

cat > case_mash_dists.txt <<'EOF'
S1.fasta.msh S2.fasta.msh 0.666666666667 0 1/1
S1.fasta.msh S3.fasta.msh 1.0 0 1/1
S2.fasta.msh S3.fasta.msh 0.666666666667 0 1/1
EOF

cat > meta.data.file.txt <<'EOF'
File Abundance Conserved Errors Lengths Picks Type
case_otu_table.txt 0 0 0 0 0 O
EOF

sh "$post" distances --output summary.tsv case_otu_table.txt
assert_lines 2 summary.tsv
sh "$post" distances --output summary.tsv case_otu_table.txt
assert_lines 2 summary.tsv
awk -F '\t' '
    NR == 2 {
        if ($1 != "case" || $2 != "O" || $3 != 3 || $4 != 3)
            exit 1
        if (($5 - 0.777777777778)^2 > 1e-18)
            exit 1
        if (($6 - 0.777777777778)^2 > 1e-18)
            exit 1
        if (($7 - 1)^2 > 1e-18)
            exit 1
    }
' summary.tsv || fail "distance summary is incorrect"

sh "$post" all_distances --output points.tsv case_otu_table.txt
assert_lines 4 points.tsv
awk -F '\t' '
    NR == 1 { next }
    $1 != "case" || $2 != "O" { exit 1 }
    $5 < 0 || $5 > 1 || $6 < 0 || $6 > 1 { exit 1 }
' points.tsv || fail "all-distances output is incorrect"

sh "$post" mantel --output mantel.tsv case_otu_table.txt case_mash_dists.txt
assert_lines 4 mantel.tsv
awk -F '\t' 'NR > 1 && (($7 - 1)^2 > 1e-18) { exit 1 }' mantel.tsv || \
    fail "Mantel/Pearson correlations are incorrect"

sh "$post" distances --jaccard binary --output binary.tsv case_otu_table.txt
awk -F '\t' '
    NR == 2 && (($5 - 0.722222222222)^2 > 1e-18) { exit 1 }
' binary.tsv || fail "binary Jaccard calculation is incorrect"

mv summary.tsv explicit-summary.tsv
sh "$post" distances --output summary.tsv
cmp summary.tsv explicit-summary.tsv >/dev/null || fail "automatic OTU discovery changed results"

cat > similarity_otu_table.txt <<'EOF'
OTU.id	S1	S2	S3
otu-1	2	1	0
otu-2	0	1	1
otu-3	0	0	1
EOF
cat > similarity_mash_dists.txt <<'EOF'
S1 S2 0.333333333333
S1 S3 0
S2 S3 0.333333333333
EOF
sh "$post" distances --mash-values similarity --meta none \
    --output similarity.tsv similarity_otu_table.txt
awk -F '\t' 'NR == 2 && (($7 - 1)^2 > 1e-18) { exit 1 }' similarity.tsv || \
    fail "similarity-valued Mash compatibility is incorrect"

if sh "$post" distances --output missing.tsv missing_otu_table.txt >/dev/null 2>&1; then
    fail "missing input was accepted"
fi

echo "All Barcode Simulator post-processing tests passed."
