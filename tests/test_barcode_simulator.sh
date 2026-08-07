#!/bin/sh

set -eu

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
simulator="$root_dir/Barcode_Simulator"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/barcode-simulator.XXXXXX")
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

sh "$simulator" --help | grep -q "total-clusters" || fail "help output is incomplete"

sh "$simulator" \
    --num-fasta 2 \
    --total-clusters 3 \
    --sequences-per-cluster 4 \
    --min-sequences-per-file 5 \
    --max-sequences-per-file 5 \
    --min-length 12 \
    --max-length 12 \
    --min-differences 2 \
    --max-differences 2 \
    --project-name simulated \
    --gen-ref references.fasta \
    --save master.txt \
    --seed 42

assert_lines 10 simulated-1.fasta
assert_lines 10 simulated-2.fasta
assert_lines 6 references.fasta
assert_lines 12 master.txt

awk '
    /^>/ { next }
    { reference[++reference_count] = $0 }
    END {
        while ((getline sequence < "master.txt") > 0) {
            cluster = int((++record - 1) / 4) + 1
            differences = 0
            for (i = 1; i <= length(sequence); i++)
                if (substr(sequence, i, 1) != substr(reference[cluster], i, 1))
                    differences++
            if (differences != 2)
                exit 1
        }
    }
' references.fasta || fail "generated variants do not have exactly two substitutions"

cat > reused.fasta <<'EOF'
>copy-one
ACGTACGT
>copy-two
TTTT
CCCC
EOF

sh "$simulator" \
    --num-fasta 1 \
    --reuse reused.fasta \
    --min-sequences-per-file 2 \
    --max-sequences-per-file 2 \
    --project-name reused \
    --add=-copy \
    --seed 9

assert_lines 4 reused-1.fasta
grep -q '^ACGTACGT$' reused-1.fasta || fail "first reused sequence is absent"
grep -q '^TTTTCCCC$' reused-1.fasta || fail "wrapped reused sequence was not joined"
awk 'NR % 2 == 1 && $0 !~ /^>FASTA-1_[12]-copy$/ { exit 1 }' reused-1.fasta || fail "headers are malformed"

if sh "$simulator" --min-length 5 --max-length 4 >/dev/null 2>&1; then
    fail "invalid length range was accepted"
fi

sh "$simulator" \
    --total-clusters 2 \
    --sequences-per-cluster 2 \
    --only \
    --save only-master.txt \
    --seed 7

assert_lines 4 only-master.txt
[ ! -e default-1.fasta ] || fail "--only unexpectedly generated a FASTA file"

echo "All Barcode Simulator tests passed."
