#!/bin/sh

set -eu

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
simulator="$root_dir/Barcode_Simulator"
experiment="$root_dir/Barcode_Experiment"
distance="$root_dir/Barcode_Distance"
post="$root_dir/Barcode_Simulator_Post"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/barcode-upgrades.XXXXXX")
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

check_json() {
    file=$1
    tool=$2
    grep -q '"tool": "'"$tool"'"' "$file" || fail "$file is missing tool marker"
    head -c 1 "$file" | grep -q '{' || fail "$file does not start with {"
    if command -v python3 >/dev/null 2>&1; then
        python3 -c 'import json,sys; json.load(open(sys.argv[1]))' "$file" \
            || fail "$file is not valid JSON"
    fi
}

cd "$test_dir"

# --loglevel is accepted everywhere and filters output.
sh "$simulator" --help | grep -q -- '--loglevel' || fail "simulator help is missing --loglevel"
sh "$experiment" --help | grep -q -- '--loglevel' || fail "experiment help is missing --loglevel"
sh "$distance" --help | grep -q -- '--loglevel' || fail "distance help is missing --loglevel"
sh "$post" --help | grep -q -- '--loglevel' || fail "post help is missing --loglevel"

sh "$simulator" --num-fasta 1 --total-clusters 2 --sequences-per-cluster 2 \
    --project-name quiet --seed 1 --loglevel error 2>quiet-stderr.txt
[ -s quiet-stderr.txt ] && fail "--loglevel error should suppress the report"
[ -f quiet-1.fasta ] || fail "quiet run did not write its FASTA"

sh "$simulator" --num-fasta 1 --total-clusters 2 --sequences-per-cluster 2 \
    --project-name loud --seed 1 --loglevel debug 2>loud-stderr.txt
grep -q "Barcode Simulator complete" loud-stderr.txt || fail "debug run lost the report"
grep -q "writing loud-1.fasta" loud-stderr.txt || fail "debug run is missing per-file progress"

sh "$experiment" --effects errors --project level --samples 1 --replicates 1 \
    --clusters 3 --variants-per-cluster 1 --fixed-depth 2 --seed 1 \
    --loglevel error 2>level-stderr.txt
[ -s level-stderr.txt ] && fail "experiment --loglevel error should suppress the report"

if sh "$simulator" --project-name bad --loglevel verbose >/dev/null 2>&1; then
    fail "invalid --loglevel was accepted"
fi
if sh "$simulator" --project-name bad --format sam >/dev/null 2>&1; then
    fail "invalid --format was accepted"
fi
if sh "$experiment" --project bad --loglevel verbose >/dev/null 2>&1; then
    fail "experiment accepted an invalid --loglevel"
fi
if sh "$distance" --dataset bad --loglevel verbose >/dev/null 2>&1; then
    fail "distance accepted an invalid --loglevel"
fi

# --dry-run estimates without writing anything.
sh "$simulator" --num-fasta 2 --total-clusters 3 --sequences-per-cluster 4 \
    --min-sequences-per-file 5 --max-sequences-per-file 6 \
    --project-name dry --seed 42 --dry-run 2>dry-stderr.txt
grep -q "dry run (no files written)" dry-stderr.txt || fail "dry run banner missing"
grep -q "master sequences: ~12" dry-stderr.txt || fail "dry run master estimate wrong"
[ -f dry-1.fasta ] && fail "--dry-run wrote files"
[ -f dry-2.fasta ] && fail "--dry-run wrote files"

sh "$experiment" --factorial --project dryexp --seed 1 --dry-run 2>dryexp-stderr.txt
grep -q "treatments: 32" dryexp-stderr.txt || fail "experiment dry run missed treatments"
[ -f dryexp-manifest.tsv ] && fail "experiment --dry-run wrote files"

# --no-disk-check is accepted and runs.
sh "$simulator" --num-fasta 1 --total-clusters 1 --sequences-per-cluster 1 \
    --project-name nochek --seed 1 --no-disk-check >/dev/null 2>&1
[ -f nochek-1.fasta ] || fail "--no-disk-check run did not write its FASTA"

# Transparent .gz reuse: compressed and plain reuse agree exactly.
sh "$simulator" --num-fasta 1 --total-clusters 2 --sequences-per-cluster 2 \
    --min-sequences-per-file 3 --max-sequences-per-file 3 \
    --min-length 20 --max-length 20 --project-name base --seed 7 >/dev/null 2>&1
gzip -c base-1.fasta > reuse.fasta.gz
sh "$simulator" --num-fasta 1 --min-sequences-per-file 2 --max-sequences-per-file 2 \
    --project-name plain --seed 7 --reuse base-1.fasta >/dev/null 2>&1
sh "$simulator" --num-fasta 1 --min-sequences-per-file 2 --max-sequences-per-file 2 \
    --project-name zipped --seed 7 --reuse reuse.fasta.gz >/dev/null 2>&1
cmp plain-1.fasta zipped-1.fasta || fail ".gz reuse differs from plain reuse"

# --format fastq: same seed gives the same headers and sequences as FASTA.
sh "$simulator" --num-fasta 1 --total-clusters 2 --sequences-per-cluster 2 \
    --min-sequences-per-file 3 --max-sequences-per-file 3 \
    --min-length 20 --max-length 20 --project-name fasta --seed 7 >/dev/null 2>&1
sh "$simulator" --num-fasta 1 --total-clusters 2 --sequences-per-cluster 2 \
    --min-sequences-per-file 3 --max-sequences-per-file 3 \
    --min-length 20 --max-length 20 --project-name fastq --seed 7 \
    --format fastq >/dev/null 2>&1
awk 'NR % 2 == 1 { sub(/^>/, "@"); headers[++h] = $0; next } \
     { seqs[h] = $0 } \
     END { print h }' fasta-1.fasta > fasta-headers.txt
awk 'NR % 4 == 1 { headers[++h] = $0; next } \
     NR % 4 == 2 { seqs[h] = $0; next } \
     NR % 4 == 3 { if ($0 != "+") exit 1; next } \
     { if (length($0) != length(seqs[h])) exit 1; \
       if ($0 !~ /^I+$/) exit 1 } \
     END { print h }' fastq-1.fastq > fastq-check.txt \
    || fail "FASTQ records are malformed"
[ "$(cat fastq-check.txt)" = "3" ] || fail "FASTQ record count is wrong"
awk 'NR % 4 == 1 { print substr($0, 2) }' fastq-1.fastq > fastq-headers.txt
awk 'NR % 4 == 2 { print }' fastq-1.fastq > fastq-seqs.txt
awk 'NR % 2 == 1 { print substr($0, 2) }' fasta-1.fasta > fasta-headers-stripped.txt
awk 'NR % 2 == 0 { print }' fasta-1.fasta > fasta-seqs.txt
cmp fasta-headers-stripped.txt fastq-headers.txt || fail "FASTQ headers differ from FASTA"
cmp fasta-seqs.txt fastq-seqs.txt || fail "FASTQ sequences differ from FASTA"

# --gzip round-trip: decompressed output matches the plain run.
sh "$simulator" --num-fasta 1 --total-clusters 2 --sequences-per-cluster 2 \
    --min-sequences-per-file 3 --max-sequences-per-file 3 \
    --min-length 20 --max-length 20 --project-name gz --seed 7 \
    --format fastq --gzip >/dev/null 2>&1
[ -f gz-1.fastq.gz ] || fail "--gzip output is missing"
gunzip -c gz-1.fastq.gz > gz-restored.fastq
cmp fastq-1.fastq gz-restored.fastq || fail "--gzip round-trip differs"

# Splitting: parts cover exactly the single-file run.
sh "$simulator" --num-fasta 1 --total-clusters 2 --sequences-per-cluster 2 \
    --min-sequences-per-file 4 --max-sequences-per-file 4 \
    --min-length 20 --max-length 20 --project-name whole --seed 7 >/dev/null 2>&1
sh "$simulator" --num-fasta 1 --total-clusters 2 --sequences-per-cluster 2 \
    --min-sequences-per-file 4 --max-sequences-per-file 4 \
    --min-length 20 --max-length 20 --project-name split --seed 7 \
    --max-bytes-per-fasta 100 >/dev/null 2>&1
[ -f split-1.part000.fasta ] || fail "split part 000 is missing"
[ -f split-1.part001.fasta ] || fail "split part 001 is missing"
[ -f split-1.fasta ] && fail "unsplit name should not exist when splitting"
cat split-1.part*.fasta > split-joined.fasta
cmp whole-1.fasta split-joined.fasta || fail "split parts do not reassemble"
for part in split-1.part*.fasta; do
    size=$(wc -c < "$part")
    [ "$size" -le 300 ] || fail "$part exceeds the split budget absurdly"
done

# --run-metadata sidecars are valid JSON with the tool marker.
sh "$simulator" --num-fasta 1 --total-clusters 1 --sequences-per-cluster 1 \
    --project-name meta --seed 1 --run-metadata meta-sim.json >/dev/null 2>&1
check_json meta-sim.json Barcode_Simulator

sh "$experiment" --effects errors --project metaexp --samples 1 --replicates 1 \
    --clusters 3 --variants-per-cluster 1 --fixed-depth 2 --seed 1 \
    --run-metadata meta-exp.json >/dev/null 2>&1
check_json meta-exp.json Barcode_Experiment

# Experiment FASTQ+gzip outputs land in the manifest under their real names.
sh "$experiment" --effects errors --project expgz --samples 1 --replicates 1 \
    --clusters 3 --variants-per-cluster 1 --fixed-depth 2 --seed 1 \
    --format fastq --gzip >/dev/null 2>&1
[ -f expgz-E-r001-s001.fastq.gz ] || fail "experiment --gzip sample is missing"
grep -q "expgz-E-r001-s001.fastq.gz" expgz-manifest.tsv \
    || fail "manifest does not list the .fastq.gz sample"
gunzip -c expgz-E-r001-s001.fastq.gz | awk 'NR % 4 == 3 { exit ($0 != "+") }' \
    || fail "experiment FASTQ sample is malformed"

# Distance reads FASTQ and .gz inputs; hand-computed fixture:
# f1 has {ACGTACGT x2, TTTT}, f3 has {ACGTACGT, TTTT}: shared 2, union 2.
printf '>s1\nACGTACGT\n>s1\nACGTACGT\n>s1\nTTTT\n' > d1.fasta
printf '@r1\nACGTACGT\n+\nIIIIIIII\n@r1\nTTTT\n+\nIIII\n' > d3.fastq
gzip -c d1.fasta > d1.fasta.gz
sh "$distance" --dataset gz d1.fasta.gz d3.fastq --no-graph --no-matrices \
    --run-metadata meta-dist.json >/dev/null 2>&1
grep -q "gz[[:space:]]*gz[[:space:]]*1[[:space:]]*d1[[:space:]]*d3[[:space:]]*3[[:space:]]*2[[:space:]]*2[[:space:]]*2[[:space:]]*0" \
    gz-distances.tsv || fail "mixed .gz/FASTQ distance row is wrong"
grep -q "sample_a.*sample_b" gz-distances.tsv || fail "distance header missing"
check_json meta-dist.json Barcode_Distance

# Malformed FASTQ is rejected, not silently misread.
printf '@r1\nACGT\n+\nIII\n' > badqual.fastq
if sh "$distance" --dataset bad badqual.fastq d3.fastq >/dev/null 2>&1; then
    fail "quality/sequence length mismatch was accepted"
fi
printf '@r1\nACGT\n' > truncated.fastq
if sh "$distance" --dataset bad truncated.fastq d3.fastq >/dev/null 2>&1; then
    fail "truncated FASTQ was accepted"
fi

# Post honors --loglevel and writes metadata.
printf '# Constructed from biom file\n#OTU ID\ts1\ts2\notu1\t2\t0\notu2\t1\t3\n' > tiny_otu_table.txt
printf 's1\ts2\t0.5\n' > tiny_mash_dists.txt
sh "$post" mantel --output post-mantel.tsv tiny_otu_table.txt tiny_mash_dists.txt \
    --run-metadata meta-post.json --loglevel error 2>post-stderr.txt
[ -s post-stderr.txt ] && fail "post --loglevel error should suppress the report"
check_json meta-post.json Barcode_Simulator_Post

echo "All upgrade tests passed."
