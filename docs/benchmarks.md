# Benchmarks

Measured numbers for representative workloads, in the spirit of ReadMixer's
performance tables (see `docs/lessons-from-readmixer.md`). All runs are
single-threaded POSIX AWK; to use more cores, shard factorial replicates
across processes (e.g. `xargs -P`), which composes with per-replicate
determinism.

## Machine

- Apple M3, 8 cores, 16 GB RAM, macOS (Darwin 25.6.0 arm64)
- BSD awk (20200816), system `gzip`, APFS SSD
- Wall time via `/usr/bin/time -l`; RSS is maximum resident set size

## Results

| Workload | Wall | Peak RSS | Output |
| --- | --- | --- | --- |
| Simulator quick start: 5 FASTAs, 20 clusters × 10 variants, 50–100 reads/file, 300–500 bp (`--seed 2026`) | 0.07 s | 2.8 MB | 176 KB |
| Same workload, `--format fastq --gzip` | 0.20 s | 2.8 MB | 28 KB |
| Experiment factorial: 32 treatments × 2 replicates × 4 samples, 12 clusters × 3 variants, depth 50 | 1.80 s | 2.6 MB | 59 MB (320 FASTAs + manifest + truth) |
| `Barcode_Distance` over that factorial manifest (64 datasets, 256 samples, 5 metrics + 320 matrices + SVG) | 6.36 s | 3.6 MB | 1.3 MB |

Memory stays flat: the tools hold one variant pool (experiment) or one
master list (simulator) plus one dataset's feature table (distance) at a
time, so peak RSS is measured in megabytes even for multi-gigabyte runs.
Runtime scales with total reads written (simulator/experiment) and with
distinct sequences × sample pairs (distance).

## Estimator accuracy

`--dry-run` for the factorial workload above predicted ~61.0 MB expected /
~241.8 MB worst case; the actual run wrote 59 MB. The worst case assumes
every sample hits maximum depth and length simultaneously, so treat it as
an upper bound for the disk preflight, not a forecast.

For the manuscript-scale `--paper` profile, `--dry-run` predicts ~56.5 GB
expected / ~200.4 GB worst case (plus ~8 GB of truth-table rows) — hence
the tens-of-gigabytes warning in the README. Run `--dry-run` first and use
`--no-truth` when read-level ground truth is not needed.

## Reproducing

```sh
sh Barcode_Simulator --num-fasta 5 --total-clusters 20 \
  --sequences-per-cluster 10 --min-sequences-per-file 50 \
  --max-sequences-per-file 100 --min-length 300 --max-length 500 \
  --project-name bench --seed 2026 --dry-run

sh Barcode_Experiment --factorial --project bench --samples 4 \
  --replicates 2 --clusters 12 --variants-per-cluster 3 \
  --fixed-depth 50 --seed 7 --dry-run
```
