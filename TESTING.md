# Testing Barcode_Simulator

All tests are POSIX shell scripts run from the repository root. No
packages need to be installed; the only requirements are a
POSIX-compatible `awk`, a POSIX shell, and (for the gzip round-trip
tests) `gzip`/`gunzip`.

```sh
sh tests/test_barcode_simulator.sh
sh tests/test_barcode_experiment.sh
sh tests/test_barcode_distance.sh
sh tests/test_barcode_post.sh
sh tests/test_paper_pipeline.sh
sh tests/test_upgrades.sh
```

Each script creates a private temporary directory with `mktemp -d`,
installs a `trap 'rm -rf' EXIT HUP INT TERM` cleanup handler, runs the
tools with fixed `--seed` values, and asserts on the outputs. A non-zero
exit means failure; passing scripts print `All ... tests passed.`

## Suite overview

| Script | What it covers |
| --- | --- |
| `test_barcode_simulator.sh` | FASTA structure, exact mutation counts, seeded/deterministic reruns, wrapped-FASTA reuse, historical header behavior, option validation, `--only` master mode |
| `test_barcode_experiment.sh` | Each A/C/E/L/N effect in isolation, the complete 32-condition factorial, manifest and truth-table contents, deterministic reruns |
| `test_barcode_distance.sh` | Quantitative and binary Jaccard, all five direct distance formulas against hand-computed fixtures, symmetric matrices, SVG graph generation, safe output handling (no input clobbering) |
| `test_barcode_post.sh` | Legacy Mash parsing, quantitative/binary Jaccard, Mantel/Pearson correlations against known values |
| `test_paper_pipeline.sh` | Regression scaffold for the paper design: all 32 A/C/E/L/N combinations × 10 replicates (320 datasets, 1,920 FASTAs) straight into `Barcode_Distance`; checks 1,600 matrices, 4,800 sample pairs, all five metrics, and the SVG. Sample counts, depth, and length are scaled down for CI |
| `test_upgrades.sh` | `--loglevel` filtering on all four tools; `--dry-run` estimates without writing files; `.gz` reuse round-trip; `--format fastq` record shape and sequence/quality agreement; `--gzip` output round-trip; `--max-bytes-per-fasta` part splitting; `--run-metadata` JSON validity; `Barcode_Distance` over FASTQ and `.gz` inputs; invalid `--loglevel`/`--format` rejection |

## Conventions for new tests

- One `fail()` helper printing `FAIL: ...` to stderr and exiting 1.
- Fixed `--seed` values everywhere; never assert on wall-clock time.
- New tool options get at least one positive test and one rejection test.
- Keep the paper scaffold scaled down (it already exercises 320 datasets);
  full `--paper` depth is for manual runs, with numbers recorded in
  `docs/benchmarks.md`.
