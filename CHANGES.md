# Changes

This changelog summarizes user-visible changes. Dates are release-merge
dates from the git history; see `git log` for the full record.

## Unreleased

- Added `--loglevel error|warning|info|debug` to all four tools
  (default `info` preserves the previous stderr report).
- Added `--dry-run` size estimates plus an automatic disk-space preflight
  (with `--no-disk-check` override) to `Barcode_Simulator` and
  `Barcode_Experiment`.
- Added transparent `.gz` input reading (`--reuse` and `Barcode_Distance`,
  via `gzip -dc`; plain text remains the default everywhere).
- Added `--run-metadata FILE` JSON sidecars recording effective parameters,
  seed, platform, and outputs written.
- Added `--format fasta|fastq` (uniform synthetic Q40 qualities) and
  `--gzip` sequence outputs to `Barcode_Simulator` and
  `Barcode_Experiment`; `Barcode_Distance` now reads FASTQ inputs too.
- Added `--max-bytes-per-fasta N` output splitting to `Barcode_Simulator`.
- Added scholarly-hygiene files: `LICENSE` (Apache-2.0), `CITATION.cff`,
  `TESTING.md`, `docs/definitions.md`, `ARCHIVAL_STRATEGY.md`,
  `Barcode_Simulator_SBOM.json`, `docs/benchmarks.md`.
- New `tests/test_upgrades.sh` covering the features above; wired into
  `.github/workflows/test.yml`.

## 2026-08-07 — Direct AWK distance analysis

- New `Barcode_Distance` tool: five distances (Jaccard, Ruzicka,
  Bray–Curtis, cosine, Hellinger) computed directly from simulated FASTAs,
  with square matrices and an AWK-authored SVG heatmap. No QIIME, DADA2,
  Mash, R, or plotting package required.
- New `tests/test_barcode_distance.sh` and `tests/test_paper_pipeline.sh`
  (32-condition, 10-replicate regression scaffold analyzed end to end).

## 2026-08-06 — AWK post-processing replaces R

- New `Barcode_Simulator_Post` AWK implementation replacing
  `Barcode_Simulator_Post.R` / `Barcode_Simulator_Post_Single.R`
  (`distances`, `all-distances`, `mantel` modes emitting tidy TSV).
- Quantitative Jaccard matching `vegan::vegdist`, `--jaccard binary` option,
  corrected Mash-distance handling (`--mash-values`).

## 2026-08-05 — Experimental-design layer restored

- New `Barcode_Experiment`: the five A/C/E/L/N effects from
  Molik, Pfrender, and Emrich (2020), single-treatment and full 2^5
  factorial generation, `--paper` manuscript-scale profile, per-treatment
  reference FASTAs, `*-manifest.tsv` design tables, and `*-truth.tsv`
  read-level ground truth with `--no-truth` escape hatch.

## 2026-08-04 — POSIX AWK rewrite of the simulator

- `Barcode_Simulator` reimplemented in POSIX AWK (historical CLI retained,
  clearer `--total-clusters` / `--sequences-per-cluster` aliases added).
- Internal Park–Miller RNG replaces `rand()`/`shuf`/`/dev/urandom` so
  `--seed` reproduces runs across AWK implementations.
- `--reuse` accepts wrapped FASTA or one-sequence-per-line input;
  `--gen-ref`, `--save`, and `--only` master/reference modes added.
- Regression suite in `tests/test_barcode_simulator.sh`; CI on
  Ubuntu and macOS.

## 2018 — Original Bash simulator

- Historical Bash implementation with GNU `getopt`, `shuf`, `sed`, and
  `/dev/urandom` dependencies, plus QIIME/Mash/cluster-submission scripts
  (retained under `scripts/` as research artifacts) supporting the 2020
  paper analysis.
