# Lessons from ReadMixer for Barcode_Simulator

Resolves the research question in
[issue #4](https://github.com/molikd/Barcode_Simulator/issues/4):
"What can we learn from <https://github.com/Geneinfosec-Inc/ReadMixer>?
Document potential Barcode_Simulator upgrades."

Analysis date: 2026-09-09. ReadMixer reviewed at `main` (v1.0.1,
Apache-2.0, Python 3.11, DOI `10.5281/zenodo.18188647`).

## 1. What ReadMixer is

ReadMixer (`Geneinfosec-Inc/ReadMixer`) scrambles/obfuscates FASTQ files so
pooled sequencing reads cannot be attributed to their source sample without
a separate key file. Its architecture is:

- `readmixer/readmixer.py` — small `argparse` CLI (`-s/-p/-o/-z/-m/-w/-c/-l`).
- `readmixer/process/scramble.py` — chunked, multiprocessed pipeline:
  find FASTQs, count reads by streaming line counts, validate `% 4 == 0`,
  shuffle each chunk with a crypto-secure RNG, write anonymized
  `@ReadMixerSecured <tag>` headers, split outputs at ~5 GB per part file,
  then merge per-input key CSVs
  (`file_name,read_description,secure_index,original_index,tag,output_file_index`).
- `readmixer/process/util.py` — `SystemCapabilities` + `ResourceManager`:
  `psutil`-based RAM/CPU/disk detection, optimal chunk/worker sizing,
  I/O buffer sizing, disk-space preflight with safety margin, memory
  monitoring with warnings above 85–90%.
- `readmixer/tests/` — `pytest` suite (19 tests) over a 10-read
  `testing.fastq.gz` fixture: FASTQ parsing, crypto RNG output shape,
  resource-manager sanity, one end-to-end scramble.
- Docs/ops layer: `README.md`, `INSTALLATION.md` (5 install paths),
  `QUICKSTART_OLD_GCC.md`, `TESTING.md`, `CHANGES.md`,
  `ARCHIVAL_STRATEGY.md` (GitHub + Zenodo DOI + Software Heritage),
  `CITATION.cff`, CycloneDX `ReadMixer_SBOM.json`, `Dockerfile` +
  `Dockerfile.conda`, `environment.yml`, `setup.sh` installer helper,
  `generate_performance_chart.py` with measured throughput
  (1.15 GB/min on a 32-core/96 GB workstation, plus scaled estimates).

## 2. Head-to-head comparison

| Concern | Barcode_Simulator (today) | ReadMixer | Takeaway |
| --- | --- | --- | --- |
| Language / deps | POSIX AWK + shell, zero packages, Linux/macOS CI | Python 3.11 + `pycryptodome`, `numpy`, `psutil`, `matplotlib` | Keep zero-dep as a feature; do not import ReadMixer's stack. |
| RNG model | Deterministic Park–Miller + `--seed`, reproducible across AWKs | Crypto RNG (`secrets`, `Crypto.Random`) for obfuscation | Deliberately different goals. Do **not** adopt crypto RNG for simulation; document the distinction. |
| Scale strategy | In-memory master list; `--paper` warns it writes tens of GB | Streaming line counts, chunked external processing, multiprocess pool, temp dir + timeout cleanup, 5 GB output splitting | Adopt the *patterns* (preflight, streaming, splitting, cleanup) without the Python machinery. |
| Input validation | DNA/IUPAC alphabet check, wrapped FASTA, min/max consistency | FASTQ `% 4` check, gzip auto-detect, permission checks | Adopt FASTQ-shape validation if FASTQ output is ever added; adopt transparent `.gz` handling. |
| Provenance | `*-manifest.tsv` (design matrix + depth) + `*-truth.tsv` (read → cluster/variant/tier/length/substitutions) | Per-input key CSV (`secure_index ↔ original_index` + tag + output part) | Already at parity conceptually. Add a run-metadata sidecar (see §3.6). |
| Logging | Fixed `report()` to stderr | `-l error/warning/info/debug/critical` with resource logging | Adopt leveled/quiet logging (§3.3). |
| Outputs | Plain FASTA only; `.gz` inputs explicitly rejected in `barcode_distance.awk:295-296` | `.fastq`/`.gz` in and out, `-z` opt-in compression | Adopt opt-in gzip + optional FASTQ-with-quality mode (§3.4). |
| Performance data | Qualitative "`--paper` writes tens of gigabytes" warning | Measured GB/min table by machine class + chart script | Adopt a small measured benchmark table (§3.7). |
| Install story | `chmod +x` + run; no packages | 5 install methods (conda env, manual conda, pip, 2× Docker) + `setup.sh` | Nothing to fix functionally, but document the zero-install path as the deliberate answer to ReadMixer's install matrix (§3.8). |
| Scholarly hygiene | No `LICENSE`, `CITATION.cff`, `CHANGES.md`, `TESTING.md`, SBOM, or archival doc | All of the above plus Zenodo DOI + Software Heritage | Adopt the cheap text-file parts (§3.8). This repo's biggest gap. |
| Glossary | Terms inline in README | `definitions.md` (base/read/tag/key-csv/fastq) | Adopt a short glossary (§3.8). |

## 3. Proposed upgrades (ordered by value/effort)

### 3.1 Pre-flight disk-space check (high value, small effort)

`Barcode_Experiment --paper` is the one command that can fill a disk. Today
there is no estimate before writing. ReadMixer's
`check_disk_space_requirements()` (input × 3.2 + CSV × 1.7, temp × 0.8,
plus 10% margin) is over-fitted to scrambling, but the *fail-fast preflight*
idea transfers directly.

Proposal: add a `--dry-run` / `--estimate` mode to `Barcode_Experiment` and
`Barcode_Simulator` that prints expected file/byte counts without writing,
plus a `df`-based free-space check in the shell launchers before invoking
AWK. Keep it POSIX (`df -P`, `awk` arithmetic); no `psutil` equivalent needed.
Calibrate the multiplier once from a scaled-down factorial run instead of
guessing.

### 3.2 Streaming + output splitting + temp hygiene (high value, medium effort)

Today `barcode_simulator.awk:206-227` builds the full master list in memory
and `generate_fastas()` reshuffles it per file; `Barcode_Experiment` holds
variant pools across replicates. That is correct for exact Hamming-distance
guarantees and reproducibility, but `--paper`-scale runs deserve:

- stream `--save`/`--gen-ref` writes instead of buffering (already mostly
  true; audit for accidental accumulation in `Barcode_Distance` feature
  tables, which grow with distinct sequences);
- a documented split threshold (e.g. new `--max-bytes-per-fasta` or a
  launcher-level split, mirroring ReadMixer's 5 GB part files:
  `NAME-1.part000.fasta`, …);
- `mktemp -d` + `trap 'rm -rf' EXIT HUP INT TERM` in any new multi-step
  launcher (the pattern already exists in `tests/`, but not in the tools).

Do not adopt multiprocessing: single-threaded deterministic AWK is the point.
Instead document `xargs -P` sharding of factorial replicates, which composes
with the existing per-replicate determinism.

### 3.3 Leveled logging (medium value, small effort)

ReadMixer's `-l info` resource summary (RAM, CPUs, chunk/worker choice, disk
check, phase timings) is genuinely useful. Mirror it with
`--loglevel error|warning|info|debug` (plus `--quiet` alias) on all four
tools, defaulting to the current terse stderr report so CI output is
unchanged. `info` should echo seed, effective design matrix, expected bytes
from §3.1, and elapsed seconds; `debug` can trace per-replicate progress.
This also gives a home for the memory/disk warnings from §3.1–§3.2.

### 3.4 Transparent gzip + optional FASTQ output (medium value, medium effort)

`barcode_distance.awk:295-296` currently fails on `*.gz` with "decompress it
first", and `--reuse` reads only plain text via `getline`. ReadMixer
auto-detects `.gz` both ways and keeps compression opt-in (`-z` default off)
— a good precedent. Proposal:

- accept `*.gz` on `--reuse` and in `Barcode_Distance` via `gzip -dc` pipe
  (`"gzip -dc -- \"...\" | getline"`), keeping plain-FASTA default;
- add opt-in `--format fastq` (with synthetic qualities, e.g. fixed `I`, or a
  configurable quality) and `--gzip` for outputs, so simulated barcodes can
  feed FASTQ-only downstream tools without a separate converter.

Add the corresponding `% 4`-style shape validation for any FASTQ path, plus
sequence/quality length agreement — the one ReadMixer validation this repo
lacks (because it has never emitted FASTQ).

### 3.5 RNG: explicitly *not* adopting crypto randomness (clarity, no code)

ReadMixer uses `secrets`/`pycryptodome` because unpredictability *is* the
product. Barcode_Simulator's Park–Miller stream exists so `--seed` reproduces
bit-identical runs across AWK implementations. Importing a crypto RNG would
destroy reproducibility and add a dependency for no scientific gain. No change
needed — but state this in the README so a future contributor does not
"upgrade" the RNG on ReadMixer-inspired advice. The correct shared lesson is
"choose the RNG for the threat/validity model", not "use crypto randomness".

### 3.6 Run-metadata sidecar (medium value, small effort)

ReadMixer's key CSV is the closest analog to `*-truth.tsv`/`*-manifest.tsv`,
and this repo already matches it. The remaining gap is machine-readable
run metadata: ReadMixer logs RAM/CPU/chunk/worker/disk decisions per run.
Add a `--run-metadata FILE` (or always-write `*-run.json`) recording command
line, seed, tool versions (`awk -W version` / `awk --version` capture is
implementation-dependent, so record `uname -srm` + tool git SHA), input
checksums, and effective parameters. That makes `*-truth.tsv` rows
interpretable years later and is a natural companion to §3.8's citation work.

### 3.7 Measured performance table (medium value, small effort)

ReadMixer's best documentation habit is its benchmark table: real GB/min on a
named machine plus labeled estimates for smaller classes, generated by
`generate_performance_chart.py`. Copy the habit, not the script: time the
scaled-down `test_paper_pipeline.sh` factorial and one full-depth single
treatment on a named machine, record wall seconds, peak RSS (`/usr/bin/time
-v` on Linux, `time -l` on macOS), and output bytes, and publish the table in
this doc. Keep charting dependency-free (the existing AWK-authored SVG path in
`Barcode_Distance` already proves the pattern).

### 3.8 Scholarly and ops hygiene (low effort, overdue)

This is the largest gap relative to ReadMixer and the cheapest to close:

- `LICENSE` — the repo currently ships **no license file**; pick one
  (ReadMixer uses Apache-2.0) so reuse is unambiguous.
- `CITATION.cff` — authors, DOI-or-URL, version; enables GitHub's "Cite this
  repository" button.
- `CHANGES.md` — move the release-note content currently buried in commit
  messages into a changelog.
- `TESTING.md` — the `tests/*.sh` suite is thorough but undiscoverable;
  document the fixture strategy and the scaled-down 320-dataset regression
  scaffold, mirroring ReadMixer's test table.
- `docs/definitions.md` — short glossary (cluster, variant, Hamming distance,
  treatment/replicate, manifest, truth table, feature); ReadMixer's
  `definitions.md` shows the format.
- `ARCHIVAL_STRATEGY.md` + SBOM — a paragraph stating intent (Zenodo DOI +
  Software Heritage, as ReadMixer does) and a minimal CycloneDX SBOM noting
  "no runtime dependencies beyond POSIX AWK + shell".
- Install story — no `Dockerfile`/`environment.yml` needed; instead document
  *why not*: zero-dependency POSIX is the portability strategy, with CI on
  `ubuntu-latest` + `macos-latest` as the proof. A one-line `Dockerfile`
  (`FROM alpine + awk`) is optional and probably unnecessary.

### 3.9 What not to adopt

- Python runtime, `numpy`/`psutil`/`matplotlib`/`pycryptodome` dependencies.
- Auto-tuning worker pools or buffer-size ladders — meaningless for
  single-threaded AWK; the §3.1 estimate + §3.3 logging covers the need.
- Conda/pip/Docker install matrix — the fix for ReadMixer's GCC pain does
  not apply to a zero-dependency AWK project.
- Anonymized-header obfuscation mode as a core feature — out of scope for a
  simulator; if privacy-preserving sharing is ever wanted, it belongs in a
  separate small tool that consumes `*-truth.tsv`, not in the simulators.

## 4. Suggested roadmap

- Small (docs only): `LICENSE`, `CITATION.cff`, `CHANGES.md`, `TESTING.md`,
  glossary, archival note — closes most of the ReadMixer gap with no code risk.
- Medium (code, POSIX-safe): `--dry-run` disk estimate, `--loglevel`,
  transparent `.gz` reads, `--run-metadata` sidecar, benchmark table.
- Large (new capability): `--format fastq` + `--gzip` outputs, FASTA output
  splitting, RSS-bounded audit of `Barcode_Distance` feature tables.

## 5. References

- ReadMixer: <https://github.com/Geneinfosec-Inc/ReadMixer> (reviewed files:
  `README.md`, `readmixer/readmixer.py`, `readmixer/process/scramble.py`,
  `readmixer/process/util.py`, `TESTING.md`, `CHANGES.md`,
  `ARCHIVAL_STRATEGY.md`, `CITATION.cff`, `ReadMixer_SBOM.json`,
  `readmixer/INSTALLATION.md`, `readmixer/definitions.md`).
- This repo: `README.md`, `barcode_simulator.awk`, `barcode_experiment.awk`,
  `barcode_distance.awk` (notably `load_fasta()` `.gz` rejection),
  `barcode_simulator_post.awk`, `tests/*.sh`, `.github/workflows/test.yml`.
