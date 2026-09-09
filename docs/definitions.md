# Definitions

Short glossary for readers who do not work with metabarcoding data every
day. Terms match their use in the README, `--help` output, manifests, and
truth tables.

| Term | Definition |
| --- | --- |
| barcode | A short DNA sequence used to identify or group biological samples or taxa. |
| cluster | A group of related simulated sequences descended from one generated reference. (Historical option names call these OTUs.) |
| variant | One mutated copy of a cluster's reference sequence; carries an exact Hamming distance from that reference. |
| Hamming distance | The number of substituted positions between a variant and its reference. Replacement bases always differ, so the count is exact. |
| reference | The unmutated per-cluster sequence written by `--gen-ref` or per-treatment reference files. |
| master list | The full generated (or reused) sequence pool from which output files are sampled; saved with `--save`, one sequence per line. |
| treatment | One experimental condition, named by its effect codes (`A`, `C`, `E`, `L`, `N`; `O` is the no-effect baseline). |
| replicate | One independent repeat of a treatment with its own derived random stream (same `--seed`, different replicate index). |
| sample | One output FASTA/FASTQ file of sequencing-like reads drawn for a treatment and replicate. |
| read | One sequence record in a sample file (with replacement in experiment mode, without replacement in `Barcode_Simulator` sampling). |
| manifest | The sample-level `*-manifest.tsv` design table: one row per sample file with treatment, replicate, depth, and effect levels. |
| truth table | The read-level `*-truth.tsv` ground-truth table mapping every read to its cluster, variant, abundance tier, length, and substitution count. |
| run metadata | The optional `--run-metadata` JSON sidecar recording effective parameters, seed, platform, and outputs for a run. |
| key file | ReadMixer term for its per-input CSV mapping scrambled reads back to originals; the analog here is the truth table plus manifest. |
| feature | In `Barcode_Distance`, one distinct complete barcode sequence; distances compare per-sample feature-count vectors. |
| dry run | An estimation-only invocation (`--dry-run`) that prints expected outputs and disk needs without writing files. |
