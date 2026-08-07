# Barcode Simulator

Barcode Simulator creates controlled sets of related DNA barcode, amplicon, or
gene-copy sequences and samples them into FASTA files. It is a POSIX-AWK rewrite
of the Bash program originally written for an OTU-versus-k-mer comparison,
retaining the historical command-line interface without its GNU `getopt`,
`shuf`, `sed`, or `/dev/urandom` dependencies.

The simulator can vary:

- the number and length of reference clusters;
- the number of sequence variants generated per cluster;
- the exact number of substitutions in each variant; and
- the number of sequences sampled, without replacement, into each output file.

## Requirements

- a POSIX-compatible `awk`;
- a POSIX shell for the `Barcode_Simulator` launcher and test suite.

Linux and macOS are tested in GitHub Actions. No packages need to be installed.

## Quick start

```sh
chmod +x Barcode_Simulator

./Barcode_Simulator \
  --num-fasta 5 \
  --total-clusters 20 \
  --sequences-per-cluster 10 \
  --min-sequences-per-file 50 \
  --max-sequences-per-file 100 \
  --min-length 300 \
  --max-length 500 \
  --min-differences 0 \
  --max-differences 8 \
  --project-name barcode \
  --gen-ref barcode-reference.fasta \
  --save barcode-master.txt \
  --seed 2026
```

This writes `barcode-1.fasta` through `barcode-5.fasta`, a reference FASTA, and
a reusable master list. The same seed and AWK implementation reproduce the same
random run.

Run `./Barcode_Simulator --help` for the complete option list. The old options,
including `--total-otus`, `--number-of-seq-per-otu`, and
`--max-number-otus-per-file`, remain supported. New aliases use the more literal
terms “cluster” and “sequence.”

## Reusing barcode or gene-copy FASTA files

`--reuse` accepts either a normal FASTA file (including wrapped sequences) or a
plain file containing one DNA sequence per line:

```sh
./Barcode_Simulator \
  --reuse gene-copies.fasta \
  --num-fasta 10 \
  --min-sequences-per-file 25 \
  --max-sequences-per-file 50 \
  --project-name gene-copy-sample \
  --seed 73
```

Sampling is without replacement within each output file. A sequence may appear
again in a different output file. If the requested sample is larger than the
master set, every available sequence is written once.

## Simulation method

For each cluster, the program draws a reference length uniformly from the
configured inclusive range and generates an equiprobable A/C/G/T sequence. Each
variant receives a uniformly drawn number of substitutions. Mutation positions
are unique and replacement bases must differ from the reference, so the stated
number of differences is the exact Hamming distance from that cluster's
reference. Output sequences are selected with a Fisher–Yates shuffle.

AWK's pseudorandom-number algorithm is implementation-defined. `--seed` makes a
run reproducible on the same AWK implementation; byte-identical results are not
promised across different AWK implementations.

## Outputs

| Option | Output |
| --- | --- |
| `--project-name NAME` | `NAME-1.fasta`, `NAME-2.fasta`, and so on |
| `--gen-ref FILE` | One reference FASTA record per generated cluster |
| `--save FILE` | All generated/reused master sequences, one per line |
| `--only` | Suppresses sampled FASTAs while allowing master/reference output |

Output files are replaced when the same command is rerun; records are not
silently appended to stale results.

## Tests

```sh
sh tests/test_barcode_simulator.sh
```

The tests cover FASTA structure, exact mutation counts, seeded generation,
wrapped FASTA reuse, historical header behavior, validation, and `--only` mode.

## Historical analysis scripts

The `scripts/` directory contains the original downstream R, QIIME, Mash, and
cluster-submission scripts used for the 2018 analysis. They are retained as
historical research artifacts and are not required by the AWK simulator.
