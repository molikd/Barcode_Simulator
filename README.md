# Barcode Simulator

Barcode Simulator creates controlled sets of related DNA barcode, amplicon, or
gene-copy sequences and samples them into FASTA files. It is a POSIX-AWK rewrite
of the Bash program originally written for an OTU-versus-k-mer comparison,
retaining the historical command-line interface without its GNU `getopt`,
`shuf`, `sed`, or `/dev/urandom` dependencies.

The lower-level simulator can vary:

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
chmod +x Barcode_Simulator Barcode_Experiment Barcode_Simulator_Post

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
a reusable master list. The same seed and options reproduce the same random run.

Run `./Barcode_Simulator --help` for the complete option list. The old options,
including `--total-otus`, `--number-of-seq-per-otu`, and
`--max-number-otus-per-file`, remain supported. New aliases use the more literal
terms “cluster” and “sequence.”

## Experimental effects and factorial designs

`Barcode_Experiment` restores the experimental-design layer used in
[Molik, Pfrender, and Emrich (2020)](https://doi.org/10.3390/mps3010022).
It can independently enable the five properties examined in that study:

| Code | Effect | Baseline | Enabled behavior |
| --- | --- | --- | --- |
| `A` | Abundance distribution | clusters sampled uniformly | high, middling, and low tiers receive equal expected read mass across cluster groups in a 1:2:3 size ratio |
| `C` | Conserved region | no shared prefix | the historical 24 bp conserved sequence is prepended to every reference and variant |
| `E` | Sequence errors | zero substitutions | each variant receives 1–10 exact substitutions |
| `L` | Sequence length | fixed at 500 bp | cluster lengths vary uniformly from 350–500 bp |
| `N` | Sample depth (“picks”) | fixed at 136 reads | each sample contains 14–1360 reads |

Generate one treatment by listing effects or using the historical switches:

```sh
./Barcode_Experiment \
  --effects abundance,errors,depth \
  --project AEN \
  --samples 68 \
  --replicates 10 \
  --seed 2020

# Equivalent effect selection:
./Barcode_Experiment -A -E -N --project AEN --seed 2020
```

Generate the baseline and all 31 non-empty combinations as a complete
\(2^5\) factorial design:

```sh
./Barcode_Experiment \
  --factorial \
  --project metabarcode-factorial \
  --samples 68 \
  --replicates 10 \
  --seed 2020
```

The ordinary defaults reproduce the levels in the historical
`Run_Simulation.sh` (`136` fixed reads or `14`–`1360` variable reads). The
paper describes a tenfold depth and a complete 32-condition, 10-replicate
design. Use the explicit paper profile to generate those 320 datasets:

```sh
./Barcode_Experiment \
  --paper \
  --project metabarcode-paper \
  --seed 2020 \
  --no-truth
```

`--paper` selects 68 samples, 68 clusters, 10 variants per cluster, 1360 fixed
reads or 140–13,600 variable reads, 500 bp baseline sequences, 350–500 bp
variable sequences, 1–10 errors, and the full factorial. This manuscript-scale
run writes tens of gigabytes of FASTA, so plan storage accordingly. Every
level remains configurable; explicit size options override the profile.
Run `./Barcode_Experiment --help` for the full set.
Sampling in experiment mode is with replacement, representing sequencing
reads. Reference templates are paired within each replicate, and abundance or
depth changes do not regenerate the underlying variant pool.

Every experiment writes:

- one FASTA per treatment, replicate, and sample;
- one reference FASTA per treatment and replicate;
- a sample-level `*-manifest.tsv` recording the design matrix and read depth;
- a read-level `*-truth.tsv` mapping every read to its cluster, variant,
  abundance tier, sequence length, and number of substitutions.

Treatment names retain the historical codes (`A`, `C`, `E`, `L`, and `N`),
with `O` representing the no-effect baseline. Use `--no-truth` when a
read-level truth table would be unnecessarily large.

## Post-processing without R

`Barcode_Simulator_Post` replaces both `Barcode_Simulator_Post.R` and
`Barcode_Simulator_Post_Single.R`. It implements their numeric analyses in
POSIX AWK and writes tidy TSV data instead of opening an R plotting device:

```sh
# Mean OTU and Mash distances plus the Mantel/Pearson correlation
./Barcode_Simulator_Post distances --output distances.tsv

# Every paired point formerly sent to all_distances.png
./Barcode_Simulator_Post all-distances --output all_distances.tsv

# Correlate every unique pair of OTU/Mash distance matrices
./Barcode_Simulator_Post mantel --output mantel.tsv
```

When file arguments are omitted, the command discovers `*_otu_table*` and
`*_mash_dists*` inputs in the current directory. An OTU table passed to a
distance mode is paired with the corresponding `_mash_dists.txt` file.

The default quantitative Jaccard calculation matches the old
`vegan::vegdist(..., method="jaccard")` behavior. Use `--jaccard binary` for
presence/absence Jaccard distances. Standard `mash dist` output is treated as
distance data directly; the R scripts incorrectly subtracted those values from
one. `--mash-values similarity` remains available for already-inverted legacy
files.

`distances.tsv` and `all_distances.tsv` are ready for plotting in R, Python,
Vega-Lite, a spreadsheet, or another visualization tool without coupling the
simulation workflow to one graphics stack.

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

Both runners use an internal Park–Miller pseudorandom-number generator rather
than AWK's implementation-defined `rand()`. A supplied `--seed` therefore makes
the sequence and sampling streams reproducible across POSIX AWK implementations.

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
sh tests/test_barcode_experiment.sh
sh tests/test_barcode_post.sh
sh tests/test_paper_pipeline.sh
```

The tests cover FASTA structure, exact mutation counts, seeded generation,
wrapped FASTA reuse, historical header behavior, validation, `--only` mode,
each experimental effect, the complete 32-condition factorial, manifests,
truth tables, deterministic reruns, quantitative and binary Jaccard distances,
Mash parsing, and Mantel/Pearson correlations.

`test_paper_pipeline.sh` is a no-R, end-to-end regression of the experiment in
Section 2.3 of the paper. It runs all 32 A/C/E/L/N combinations with all ten
replicates (320 datasets), derives cluster-count OTU proxies and exact-sequence
ASV proxies, computes exact k-mer Jaccard distances in POSIX AWK, and sends all
640 method/dataset pairs through `Barcode_Simulator_Post`. It keeps the full
factorial and replication structure while scaling samples, read depth, and
sequence length down for CI. It checks that the error effect increases ASV
distance and that variable depth increases OTU distance.

This regression validates the AWK workflow and the qualitative effects; it is
not a claim of bit-for-bit reproduction of the published QIIME, DADA2, and Mash
outputs. Those tools use clustering, error models, and MinHash sketches rather
than the deterministic truth-based proxies used by the dependency-free test.

## Historical analysis scripts

The `scripts/` directory retains the original QIIME, Mash, and
cluster-submission scripts used for the 2018 analysis as historical research
artifacts. `Barcode_Experiment` and `Barcode_Simulator_Post` replace their
simulation and R orchestration without GNU utilities, a cluster scheduler, R,
`vegan`, `expss`, or `parallel`.
