# parallel-hicpro

Micro-C-oriented parallel post-alignment processing built on
[HiC-Pro 3.1.0](https://github.com/nservant/HiC-Pro). The slow Python pairing
and DNase/Micro-C valid-pair classification stages are replaced with Rust,
technical replicates retain HiC-Pro's sample-level merge/deduplication
semantics, and the final deliverable is a balanced multi-resolution `.mcool`
instead of an ICE-normalized matrix.

## Primary scope

> **This fork is a focused Micro-C pipeline, not a generally validated replacement for every HiC-Pro protocol or reference.**

The optimized and validated target is **Micro-C aligned to the GRCh38 no-alt plus hs38d1 decoy analysis set**, specifically the reference distributed as `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna` and used here with the Bowtie2 prefix `GRCh38_noalt_decoy_as`. The supplied centromeric N-run breakpoints, real-data equivalence tests, and default GRCh38 example configuration all refer to this exact reference build.

The Rust Micro-C classifier assumes the no-restriction-fragment HiC-Pro path (`GENOME_FRAGMENT` empty). Other assemblies, alternate GRCh38 references, and restriction-enzyme Hi-C may require new chromosome sizes, breakpoints, configuration, and independent output validation. Do not reuse the bundled breakpoint table on another FASTA.

The original HiC-Pro README is preserved as
[README_UPSTREAM.md](README_UPSTREAM.md). Please cite HiC-Pro when using this
pipeline.

## What changed

- `mergeSAM.py` replacement with parallel batching and threaded BAM I/O.
- `mapped_2hic_dnase.py` replacement for Micro-C, validated against its output.
- Parallel GNU sort after valid-pair generation.
- Optional lossless sharding at one unalignable reference coordinate inside
  obvious centromeric N-runs. The shards are intermediate only; final output is
  merged genome-wide.
- Sample-level technical-replicate merge and cross-replicate duplicate removal.
- `allValidPairs -> cool -> mcool -> cooler balance` using the configured
  chromosome set and resolutions.

ICE is intentionally not part of the optimized workflow.

## Installation

Requirements: Linux, mamba, GNU coreutils, Bowtie2, and enough temporary disk
space for BAM sorting and Cooler chunk files.

```bash
git clone https://github.com/hahaNOCLUE/parallel-hicpro.git
cd parallel-hicpro
mamba env create -f environment-rust.yml
mamba run -n hicpro-rust cargo build --release \
  --manifest-path rust/hicpro-rs/Cargo.toml
cp config-system.example.txt config-system.txt
```

Edit `config-system.txt` and replace `/path/to/parallel-hicpro` with the clone
path. `config-system.txt` is ignored because it is machine-specific.

## Micro-C configuration

Copy [config-microc-grch38.example.txt](config-microc-grch38.example.txt) and
set at least:

- `BOWTIE2_IDX_PATH`
- `REFERENCE_GENOME`
- `GENOME_SIZE`
- read suffixes `PAIR1_EXT` and `PAIR2_EXT`

`GENOME_FRAGMENT` must remain empty for Micro-C. `GENOME_SIZE` controls which
chromosomes enter the final `.mcool`; alignments to decoys may be retained in
validPairs while the standard contact map uses primary chromosomes only.

The default resolutions are:

```text
500 1000 2500 5000 10000 25000 50000 100000 500000
```

## Input layout and technical replicates

Put every technical replicate of one biological sample in the same sample
directory. File prefixes must remain distinct.

```text
input/
└── sample1/
    ├── runA_sample1_1.fastp.fq.gz
    ├── runA_sample1_2.fastp.fq.gz
    ├── runB_sample1_1.fastp.fq.gz
    └── runB_sample1_2.fastp.fq.gz
```

Each replicate is mapped and classified independently. HiC-Pro then merges all
replicate `.validPairs` files and removes duplicates across the complete
sample. Chromosome/centromere shards, when enabled, are also intermediate and
must never be treated as separate final maps.

## Run

Activate the environment, then run only the stages used by this fork:

```bash
mamba activate hicpro-rust
export USE_RUST=1
export COOLER_BIN="$CONDA_PREFIX/bin/cooler"

bin/HiC-Pro \
  -i /path/to/input \
  -o /path/to/output \
  -c /path/to/config-microc.txt \
  -s mapping \
  -s proc_hic \
  -s merge_persample \
  -s build_contact_maps
```

The final file is:

```text
output/hic_results/mcool/<sample>.mcool
```

All configured resolutions are stored in one genome-wide file. `cooler
balance` writes a `weight` column at every resolution. Very sparse resolutions
can legitimately fail to converge; inspect the per-sample `.mcool.log` and the
Cooler `converged` metadata before quantitative use.

## Validation

See [VALIDATION.md](VALIDATION.md). On the full validation dataset, Rust
classification reproduced Python `validPairs`, interaction-class files, and
statistics byte-for-byte. The corrected 1-based HiC-Pro-to-Cooler binning was
also checked pixel-by-pixel against HiC-Pro's raw 500 bp matrix.

## License and attribution

This repository retains the upstream HiC-Pro license and history. The Rust
extensions are distributed under the same repository license. HiC-Pro citation:

Servant N. et al. *HiC-Pro: an optimized and flexible pipeline for Hi-C data
processing.* Genome Biology 16, 259 (2015).
