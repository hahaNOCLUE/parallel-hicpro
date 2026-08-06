# parallel-hicpro

Parallel Micro-C and Hi-C processing built on
[HiC-Pro 3.1.0](https://github.com/nservant/HiC-Pro). The slow Python pairing
and DNase/Micro-C valid-pair classification stages are replaced with Rust,
technical replicates retain HiC-Pro's sample-level merge/deduplication
semantics, and the final deliverable is a balanced multi-resolution `.mcool`
instead of an ICE-normalized matrix.

## Primary scope

> **Validated protocols: Micro-C and single-enzyme Hi-C. Multi-enzyme Hi-C has not been tested.**

The primary target is **GRCh38 no-alt analysis set**, specifically
`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`, used here with the Bowtie2
prefix `GRCh38_noalt_as`. Other assemblies or GRCh38 representations require
matching Bowtie2 indexes, chromosome sizes, restriction-fragment annotations
(for Hi-C), and independent output validation.

Micro-C uses the Rust no-restriction-fragment path (`GENOME_FRAGMENT` empty).
Single-enzyme Hi-C uses HiC-Pro's restriction-fragment classifier and requires
a BED file generated for the same FASTA and enzyme. Multi-enzyme digestion is
supported upstream by HiC-Pro but has not been tested in this fork.

The bundled `GRCh38_noalt_decoy_as` breakpoint table belongs to the distinct
no-alt-plus-hs38d1-decoy reference. It is optional, is not used by default, and
must not be used with `GRCh38_noalt_as` or another FASTA.

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

## Configuration

Copy [config-microc-grch38.example.txt](config-microc-grch38.example.txt) and
set at least:

- `BOWTIE2_IDX_PATH`
- `REFERENCE_GENOME`
- `GENOME_SIZE`
- read suffixes `PAIR1_EXT` and `PAIR2_EXT`

The example is configured for single-enzyme MboI Hi-C on `GRCh38_noalt_as`.
Set `GENOME_FRAGMENT` to the matching digest BED and `LIGATION_SITE` to the
appropriate ligation junction.

For Micro-C, clear both values:

```text
GENOME_FRAGMENT =
LIGATION_SITE =
```

`GENOME_SIZE` controls which chromosomes enter the final `.mcool`.

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
