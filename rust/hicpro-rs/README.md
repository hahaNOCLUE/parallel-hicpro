# hicpro-rs

Parallel Rust implementations of the CPU-bound HiC-Pro post-alignment stages.

Build in the reproducible environment:

```bash
mamba env create -f environment-rust.yml
mamba run -n hicpro-rust cargo build --release --manifest-path rust/hicpro-rs/Cargo.toml
```

`hicpro-rs pair` replaces `mergeSAM.py`. It preserves input order while filtering
read pairs in parallel batches and uses HTSlib threads for BAM I/O.

`hicpro-microc-rs` replaces `mapped_2hic_dnase.py` for Micro-C:

```bash
rust/hicpro-rs/target/release/hicpro-microc-rs \
  --mappedReadsFile sample.bwt2pairs.bam \
  --outputDir output \
  --minCisDist 1000 \
  --threads 16 \
  --breakpoints annotations/GRCh38_noalt_decoy_as.breakpoints.high_confidence.tsv \
  --shard-dir output/shards
```

Breakpoints are single zero-based coordinates. A pair is assigned from the 5'
coordinates using `< breakpoint` and `>= breakpoint`; no centromere interval is
excluded, so every valid pair is written exactly once to the main output and one
optional shard. Only chromosomes with a centromeric N run of at least 1 kb are
split in the high-confidence file.
