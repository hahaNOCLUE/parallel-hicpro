# Rust post-alignment validation

## Scope

The pipeline has been tested with Micro-C and single-enzyme Hi-C. Multi-enzyme
Hi-C has not been tested. The primary supported reference is the GRCh38 no-alt
analysis set (`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`; Bowtie2 prefix
`GRCh38_noalt_as`). Reference-specific indexes, chromosome sizes, and digest
BED files must not be mixed between assemblies.

The detailed equivalence test below predates that primary-reference choice and
used Micro-C aligned to the GRCh38 no-alt plus hs38d1 decoy analysis set. It
validates the Rust Micro-C classifier, not the bundled single-enzyme example.
ICE normalization is intentionally outside the optimized path; create a
COOL/MCOOL file and run `cooler balance` only after sample-level valid-pair
merge and duplicate removal.

Technical replicates remain separate through pairing and interaction
classification. HiC-Pro's `merge_valid_interactions.sh` still performs the
sample-level merge and removes duplicates across all replicates.

## O33_1 validation

Input:

an internal GRCh38 no-alt-plus-hs38d1-decoy Micro-C sample (`O33_1`)

The Rust Micro-C classifier processed 48,412,320 pairs. Its sorted
`validPairs`, `FiltPairs`, `SinglePairs`, `DumpPairs`, and `RSstat` outputs
matched the existing Python outputs byte for byte.

| Class | Pairs |
| --- | ---: |
| Valid | 13,013,504 |
| Filtered by cis distance | 35,398,816 |
| Single-end | 0 |
| Dumped | 0 |

Valid orientations were FF 3,249,277; RR 3,248,902; RF 3,185,023; and FR
3,330,302. The Rust classification run used 64 configured threads and took
1:49.87 wall time with about 622 MiB maximum RSS. The following parallel sort
took 16.85 seconds. No Python baseline was timed in the same run, so these
numbers must not be interpreted as a speedup ratio.

## Breakpoints and shards

`annotations/GRCh38_noalt_decoy_as.breakpoints.high_confidence.tsv` contains
only chromosomes with an obvious centromeric reference N-run of at least 1 kb.
Each entry is one coordinate inside the N-run. It is a partition key, not an
excluded interval: positions below it go left and all other positions go
right, so no valid pair is discarded.

Sharding is optional (`RUST_SHARD_VALIDPAIRS = 0` by default). It is not yet
used by HiC-Pro's sample-level technical-replicate merge; enabling it currently
creates auxiliary files in addition to the canonical `validPairs` output.
The table is specific to `GRCh38_noalt_decoy_as`; do not use it with the
primary `GRCh38_noalt_as` reference. Set `HICPRO_BREAKPOINTS` explicitly only
when the table matches the alignment reference.

## Cooler equivalence

For the 500 bp validation map, the corrected HiC-Pro 1-based pair coordinates produced the same 25-chromosome axis, 6,176,584 bins, 11,987,614 nonzero pixels, and 12,209,953 total contacts as conversion of the HiC-Pro raw matrix. Chromosome, bin, pixel-count, and index datasets matched element-by-element.
