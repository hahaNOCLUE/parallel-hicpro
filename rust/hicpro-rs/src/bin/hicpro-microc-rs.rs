use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{Read, Record};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    version,
    about = "Parallel HiC-Pro interaction classification for Micro-C"
)]
struct Args {
    #[arg(short = 'r', long = "mappedReadsFile")]
    mapped_reads: PathBuf,

    #[arg(short = 'o', long = "outputDir", default_value = ".")]
    output_dir: PathBuf,

    #[arg(short = 'd', long = "minCisDist")]
    min_cis_dist: Option<i64>,

    #[arg(short = 'g', long = "gtag")]
    genotype_tag: Option<String>,

    #[arg(short = 'a', long = "all")]
    all_outputs: bool,

    #[arg(short = 'v', long = "verbose")]
    verbose: bool,

    #[arg(long, default_value_t = 1)]
    threads: usize,

    #[arg(long, default_value_t = 65_536)]
    batch_size: usize,

    /// Two-column TSV: chromosome and a zero-based breakpoint.
    #[arg(long)]
    breakpoints: Option<PathBuf>,

    /// Additionally write valid pairs into chromosome-side pair shards.
    #[arg(long, requires = "breakpoints")]
    shard_dir: Option<PathBuf>,
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum Class {
    Valid,
    Singleton,
    Filtered,
    Dump,
}

impl Class {
    fn suffix(self) -> &'static str {
        match self {
            Self::Valid => "validPairs",
            Self::Singleton => "SinglePairs",
            Self::Filtered => "FiltPairs",
            Self::Dump => "DumpPairs",
        }
    }
}

struct ResultRow {
    class: Class,
    line: Option<String>,
    orientation: Option<&'static str>,
    genotype: Option<(i64, i64)>,
    shard: Option<String>,
}

#[derive(Default)]
struct Stats {
    valid: u64,
    ff: u64,
    rr: u64,
    rf: u64,
    fr: u64,
    singleton: u64,
    filtered: u64,
    dumped: u64,
    g1g1: u64,
    g2g2: u64,
    g1u: u64,
    ug1: u64,
    g2u: u64,
    ug2: u64,
    g1g2: u64,
    g2g1: u64,
    uu: u64,
    conflict: u64,
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.threads == 0 || args.batch_size == 0 {
        bail!("--threads and --batch-size must be at least 1");
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;
    run(args)
}

fn run(args: Args) -> Result<()> {
    fs::create_dir_all(&args.output_dir)?;
    if let Some(dir) = &args.shard_dir {
        fs::create_dir_all(dir)?;
    }
    let breakpoints = load_breakpoints(args.breakpoints.as_deref())?;
    let tag = parse_tag(args.genotype_tag.as_deref())?;

    let mut reader = bam::Reader::from_path(&args.mapped_reads)
        .with_context(|| format!("cannot open {}", args.mapped_reads.display()))?;
    reader.set_threads(io_threads(args.threads))?;
    let names: Vec<String> = reader
        .header()
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect();

    let base = base_name(&args.mapped_reads)?;
    let mut outputs = open_outputs(&args.output_dir, &base, args.all_outputs)?;
    let mut shards: HashMap<String, BufWriter<File>> = HashMap::new();
    let mut records = reader.records();
    let mut stats = Stats::default();
    let mut processed = 0_u64;

    loop {
        let mut batch = Vec::with_capacity(args.batch_size);
        for _ in 0..args.batch_size {
            let Some(r1) = records.next() else { break };
            let r2 = records
                .next()
                .context("paired BAM ends with an incomplete read pair")??;
            batch.push((r1?, r2));
        }
        if batch.is_empty() {
            break;
        }

        let rows: Result<Vec<_>> = batch
            .into_par_iter()
            .map(|(r1, r2)| classify(r1, r2, &names, &breakpoints, tag.as_ref(), &args))
            .collect();
        for row in rows? {
            update_stats(&mut stats, &row);
            if let Some(line) = row.line {
                outputs
                    .get_mut(row.class.suffix())
                    .unwrap()
                    .write_all(line.as_bytes())?;
                if let (Some(dir), Some(key)) = (&args.shard_dir, row.shard) {
                    if !shards.contains_key(&key) {
                        let path = dir.join(format!("{base}.{key}.validPairs"));
                        shards.insert(key.clone(), BufWriter::new(File::create(path)?));
                    }
                    shards.get_mut(&key).unwrap().write_all(line.as_bytes())?;
                }
            }
            processed += 1;
        }
        if args.verbose && processed % 100_000 < args.batch_size as u64 {
            eprintln!("## processed {processed} read pairs");
        }
    }
    drop(outputs);
    drop(shards);
    write_stats(
        &args.output_dir.join(format!("{base}.RSstat")),
        &stats,
        tag.is_some(),
    )?;
    Ok(())
}

fn classify(
    r1: Record,
    r2: Record,
    names: &[String],
    breakpoints: &HashMap<String, i64>,
    tag: Option<&[u8; 2]>,
    args: &Args,
) -> Result<ResultRow> {
    if !r1.is_first_in_template() || !r2.is_last_in_template() {
        bail!(
            "expected adjacent read1/read2 records at {}",
            String::from_utf8_lossy(r1.qname())
        );
    }
    if normalized_name(r1.qname()) != normalized_name(r2.qname()) {
        bail!(
            "adjacent records have different names: {} and {}",
            String::from_utf8_lossy(r1.qname()),
            String::from_utf8_lossy(r2.qname())
        );
    }

    let distance = if !r1.is_unmapped() && !r2.is_unmapped() && r1.tid() == r2.tid() {
        Some((read_start(&r1) - read_start(&r2)).abs())
    } else {
        None
    };
    let class = if r1.is_unmapped() || r2.is_unmapped() {
        Class::Singleton
    } else if args
        .min_cis_dist
        .is_some_and(|min| distance.is_some_and(|value| value < min))
    {
        Class::Filtered
    } else {
        Class::Valid
    };

    let orientation = (class == Class::Valid).then(|| orientation(&r1, &r2));
    let genotype = tag.map(|t| {
        (
            aux_integer(&r1, t).unwrap_or(-1),
            aux_integer(&r2, t).unwrap_or(-1),
        )
    });
    let should_write = class == Class::Valid || args.all_outputs;
    let line = should_write.then(|| valid_pairs_line(&r1, &r2, names, tag));
    let shard = if class == Class::Valid && !breakpoints.is_empty() {
        let (left, _) = ordered_pair(&r1, &r2);
        Some(side_name(
            &names[left.tid() as usize],
            read_start(left),
            breakpoints,
        ))
    } else {
        None
    };
    Ok(ResultRow {
        class,
        line,
        orientation,
        genotype,
        shard,
    })
}

fn normalized_name(name: &[u8]) -> &[u8] {
    &name[..name
        .iter()
        .position(|b| *b == b'/' || *b == b' ')
        .unwrap_or(name.len())]
}

fn aligned_reference_len(record: &Record) -> i64 {
    record
        .cigar()
        .iter()
        .map(|op| match op {
            Cigar::Match(n)
            | Cigar::Del(n)
            | Cigar::RefSkip(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n) => *n as i64,
            _ => 0,
        })
        .sum()
}

fn read_start(record: &Record) -> i64 {
    if record.is_reverse() {
        record.pos() + aligned_reference_len(record) - 1
    } else {
        record.pos()
    }
}

fn ordered_pair<'a>(r1: &'a Record, r2: &'a Record) -> (&'a Record, &'a Record) {
    let first_is_r1 = if r1.tid() == r2.tid() {
        read_start(r1) < read_start(r2)
    } else {
        r1.tid() < r2.tid()
    };
    if first_is_r1 {
        (r1, r2)
    } else {
        (r2, r1)
    }
}

fn orientation(r1: &Record, r2: &Record) -> &'static str {
    let (left, right) = ordered_pair(r1, r2);
    match (left.is_reverse(), right.is_reverse()) {
        (false, false) => "FF",
        (true, true) => "RR",
        (false, true) => "FR",
        (true, false) => "RF",
    }
}

fn strand(record: &Record) -> char {
    if record.is_reverse() {
        '-'
    } else {
        '+'
    }
}

fn valid_pairs_line(r1: &Record, r2: &Record, names: &[String], tag: Option<&[u8; 2]>) -> String {
    if !r1.is_unmapped() && !r2.is_unmapped() {
        let (left, right) = ordered_pair(r1, r2);
        let htag = tag
            .map(|t| {
                format!(
                    "{}-{}",
                    aux_integer(left, t).map_or("None".to_string(), |v| v.to_string()),
                    aux_integer(right, t).map_or("None".to_string(), |v| v.to_string())
                )
            })
            .unwrap_or_default();
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\tNA\tNA\tNA\t{}\t{}\t{}\n",
            String::from_utf8_lossy(left.qname()),
            names[left.tid() as usize],
            read_start(left) + 1,
            strand(left),
            names[right.tid() as usize],
            read_start(right) + 1,
            strand(right),
            left.mapq(),
            right.mapq(),
            htag
        )
    } else if !r1.is_unmapped() {
        format!(
            "{}\t{}\t{}\t{}\t*\t*\t*\t*\t*\t*\t{}\t*\n",
            String::from_utf8_lossy(r1.qname()),
            names[r1.tid() as usize],
            read_start(r1) + 1,
            strand(r1),
            r1.mapq()
        )
    } else {
        format!(
            "{}\t*\t*\t*\t{}\t{}\t{}\t*\t*\t*\t*\t{}\n",
            String::from_utf8_lossy(r2.qname()),
            names[r2.tid() as usize],
            read_start(r2) + 1,
            strand(r2),
            r2.mapq()
        )
    }
}

fn parse_tag(tag: Option<&str>) -> Result<Option<[u8; 2]>> {
    match tag {
        None => Ok(None),
        Some(value) if value.is_ascii() && value.len() == 2 => {
            Ok(Some([value.as_bytes()[0], value.as_bytes()[1]]))
        }
        Some(_) => bail!("BAM auxiliary tag must be two ASCII characters"),
    }
}

fn aux_integer(record: &Record, tag: &[u8; 2]) -> Option<i64> {
    match record.aux(tag).ok()? {
        Aux::I8(v) => Some(v as i64),
        Aux::U8(v) => Some(v as i64),
        Aux::I16(v) => Some(v as i64),
        Aux::U16(v) => Some(v as i64),
        Aux::I32(v) => Some(v as i64),
        Aux::U32(v) => Some(v as i64),
        _ => None,
    }
}

fn load_breakpoints(path: Option<&Path>) -> Result<HashMap<String, i64>> {
    let Some(path) = path else {
        return Ok(HashMap::new());
    };
    let mut points = HashMap::new();
    for (line_no, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let chrom = fields.next().unwrap();
        let point = fields
            .next()
            .with_context(|| format!("missing breakpoint at line {}", line_no + 1))?
            .parse()?;
        if points.insert(chrom.to_string(), point).is_some() {
            bail!("duplicate breakpoint for {chrom}");
        }
    }
    Ok(points)
}

fn side_name(chrom: &str, pos: i64, points: &HashMap<String, i64>) -> String {
    match points.get(chrom) {
        Some(point) if pos < *point => format!("{chrom}L"),
        Some(_) => format!("{chrom}R"),
        None if is_primary_chrom(chrom) => chrom.to_string(),
        None => "other".to_string(),
    }
}

fn is_primary_chrom(chrom: &str) -> bool {
    if matches!(chrom, "chrX" | "chrY" | "chrM") {
        return true;
    }
    chrom
        .strip_prefix("chr")
        .and_then(|value| value.parse::<u8>().ok())
        .is_some_and(|value| (1..=22).contains(&value))
}

fn open_outputs(
    dir: &Path,
    base: &str,
    all: bool,
) -> Result<HashMap<&'static str, BufWriter<File>>> {
    let classes = if all {
        vec![Class::Valid, Class::Singleton, Class::Filtered, Class::Dump]
    } else {
        vec![Class::Valid]
    };
    let mut outputs = HashMap::new();
    for class in classes {
        outputs.insert(
            class.suffix(),
            BufWriter::new(File::create(
                dir.join(format!("{base}.{}", class.suffix())),
            )?),
        );
    }
    Ok(outputs)
}

fn update_stats(stats: &mut Stats, row: &ResultRow) {
    match row.class {
        Class::Valid => {
            stats.valid += 1;
            match row.orientation {
                Some("FF") => stats.ff += 1,
                Some("RR") => stats.rr += 1,
                Some("RF") => stats.rf += 1,
                Some("FR") => stats.fr += 1,
                _ => stats.dumped += 1,
            }
        }
        Class::Singleton => stats.singleton += 1,
        Class::Filtered => stats.filtered += 1,
        Class::Dump => stats.dumped += 1,
    }
    if let Some(pair) = row.genotype {
        match pair {
            (1, 1) => stats.g1g1 += 1,
            (2, 2) => stats.g2g2 += 1,
            (1, 0) => stats.g1u += 1,
            (0, 1) => stats.ug1 += 1,
            (2, 0) => stats.g2u += 1,
            (0, 2) => stats.ug2 += 1,
            (1, 2) => stats.g1g2 += 1,
            (2, 1) => stats.g2g1 += 1,
            (3, _) | (_, 3) => stats.conflict += 1,
            _ => stats.uu += 1,
        }
    }
}

fn write_stats(path: &Path, stats: &Stats, allele: bool) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "## Hi-C processing - no restriction fragments")?;
    writeln!(out, "Valid_interaction_pairs\t{}", stats.valid)?;
    writeln!(out, "Valid_interaction_pairs_FF\t{}", stats.ff)?;
    writeln!(out, "Valid_interaction_pairs_RR\t{}", stats.rr)?;
    writeln!(out, "Valid_interaction_pairs_RF\t{}", stats.rf)?;
    writeln!(out, "Valid_interaction_pairs_FR\t{}", stats.fr)?;
    writeln!(out, "Single-end_pairs\t{}", stats.singleton)?;
    writeln!(out, "Filtered_pairs\t{}", stats.filtered)?;
    writeln!(out, "Dumped_pairs\t{}", stats.dumped)?;
    if allele {
        writeln!(out, "## ======================================")?;
        writeln!(out, "## Allele specific information")?;
        writeln!(out, "Valid_pairs_from_ref_genome_(1-1)\t{}", stats.g1g1)?;
        writeln!(
            out,
            "Valid_pairs_from_ref_genome_with_one_unassigned_mate_(0-1/1-0)\t{}",
            stats.g1u + stats.ug1
        )?;
        writeln!(out, "Valid_pairs_from_alt_genome_(2-2)\t{}", stats.g2g2)?;
        writeln!(
            out,
            "Valid_pairs_from_alt_genome_with_one_unassigned_mate_(0-2/2-0)\t{}",
            stats.g2u + stats.ug2
        )?;
        writeln!(
            out,
            "Valid_pairs_from_alt_and_ref_genome_(1-2/2-1)\t{}",
            stats.g1g2 + stats.g2g1
        )?;
        writeln!(
            out,
            "Valid_pairs_with_both_unassigned_mated_(0-0)\t{}",
            stats.uu
        )?;
        writeln!(
            out,
            "Valid_pairs_with_at_least_one_conflicting_mate_(3-)\t{}",
            stats.conflict
        )?;
    }
    Ok(())
}

fn base_name(path: &Path) -> Result<String> {
    let name = path
        .file_name()
        .context("input path has no file name")?
        .to_string_lossy();
    Ok(name
        .strip_suffix(".bam")
        .or_else(|| name.strip_suffix(".sam"))
        .unwrap_or(&name)
        .to_string())
}

fn io_threads(total: usize) -> usize {
    total.saturating_div(4).clamp(1, 4)
}
