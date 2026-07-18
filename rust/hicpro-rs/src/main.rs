use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{Read, Record};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

const DEFAULT_BATCH_SIZE: usize = 65_536;

#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Pair name-sorted R1 and R2 BAM files and apply HiC-Pro mapping filters.
    Pair(PairArgs),
    /// Classify paired alignments against restriction fragments.
    Classify(ClassifyArgs),
}

#[derive(Args)]
struct CommonParallelArgs {
    /// Total worker count, including BAM compression/decompression workers.
    #[arg(long, default_value_t = 1)]
    threads: usize,

    /// Number of read pairs processed in each parallel batch.
    #[arg(long, default_value_t = DEFAULT_BATCH_SIZE)]
    batch_size: usize,
}

#[derive(Args)]
struct PairArgs {
    #[arg(short = 'f', long = "forward")]
    forward: PathBuf,

    #[arg(short = 'r', long = "reverse")]
    reverse: PathBuf,

    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    #[arg(short = 'q', long = "qual")]
    min_mapq: Option<u8>,

    #[arg(short = 's', long = "single")]
    report_single: bool,

    #[arg(short = 'm', long = "multi")]
    report_multi: bool,

    #[arg(short = 't', long = "stat")]
    stat: bool,

    #[arg(short = 'v', long = "verbose")]
    verbose: bool,

    #[command(flatten)]
    parallel: CommonParallelArgs,
}

#[derive(Args)]
struct ClassifyArgs {
    #[arg(short = 'f', long = "fragmentFile")]
    fragment_file: PathBuf,

    #[arg(short = 'r', long = "mappedReadsFile")]
    mapped_reads: PathBuf,

    #[arg(short = 'o', long = "outputDir", default_value = ".")]
    output_dir: PathBuf,

    #[arg(short = 's', long = "shortestInsertSize")]
    min_insert_size: Option<i64>,

    #[arg(short = 'l', long = "longestInsertSize")]
    max_insert_size: Option<i64>,

    #[arg(short = 't', long = "shortestFragmentLength")]
    min_frag_size: Option<i64>,

    #[arg(short = 'm', long = "longestFragmentLength")]
    max_frag_size: Option<i64>,

    #[arg(short = 'd', long = "minCisDist")]
    min_cis_dist: Option<i64>,

    #[arg(short = 'g', long = "gtag")]
    genotype_tag: Option<String>,

    #[arg(short = 'a', long = "all")]
    all_outputs: bool,

    #[arg(short = 'S', long = "sam")]
    interaction_bam: bool,

    #[arg(short = 'v', long = "verbose")]
    verbose: bool,

    #[command(flatten)]
    parallel: CommonParallelArgs,
}

#[derive(Default)]
struct PairStats {
    processed: u64,
    unmapped: u64,
    lowq_pairs: u64,
    unique_pairs: u64,
    multi_pairs: u64,
    singletons: u64,
    lowq_singletons: u64,
    unique_singletons: u64,
    multi_singletons: u64,
    reported: u64,
}

enum PairStatus {
    Unmapped,
    LowQualityPair,
    UniquePair,
    MultiPair,
    Singleton,
    LowQualitySingleton,
    UniqueSingleton,
    MultiSingleton,
}

struct PairResult {
    status: PairStatus,
    records: Option<(Record, Record)>,
}

#[derive(Clone)]
struct Fragment {
    start: i64,
    end: i64,
    name: Arc<str>,
    filtered: bool,
}

type FragmentMap = HashMap<String, Vec<Fragment>>;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Interaction {
    Valid,
    DanglingEnd,
    Religation,
    SelfCircle,
    Singleton,
    Filtered,
    Dump,
}

impl Interaction {
    fn code(self) -> &'static str {
        match self {
            Self::Valid => "VI",
            Self::DanglingEnd => "DE",
            Self::Religation => "RE",
            Self::SelfCircle => "SC",
            Self::Singleton => "SI",
            Self::Filtered => "FILT",
            Self::Dump => "DUMP",
        }
    }

    fn suffix(self) -> &'static str {
        match self {
            Self::Valid => "validPairs",
            Self::DanglingEnd => "DEPairs",
            Self::Religation => "REPairs",
            Self::SelfCircle => "SCPairs",
            Self::Singleton => "SinglePairs",
            Self::Filtered => "FiltPairs",
            Self::Dump => "DumpPairs",
        }
    }
}

#[derive(Default)]
struct ClassStats {
    valid: u64,
    ff: u64,
    rr: u64,
    rf: u64,
    fr: u64,
    dangling: u64,
    religation: u64,
    self_circle: u64,
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

struct ClassResult {
    interaction: Interaction,
    line: Option<String>,
    orientation: Option<&'static str>,
    genotype_pair: Option<(i64, i64)>,
    records: Option<(Record, Record)>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Pair(args) => run_pair(args),
        Command::Classify(args) => run_classify(args),
    }
}

fn validate_parallel(args: &CommonParallelArgs) -> Result<()> {
    if args.threads == 0 {
        bail!("--threads must be at least 1");
    }
    if args.batch_size == 0 {
        bail!("--batch-size must be at least 1");
    }
    Ok(())
}

fn configure_pool(threads: usize) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .context("failed to initialize Rayon thread pool")
}

fn io_threads(total: usize) -> usize {
    total.saturating_div(4).clamp(1, 4)
}

fn normalized_qname(name: &[u8]) -> &[u8] {
    let end = name
        .iter()
        .position(|b| *b == b'/' || *b == b' ')
        .unwrap_or(name.len());
    &name[..end]
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

fn is_unique_bowtie2(record: &Record) -> bool {
    if record.is_unmapped() {
        return false;
    }
    let Some(primary) = aux_integer(record, b"AS") else {
        return false;
    };
    match aux_integer(record, b"XS") {
        Some(secondary) => primary > secondary,
        None => true,
    }
}

fn set_pair_fields(r1: &mut Record, r2: &mut Record) {
    let mut f1 = r1.flags();
    let mut f2 = r2.flags();

    if f1 & 0x4 != 0 {
        f1 |= 0x8;
    }
    if f2 & 0x4 != 0 {
        f2 |= 0x8;
    }
    if f1 & 0x4 == 0 && f2 & 0x4 == 0 {
        f1 |= 0x1 | 0x2;
        f2 |= 0x1 | 0x2;
    }
    if f1 & 0x10 != 0 {
        f2 |= 0x20;
    }
    if f2 & 0x10 != 0 {
        f1 |= 0x20;
    }
    f1 |= 0x40;
    f2 |= 0x80;
    r1.set_flags(f1);
    r2.set_flags(f2);

    if r1.tid() == r2.tid() {
        r1.set_mtid(r1.tid());
        r2.set_mtid(r1.tid());
    } else {
        r1.set_mtid(r2.tid());
        r2.set_mtid(r1.tid());
    }
    r1.set_mpos(r2.pos());
    r2.set_mpos(r1.pos());
}

fn process_mapping_pair(mut r1: Record, mut r2: Record, args: &PairArgs) -> Result<PairResult> {
    if normalized_qname(r1.qname()) != normalized_qname(r2.qname()) {
        bail!(
            "forward and reverse reads are not paired: {} vs {}",
            String::from_utf8_lossy(r1.qname()),
            String::from_utf8_lossy(r2.qname())
        );
    }

    let status;
    let keep;
    if r1.is_unmapped() && r2.is_unmapped() {
        status = PairStatus::Unmapped;
        keep = false;
    } else if !r1.is_unmapped() && !r2.is_unmapped() {
        if args
            .min_mapq
            .is_some_and(|q| r1.mapq() < q || r2.mapq() < q)
        {
            status = PairStatus::LowQualityPair;
            keep = false;
        } else if is_unique_bowtie2(&r1) && is_unique_bowtie2(&r2) {
            status = PairStatus::UniquePair;
            keep = true;
        } else {
            status = PairStatus::MultiPair;
            keep = args.report_multi;
        }
    } else {
        if !args.report_single {
            return Ok(PairResult {
                status: PairStatus::Singleton,
                records: None,
            });
        }
        let mapped = if r1.is_unmapped() { &r2 } else { &r1 };
        if args.min_mapq.is_some_and(|q| mapped.mapq() < q) {
            status = PairStatus::LowQualitySingleton;
            keep = false;
        } else if is_unique_bowtie2(mapped) {
            status = PairStatus::UniqueSingleton;
            keep = true;
        } else {
            status = PairStatus::MultiSingleton;
            keep = args.report_multi;
        }
    }

    if keep {
        set_pair_fields(&mut r1, &mut r2);
        Ok(PairResult {
            status,
            records: Some((r1, r2)),
        })
    } else {
        Ok(PairResult {
            status,
            records: None,
        })
    }
}

fn update_pair_stats(stats: &mut PairStats, result: &PairResult) {
    stats.processed += 1;
    match result.status {
        PairStatus::Unmapped => stats.unmapped += 1,
        PairStatus::LowQualityPair => stats.lowq_pairs += 1,
        PairStatus::UniquePair => stats.unique_pairs += 1,
        PairStatus::MultiPair => stats.multi_pairs += 1,
        PairStatus::Singleton => stats.singletons += 1,
        PairStatus::LowQualitySingleton => {
            stats.singletons += 1;
            stats.lowq_singletons += 1;
        }
        PairStatus::UniqueSingleton => {
            stats.singletons += 1;
            stats.unique_singletons += 1;
        }
        PairStatus::MultiSingleton => {
            stats.singletons += 1;
            stats.multi_singletons += 1;
        }
    }
    if result.records.is_some() {
        stats.reported += 1;
    }
}

fn run_pair(args: PairArgs) -> Result<()> {
    validate_parallel(&args.parallel)?;
    configure_pool(args.parallel.threads)?;

    let mut r1_reader = bam::Reader::from_path(&args.forward)
        .with_context(|| format!("cannot open {}", args.forward.display()))?;
    let mut r2_reader = bam::Reader::from_path(&args.reverse)
        .with_context(|| format!("cannot open {}", args.reverse.display()))?;
    let bgzf_threads = io_threads(args.parallel.threads);
    r1_reader.set_threads(bgzf_threads)?;
    r2_reader.set_threads(bgzf_threads)?;

    let header = bam::Header::from_template(r1_reader.header());
    let mut writer = bam::Writer::from_path(&args.output, &header, bam::Format::Bam)
        .with_context(|| format!("cannot create {}", args.output.display()))?;
    writer.set_threads(bgzf_threads)?;

    let mut it1 = r1_reader.records().fuse();
    let mut it2 = r2_reader.records().fuse();
    let mut stats = PairStats::default();

    loop {
        let mut batch = Vec::with_capacity(args.parallel.batch_size);
        for _ in 0..args.parallel.batch_size {
            match (it1.next(), it2.next()) {
                (None, None) => break,
                (Some(a), Some(b)) => batch.push((a?, b?)),
                _ => bail!("R1 and R2 BAM files contain different record counts"),
            }
        }
        if batch.is_empty() {
            break;
        }

        let results: Result<Vec<_>> = batch
            .into_par_iter()
            .map(|(r1, r2)| process_mapping_pair(r1, r2, &args))
            .collect();
        for result in results? {
            update_pair_stats(&mut stats, &result);
            if let Some((r1, r2)) = result.records {
                writer.write(&r1)?;
                writer.write(&r2)?;
            }
        }
        if args.verbose && stats.processed % 1_000_000 < args.parallel.batch_size as u64 {
            eprintln!("## {} read pairs", stats.processed);
        }
    }
    drop(writer);

    if args.stat {
        let stat_path = args.output.with_extension("pairstat");
        write_pair_stats(&stat_path, &stats)?;
    }
    Ok(())
}

fn python_percent(value: u64, total: u64) -> String {
    if total == 0 {
        return "0.0".to_string();
    }
    let mut text = format!("{:.3}", value as f64 * 100.0 / total as f64);
    while text.ends_with('0') && text.contains('.') {
        text.pop();
    }
    if text.ends_with('.') {
        text.push('0');
    }
    text
}

fn write_pair_stats(path: &Path, s: &PairStats) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    let rows = [
        ("Total_pairs_processed", s.processed),
        ("Unmapped_pairs", s.unmapped),
        ("Low_qual_pairs", s.lowq_pairs),
        ("Unique_paired_alignments", s.unique_pairs),
        ("Multiple_pairs_alignments", s.multi_pairs),
        ("Pairs_with_singleton", s.singletons),
        ("Low_qual_singleton", s.lowq_singletons),
        ("Unique_singleton_alignments", s.unique_singletons),
        ("Multiple_singleton_alignments", s.multi_singletons),
        ("Reported_pairs", s.reported),
    ];
    for (name, value) in rows {
        writeln!(
            out,
            "{name}\t{value}\t{}",
            python_percent(value, s.processed)
        )?;
    }
    Ok(())
}

fn load_fragments(args: &ClassifyArgs) -> Result<FragmentMap> {
    let input = BufReader::new(
        File::open(&args.fragment_file)
            .with_context(|| format!("cannot open {}", args.fragment_file.display()))?,
    );
    let mut fragments: FragmentMap = HashMap::new();
    for (line_no, line) in input.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            bail!("invalid fragment BED at line {}", line_no + 1);
        }
        let start: i64 = fields[1].parse()?;
        let end: i64 = fields[2].parse()?;
        let length = (end - start).abs();
        let filtered = args.min_frag_size.is_some_and(|v| length < v)
            || args.max_frag_size.is_some_and(|v| length > v);
        fragments
            .entry(fields[0].to_string())
            .or_default()
            .push(Fragment {
                start,
                end,
                name: Arc::from(fields[3]),
                filtered,
            });
    }
    for list in fragments.values_mut() {
        list.sort_unstable_by_key(|f| f.start);
        for adjacent in list.windows(2) {
            if adjacent[0].end > adjacent[1].start {
                bail!("overlapping restriction fragments are not supported");
            }
        }
    }
    Ok(fragments)
}

fn reference_names(header: &bam::HeaderView) -> Vec<String> {
    header
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect()
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

fn read_middle(record: &Record) -> i64 {
    record.pos() + aligned_reference_len(record) / 2
}

fn ordered_pair<'a>(r1: &'a Record, r2: &'a Record) -> (&'a Record, &'a Record, bool) {
    let first_is_r1 = if r1.tid() == r2.tid() {
        read_start(r1) < read_start(r2)
    } else {
        r1.tid() < r2.tid()
    };
    if first_is_r1 {
        (r1, r2, true)
    } else {
        (r2, r1, false)
    }
}

fn find_fragment<'a>(
    fragments: &'a FragmentMap,
    chrom: &str,
    record: &Record,
) -> Option<(usize, &'a Fragment)> {
    let list = fragments.get(chrom)?;
    let pos = read_middle(record);
    let idx = list.partition_point(|f| f.start <= pos).checked_sub(1)?;
    let fragment = &list[idx];
    (fragment.start <= pos && pos < fragment.end).then_some((idx, fragment))
}

fn genotype(record: &Record, tag: Option<&[u8; 2]>) -> Option<i64> {
    tag.and_then(|t| aux_integer(record, t))
}

fn aux_tag(tag: Option<&String>) -> Result<Option<[u8; 2]>> {
    match tag {
        None => Ok(None),
        Some(tag) if tag.as_bytes().len() == 2 => Ok(Some([tag.as_bytes()[0], tag.as_bytes()[1]])),
        Some(_) => bail!("BAM auxiliary tag must contain exactly two ASCII characters"),
    }
}

fn classify_pair(
    mut r1: Record,
    mut r2: Record,
    names: &[String],
    fragments: &FragmentMap,
    args: &ClassifyArgs,
    gtag: Option<&[u8; 2]>,
) -> Result<ClassResult> {
    if !r1.is_first_in_template() || !r2.is_last_in_template() {
        bail!(
            "paired BAM must contain adjacent read1/read2 records: {}",
            String::from_utf8_lossy(r1.qname())
        );
    }
    if normalized_qname(r1.qname()) != normalized_qname(r2.qname()) {
        bail!("adjacent paired BAM records have different names");
    }

    let chrom1 = (!r1.is_unmapped()).then(|| &names[r1.tid() as usize]);
    let chrom2 = (!r2.is_unmapped()).then(|| &names[r2.tid() as usize]);
    let frag1 = chrom1.and_then(|c| find_fragment(fragments, c, &r1));
    let frag2 = chrom2.and_then(|c| find_fragment(fragments, c, &r2));

    let mut interaction =
        if !r1.is_unmapped() && !r2.is_unmapped() && frag1.is_some() && frag2.is_some() {
            let (idx1, f1) = frag1.unwrap();
            let (idx2, f2) = frag2.unwrap();
            if r1.tid() == r2.tid() && idx1 == idx2 {
                let (left, right, _) = ordered_pair(&r1, &r2);
                match (left.is_reverse(), right.is_reverse()) {
                    (true, false) => Interaction::SelfCircle,
                    (false, true) => Interaction::DanglingEnd,
                    _ => Interaction::Dump,
                }
            } else if r1.tid() == r2.tid()
                && if f1.start < f2.start {
                    f2.start - f1.end == 0
                } else {
                    f1.start - f2.end == 0
                }
            {
                Interaction::Religation
            } else {
                Interaction::Valid
            }
        } else if r1.is_unmapped() || r2.is_unmapped() {
            Interaction::Singleton
        } else {
            Interaction::Dump
        };

    let distance = fragment_size(
        &r1,
        &r2,
        frag1.map(|x| x.1),
        frag2.map(|x| x.1),
        interaction,
    );
    let cis_distance = if !r1.is_unmapped() && !r2.is_unmapped() && r1.tid() == r2.tid() {
        Some((read_start(&r1) - read_start(&r2)).abs())
    } else {
        None
    };

    if frag1.is_some_and(|(_, f)| f.filtered) || frag2.is_some_and(|(_, f)| f.filtered) {
        interaction = Interaction::Filtered;
    }
    if args
        .min_insert_size
        .is_some_and(|v| distance.is_some_and(|d| d < v))
        || args
            .max_insert_size
            .is_some_and(|v| distance.is_some_and(|d| d > v))
    {
        interaction = Interaction::Filtered;
    }
    if interaction == Interaction::Valid
        && args
            .min_cis_dist
            .is_some_and(|v| cis_distance.is_some_and(|d| d < v))
    {
        interaction = Interaction::Filtered;
    }

    let orientation = (interaction == Interaction::Valid).then(|| {
        let (left, right, _) = ordered_pair(&r1, &r2);
        match (left.is_reverse(), right.is_reverse()) {
            (false, false) => "FF",
            (true, true) => "RR",
            (false, true) => "FR",
            (true, false) => "RF",
        }
    });
    let genotype_pair = if interaction == Interaction::Valid && gtag.is_some() {
        Some((
            genotype(&r1, gtag).unwrap_or(-1),
            genotype(&r2, gtag).unwrap_or(-1),
        ))
    } else {
        None
    };

    let should_write = interaction == Interaction::Valid || args.all_outputs;
    let line = should_write.then(|| {
        format_pair_line(
            &r1,
            &r2,
            names,
            frag1.map(|x| x.1),
            frag2.map(|x| x.1),
            distance,
            gtag,
        )
    });

    let records = if args.interaction_bam && should_write {
        r1.push_aux(b"CT", Aux::String(interaction.code()))?;
        r2.push_aux(b"CT", Aux::String(interaction.code()))?;
        Some((r1, r2))
    } else {
        None
    };

    Ok(ClassResult {
        interaction,
        line,
        orientation,
        genotype_pair,
        records,
    })
}

fn fragment_size(
    r1: &Record,
    r2: &Record,
    frag1: Option<&Fragment>,
    frag2: Option<&Fragment>,
    interaction: Interaction,
) -> Option<i64> {
    if r1.is_unmapped() || r2.is_unmapped() {
        return None;
    }
    let (left, right, first_is_r1) = ordered_pair(r1, r2);
    let (left_frag, right_frag) = if first_is_r1 {
        (frag1, frag2)
    } else {
        (frag2, frag1)
    };
    let left_pos = read_start(left);
    let right_pos = read_start(right);
    match interaction {
        Interaction::DanglingEnd | Interaction::Religation => Some(right_pos - left_pos),
        Interaction::SelfCircle => Some(left_pos - left_frag?.start + right_frag?.end - right_pos),
        Interaction::Valid => {
            let lf = left_frag?;
            let rf = right_frag?;
            let d1 = if left.is_reverse() {
                left_pos - lf.start
            } else {
                lf.end - left_pos
            };
            let d2 = if right.is_reverse() {
                right_pos - rf.start
            } else {
                rf.end - right_pos
            };
            Some(d1 + d2)
        }
        _ => None,
    }
}

fn record_strand(record: &Record) -> char {
    if record.is_reverse() {
        '-'
    } else {
        '+'
    }
}

fn format_pair_line(
    r1: &Record,
    r2: &Record,
    names: &[String],
    frag1: Option<&Fragment>,
    frag2: Option<&Fragment>,
    distance: Option<i64>,
    gtag: Option<&[u8; 2]>,
) -> String {
    if !r1.is_unmapped() && !r2.is_unmapped() {
        let (left, right, first_is_r1) = ordered_pair(r1, r2);
        let (lf, rf) = if first_is_r1 {
            (frag1, frag2)
        } else {
            (frag2, frag1)
        };
        let htag = gtag
            .map(|tag| {
                format!(
                    "{}-{}",
                    genotype(left, Some(tag)).map_or("None".to_string(), |v| v.to_string()),
                    genotype(right, Some(tag)).map_or("None".to_string(), |v| v.to_string())
                )
            })
            .unwrap_or_default();
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            String::from_utf8_lossy(left.qname()),
            names[left.tid() as usize],
            read_start(left) + 1,
            record_strand(left),
            names[right.tid() as usize],
            read_start(right) + 1,
            record_strand(right),
            distance.map_or("None".to_string(), |v| v.to_string()),
            lf.map_or("None", |f| f.name.as_ref()),
            rf.map_or("None", |f| f.name.as_ref()),
            left.mapq(),
            right.mapq(),
            htag
        )
    } else if !r1.is_unmapped() {
        format!(
            "{}\t{}\t{}\t{}\t*\t*\t*\t*\t{}\t*\t{}\t*\n",
            String::from_utf8_lossy(r1.qname()),
            names[r1.tid() as usize],
            read_start(r1) + 1,
            record_strand(r1),
            frag1.map_or("None", |f| f.name.as_ref()),
            r1.mapq()
        )
    } else {
        format!(
            "{}\t*\t*\t*\t{}\t{}\t{}\t*\t*\t{}\t*\t{}\n",
            String::from_utf8_lossy(r2.qname()),
            names[r2.tid() as usize],
            read_start(r2) + 1,
            record_strand(r2),
            frag2.map_or("None", |f| f.name.as_ref()),
            r2.mapq()
        )
    }
}

fn base_name(path: &Path) -> Result<String> {
    let name = path
        .file_name()
        .context("mapped reads path has no file name")?
        .to_string_lossy();
    Ok(name
        .strip_suffix(".bam")
        .or_else(|| name.strip_suffix(".sam"))
        .unwrap_or(&name)
        .to_string())
}

fn open_pair_outputs(
    output_dir: &Path,
    base: &str,
    all: bool,
) -> Result<HashMap<&'static str, BufWriter<File>>> {
    let mut outputs = HashMap::new();
    let classes = if all {
        vec![
            Interaction::Valid,
            Interaction::DanglingEnd,
            Interaction::Religation,
            Interaction::SelfCircle,
            Interaction::Singleton,
            Interaction::Filtered,
            Interaction::Dump,
        ]
    } else {
        vec![Interaction::Valid]
    };
    for class in classes {
        let path = output_dir.join(format!("{base}.{}", class.suffix()));
        outputs.insert(class.code(), BufWriter::new(File::create(path)?));
    }
    Ok(outputs)
}

fn update_class_stats(stats: &mut ClassStats, result: &ClassResult) {
    match result.interaction {
        Interaction::Valid => {
            stats.valid += 1;
            match result.orientation {
                Some("FF") => stats.ff += 1,
                Some("RR") => stats.rr += 1,
                Some("RF") => stats.rf += 1,
                Some("FR") => stats.fr += 1,
                _ => {}
            }
            if let Some((a, b)) = result.genotype_pair {
                match (a, b) {
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
        Interaction::DanglingEnd => stats.dangling += 1,
        Interaction::Religation => stats.religation += 1,
        Interaction::SelfCircle => stats.self_circle += 1,
        Interaction::Singleton => stats.singleton += 1,
        Interaction::Filtered => stats.filtered += 1,
        Interaction::Dump => stats.dumped += 1,
    }
}

fn run_classify(args: ClassifyArgs) -> Result<()> {
    validate_parallel(&args.parallel)?;
    configure_pool(args.parallel.threads)?;
    fs::create_dir_all(&args.output_dir)?;
    let fragments = load_fragments(&args)?;
    let gtag = aux_tag(args.genotype_tag.as_ref())?;

    let mut reader = bam::Reader::from_path(&args.mapped_reads)
        .with_context(|| format!("cannot open {}", args.mapped_reads.display()))?;
    reader.set_threads(io_threads(args.parallel.threads))?;
    let names = reference_names(reader.header());
    let header = bam::Header::from_template(reader.header());
    let base = base_name(&args.mapped_reads)?;
    let mut outputs = open_pair_outputs(&args.output_dir, &base, args.all_outputs)?;
    let mut bam_writer = if args.interaction_bam {
        let path = args.output_dir.join(format!("{base}_interaction.bam"));
        let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam)?;
        writer.set_threads(io_threads(args.parallel.threads))?;
        Some(writer)
    } else {
        None
    };

    let mut records = reader.records();
    let mut stats = ClassStats::default();
    let mut processed = 0_u64;
    loop {
        let mut batch = Vec::with_capacity(args.parallel.batch_size);
        for _ in 0..args.parallel.batch_size {
            match records.next() {
                None => break,
                Some(first) => {
                    let first = first?;
                    let second = records
                        .next()
                        .context("paired BAM ends with an incomplete read pair")??;
                    batch.push((first, second));
                }
            }
        }
        if batch.is_empty() {
            break;
        }
        let results: Result<Vec<_>> = batch
            .into_par_iter()
            .map(|(r1, r2)| classify_pair(r1, r2, &names, &fragments, &args, gtag.as_ref()))
            .collect();
        for result in results? {
            update_class_stats(&mut stats, &result);
            if let Some(line) = result.line {
                if let Some(out) = outputs.get_mut(result.interaction.code()) {
                    out.write_all(line.as_bytes())?;
                }
            }
            if let (Some(writer), Some((r1, r2))) = (&mut bam_writer, result.records) {
                writer.write(&r1)?;
                writer.write(&r2)?;
            }
            processed += 1;
        }
        if args.verbose && processed % 100_000 < args.parallel.batch_size as u64 {
            eprintln!("## {processed} read pairs");
        }
    }
    drop(outputs);
    drop(bam_writer);
    write_class_stats(
        &args.output_dir.join(format!("{base}.RSstat")),
        &stats,
        gtag.is_some(),
    )?;
    Ok(())
}

fn write_class_stats(path: &Path, s: &ClassStats, allele_specific: bool) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "## Hi-C processing")?;
    writeln!(out, "Valid_interaction_pairs\t{}", s.valid)?;
    writeln!(out, "Valid_interaction_pairs_FF\t{}", s.ff)?;
    writeln!(out, "Valid_interaction_pairs_RR\t{}", s.rr)?;
    writeln!(out, "Valid_interaction_pairs_RF\t{}", s.rf)?;
    writeln!(out, "Valid_interaction_pairs_FR\t{}", s.fr)?;
    writeln!(out, "Dangling_end_pairs\t{}", s.dangling)?;
    writeln!(out, "Religation_pairs\t{}", s.religation)?;
    writeln!(out, "Self_Cycle_pairs\t{}", s.self_circle)?;
    writeln!(out, "Single-end_pairs\t{}", s.singleton)?;
    writeln!(out, "Filtered_pairs\t{}", s.filtered)?;
    writeln!(out, "Dumped_pairs\t{}", s.dumped)?;
    if allele_specific {
        writeln!(out, "## ======================================")?;
        writeln!(out, "## Allele specific information")?;
        writeln!(out, "Valid_pairs_from_ref_genome_(1-1)\t{}", s.g1g1)?;
        writeln!(
            out,
            "Valid_pairs_from_ref_genome_with_one_unassigned_mate_(0-1/1-0)\t{}",
            s.ug1 + s.g1u
        )?;
        writeln!(out, "Valid_pairs_from_alt_genome_(2-2)\t{}", s.g2g2)?;
        writeln!(
            out,
            "Valid_pairs_from_alt_genome_with_one_unassigned_mate_(0-2/2-0)\t{}",
            s.ug2 + s.g2u
        )?;
        writeln!(
            out,
            "Valid_pairs_from_alt_and_ref_genome_(1-2/2-1)\t{}",
            s.g1g2 + s.g2g1
        )?;
        writeln!(
            out,
            "Valid_pairs_with_both_unassigned_mated_(0-0)\t{}",
            s.uu
        )?;
        writeln!(
            out,
            "Valid_pairs_with_at_least_one_conflicting_mate_(3-)\t{}",
            s.conflict
        )?;
    }
    Ok(())
}
