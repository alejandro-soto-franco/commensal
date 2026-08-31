//! Dereplicate a FASTQ and flag its bimeras, reporting what each step cost.

use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use commensal::amplicon::{BimeraOptions, dereplicate, flag_bimeras, read_fastq};

pub fn run(args: Vec<String>) {
    if args.is_empty() {
        eprintln!("usage: commensal amplicon <reads.fastq[.gz]> [more ...] [--out table.tsv]");
        std::process::exit(2);
    }
    let mut inputs: Vec<PathBuf> = Vec::new();
    let mut out: Option<PathBuf> = None;
    let mut it = args.into_iter();
    while let Some(a) = it.next() {
        if a == "--out" {
            out = it.next().map(PathBuf::from);
        } else {
            inputs.push(PathBuf::from(a));
        }
    }

    let t = Instant::now();
    let mut reads: Vec<Vec<u8>> = Vec::new();
    for path in &inputs {
        let file = File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let parsed = if path.extension().is_some_and(|e| e == "gz") {
            read_fastq(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
        } else {
            read_fastq(BufReader::new(file))
        };
        reads.extend(parsed.expect("fastq"));
    }
    let read_s = t.elapsed().as_secs_f64();

    let t = Instant::now();
    let uniques = dereplicate(&reads);
    let derep_s = t.elapsed().as_secs_f64();

    let t = Instant::now();
    let flags = flag_bimeras(&uniques, &BimeraOptions::default());
    let bimera_s = t.elapsed().as_secs_f64();

    let chimeric: usize = flags.iter().filter(|f| **f).count();
    let kept_reads: u32 = uniques
        .iter()
        .zip(&flags)
        .filter(|(_, f)| !**f)
        .map(|(u, _)| u.abundance)
        .sum();

    println!("reads          {}", reads.len());
    println!("unique         {}", uniques.len());
    println!("bimeras        {chimeric}");
    println!("reads kept     {kept_reads}");
    println!("read   {read_s:.3} s");
    println!("derep  {derep_s:.3} s");
    println!("bimera {bimera_s:.3} s");
    println!("total  {:.3} s", read_s + derep_s + bimera_s);

    if let Some(path) = out {
        let mut w = BufWriter::new(File::create(&path).expect("output"));
        writeln!(w, "sequence\tabundance\tbimera").expect("write");
        for (u, f) in uniques.iter().zip(&flags) {
            writeln!(
                w,
                "{}\t{}\t{}",
                String::from_utf8_lossy(&u.sequence),
                u.abundance,
                if *f { 1 } else { 0 }
            )
            .expect("write");
        }
        println!("wrote {}", path.display());
    }
}
