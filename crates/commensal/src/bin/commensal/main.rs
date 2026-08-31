//! Dispatch to one of the crate's three tables: amplicon, consortium, receptor.

mod amplicon;
mod consortium;
mod receptor;

fn usage() -> ! {
    eprintln!("usage: commensal <amplicon|consortium|receptor> [args ...]");
    eprintln!();
    eprintln!("  amplicon <reads.fastq[.gz]> [more ...] [--out table.tsv]");
    eprintln!("      dereplicate reads and flag bimeras");
    eprintln!("  consortium");
    eprintln!("      colonisation and butyrate production for the reference synbiotic");
    eprintln!("  receptor");
    eprintln!("      binding against activation for flagellin at TLR5");
    std::process::exit(2);
}

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        usage();
    }
    let subcommand = args.remove(0);
    match subcommand.as_str() {
        "amplicon" => amplicon::run(args),
        "consortium" => consortium::run(),
        "receptor" => receptor::run(),
        _ => usage(),
    }
}
