use check_build::{detect_build, Reference, Verifier};
use clap::Parser;
use std::path::Path;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Verifies a plain-text VCF against hg19 and hg38 reference genomes."
)]
struct Args {
    /// Path to a non-gzipped VCF file
    vcf: String,

    /// Just detect the build (simple output)
    #[arg(short, long)]
    detect: bool,

    /// Suppress progress bars
    #[arg(short, long)]
    quiet: bool,

    /// Suppress mismatch details (summary only)
    #[arg(short = 's', long)]
    summary_only: bool,

    /// Only verify against hg19
    #[arg(long, conflicts_with = "hg38_only")]
    hg19_only: bool,

    /// Only verify against hg38
    #[arg(long, conflicts_with = "hg19_only")]
    hg38_only: bool,

    /// Custom path to hg19.fa reference
    #[arg(long)]
    hg19_path: Option<String>,

    /// Custom path to hg38.fa reference
    #[arg(long)]
    hg38_path: Option<String>,
}

fn main() {
    let args = Args::parse();

    if !args.vcf.ends_with(".vcf") {
        eprintln!("ERROR: Only accepts a .vcf file (not .gz).");
        std::process::exit(1);
    }
    if !Path::new(&args.vcf).exists() {
        eprintln!("ERROR: VCF file '{}' not found.", args.vcf);
        std::process::exit(1);
    }

    let mut verifier = Verifier::new(&args.vcf);

    if args.quiet {
        verifier = verifier.quiet();
    }
    if args.summary_only {
        verifier = verifier.silent();
    }
    if let Some(path) = args.hg19_path {
        verifier = verifier.hg19_path(path);
    }
    if let Some(path) = args.hg38_path.clone() {
        verifier = verifier.hg38_path(path);
    }

    // Simple detect mode - just print match rates
    if args.detect {
        match detect_build(&args.vcf) {
            Ok(result) => {
                println!("{}", result);
                if let Some(build) = result.better_match() {
                    if result.is_clear_match() {
                        println!("Best match: {:?}", build);
                    }
                }
            }
            Err(e) => {
                eprintln!("Error: {}", e);
                std::process::exit(1);
            }
        }
        return;
    }
    if args.hg19_only {
        println!("Verifying against hg19 only...");
        match verifier.verify_single(Reference::Hg19) {
            Ok((lines, mismatches)) => {
                println!("\n===================");
                println!("hg19: {} lines, {} mismatches", lines, mismatches);
                if mismatches == 0 {
                    println!("All checks passed.");
                } else {
                    std::process::exit(1);
                }
            }
            Err(e) => {
                eprintln!("ERROR: {}", e);
                std::process::exit(1);
            }
        }
        return;
    }

    if args.hg38_only {
        println!("Verifying against hg38 only...");
        match verifier.verify_single(Reference::Hg38) {
            Ok((lines, mismatches)) => {
                println!("\n===================");
                println!("hg38: {} lines, {} mismatches", lines, mismatches);
                if mismatches == 0 {
                    println!("All checks passed.");
                } else {
                    std::process::exit(1);
                }
            }
            Err(e) => {
                eprintln!("ERROR: {}", e);
                std::process::exit(1);
            }
        }
        return;
    }

    // Default: verify both
    println!("Verifying against hg19 and hg38...");
    match verifier.verify_both() {
        Ok(result) => {
            println!("\n===================");
            println!("Verification Summary:");
            println!(
                "  hg19: {} lines, {} mismatches ({:.1}% match)",
                result.hg19_lines,
                result.hg19_mismatches,
                result.match_rate(Reference::Hg19)
            );
            println!(
                "  hg38: {} lines, {} mismatches ({:.1}% match)",
                result.hg38_lines,
                result.hg38_mismatches,
                result.match_rate(Reference::Hg38)
            );

            if let Some(build) = result.likely_build() {
                println!("\nLikely reference build: {:?}", build);
            }

            if result.all_passed() {
                println!("All checks passed.");
            } else {
                std::process::exit(1);
            }
        }
        Err(e) => {
            eprintln!("ERROR: {}", e);
            std::process::exit(1);
        }
    }
}
