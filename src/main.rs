use clap::Parser;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use reqwest::blocking::Client;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path;
use tempfile::NamedTempFile;

/// URLs for uncompressed FASTA files
const HG19_URL: &str = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa";
const HG38_URL: &str = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa";

/// Names for the downloaded FASTA files
const HG19_FA: &str = "hg19.fa";
const HG38_FA: &str = "hg38.fa";

/// Holds the start offset and length of each chromosome in our flattened file.
#[derive(Debug)]
struct ChromOffset {
    start: u64,
    length: u64,
}

/// Mapping: chromosome name -> ChromOffset
type ChromMap = HashMap<String, ChromOffset>;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Verifies a plain-text .vcf (no gz) against both hg19 and hg38 with minimal memory usage."
)]
struct Args {
    /// Path to a non-gzipped VCF file.
    vcf: String,
}

fn main() {
    let args = Args::parse();

    // 1) Validate the input is a plain-text .vcf file.
    if !args.vcf.ends_with(".vcf") {
        eprintln!("ERROR: Only accepts a .vcf file (not .gz).");
        std::process::exit(1);
    }
    if !Path::new(&args.vcf).exists() {
        eprintln!("ERROR: VCF file '{}' not found.", args.vcf);
        std::process::exit(1);
    }

    // 2) check FASTA files are downloaded.
    check_file(HG19_FA, HG19_URL);
    check_file(HG38_FA, HG38_URL);

    // 3) Flatten each FASTA to a temporary file and build an offset map.
    println!("Flattening and indexing {}...", HG19_FA);
    let (hg19_flat, hg19_map) = match flatten_fasta(HG19_FA, "hg19") {
        Ok(val) => val,
        Err(e) => {
            eprintln!("ERROR flattening {}: {}", HG19_FA, e);
            std::process::exit(1);
        }
    };

    println!("Flattening and indexing {}...", HG38_FA);
    let (hg38_flat, hg38_map) = match flatten_fasta(HG38_FA, "hg38") {
        Ok(val) => val,
        Err(e) => {
            eprintln!("ERROR flattening {}: {}", HG38_FA, e);
            std::process::exit(1);
        }
    };

    // 4) Count the number of non-header records in the VCF for progress display.
    println!("Counting records in VCF {}...", args.vcf);
    let total_vcf_records = match count_vcf_records(&args.vcf) {
        Ok(count) => count,
        Err(e) => {
            eprintln!("ERROR counting VCF lines: {}", e);
            std::process::exit(1);
        }
    };

    // 5) Verify the VCF against each reference.
    println!("\n===================");
    println!("Verifying against hg19:");
    let (lines_hg19, mismatches_hg19) =
        verify_vcf(&hg19_flat, &hg19_map, &args.vcf, total_vcf_records, "hg19");
    println!("Finished checking hg19.");
    println!("Lines processed: {}", lines_hg19);
    println!("Mismatches: {}", mismatches_hg19);

    println!("\n===================");
    println!("Verifying against hg38:");
    let (lines_hg38, mismatches_hg38) =
        verify_vcf(&hg38_flat, &hg38_map, &args.vcf, total_vcf_records, "hg38");
    println!("Finished checking hg38.");
    println!("Lines processed: {}", lines_hg38);
    println!("Mismatches: {}", mismatches_hg38);

    // 6) Final summary and comparison.
    println!("\n===================");
    println!("Verification Summary:");
    println!("  - hg19 => {} lines, {} mismatches", lines_hg19, mismatches_hg19);
    println!("  - hg38 => {} lines, {} mismatches", lines_hg38, mismatches_hg38);
    if mismatches_hg19 + mismatches_hg38 == 0 {
        println!("All checks passed (no mismatches on either reference).");
    } else {
        println!("Some mismatches occurred. See details above.");
    }

    // Temporary files are cleaned up automatically.
    if mismatches_hg19 + mismatches_hg38 > 0 {
        std::process::exit(1);
    }
}

///  that `dest_path` exists locally; if not, downloads it from `url` with a progress bar.
fn check_file(dest_path: &str, url: &str) {
    if Path::new(dest_path).exists() {
        println!("{} already exists, skipping download.", dest_path);
    } else {
        println!("{} not found; downloading from {}...", dest_path, url);
        if let Err(e) = download_file(url, dest_path) {
            eprintln!("ERROR downloading {}: {}", dest_path, e);
            std::process::exit(1);
        }
        println!("Download complete: {}", dest_path);
    }
}

/// Downloads `url` to `dest_path` using a blocking reqwest call with an indicatif progress bar.
fn download_file(url: &str, dest_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let client = Client::new();
    let mut response = client.get(url).send()?;
    if !response.status().is_success() {
        return Err(format!("HTTP error, status code: {}", response.status()).into());
    }
    let total_size = response
        .content_length()
        .ok_or_else(|| "Could not get Content-Length from server.".to_string())?;
    let pb = ProgressBar::new(total_size);
    pb.set_draw_target(ProgressDrawTarget::stderr());
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{msg} [{bar:40.cyan-blue}] {bytes} of {total_bytes} bytes ({eta})")
            .progress_chars("=>-")
    );
    pb.set_message(format!("Downloading {}", dest_path));
    let mut out_file = File::create(dest_path)?;
    let mut downloaded = 0u64;
    let mut buffer = [0u8; 8192];
    loop {
        let bytes_read = response.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        out_file.write_all(&buffer[..bytes_read])?;
        downloaded += bytes_read as u64;
        pb.set_position(downloaded);
    }
    pb.finish_with_message(format!("Downloaded {}", dest_path));
    Ok(())
}

/// Reads a FASTA file line-by-line, removes newlines, and writes the raw bases to a temporary flattened file.
/// It also builds an in-memory map of chromosome offsets. Sequences are stored in uppercase.
fn flatten_fasta(
    fasta_path: &str,
    tag: &str,
) -> Result<(NamedTempFile, ChromMap), Box<dyn std::error::Error>> {
    let fasta_in = File::open(fasta_path)?;
    let mut reader = BufReader::new(fasta_in);
    let mut flat_file = NamedTempFile::new()?;
    let mut chrom_map: ChromMap = HashMap::new();
    let mut current_chrom = String::new();
    let mut current_chrom_start = 0u64;
    let mut current_chrom_len = 0u64;
    let mut line_buf = String::new();
    let mut total_bytes_written = 0u64;
    let fasta_size = std::fs::metadata(fasta_path)
        .ok()
        .map(|m| m.len())
        .unwrap_or(0);
    let pb = ProgressBar::new(fasta_size);
    pb.set_draw_target(ProgressDrawTarget::stderr());
    pb.set_style(
        ProgressStyle::default_bar()
            .template(&format!(
                "Flattening {} {{spinner}} [{{bar:40.cyan-blue}}] {{bytes}} of {{total_bytes}} bytes ({{eta}})",
                tag
            )),
    );
    pb.set_message(format!("Flattening {}", fasta_path));
    loop {
        let bytes_read = reader.read_line(&mut line_buf)?;
        if bytes_read == 0 {
            if !current_chrom.is_empty() {
                chrom_map.insert(
                    current_chrom.clone(),
                    ChromOffset {
                        start: current_chrom_start,
                        length: current_chrom_len,
                    },
                );
            }
            break;
        }
        pb.inc(bytes_read as u64);
        if line_buf.starts_with('>') {
            if !current_chrom.is_empty() {
                chrom_map.insert(
                    current_chrom.clone(),
                    ChromOffset {
                        start: current_chrom_start,
                        length: current_chrom_len,
                    },
                );
            }
            current_chrom.clear();
            current_chrom_start = total_bytes_written;
            current_chrom_len = 0;
            let header = line_buf
                .trim_start_matches('>')
                .split_whitespace()
                .next()
                .unwrap_or("");
            current_chrom.push_str(header);
        } else {
            let seq_line = line_buf.trim_end();
            let seq_line_upper = seq_line.to_uppercase();
            flat_file.write_all(seq_line_upper.as_bytes())?;
            let len_written = seq_line_upper.len() as u64;
            total_bytes_written += len_written;
            current_chrom_len += len_written;
        }
        line_buf.clear();
    }
    pb.finish_with_message(format!("Flattened {} into temporary file", fasta_path));
    Ok((flat_file, chrom_map))
}

/// Quickly counts the number of data lines (non-header) in the VCF.
fn count_vcf_records(vcf_path: &str) -> Result<u64, Box<dyn std::error::Error>> {
    let file = File::open(vcf_path)?;
    let reader = BufReader::new(file);
    let mut count = 0u64;
    for line_res in reader.lines() {
        let line = line_res?;
        if !line.starts_with('#') && !line.trim().is_empty() {
            count += 1;
        }
    }
    Ok(count)
}

/// Opens the flattened FASTA file, streams through the VCF, and checks that each REF allele matches
/// the corresponding substring from the reference (using the chromosome offset map). A progress bar is shown.
fn verify_vcf(
    flat_fasta: &NamedTempFile,
    chrom_map: &ChromMap,
    vcf_path: &str,
    total_lines: u64,
    label: &str,
) -> (u64, u64) {
    let mut flat_reader = File::open(flat_fasta.path())
        .unwrap_or_else(|_| panic!("ERROR: Could not open flattened file for {}.", label));
    let pb = ProgressBar::new(total_lines);
    pb.set_draw_target(ProgressDrawTarget::stderr());
    pb.set_style(
        ProgressStyle::default_bar()
            .template(&format!(
                "{}: Checking VCF records {{spinner}} [{{bar:40.cyan-blue}}] {{pos}} of {{len}} ({{eta}})",
                label
            )),
    );
    pb.set_message(format!("Verifying {} vs {}", vcf_path, label));
    let file = File::open(vcf_path)
        .unwrap_or_else(|_| panic!("ERROR: Could not open VCF file: {}", vcf_path));
    let reader = BufReader::new(file);
    let mut lineno = 0u64;
    let mut mismatches = 0u64;
    for line_res in reader.lines() {
        let line = match line_res {
            Ok(x) => x,
            Err(e) => {
                eprintln!("Error reading VCF line {}: {}", lineno, e);
                break;
            }
        };
        lineno += 1;
        if line.starts_with('#') {
            continue;
        }
        pb.inc(1);
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 4 {
            eprintln!("WARNING: Malformed line {}: {}", lineno, line);
            continue;
        }
        let chrom = cols[0];
        let pos_str = cols[1];
        let ref_allele = cols[3];
        let pos_1based = match pos_str.parse::<u64>() {
            Ok(p) => p,
            Err(_) => {
                eprintln!("WARNING: Invalid position on line {}: {}", lineno, line);
                continue;
            }
        };
        let pos_0based = pos_1based.saturating_sub(1);
        let (start, length) = if let Some(offset) = chrom_map.get(chrom) {
            (offset.start, offset.length)
        } else {
            if chrom.starts_with("chr") {
                let alt = chrom.trim_start_matches("chr");
                if let Some(off) = chrom_map.get(alt) {
                    (off.start, off.length)
                } else {
                    eprintln!("WARNING: Chrom {} not found in FASTA (line {})", chrom, lineno);
                    continue;
                }
            } else {
                let alt = format!("chr{}", chrom);
                if let Some(off) = chrom_map.get(&alt) {
                    (off.start, off.length)
                } else {
                    eprintln!("WARNING: Chrom {} not found in FASTA (line {})", chrom, lineno);
                    continue;
                }
            }
        };
        let end_pos = pos_0based + ref_allele.len() as u64;
        if end_pos > length {
            eprintln!(
                "WARNING: Out-of-bounds for chrom {} at line {} (pos={}, ref_len={}, chrom_len={})",
                chrom, lineno, pos_1based, ref_allele.len(), length
            );
            continue;
        }
        let absolute_offset = start + pos_0based;
        if let Err(e) = flat_reader.seek(SeekFrom::Start(absolute_offset)) {
            eprintln!("ERROR seeking in flattened FASTA at line {}: {}", lineno, e);
            break;
        }
        let mut buffer = vec![0u8; ref_allele.len()];
        if let Err(e) = flat_reader.read_exact(&mut buffer) {
            eprintln!("ERROR reading substring at line {}: {}", lineno, e);
            break;
        }
        let extracted_str = String::from_utf8_lossy(&buffer);
        if !equals_ignore_case(&extracted_str, ref_allele) {
            eprintln!(
                "Mismatch line {}:\n  Chrom: {} Pos: {} REF='{}'\n  FASTA => '{}'\n",
                lineno, chrom, pos_str, ref_allele, extracted_str
            );
            mismatches += 1;
        }
    }
    pb.finish_with_message(format!("Done checking {}", label));
    (lineno, mismatches)
}

fn equals_ignore_case(a: &str, b: &str) -> bool {
    a.len() == b.len() && a.eq_ignore_ascii_case(b)
}
