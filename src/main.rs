use clap::Parser;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use reqwest::blocking::Client;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use tempfile::tempdir;

/// URLs for uncompressed FASTA files
const HG19_URL: &str = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa";
const HG38_URL: &str = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa";

/// Local filenames for references
const HG19_FA: &str = "hg19.fa";
const HG38_FA: &str = "hg38.fa";

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Verifies a plain-text .vcf (no gz) against both hg19 and hg38 in a streaming manner without large memory usage."
)]
struct Args {
    /// Path to a non-gzipped VCF file.
    vcf: String,
}

/// This struct tracks the current contig name/sequence and global verification counters.
/// By using a struct, we avoid multiple mutable borrows to the same data.
struct ContigState {
    current_contig: String,
    current_seq: Vec<u8>,
    lines_checked: u64,
    mismatches: u64,
}

impl ContigState {
    fn new() -> Self {
        ContigState {
            current_contig: String::new(),
            current_seq: Vec::new(),
            lines_checked: 0,
            mismatches: 0,
        }
    }

    /// Flushes the current contig by verifying it against the splitted VCF lines, then clears name/sequence.
    fn flush_contig(&mut self, contig_file_map: &HashMap<String, PathBuf>) {
        if self.current_contig.is_empty() {
            return;
        }
        let cnt = verify_contig(
            &self.current_contig,
            &self.current_seq,
            contig_file_map,
            &mut self.mismatches,
        );
        self.lines_checked += cnt;

        // Reset
        self.current_contig.clear();
        self.current_seq.clear();
    }
}

fn main() {
    let args = Args::parse();

    // 1) Validate VCF input
    if !args.vcf.ends_with(".vcf") {
        eprintln!("ERROR: Only accepts a .vcf file (not .gz).");
        std::process::exit(1);
    }
    if !Path::new(&args.vcf).exists() {
        eprintln!("ERROR: VCF file '{}' not found.", args.vcf);
        std::process::exit(1);
    }

    // 2) reference FASTA is downloaded
    check_file(HG19_FA, HG19_URL);
    check_file(HG38_FA, HG38_URL);

    // 3) Split the VCF into per-contig files (in a temp directory),
    //    so we can retrieve lines for each contig on demand with minimal memory usage.
    println!("Splitting VCF by contig into temporary files...");
    let split_dir = match tempdir() {
        Ok(d) => d,
        Err(e) => {
            eprintln!("ERROR creating temp dir for splitted VCF: {}", e);
            std::process::exit(1);
        }
    };
    let contig_file_map = match split_vcf_by_contig(&args.vcf, &split_dir) {
        Ok(map) => map,
        Err(e) => {
            eprintln!("ERROR splitting VCF: {}", e);
            std::process::exit(1);
        }
    };
    println!(
        "VCF split complete. Found {} contigs in the VCF.",
        contig_file_map.len()
    );

    // 4) Verify the VCF against each reference (hg19 & hg38) using streaming.
    println!("\n===================");
    println!("Verifying against hg19 (streaming)...");
    let (lines_hg19, mismatches_hg19) = verify_reference_streaming(HG19_FA, &contig_file_map);
    println!("Finished checking hg19.");
    println!("Lines processed: {}", lines_hg19);
    println!("Mismatches: {}", mismatches_hg19);

    println!("\n===================");
    println!("Verifying against hg38 (streaming)...");
    let (lines_hg38, mismatches_hg38) = verify_reference_streaming(HG38_FA, &contig_file_map);
    println!("Finished checking hg38.");
    println!("Lines processed: {}", lines_hg38);
    println!("Mismatches: {}", mismatches_hg38);

    // 5) Final summary
    println!("\n===================");
    println!("Verification Summary:");
    println!(
        "  - hg19 => {} lines, {} mismatches",
        lines_hg19, mismatches_hg19
    );
    println!(
        "  - hg38 => {} lines, {} mismatches",
        lines_hg38, mismatches_hg38
    );
    if mismatches_hg19 + mismatches_hg38 == 0 {
        println!("All checks passed (no mismatches on either reference).");
    } else {
        println!("Some mismatches occurred. See details above.");
        std::process::exit(1);
    }

    // Temp files auto-cleaned when 'split_dir' goes out of scope.
}

/// `dest_path` exists locally; if not, download from `url` with a progress bar.
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

/// Download `url` to `dest_path` using reqwest + indicatif progress bar.
/// We read in larger chunks, only updating the progress bar after each chunk.
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
            .progress_chars("=>-"),
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

/// Splits all VCF lines (non-header) by contig into separate temporary files,
/// stored in `split_dir`. Returns a map contig -> PathBuf (temp file location).
///
/// Memory usage is minimal: we stream line by line. If the VCF is large,
/// we never hold more than one line in memory at a time.
fn split_vcf_by_contig(
    vcf_path: &str,
    split_dir: &tempfile::TempDir,
) -> Result<HashMap<String, PathBuf>, Box<dyn std::error::Error>> {
    let f = File::open(vcf_path)?;
    let reader = BufReader::new(f);

    let mut handles: HashMap<String, File> = HashMap::new();
    let mut file_paths: HashMap<String, PathBuf> = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res?;
        if line.starts_with('#') || line.trim().is_empty() {
            // skip headers / empties
            continue;
        }
        // columns: contig, pos, ...
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.is_empty() {
            continue;
        }
        let contig = cols[0].to_string();

        // If not already have a file for this contig, open/create one
        if !handles.contains_key(&contig) {
            let contig_file = split_dir.path().join(format!("{}.tmp", contig));
            let f = OpenOptions::new()
                .create(true)
                .append(true)
                .write(true)
                .open(&contig_file)?;
            handles.insert(contig.clone(), f);
            file_paths.insert(contig.clone(), contig_file);
        }
        // Write this line to the contig's file
        if let Some(h) = handles.get_mut(&contig) {
            writeln!(h, "{}", line)?;
        }
    }

    Ok(file_paths)
}

/// Streams the given reference FASTA contig by contig, verifying the splitted VCF lines for each.
fn verify_reference_streaming(
    fasta_path: &str,
    contig_file_map: &HashMap<String, PathBuf>,
) -> (u64, u64) {
    // We'll parse the FASTA in chunked reads, building lines, then detect contigs with '>'.
    // We'll keep track in a struct `ContigState`.

    let meta = match std::fs::metadata(fasta_path) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("ERROR reading metadata for {}: {}", fasta_path, e);
            return (0, 0);
        }
    };
    let total_len = meta.len();

    let pb = ProgressBar::new(total_len);
    pb.set_draw_target(ProgressDrawTarget::stderr());
    pb.set_style(
        ProgressStyle::default_bar()
            .template(&format!(
                "Streaming reference {} {{spinner}} [{{bar:40.cyan-blue}}] {{bytes}} of {{total_bytes}} ({{eta}})",
                fasta_path
            ))
            .progress_chars("=>-"),
    );

    // Open FASTA
    let f = match File::open(fasta_path) {
        Ok(x) => x,
        Err(e) => {
            eprintln!("ERROR opening {}: {}", fasta_path, e);
            return (0, 0);
        }
    };
    let mut reader = BufReader::new(f);

    let mut chunk = [0u8; 65536]; // 64 KB
    let mut read_bytes = 0u64;

    let mut line_buf = Vec::<u8>::new();
    let mut state = ContigState::new();

    // Read in chunks, gather lines ourselves
    loop {
        let n = match reader.read(&mut chunk) {
            Ok(x) => x,
            Err(e) => {
                eprintln!("ERROR reading {}: {}", fasta_path, e);
                break;
            }
        };
        if n == 0 {
            // EOF
            break;
        }
        read_bytes += n as u64;
        pb.set_position(read_bytes);

        let slice = &chunk[..n];
        for &b in slice {
            if b == b'\n' || b == b'\r' {
                // We've hit a line boundary
                process_fasta_line(&mut state, &mut line_buf, contig_file_map);
                line_buf.clear();
            } else {
                line_buf.push(b);
            }
        }
    }
    // Process leftover line
    if !line_buf.is_empty() {
        process_fasta_line(&mut state, &mut line_buf, contig_file_map);
        line_buf.clear();
    }
    // Flush final contig
    state.flush_contig(contig_file_map);

    pb.finish_with_message(format!("Done streaming {}.", fasta_path));
    (state.lines_checked, state.mismatches)
}

/// Process one line from the FASTA. If it starts with '>', it's a new contig. Otherwise, sequence data.
fn process_fasta_line(
    state: &mut ContigState,
    line_buf: &mut Vec<u8>,
    contig_file_map: &HashMap<String, PathBuf>,
) {
    if line_buf.is_empty() {
        return;
    }
    if line_buf[0] == b'>' {
        // We have a new contig
        // Flush the old one
        state.flush_contig(contig_file_map);

        // Parse contig name
        let line_str = String::from_utf8_lossy(&line_buf[1..]).to_string();
        let header = line_str.split_whitespace().next().unwrap_or("").to_string();
        state.current_contig = header;
    } else {
        // Sequence line => uppercase + append
        for b in line_buf.iter_mut() {
            b.make_ascii_uppercase();
        }
        state.current_seq.extend_from_slice(line_buf);
    }
}

/// Verifies splitted VCF lines for `contig` against the in-memory `seq`.
/// Updates `mismatch_count` in place, returns how many lines processed for that contig.
fn verify_contig(
    contig: &str,
    seq: &[u8],
    contig_file_map: &HashMap<String, PathBuf>,
    mismatch_count: &mut u64,
) -> u64 {
    // The VCF might reference contig as "chrXYZ", "XYZ", etc.
    let candidates = [
        contig.to_string(),
        contig.trim_start_matches("chr").to_string(),
        format!("chr{}", contig),
    ];
    let mut lines_checked = 0;

    // Find splitted file path
    let mut file_path_opt = None;
    for c in &candidates {
        if let Some(p) = contig_file_map.get(c) {
            file_path_opt = Some(p);
            break;
        }
    }

    let contig_file = match file_path_opt {
        Some(p) => p,
        None => {
            // no lines for this contig
            return 0;
        }
    };

    let f = match File::open(contig_file) {
        Ok(x) => x,
        Err(e) => {
            eprintln!(
                "WARNING: Could not open splitted VCF file for contig {}: {}",
                contig, e
            );
            return 0;
        }
    };
    let reader = BufReader::new(f);

    for line_res in reader.lines() {
        let line = match line_res {
            Ok(x) => x,
            Err(e) => {
                eprintln!(
                    "WARNING: Error reading splitted line for contig {}: {}",
                    contig, e
                );
                break;
            }
        };
        lines_checked += 1;

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 4 {
            eprintln!("WARNING: Malformed line in splitted file: {}", line);
            continue;
        }
        let pos_str = cols[1];
        let ref_allele = cols[3];

        let pos_1based = match pos_str.parse::<u64>() {
            Ok(p) => p,
            Err(_) => {
                eprintln!("WARNING: Invalid position in splitted file line: {}", line);
                continue;
            }
        };
        let pos_0based = pos_1based.saturating_sub(1);
        let end_pos = pos_0based + ref_allele.len() as u64;
        if end_pos > seq.len() as u64 {
            eprintln!(
                "WARNING: Out-of-bounds: contig={}, pos={}, ref_len={}, seq_len={}",
                contig,
                pos_1based,
                ref_allele.len(),
                seq.len()
            );
            continue;
        }

        let slice = &seq[pos_0based as usize..end_pos as usize];
        if !equals_ignore_case(slice, ref_allele.as_bytes()) {
            eprintln!(
                "Mismatch:\n  Contig: {} Pos: {} REF='{}'\n  FASTA => '{}'\n",
                contig,
                pos_str,
                ref_allele,
                String::from_utf8_lossy(slice)
            );
            *mismatch_count += 1;
        }
    }
    lines_checked
}

/// Compare two byte slices ignoring ASCII case (length must match).
fn equals_ignore_case(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    for (x, y) in a.iter().zip(b) {
        if x.to_ascii_uppercase() != y.to_ascii_uppercase() {
            return false;
        }
    }
    true
}
