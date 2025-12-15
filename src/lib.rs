//! # check_build
//!
//! A library and CLI tool to verify VCF files against hg19 and hg38 reference genomes.
//!
//! ## Simplest Usage: "What build is my file?"
//!
//! ```rust,no_run
//! use check_build::detect_build;
//!
//! let result = detect_build("my_file.vcf").unwrap();
//! println!("{}", result);  // e.g., "hg38 (100.0% match)"
//! ```
//!
//! ## Standard Usage
//!
//! ```rust,no_run
//! use check_build::Verifier;
//!
//! let result = Verifier::new("file.vcf")
//!     .verify_both()
//!     .unwrap();
//!
//! println!("hg19: {} mismatches", result.hg19_mismatches);
//! println!("hg38: {} mismatches", result.hg38_mismatches);
//!
//! if let Some(build) = result.likely_build() {
//!     println!("Detected build: {:?}", build);
//! }
//! ```
//!
//! ## Advanced Usage
//!
//! ```rust,no_run
//! use check_build::{Verifier, Reference};
//!
//! // Verify against just one reference
//! let (lines, mismatches) = Verifier::new("file.vcf")
//!     .verify_single(Reference::Hg38)
//!     .unwrap();
//!
//! // Use custom reference paths, no progress output
//! let result = Verifier::new("file.vcf")
//!     .hg19_path("/data/refs/hg19.fa")
//!     .hg38_path("/data/refs/hg38.fa")
//!     .quiet()
//!     .silent()
//!     .verify_both()
//!     .unwrap();
//! ```

use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use memchr::memchr;
use rayon::prelude::*;
use reqwest::blocking::Client;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tempfile::TempDir;

// ============================================================================
// Public Constants
// ============================================================================

/// URL for hg19 reference FASTA
pub const HG19_URL: &str = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa";
/// URL for hg38 reference FASTA
pub const HG38_URL: &str = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa";

/// Default local filename for hg19
pub const HG19_DEFAULT_PATH: &str = "hg19.fa";
/// Default local filename for hg38
pub const HG38_DEFAULT_PATH: &str = "hg38.fa";

// ============================================================================
// Internal Constants
// ============================================================================

const EXPECTED_MAX_CONTIG_SIZE: usize = 300_000_000;
const DOWNLOAD_BUFFER_SIZE: usize = 131_072;
const FASTA_READ_BUFFER_SIZE: usize = 262_144;

// ============================================================================
// Public Types
// ============================================================================

/// Reference genome to verify against
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reference {
    /// GRCh37/hg19
    Hg19,
    /// GRCh38/hg38
    Hg38,
}

impl Reference {
    /// Get the default download URL for this reference
    pub fn url(&self) -> &'static str {
        match self {
            Reference::Hg19 => HG19_URL,
            Reference::Hg38 => HG38_URL,
        }
    }

    /// Get the default local filename for this reference
    pub fn default_path(&self) -> &'static str {
        match self {
            Reference::Hg19 => HG19_DEFAULT_PATH,
            Reference::Hg38 => HG38_DEFAULT_PATH,
        }
    }
}

/// Result of verifying a VCF against both references
#[derive(Debug, Clone, Default)]
pub struct VerificationResult {
    /// Lines verified against hg19
    pub hg19_lines: u64,
    /// Mismatches found against hg19
    pub hg19_mismatches: u64,
    /// Lines verified against hg38
    pub hg38_lines: u64,
    /// Mismatches found against hg38
    pub hg38_mismatches: u64,
}

impl VerificationResult {
    /// Returns true if no mismatches on either reference
    pub fn all_passed(&self) -> bool {
        self.hg19_mismatches == 0 && self.hg38_mismatches == 0
    }

    /// Infer which build the VCF is aligned to based on mismatch patterns
    ///
    /// Returns `Some(Reference::Hg19)` if hg19 has 0 mismatches but hg38 doesn't,
    /// `Some(Reference::Hg38)` if vice versa, or `None` if ambiguous.
    pub fn likely_build(&self) -> Option<Reference> {
        match (self.hg19_mismatches, self.hg38_mismatches) {
            (0, m) if m > 0 => Some(Reference::Hg19),
            (m, 0) if m > 0 => Some(Reference::Hg38),
            _ => None,
        }
    }

    /// Get percentage match rate for a reference
    pub fn match_rate(&self, reference: Reference) -> f64 {
        let (lines, mismatches) = match reference {
            Reference::Hg19 => (self.hg19_lines, self.hg19_mismatches),
            Reference::Hg38 => (self.hg38_lines, self.hg38_mismatches),
        };
        if lines == 0 {
            return 0.0;
        }
        ((lines - mismatches) as f64 / lines as f64) * 100.0
    }
}

/// Result of build detection - just the data, caller interprets
#[derive(Debug, Clone)]
pub struct BuildResult {
    /// Match rate against hg19 (0.0 to 100.0)
    pub hg19_match_rate: f64,
    /// Match rate against hg38 (0.0 to 100.0)
    pub hg38_match_rate: f64,
    /// Lines checked against hg19
    pub hg19_lines: u64,
    /// Lines checked against hg38
    pub hg38_lines: u64,
    /// Mismatches against hg19
    pub hg19_mismatches: u64,
    /// Mismatches against hg38
    pub hg38_mismatches: u64,
}

impl BuildResult {
    /// Returns the build with higher match rate, or None if no data
    pub fn better_match(&self) -> Option<Reference> {
        if self.hg19_lines == 0 && self.hg38_lines == 0 {
            return None;
        }
        if self.hg19_match_rate > self.hg38_match_rate {
            Some(Reference::Hg19)
        } else {
            Some(Reference::Hg38)
        }
    }

    /// Returns true if one build has 0 mismatches and the other doesn't
    pub fn is_clear_match(&self) -> bool {
        (self.hg19_mismatches == 0) != (self.hg38_mismatches == 0)
    }
}

impl std::fmt::Display for BuildResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "hg19: {:.1}% ({}/{} matched), hg38: {:.1}% ({}/{} matched)",
            self.hg19_match_rate,
            self.hg19_lines - self.hg19_mismatches,
            self.hg19_lines,
            self.hg38_match_rate,
            self.hg38_lines - self.hg38_mismatches,
            self.hg38_lines
        )
    }
}

// ============================================================================
// Convenience Functions
// ============================================================================

/// Simplest way to check a VCF file against both references
///
/// Downloads references if needed, verifies against both, returns raw match rates.
///
/// # Example
/// ```rust,no_run
/// use check_build::detect_build;
///
/// let result = detect_build("sample.vcf").unwrap();
/// println!("hg19: {:.1}%", result.hg19_match_rate);
/// println!("hg38: {:.1}%", result.hg38_match_rate);
///
/// if let Some(better) = result.better_match() {
///     println!("Better match: {:?}", better);
/// }
/// ```
pub fn detect_build(vcf_path: impl Into<String>) -> Result<BuildResult, VerifyError> {
    let r = Verifier::new(vcf_path).silent().verify_both()?;
    Ok(BuildResult {
        hg19_match_rate: r.match_rate(Reference::Hg19),
        hg38_match_rate: r.match_rate(Reference::Hg38),
        hg19_lines: r.hg19_lines,
        hg38_lines: r.hg38_lines,
        hg19_mismatches: r.hg19_mismatches,
        hg38_mismatches: r.hg38_mismatches,
    })
}

/// Error type for verification operations
#[derive(Debug)]
pub enum VerifyError {
    /// VCF file not found or invalid format
    InvalidVcf(String),
    /// Reference file not found
    ReferenceNotFound(String),
    /// I/O error
    Io(std::io::Error),
    /// Network/download error
    Download(String),
}

impl std::fmt::Display for VerifyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VerifyError::InvalidVcf(msg) => write!(f, "Invalid VCF: {}", msg),
            VerifyError::ReferenceNotFound(path) => write!(f, "Reference not found: {}", path),
            VerifyError::Io(e) => write!(f, "I/O error: {}", e),
            VerifyError::Download(msg) => write!(f, "Download failed: {}", msg),
        }
    }
}

impl std::error::Error for VerifyError {}

impl From<std::io::Error> for VerifyError {
    fn from(e: std::io::Error) -> Self {
        VerifyError::Io(e)
    }
}

// ============================================================================
// Verifier (Builder Pattern)
// ============================================================================

/// Builder for VCF verification with fluent API
///
/// # Example
/// ```rust,no_run
/// use check_build::Verifier;
///
/// let result = Verifier::new("sample.vcf")
///     .quiet()           // no progress bars
///     .silent()          // no mismatch output
///     .verify_both()
///     .unwrap();
/// ```
#[derive(Debug, Clone)]
pub struct Verifier {
    vcf_path: String,
    hg19_path: String,
    hg38_path: String,
    show_progress: bool,
    verbose: bool,
    auto_download: bool,
}

impl Verifier {
    /// Create a new verifier for the given VCF file
    pub fn new(vcf_path: impl Into<String>) -> Self {
        Verifier {
            vcf_path: vcf_path.into(),
            hg19_path: HG19_DEFAULT_PATH.to_string(),
            hg38_path: HG38_DEFAULT_PATH.to_string(),
            show_progress: true,
            verbose: true,
            auto_download: true,
        }
    }

    /// Set custom path for hg19 reference
    pub fn hg19_path(mut self, path: impl Into<String>) -> Self {
        self.hg19_path = path.into();
        self
    }

    /// Set custom path for hg38 reference
    pub fn hg38_path(mut self, path: impl Into<String>) -> Self {
        self.hg38_path = path.into();
        self
    }

    /// Disable progress bars
    pub fn quiet(mut self) -> Self {
        self.show_progress = false;
        self
    }

    /// Disable mismatch detail output (summary only)
    pub fn silent(mut self) -> Self {
        self.verbose = false;
        self
    }

    /// Disable automatic download of missing references
    pub fn no_download(mut self) -> Self {
        self.auto_download = false;
        self
    }

    /// Verify against both hg19 and hg38 in parallel
    pub fn verify_both(&self) -> Result<VerificationResult, VerifyError> {
        self.validate_vcf()?;
        self.ensure_references()?;

        let split_dir = tempfile::tempdir()?;
        let contig_file_map = split_vcf_by_contig(&self.vcf_path, &split_dir)?;
        let contig_file_map = Arc::new(contig_file_map);

        let configs = [
            (&self.hg19_path, Reference::Hg19),
            (&self.hg38_path, Reference::Hg38),
        ];

        let results: Vec<(Reference, u64, u64)> = configs
            .par_iter()
            .map(|(path, reference)| {
                let (lines, mismatches) = verify_reference_streaming(
                    path,
                    &contig_file_map,
                    self.show_progress,
                    self.verbose,
                );
                (*reference, lines, mismatches)
            })
            .collect();

        let mut result = VerificationResult::default();
        for (reference, lines, mismatches) in results {
            match reference {
                Reference::Hg19 => {
                    result.hg19_lines = lines;
                    result.hg19_mismatches = mismatches;
                }
                Reference::Hg38 => {
                    result.hg38_lines = lines;
                    result.hg38_mismatches = mismatches;
                }
            }
        }

        Ok(result)
    }

    /// Verify against a single reference
    ///
    /// Returns `(lines_checked, mismatches)`
    pub fn verify_single(&self, reference: Reference) -> Result<(u64, u64), VerifyError> {
        self.validate_vcf()?;

        let ref_path = match reference {
            Reference::Hg19 => &self.hg19_path,
            Reference::Hg38 => &self.hg38_path,
        };

        if self.auto_download {
            ensure_reference(ref_path, reference.url(), self.show_progress)?;
        } else if !Path::new(ref_path).exists() {
            return Err(VerifyError::ReferenceNotFound(ref_path.clone()));
        }

        let split_dir = tempfile::tempdir()?;
        let contig_file_map = split_vcf_by_contig(&self.vcf_path, &split_dir)?;

        Ok(verify_reference_streaming(
            ref_path,
            &contig_file_map,
            self.show_progress,
            self.verbose,
        ))
    }

    fn validate_vcf(&self) -> Result<(), VerifyError> {
        if !self.vcf_path.ends_with(".vcf") {
            return Err(VerifyError::InvalidVcf(
                "File must have .vcf extension (not .gz)".to_string(),
            ));
        }
        if !Path::new(&self.vcf_path).exists() {
            return Err(VerifyError::InvalidVcf(format!(
                "File not found: {}",
                self.vcf_path
            )));
        }
        Ok(())
    }

    fn ensure_references(&self) -> Result<(), VerifyError> {
        if self.auto_download {
            ensure_reference(&self.hg19_path, HG19_URL, self.show_progress)?;
            ensure_reference(&self.hg38_path, HG38_URL, self.show_progress)?;
        } else {
            if !Path::new(&self.hg19_path).exists() {
                return Err(VerifyError::ReferenceNotFound(self.hg19_path.clone()));
            }
            if !Path::new(&self.hg38_path).exists() {
                return Err(VerifyError::ReferenceNotFound(self.hg38_path.clone()));
            }
        }
        Ok(())
    }
}

// ============================================================================
// Standalone Functions (for advanced users)
// ============================================================================

/// Ensure a reference file exists, downloading if necessary
pub fn ensure_reference(path: &str, url: &str, show_progress: bool) -> Result<(), VerifyError> {
    if Path::new(path).exists() {
        return Ok(());
    }

    download_file(url, path, show_progress)?;

    // Check for gzip and decompress
    let mut file = File::open(path)?;
    let mut buffer = [0; 2];
    if file.read_exact(&mut buffer).is_ok() && buffer == [0x1f, 0x8b] {
        decompress_gzip_in_place(path)?;
    }

    Ok(())
}

/// Generate contig name candidates for flexible matching
///
/// Handles both "chr1" -> "1" and "1" -> "chr1" conversions
pub fn get_contig_candidates(contig: &str) -> Vec<String> {
    let mut candidates = Vec::with_capacity(2);
    candidates.push(contig.to_string());

    if let Some(stripped) = contig.strip_prefix("chr") {
        if !stripped.is_empty() {
            candidates.push(stripped.to_string());
        }
    } else {
        candidates.push(format!("chr{}", contig));
    }

    candidates
}

/// Compare byte slices ignoring ASCII case
#[inline]
pub fn equals_ignore_case(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|(x, y)| x.to_ascii_uppercase() == y.to_ascii_uppercase())
}

// ============================================================================
// Internal Implementation
// ============================================================================

fn download_file(url: &str, dest_path: &str, show_progress: bool) -> Result<(), VerifyError> {
    if show_progress {
        eprintln!("Downloading {} to {}...", url, dest_path);
    }

    let client = Client::new();
    let mut response = client
        .get(url)
        .send()
        .map_err(|e| VerifyError::Download(e.to_string()))?;

    if !response.status().is_success() {
        return Err(VerifyError::Download(format!(
            "HTTP {}",
            response.status()
        )));
    }

    let total_size = response.content_length();

    let pb = if show_progress {
        let pb = ProgressBar::new(total_size.unwrap_or(0));
        pb.set_draw_target(ProgressDrawTarget::stderr());
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{msg} [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")
                .expect("template")
                .progress_chars("=>-"),
        );
        pb.set_message(dest_path.to_string());
        Some(pb)
    } else {
        None
    };

    let out_file = File::create(dest_path)?;
    let mut writer = BufWriter::with_capacity(DOWNLOAD_BUFFER_SIZE, out_file);
    let mut downloaded = 0u64;
    let mut buffer = [0u8; DOWNLOAD_BUFFER_SIZE];

    loop {
        let n = response.read(&mut buffer).map_err(VerifyError::Io)?;
        if n == 0 {
            break;
        }
        writer.write_all(&buffer[..n])?;
        downloaded += n as u64;
        if let Some(ref pb) = pb {
            pb.set_position(downloaded);
        }
    }
    writer.flush()?;

    if let Some(pb) = pb {
        pb.finish_with_message(format!("Downloaded {}", dest_path));
    }

    Ok(())
}

fn decompress_gzip_in_place(path: &str) -> Result<(), VerifyError> {
    let compressed = std::fs::read(path)?;
    let mut decoder = GzDecoder::new(&compressed[..]);
    let mut decompressed = Vec::new();
    decoder.read_to_end(&mut decompressed).map_err(VerifyError::Io)?;

    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(DOWNLOAD_BUFFER_SIZE, file);
    writer.write_all(&decompressed)?;
    writer.flush()?;

    Ok(())
}

/// Split VCF by contig into temporary files
pub fn split_vcf_by_contig(
    vcf_path: &str,
    split_dir: &TempDir,
) -> Result<HashMap<String, PathBuf>, VerifyError> {
    let f = File::open(vcf_path)?;
    let reader = BufReader::with_capacity(FASTA_READ_BUFFER_SIZE, f);

    let mut handles: HashMap<String, BufWriter<File>> = HashMap::new();
    let mut file_paths: HashMap<String, PathBuf> = HashMap::new();

    for line_res in reader.lines() {
        let line = line_res?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let mut parts = line.split('\t');
        let contig = match parts.next() {
            Some(c) if !c.is_empty() => c,
            _ => continue,
        };

        // Validate: need POS and REF columns
        let pos = parts.next();
        let _id = parts.next();
        let ref_allele = parts.next();

        if pos.is_none() || ref_allele.is_none() || pos.unwrap().parse::<u64>().is_err() {
            continue;
        }

        let contig_owned = contig.to_string();

        if !handles.contains_key(&contig_owned) {
            let contig_file = split_dir.path().join(format!("{}.tmp", contig_owned));
            let f = OpenOptions::new()
                .create(true)
                .append(true)
                .write(true)
                .open(&contig_file)?;
            handles.insert(contig_owned.clone(), BufWriter::with_capacity(65536, f));
            file_paths.insert(contig_owned.clone(), contig_file);
        }

        if let Some(h) = handles.get_mut(&contig_owned) {
            writeln!(h, "{}", line)?;
        }
    }

    for (_, mut w) in handles {
        w.flush()?;
    }

    Ok(file_paths)
}

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
            current_seq: Vec::with_capacity(EXPECTED_MAX_CONTIG_SIZE),
            lines_checked: 0,
            mismatches: 0,
        }
    }

    fn flush(&mut self, contig_file_map: &HashMap<String, PathBuf>, verbose: bool) {
        if self.current_contig.is_empty() {
            return;
        }
        let cnt = verify_contig(
            &self.current_contig,
            &self.current_seq,
            contig_file_map,
            &mut self.mismatches,
            verbose,
        );
        self.lines_checked += cnt;
        self.current_contig.clear();
        self.current_seq.clear();
    }
}

fn verify_reference_streaming(
    fasta_path: &str,
    contig_file_map: &HashMap<String, PathBuf>,
    show_progress: bool,
    verbose: bool,
) -> (u64, u64) {
    let meta = match std::fs::metadata(fasta_path) {
        Ok(m) => m,
        Err(_) => return (0, 0),
    };

    let pb = if show_progress {
        let pb = ProgressBar::new(meta.len());
        pb.set_draw_target(ProgressDrawTarget::stderr());
        pb.set_style(
            ProgressStyle::default_bar()
                .template(&format!(
                    "{} {{spinner}} [{{bar:40.cyan/blue}}] {{bytes}}/{{total_bytes}} ({{eta}})",
                    fasta_path
                ))
                .expect("template")
                .progress_chars("=>-"),
        );
        Some(pb)
    } else {
        None
    };

    let f = match File::open(fasta_path) {
        Ok(x) => x,
        Err(_) => return (0, 0),
    };
    let mut reader = BufReader::with_capacity(FASTA_READ_BUFFER_SIZE, f);

    let mut chunk = vec![0u8; FASTA_READ_BUFFER_SIZE];
    let mut read_bytes = 0u64;
    let mut line_buf = Vec::<u8>::with_capacity(1024);
    let mut state = ContigState::new();

    loop {
        let n = match reader.read(&mut chunk) {
            Ok(x) => x,
            Err(_) => break,
        };
        if n == 0 {
            break;
        }
        read_bytes += n as u64;
        if let Some(ref pb) = pb {
            pb.set_position(read_bytes);
        }
        process_chunk(&mut state, &mut line_buf, &chunk[..n], contig_file_map, verbose);
    }

    if !line_buf.is_empty() {
        process_line(&mut state, &mut line_buf, contig_file_map, verbose);
    }
    state.flush(contig_file_map, verbose);

    if let Some(pb) = pb {
        pb.finish_with_message(format!("Done {}", fasta_path));
    }

    (state.lines_checked, state.mismatches)
}

fn process_chunk(
    state: &mut ContigState,
    line_buf: &mut Vec<u8>,
    chunk: &[u8],
    contig_file_map: &HashMap<String, PathBuf>,
    verbose: bool,
) {
    let mut start = 0;
    while start < chunk.len() {
        match memchr(b'\n', &chunk[start..]) {
            Some(pos) => {
                let end = start + pos;
                let line_end = if end > start && chunk[end - 1] == b'\r' {
                    end - 1
                } else {
                    end
                };
                line_buf.extend_from_slice(&chunk[start..line_end]);
                process_line(state, line_buf, contig_file_map, verbose);
                line_buf.clear();
                start = end + 1;
            }
            None => {
                for &b in &chunk[start..] {
                    if b != b'\r' {
                        line_buf.push(b);
                    }
                }
                break;
            }
        }
    }
}

fn process_line(
    state: &mut ContigState,
    line_buf: &mut Vec<u8>,
    contig_file_map: &HashMap<String, PathBuf>,
    verbose: bool,
) {
    if line_buf.is_empty() {
        return;
    }
    if line_buf[0] == b'>' {
        state.flush(contig_file_map, verbose);
        let line_str = String::from_utf8_lossy(&line_buf[1..]);
        state.current_contig = line_str.split_whitespace().next().unwrap_or("").to_string();
    } else {
        for b in line_buf.iter_mut() {
            b.make_ascii_uppercase();
        }
        state.current_seq.extend_from_slice(line_buf);
    }
}

/// Verify VCF lines for a contig against reference sequence
pub fn verify_contig(
    contig: &str,
    seq: &[u8],
    contig_file_map: &HashMap<String, PathBuf>,
    mismatch_count: &mut u64,
    verbose: bool,
) -> u64 {
    let candidates = get_contig_candidates(contig);
    let file_path = candidates
        .iter()
        .find_map(|c| contig_file_map.get(c));

    let contig_file = match file_path {
        Some(p) => p,
        None => return 0,
    };

    let f = match File::open(contig_file) {
        Ok(x) => x,
        Err(_) => return 0,
    };

    let mut lines_checked = 0u64;
    let reader = BufReader::with_capacity(65536, f);

    for line_res in reader.lines() {
        let line = match line_res {
            Ok(x) => x,
            Err(_) => break,
        };
        lines_checked += 1;

        let mut parts = line.split('\t');
        let _contig = parts.next();
        let pos_str = match parts.next() {
            Some(p) => p,
            None => continue,
        };
        let _id = parts.next();
        let ref_allele = match parts.next() {
            Some(r) => r,
            None => continue,
        };

        let pos_1based = match pos_str.parse::<u64>() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let pos_0based = pos_1based.saturating_sub(1);
        let end_pos = pos_0based + ref_allele.len() as u64;

        if end_pos > seq.len() as u64 {
            if verbose {
                eprintln!(
                    "WARNING: Out-of-bounds: {}:{} ref_len={} seq_len={}",
                    contig, pos_1based, ref_allele.len(), seq.len()
                );
            }
            continue;
        }

        let slice = &seq[pos_0based as usize..end_pos as usize];
        if !equals_ignore_case(slice, ref_allele.as_bytes()) {
            if verbose {
                eprintln!(
                    "Mismatch: {}:{} REF='{}' FASTA='{}'",
                    contig, pos_str, ref_allele, String::from_utf8_lossy(slice)
                );
            }
            *mismatch_count += 1;
        }
    }

    lines_checked
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_equals_ignore_case() {
        assert!(equals_ignore_case(b"ACGT", b"ACGT"));
        assert!(equals_ignore_case(b"acgt", b"ACGT"));
        assert!(equals_ignore_case(b"AcGt", b"aCgT"));
        assert!(!equals_ignore_case(b"ACGT", b"ACG"));
        assert!(!equals_ignore_case(b"ACGT", b"ACGA"));
        assert!(equals_ignore_case(b"", b""));
    }

    #[test]
    fn test_get_contig_candidates() {
        let c = get_contig_candidates("chr1");
        assert!(c.contains(&"chr1".to_string()));
        assert!(c.contains(&"1".to_string()));
        assert!(!c.iter().any(|x| x.contains("chrchr")));

        let c = get_contig_candidates("1");
        assert!(c.contains(&"1".to_string()));
        assert!(c.contains(&"chr1".to_string()));
    }

    #[test]
    fn test_verification_result_likely_build() {
        let r1 = VerificationResult {
            hg19_lines: 100,
            hg19_mismatches: 0,
            hg38_lines: 100,
            hg38_mismatches: 50,
        };
        assert_eq!(r1.likely_build(), Some(Reference::Hg19));

        let r2 = VerificationResult {
            hg19_lines: 100,
            hg19_mismatches: 50,
            hg38_lines: 100,
            hg38_mismatches: 0,
        };
        assert_eq!(r2.likely_build(), Some(Reference::Hg38));

        let r3 = VerificationResult {
            hg19_lines: 100,
            hg19_mismatches: 0,
            hg38_lines: 100,
            hg38_mismatches: 0,
        };
        assert_eq!(r3.likely_build(), None);
    }

    #[test]
    fn test_vcf_split_and_verify_roundtrip() {
        use tempfile::tempdir;

        let test_dir = tempdir().unwrap();
        let vcf_path = test_dir.path().join("test.vcf");

        {
            let mut f = File::create(&vcf_path).unwrap();
            writeln!(f, "##fileformat=VCFv4.2").unwrap();
            writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
            writeln!(f, "chr1\t100\t.\tA\tG\t.\t.\t.").unwrap();
            writeln!(f, "chr1\t200\t.\tCG\tTA\t.\t.\t.").unwrap();
            writeln!(f, "chr2\t50\t.\tT\tC\t.\t.\t.").unwrap();
        }

        let split_dir = tempdir().unwrap();
        let map = split_vcf_by_contig(vcf_path.to_str().unwrap(), &split_dir).unwrap();

        assert_eq!(map.len(), 2);
        assert!(map.contains_key("chr1"));
        assert!(map.contains_key("chr2"));

        // Test matching sequence
        let mut seq = vec![b'N'; 300];
        seq[99] = b'A';
        seq[199] = b'C';
        seq[200] = b'G';

        let mut m = 0;
        let lines = verify_contig("chr1", &seq, &map, &mut m, false);
        assert_eq!(lines, 2);
        assert_eq!(m, 0);

        // Test mismatch
        let mut bad = vec![b'N'; 300];
        bad[99] = b'T'; // wrong
        bad[199] = b'C';
        bad[200] = b'G';

        let mut m2 = 0;
        verify_contig("chr1", &bad, &map, &mut m2, false);
        assert_eq!(m2, 1);
    }
}
