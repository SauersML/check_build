# check_build

`check_build` is a command-line tool written in Rust that verifies a plain-text (non-gzipped) VCF file against two reference genomes (hg19 and hg38) using a streaming, low-memory approach. It compares the REF alleles in your VCF with the corresponding bases in the reference FASTA files while only loading one contig at a time.

## Features

- **Streaming Verification:** Reads large reference FASTA files in chunks to avoid high memory usage.
- **VCF Splitting:** Splits the VCF file by contig into temporary files so that only relevant records are processed at a time.
- **Automatic Download:** Automatically downloads reference FASTA files (if not present locally).
- **Dual Reference Comparison:** Verifies VCF records against both hg19 and hg38, helping you determine which build the VCF is aligned to.
- **Minimal Dependencies:** Uses `clap` for argument parsing, `indicatif` for progress indication, `reqwest` for HTTP downloads, and `tempfile` for managing temporary files.

## Installation

Install via Cargo:

```bash
cargo install check_build
```

Alternatively, clone the repository and build from source:

```bash
git clone https://github.com/SauersML/check_build.git
cd check_build
cargo build --release
```

## Usage

`check_build` expects a plain-text VCF file (non-gzipped). To run the tool, simply execute:

```bash
check_build genome.vcf
```

During execution, the tool will:
1. Check for the existence of the `hg19.fa` and `hg38.fa` files in the working directory. If not found, it downloads them automatically.
2. Split the VCF file into temporary files by contig.
3. Stream each reference FASTA file contig-by-contig and verify each VCF record by comparing its REF allele with the corresponding reference bases.
4. Print a final summary with the total number of lines processed and mismatches for each reference.

### Example Output

```
Verification Summary:
  - hg19 => 4357415 lines, 3298348 mismatches
  - hg38 => 4728611 lines, 0 mismatches
```

This indicates that the VCF records match hg38 perfectly. If a large number of records mismatch (or are out-of-bounds) on hg19, but are aligned well to hg38, it suggests the VCF is aligned to hg38.

## How It Works

- **VCF Splitting:** The VCF file is streamed line-by-line and split into temporary files for each contig. This minimizes memory usage since only one line is processed at a time.
- **Streaming Reference Verification:** Instead of loading the entire reference genome into memory, the FASTA file is read in 64 KB chunks. Contigs are processed sequentially so that only the sequence for the current contig is held in memory. Once a contig is processed (i.e. when a new contig header is encountered), the corresponding VCF records are verified, and the memory is cleared before moving to the next contig.
- **Mismatch Reporting:** If a VCF line references a genomic position that is out-of-bounds (or if the REF allele does not match the reference), a warning is printed (a lot of these are expected). The final summary shows the total number of mismatches per reference build.
