# check_build

A fast, memory-efficient tool to verify VCF files against hg19 and hg38 reference genomes.

## Quick Start

**What build is my file?**

```bash
check_build --detect my_variants.vcf
# Output: Hg38 (100.0% match, high confidence)
```

**Full verification:**

```bash
check_build my_variants.vcf
```

## Installation

```bash
cargo install check_build
```

Or from source:

```bash
git clone https://github.com/SauersML/check_build.git
cd check_build
cargo build --release
```

## Usage

### CLI

```bash
# Simple build detection
check_build --detect sample.vcf

# Full verification with summary
check_build sample.vcf

# Quiet mode (no progress bars)
check_build -q sample.vcf

# Summary only (no mismatch details)
check_build -s sample.vcf

# Single reference
check_build --hg38-only sample.vcf

# Custom reference paths
check_build --hg19-path /data/hg19.fa --hg38-path /data/hg38.fa sample.vcf
```

### Library

Add to `Cargo.toml`:

```toml
[dependencies]
check_build = { git = "https://github.com/SauersML/check_build" }
```

**Simple usage:**

```rust
use check_build::detect_build;

let result = detect_build("sample.vcf")?;
println!("{}", result);  // "Hg38 (100.0% match, high confidence)"
```

**Full control:**

```rust
use check_build::{Verifier, Reference};

let result = Verifier::new("sample.vcf")
    .quiet()
    .verify_both()?;

println!("hg19: {:.1}% match", result.match_rate(Reference::Hg19));
println!("hg38: {:.1}% match", result.match_rate(Reference::Hg38));

// Detailed detection with edge case handling
match result.detect() {
    BuildDetection::Detected { build, confidence, .. } => {
        println!("Build: {:?} ({} confidence)", build, confidence);
    }
    BuildDetection::Ambiguous { reason, .. } => {
        println!("Cannot determine: {}", reason);
    }
    BuildDetection::Unknown { reason, .. } => {
        println!("Problem with file: {}", reason);
    }
    BuildDetection::NoData => {
        println!("No valid variants found");
    }
}
```

## Features

- **Fast**: Parallel verification of hg19/hg38 using rayon
- **Memory-efficient**: Streams references, processes one contig at a time
- **Auto-download**: Fetches reference FASTAs if not present
- **Edge case handling**: Detects ambiguous, unknown, or corrupt files
- **Dual interface**: Both CLI and library

## How It Works

1. Splits VCF by contig into temp files
2. Streams each reference FASTA (never loads full genome)
3. Verifies REF alleles match reference bases
4. Reports match rates and infers build

## Exit Codes

| Code | Meaning |
|------|---------|
| 0    | Success (build detected or verification passed) |
| 1    | Error (file not found, download failed, etc.) |
| 2    | Ambiguous (matches both builds similarly) |
| 3    | Unknown (low match on both, possibly corrupt) |
| 4    | No data (VCF had no valid variants) |

## License

MIT
