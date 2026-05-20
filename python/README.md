# check_build (Python)

Python wrapper for the [SauersML/check_build](https://github.com/SauersML/check_build)
VCF reference-build verifier.

The wrapper invokes the upstream Rust binary via subprocess and parses
its output into typed dataclasses. Process isolation means a buggy VCF
or reference path produces a Python exception, not a segfault.

## Install

```bash
pip install check_build
# also install the Rust binary:
cargo install check_build
```

The wrapper finds the binary on `PATH`, or via `$CHECK_BUILD_BIN`, or
through the `binary=` keyword argument.

## Quick start

```python
from check_build import detect_build, verify, Reference

# "Which build is this VCF?"
detection = detect_build("sample.vcf")
detection.likely_build         # Reference.HG38
detection.status               # BuildDetectionStatus.CLEAR_HG38
detection.confidence           # "high" / "medium" / "low" / "none"
detection.raw.hg38_match_rate  # 100.0

# Full verification.
result = verify("sample.vcf")
result.all_passed              # bool
result.hg19_mismatches         # int
result.match_rate(Reference.HG19)
```

## Builder

```python
from check_build import Verifier

result = (
    Verifier("sample.vcf")
        .with_reference_paths("/cache/hg19.fa", "/cache/hg38.fa")
        .quiet()
        .verify_both()
)
```

Builder methods return new `Verifier` instances — branching is safe.

## Errors

* `CheckBuildBinaryNotFound` — the CLI isn't installed / on PATH.
* `InvalidVcfError` — the VCF is missing, gzipped, or malformed.
* `CheckBuildExecutionFailed` — the binary ran but produced output we
  couldn't parse (or the reference file was missing).

All subclass `CheckBuildError`.
