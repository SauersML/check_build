# check_build (Python)

Python wrapper around the
[`SauersML/check_build`](https://github.com/SauersML/check_build) Rust CLI.
Ask "which reference build is this VCF aligned to?" and get a typed,
ergonomic answer.

```python
from check_build import detect_build

detection = detect_build("sample.vcf")
detection.likely_build        # Reference.HG38 / .HG19 / None
detection.status              # BuildDetectionStatus.CLEAR_HG38 / ...
detection.confidence          # 'high' / 'medium' / 'low' / 'none'
detection.raw.hg38_match_rate # 100.0
```

The wrapper invokes the Rust binary via subprocess and parses its
stdout into frozen dataclasses. Process isolation means malformed VCFs
or missing references produce a clean Python exception, not a segfault.

## Install

```bash
pip install check_build
# the Rust binary:
cargo install check_build
```

The wrapper finds the binary in this order:
1. `binary=` keyword argument (explicit path).
2. `$CHECK_BUILD_BIN` env var.
3. `check_build` on `PATH`.

If none resolve, you get `CheckBuildBinaryNotFound` with the suggested
`cargo install` command.

## Top-level helpers

```python
check_build.detect_build(vcf, *, hg19_path=None, hg38_path=None) -> BuildDetection
check_build.verify(vcf, *, hg19_path=None, hg38_path=None) -> VerificationResult
check_build.verify_against(vcf, Reference.HG38, *, ...) -> (lines, mismatches)
```

`detect_build` is the most common entry point. It runs `check_build
--detect` (fast; samples variants) and interprets the result.

`verify` runs the full hg19+hg38 sweep. Returns mismatch counts and
keeps the raw stderr/stdout for debugging.

## Shortcuts: skip downloads + MD5

`check_build` will download hg19/hg38 references to a cache directory
on first run. Pass them in directly and the binary skips both the
download and the MD5 check:

```python
detection = detect_build(
    "sample.vcf",
    hg19_path="/cache/hg19.fa",
    hg38_path="/cache/hg38.fa",
)
```

`Verifier` exposes the same controls via `with_reference_paths`,
`with_hg19_path`, `with_hg38_path`.

## Builder

For full control, use the `Verifier` builder. Each `with_*` method
returns a *new* instance — branching off a partially-configured
verifier is safe:

```python
from check_build import Verifier, Reference

base = Verifier("sample.vcf").with_reference_paths("/cache/hg19.fa", "/cache/hg38.fa")
quick = base.quiet()
full  = base.summary_only(False)

quick_detection = quick.detect_with_interpretation()
full_result = full.verify_both()
single = full.verify_single(Reference.HG38)
```

## Results

```python
# detect_build -> BuildDetection
detection.is_detected                    # bool — likely_build is non-None and high-confidence
detection.likely_build                   # Reference.HG19 / .HG38 / None
detection.status                         # BuildDetectionStatus enum
detection.confidence                     # 'high' / 'medium' / 'low' / 'none'
detection.raw                            # BuildResult — raw match rates + counts

# verify -> VerificationResult
result.all_passed                        # bool — zero mismatches on both refs
result.likely_build                      # Reference / None — only when one side is clean
result.hg19_lines / hg19_mismatches      # int
result.hg38_lines / hg38_mismatches      # int
result.match_rate(Reference.HG38)        # float 0..100
result.mismatch_details                  # captured mismatch lines (when summary_only=False)
```

## Errors

* `CheckBuildBinaryNotFound` — CLI not installed / not on PATH.
* `InvalidVcfError` — VCF missing, gzipped, or wrong extension.
* `CheckBuildExecutionFailed` — binary ran but produced output we
  couldn't parse, or the reference file was missing on disk.

All subclass `CheckBuildError`.
