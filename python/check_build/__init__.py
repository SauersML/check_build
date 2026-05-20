"""check_build — Python bindings for the SauersML/check_build CLI.

Verify a VCF against hg19 and hg38 reference genomes, or just ask "which
build is this?" and get a structured answer.

This package wraps the `check_build` Rust binary via subprocess. The
benefit over a native FFI binding is process isolation: if the underlying
tool ever crashes on bad input, Python sees a clean exception instead of
a segfault.

Quick start
-----------

>>> from check_build import detect_build
>>> result = detect_build("sample.vcf")
>>> result.likely_build
<Reference.HG38: 'hg38'>
>>> result.hg38_match_rate
99.7

Full verification with a builder:

>>> from check_build import Verifier
>>> v = Verifier("sample.vcf").quiet().verify_both()
>>> v.all_passed
False
>>> v.hg19_mismatches
142
"""

from ._api import (
    Reference,
    BuildResult,
    VerificationResult,
    BuildDetection,
    BuildDetectionStatus,
    Verifier,
    CheckBuildError,
    CheckBuildBinaryNotFound,
    CheckBuildExecutionFailed,
    InvalidVcfError,
    detect_build,
    verify,
    verify_against,
    locate_binary,
)

__all__ = [
    "Reference",
    "BuildResult",
    "VerificationResult",
    "BuildDetection",
    "BuildDetectionStatus",
    "Verifier",
    "CheckBuildError",
    "CheckBuildBinaryNotFound",
    "CheckBuildExecutionFailed",
    "InvalidVcfError",
    "detect_build",
    "verify",
    "verify_against",
    "locate_binary",
]

__version__ = "0.1.0"
