"""Tests for the check_build Python wrapper.

We don't ship the upstream Rust binary in CI; tests use a hand-rolled fake
binary written as a tiny Python script. Its stdout/stderr/exit-code
behaviour matches the real check_build CLI for the cases we care about.
"""

from __future__ import annotations

import stat
import textwrap
from pathlib import Path

import pytest

from check_build import (
    BuildDetectionStatus,
    CheckBuildBinaryNotFound,
    CheckBuildError,
    CheckBuildExecutionFailed,
    InvalidVcfError,
    Reference,
    Verifier,
    detect_build,
    locate_binary,
    verify,
    verify_against,
)


# ---------------------------------------------------------------------------
# Fake binary factory
# ---------------------------------------------------------------------------


_FAKE_TEMPLATE = textwrap.dedent(
    """\
    #!/usr/bin/env python3
    import sys

    args = sys.argv[1:]
    detect = "--detect" in args
    hg19_only = "--hg19-only" in args
    hg38_only = "--hg38-only" in args

    {behaviour}
    """
)


def _make_fake(tmp_path: Path, behaviour: str) -> Path:
    p = tmp_path / "fake_check_build"
    p.write_text(_FAKE_TEMPLATE.format(behaviour=behaviour))
    p.chmod(p.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return p


def _make_vcf(tmp_path: Path, name: str = "sample.vcf") -> Path:
    vcf = tmp_path / name
    vcf.write_text(
        "##fileformat=VCFv4.5\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t12345\trs1\tA\tG\t.\t.\t.\n"
    )
    return vcf


# ---------------------------------------------------------------------------
# Locator
# ---------------------------------------------------------------------------


def test_locate_binary_explicit(tmp_path):
    fake = _make_fake(tmp_path, "print('hi')")
    assert locate_binary(fake) == fake


def test_locate_binary_env(tmp_path, monkeypatch):
    fake = _make_fake(tmp_path, "print('hi')")
    monkeypatch.setenv("CHECK_BUILD_BIN", str(fake))
    assert locate_binary() == fake


def test_locate_binary_missing(monkeypatch, tmp_path):
    monkeypatch.delenv("CHECK_BUILD_BIN", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))  # empty
    with pytest.raises(CheckBuildBinaryNotFound):
        locate_binary()


def test_locate_binary_env_pointing_nowhere(monkeypatch, tmp_path):
    monkeypatch.setenv("CHECK_BUILD_BIN", str(tmp_path / "nope"))
    with pytest.raises(CheckBuildBinaryNotFound):
        locate_binary()


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def test_verifier_rejects_missing_vcf(tmp_path):
    with pytest.raises(InvalidVcfError):
        Verifier(tmp_path / "missing.vcf")


def test_verifier_rejects_gz(tmp_path):
    p = tmp_path / "sample.vcf.gz"
    p.write_bytes(b"")
    with pytest.raises(InvalidVcfError):
        Verifier(p)


# ---------------------------------------------------------------------------
# Detect parsing
# ---------------------------------------------------------------------------


def test_detect_clear_hg38(tmp_path):
    fake = _make_fake(
        tmp_path,
        "print('hg19: 92.4% (924/1000 matched), hg38: 100.0% (1000/1000 matched)')\n"
        "print('Best match: Hg38')",
    )
    vcf = _make_vcf(tmp_path)
    detection = detect_build(vcf, binary=fake)
    assert detection.status is BuildDetectionStatus.CLEAR_HG38
    assert detection.likely_build is Reference.HG38
    assert detection.confidence == "high"
    assert detection.raw.hg38_match_rate == 100.0
    assert detection.raw.hg38_mismatches == 0
    assert detection.raw.hg19_mismatches == 76


def test_detect_no_data(tmp_path):
    fake = _make_fake(
        tmp_path,
        "print('hg19: 0.0% (0/0 matched), hg38: 0.0% (0/0 matched)')",
    )
    vcf = _make_vcf(tmp_path)
    detection = detect_build(vcf, binary=fake)
    assert detection.status is BuildDetectionStatus.NO_DATA
    assert detection.likely_build is None


def test_detect_ambiguous(tmp_path):
    fake = _make_fake(
        tmp_path,
        "print('hg19: 99.7% (997/1000 matched), hg38: 99.9% (999/1000 matched)')",
    )
    vcf = _make_vcf(tmp_path)
    detection = detect_build(vcf, binary=fake)
    assert detection.status is BuildDetectionStatus.AMBIGUOUS
    assert detection.confidence == "low"
    # Better match still reported even when low-confidence.
    assert detection.likely_build is Reference.HG38


def test_detect_propagates_invalid_vcf(tmp_path):
    fake = _make_fake(
        tmp_path,
        "import sys; sys.stderr.write('ERROR: Invalid VCF: bogus line\\n'); sys.exit(1)",
    )
    vcf = _make_vcf(tmp_path)
    with pytest.raises(InvalidVcfError):
        detect_build(vcf, binary=fake)


def test_detect_raises_on_unparseable_stdout(tmp_path):
    fake = _make_fake(tmp_path, "print('unexpected garbage')")
    vcf = _make_vcf(tmp_path)
    with pytest.raises(CheckBuildExecutionFailed):
        detect_build(vcf, binary=fake)


# ---------------------------------------------------------------------------
# verify_both parsing
# ---------------------------------------------------------------------------


_VERIFY_OUT = textwrap.dedent(
    """\
    Verifying against hg19 and hg38...
    chr1:12345 expected A got C
    ===================
    Verification Summary:
      hg19: 1000 lines, 142 mismatches (85.8% match)
      hg38: 1000 lines, 0 mismatches (100.0% match)
    Likely reference build: Hg38
    """
)


def test_verify_both(tmp_path):
    fake = _make_fake(
        tmp_path,
        f"import sys; sys.stdout.write({_VERIFY_OUT!r}); sys.exit(1)",
    )
    vcf = _make_vcf(tmp_path)
    result = verify(vcf, binary=fake)
    assert result.hg19_lines == 1000
    assert result.hg19_mismatches == 142
    assert result.hg38_lines == 1000
    assert result.hg38_mismatches == 0
    assert not result.all_passed
    assert result.likely_build is Reference.HG38
    assert result.match_rate(Reference.HG38) == 100.0
    assert "chr1:12345 expected A got C" in result.mismatch_details


def test_verify_both_all_passed(tmp_path):
    out = textwrap.dedent(
        """\
        Verifying against hg19 and hg38...
        ===================
        Verification Summary:
          hg19: 1000 lines, 0 mismatches (100.0% match)
          hg38: 1000 lines, 0 mismatches (100.0% match)
        """
    )
    fake = _make_fake(tmp_path, f"import sys; sys.stdout.write({out!r}); sys.exit(0)")
    vcf = _make_vcf(tmp_path)
    result = verify(vcf, binary=fake)
    assert result.all_passed
    assert result.likely_build is None


def test_verify_single(tmp_path):
    out = textwrap.dedent(
        """\
        Verifying against hg19 only...
        ===================
        hg19: 500 lines, 3 mismatches
        """
    )
    fake = _make_fake(tmp_path, f"import sys; sys.stdout.write({out!r}); sys.exit(1)")
    vcf = _make_vcf(tmp_path)
    lines, mm = verify_against(vcf, Reference.HG19, binary=fake)
    assert (lines, mm) == (500, 3)


def test_verify_propagates_reference_missing(tmp_path):
    fake = _make_fake(
        tmp_path,
        "import sys; sys.stderr.write('ERROR: Reference not found: /no/such/path\\n'); sys.exit(1)",
    )
    vcf = _make_vcf(tmp_path)
    with pytest.raises(CheckBuildExecutionFailed):
        verify(vcf, binary=fake)


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------


def test_verifier_passes_paths_through(tmp_path):
    # Capture invocation by writing argv to a file.
    log = tmp_path / "argv.txt"
    fake = _make_fake(
        tmp_path,
        f"import sys, json; open({str(log)!r}, 'w').write(json.dumps(sys.argv));"
        " print('hg19: 50.0% (5/10 matched), hg38: 90.0% (9/10 matched)'); sys.exit(0)",
    )
    vcf = _make_vcf(tmp_path)

    (
        Verifier(vcf, binary=fake)
        .quiet()
        .with_reference_paths(tmp_path / "hg19.fa", tmp_path / "hg38.fa")
        .detect()
    )

    import json

    argv = json.loads(log.read_text())
    assert "--quiet" in argv
    assert argv[argv.index("--hg19-path") + 1] == str(tmp_path / "hg19.fa")
    assert argv[argv.index("--hg38-path") + 1] == str(tmp_path / "hg38.fa")
    assert "--detect" in argv
    # VCF path is last.
    assert argv[-1] == str(vcf)


def test_verifier_replace_does_not_mutate(tmp_path):
    fake = _make_fake(tmp_path, "")
    vcf = _make_vcf(tmp_path)
    base = Verifier(vcf, binary=fake, quiet=False)
    new = base.quiet()
    assert base._quiet is False  # noqa: SLF001 (test introspection)
    assert new._quiet is True  # noqa: SLF001
    assert new is not base


def test_verifier_with_timeout_kills_runaway(tmp_path):
    import subprocess

    fake = _make_fake(tmp_path, "import time; time.sleep(30); print('hi')")
    vcf = _make_vcf(tmp_path)
    with pytest.raises(subprocess.TimeoutExpired):
        Verifier(vcf, binary=fake, timeout=0.3).detect()


# ---------------------------------------------------------------------------
# BuildResult helpers
# ---------------------------------------------------------------------------


def test_build_result_str_and_helpers(tmp_path):
    fake = _make_fake(
        tmp_path,
        "print('hg19: 70.0% (70/100 matched), hg38: 100.0% (100/100 matched)')",
    )
    vcf = _make_vcf(tmp_path)
    detection = detect_build(vcf, binary=fake)
    raw = detection.raw
    assert raw.is_clear_match
    assert raw.better_match is Reference.HG38
    assert "hg19" in str(raw)
    assert "hg38" in str(raw)


def test_error_hierarchy():
    assert issubclass(InvalidVcfError, CheckBuildError)
    assert issubclass(CheckBuildBinaryNotFound, CheckBuildError)
    assert issubclass(CheckBuildExecutionFailed, CheckBuildError)
