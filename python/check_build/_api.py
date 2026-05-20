"""Python wrapper around the `check_build` CLI.

The wrapper:
* Locates the `check_build` binary (PATH, env var, explicit constructor arg).
* Translates kwargs to CLI flags via a Verifier builder.
* Parses the CLI's stdout into typed dataclasses.
* Maps non-zero exits and stderr patterns onto a small exception hierarchy.

Design notes
------------
The CLI prints human-readable summaries, not machine-readable JSON, so we
parse stdout with a couple of small regexes — but the contract is strict
enough that any future incompatible CLI change will surface as a clear
`CheckBuildExecutionFailed` rather than a silent miscalculation.

We treat exit code 1 with mismatches as a *result*, not an error. The
upstream CLI exits 1 when mismatches are found, but the Python wrapper
should return the VerificationResult regardless — callers can inspect
`all_passed` themselves.
"""

from __future__ import annotations

import enum
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple, Union

PathLike = Union[str, os.PathLike]


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


class CheckBuildError(Exception):
    """Base class for all check_build wrapper errors."""


class CheckBuildBinaryNotFound(CheckBuildError, FileNotFoundError):
    """The `check_build` binary could not be located.

    Install it with::

        cargo install check_build

    Or set the ``CHECK_BUILD_BIN`` env var (or pass ``binary=`` explicitly).
    """


class CheckBuildExecutionFailed(CheckBuildError, RuntimeError):
    """The `check_build` binary ran but its output couldn't be parsed."""

    def __init__(self, message: str, *, stdout: str = "", stderr: str = "", returncode: int = 0):
        super().__init__(message)
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode


class InvalidVcfError(CheckBuildError, ValueError):
    """The input VCF is missing, malformed, or has the wrong extension."""


# ---------------------------------------------------------------------------
# Enums and dataclasses
# ---------------------------------------------------------------------------


class Reference(str, enum.Enum):
    HG19 = "hg19"
    HG38 = "hg38"

    @property
    def cli_flag(self) -> str:
        return f"--{self.value}-only"

    @property
    def path_flag(self) -> str:
        return f"--{self.value}-path"


class BuildDetectionStatus(str, enum.Enum):
    """Coarse outcome of a `detect_build` call."""

    CLEAR_HG19 = "clear_hg19"
    CLEAR_HG38 = "clear_hg38"
    AMBIGUOUS = "ambiguous"
    NO_DATA = "no_data"


@dataclass(frozen=True)
class BuildResult:
    """Raw match-rate output from `check_build --detect`."""

    hg19_match_rate: float
    hg38_match_rate: float
    hg19_lines: int = 0
    hg38_lines: int = 0
    hg19_mismatches: int = 0
    hg38_mismatches: int = 0

    @property
    def better_match(self) -> Optional[Reference]:
        if self.hg19_lines == 0 and self.hg38_lines == 0:
            return None
        return Reference.HG19 if self.hg19_match_rate > self.hg38_match_rate else Reference.HG38

    @property
    def is_clear_match(self) -> bool:
        """True iff exactly one build has zero mismatches."""
        return (self.hg19_mismatches == 0) ^ (self.hg38_mismatches == 0)

    def __str__(self) -> str:
        return (
            f"hg19: {self.hg19_match_rate:.1f}% ({self.hg19_lines - self.hg19_mismatches}/"
            f"{self.hg19_lines} matched), hg38: {self.hg38_match_rate:.1f}% "
            f"({self.hg38_lines - self.hg38_mismatches}/{self.hg38_lines} matched)"
        )


@dataclass(frozen=True)
class BuildDetection:
    """Higher-level interpretation of a `BuildResult`."""

    status: BuildDetectionStatus
    likely_build: Optional[Reference]
    confidence: str
    raw: BuildResult

    @property
    def is_detected(self) -> bool:
        return self.status in (BuildDetectionStatus.CLEAR_HG19, BuildDetectionStatus.CLEAR_HG38)


@dataclass(frozen=True)
class VerificationResult:
    """Full `verify_both` output."""

    hg19_lines: int
    hg19_mismatches: int
    hg38_lines: int
    hg38_mismatches: int
    mismatch_details: Tuple[str, ...] = field(default_factory=tuple)

    @property
    def all_passed(self) -> bool:
        return self.hg19_mismatches == 0 and self.hg38_mismatches == 0

    def match_rate(self, reference: Reference) -> float:
        if reference is Reference.HG19:
            lines, mm = self.hg19_lines, self.hg19_mismatches
        else:
            lines, mm = self.hg38_lines, self.hg38_mismatches
        if lines == 0:
            return 0.0
        return (lines - mm) / lines * 100.0

    @property
    def likely_build(self) -> Optional[Reference]:
        if self.hg19_mismatches == 0 and self.hg38_mismatches > 0:
            return Reference.HG19
        if self.hg38_mismatches == 0 and self.hg19_mismatches > 0:
            return Reference.HG38
        return None


# ---------------------------------------------------------------------------
# Binary location
# ---------------------------------------------------------------------------


def locate_binary(override: Optional[PathLike] = None) -> Path:
    """Locate the `check_build` binary or raise `CheckBuildBinaryNotFound`.

    Resolution order:
      1. ``override`` argument, if given.
      2. ``CHECK_BUILD_BIN`` environment variable.
      3. ``check_build`` on ``PATH``.
    """
    if override is not None:
        p = Path(override)
        if not p.exists():
            raise CheckBuildBinaryNotFound(f"check_build binary not found at {p}")
        return p
    env = os.environ.get("CHECK_BUILD_BIN")
    if env:
        p = Path(env)
        if not p.exists():
            raise CheckBuildBinaryNotFound(f"$CHECK_BUILD_BIN points to nonexistent {p}")
        return p
    which = shutil.which("check_build")
    if which:
        return Path(which)
    raise CheckBuildBinaryNotFound(
        "Could not locate `check_build`. Install it with `cargo install check_build`, "
        "or set CHECK_BUILD_BIN, or pass binary=... to the Verifier."
    )


# ---------------------------------------------------------------------------
# Verifier
# ---------------------------------------------------------------------------


class Verifier:
    """Builder for VCF verification calls.

    Example::

        result = (
            Verifier("sample.vcf")
                .quiet()
                .with_reference_paths("/cache/hg19.fa", "/cache/hg38.fa")
                .verify_both()
        )
    """

    __slots__ = (
        "_vcf_path",
        "_binary",
        "_quiet",
        "_summary_only",
        "_hg19_path",
        "_hg38_path",
        "_timeout",
    )

    def __init__(
        self,
        vcf_path: PathLike,
        *,
        binary: Optional[PathLike] = None,
        quiet: bool = True,
        summary_only: bool = False,
        hg19_path: Optional[PathLike] = None,
        hg38_path: Optional[PathLike] = None,
        timeout: Optional[float] = None,
    ) -> None:
        p = Path(vcf_path)
        if not p.exists():
            raise InvalidVcfError(f"VCF file not found: {p}")
        if p.suffix.lower() != ".vcf":
            raise InvalidVcfError(
                f"check_build only accepts plain .vcf (got {p.suffix!r}); decompress .gz first."
            )
        self._vcf_path = p
        self._binary = binary
        self._quiet = quiet
        self._summary_only = summary_only
        self._hg19_path = Path(hg19_path) if hg19_path else None
        self._hg38_path = Path(hg38_path) if hg38_path else None
        self._timeout = timeout

    # Builder methods return *new* instances so users can branch off of a
    # partially-configured Verifier without surprising mutation.

    def _replace(self, **kwargs) -> "Verifier":
        new = Verifier.__new__(Verifier)
        for slot in self.__slots__:
            object.__setattr__(new, slot, kwargs.get(slot[1:], getattr(self, slot)))
        return new

    def quiet(self, value: bool = True) -> "Verifier":
        return self._replace(quiet=value)

    def summary_only(self, value: bool = True) -> "Verifier":
        return self._replace(summary_only=value)

    def with_reference_paths(self, hg19: PathLike, hg38: PathLike) -> "Verifier":
        return self._replace(hg19_path=Path(hg19), hg38_path=Path(hg38))

    def with_hg19_path(self, path: PathLike) -> "Verifier":
        return self._replace(hg19_path=Path(path))

    def with_hg38_path(self, path: PathLike) -> "Verifier":
        return self._replace(hg38_path=Path(path))

    def with_binary(self, path: PathLike) -> "Verifier":
        return self._replace(binary=Path(path))

    def with_timeout(self, seconds: Optional[float]) -> "Verifier":
        return self._replace(timeout=seconds)

    # --- Execution ---------------------------------------------------------

    def detect(self) -> BuildResult:
        """Run ``check_build --detect`` and return the raw match rates."""
        completed = self._run(["--detect"])
        return _parse_detect_output(completed.stdout, completed.stderr, completed.returncode)

    def detect_with_interpretation(self) -> BuildDetection:
        """Like ``detect`` but also returns a coarse interpretation."""
        raw = self.detect()
        if raw.hg19_lines == 0 and raw.hg38_lines == 0:
            return BuildDetection(
                status=BuildDetectionStatus.NO_DATA,
                likely_build=None,
                confidence="none",
                raw=raw,
            )
        if raw.hg19_mismatches == 0 and raw.hg38_mismatches > 0:
            return BuildDetection(
                status=BuildDetectionStatus.CLEAR_HG19,
                likely_build=Reference.HG19,
                confidence="high",
                raw=raw,
            )
        if raw.hg38_mismatches == 0 and raw.hg19_mismatches > 0:
            return BuildDetection(
                status=BuildDetectionStatus.CLEAR_HG38,
                likely_build=Reference.HG38,
                confidence="high",
                raw=raw,
            )
        better = raw.better_match
        diff = abs(raw.hg19_match_rate - raw.hg38_match_rate)
        if diff < 0.5:
            return BuildDetection(
                status=BuildDetectionStatus.AMBIGUOUS,
                likely_build=better,
                confidence="low",
                raw=raw,
            )
        return BuildDetection(
            status=(
                BuildDetectionStatus.CLEAR_HG19 if better is Reference.HG19
                else BuildDetectionStatus.CLEAR_HG38
            ),
            likely_build=better,
            confidence="medium",
            raw=raw,
        )

    def verify_both(self) -> VerificationResult:
        """Run a full hg19+hg38 verification."""
        completed = self._run([])
        return _parse_verify_output(completed.stdout, completed.stderr, completed.returncode)

    def verify_single(self, reference: Reference) -> Tuple[int, int]:
        """Run verification against a single reference. Returns ``(lines, mismatches)``."""
        completed = self._run([reference.cli_flag])
        return _parse_single_output(reference, completed.stdout, completed.stderr, completed.returncode)

    # --- Plumbing ----------------------------------------------------------

    def _run(self, extra_args: List[str]) -> subprocess.CompletedProcess:
        binary = locate_binary(self._binary)
        argv = [str(binary)]
        if self._quiet:
            argv.append("--quiet")
        if self._summary_only:
            argv.append("--summary-only")
        if self._hg19_path is not None:
            argv += ["--hg19-path", str(self._hg19_path)]
        if self._hg38_path is not None:
            argv += ["--hg38-path", str(self._hg38_path)]
        argv += list(extra_args)
        argv.append(str(self._vcf_path))
        try:
            completed = subprocess.run(
                argv,
                capture_output=True,
                text=True,
                timeout=self._timeout,
                check=False,
            )
        except FileNotFoundError as e:
            raise CheckBuildBinaryNotFound(str(e)) from e
        # Catastrophic failures (signals etc.) get raised, but exit code 1
        # with valid stdout is treated as a result.
        if completed.returncode < 0:
            raise CheckBuildExecutionFailed(
                f"check_build was killed by signal {-completed.returncode}",
                stdout=completed.stdout,
                stderr=completed.stderr,
                returncode=completed.returncode,
            )
        return completed


# ---------------------------------------------------------------------------
# Top-level convenience
# ---------------------------------------------------------------------------


def detect_build(
    vcf_path: PathLike,
    *,
    binary: Optional[PathLike] = None,
    hg19_path: Optional[PathLike] = None,
    hg38_path: Optional[PathLike] = None,
    timeout: Optional[float] = None,
) -> BuildDetection:
    """One-shot helper: run ``--detect`` and return an interpreted result."""
    return Verifier(
        vcf_path,
        binary=binary,
        quiet=True,
        hg19_path=hg19_path,
        hg38_path=hg38_path,
        timeout=timeout,
    ).detect_with_interpretation()


def verify(
    vcf_path: PathLike,
    *,
    binary: Optional[PathLike] = None,
    hg19_path: Optional[PathLike] = None,
    hg38_path: Optional[PathLike] = None,
    timeout: Optional[float] = None,
) -> VerificationResult:
    """One-shot helper: run a full both-references verification."""
    return Verifier(
        vcf_path,
        binary=binary,
        quiet=True,
        summary_only=True,
        hg19_path=hg19_path,
        hg38_path=hg38_path,
        timeout=timeout,
    ).verify_both()


def verify_against(
    vcf_path: PathLike,
    reference: Reference,
    *,
    binary: Optional[PathLike] = None,
    hg19_path: Optional[PathLike] = None,
    hg38_path: Optional[PathLike] = None,
    timeout: Optional[float] = None,
) -> Tuple[int, int]:
    return Verifier(
        vcf_path,
        binary=binary,
        quiet=True,
        summary_only=True,
        hg19_path=hg19_path,
        hg38_path=hg38_path,
        timeout=timeout,
    ).verify_single(reference)


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------


_DETECT_RE = re.compile(
    r"hg19:\s*([\d.]+)%\s*\((\d+)/(\d+)\s*matched\),\s*"
    r"hg38:\s*([\d.]+)%\s*\((\d+)/(\d+)\s*matched\)",
    re.IGNORECASE,
)

_VERIFY_HG_RE = re.compile(
    r"(hg19|hg38):\s*(\d+)\s*lines,\s*(\d+)\s*mismatches",
    re.IGNORECASE,
)


def _expected_failure(stdout: str, stderr: str, returncode: int, *, allow_nonzero: bool = True):
    """Common preamble: raise on errors that aren't just 'mismatches found'."""
    blob = f"{stderr}\n{stdout}"
    lower = blob.lower()
    if "error: invalid vcf" in lower or "invalid vcf" in lower:
        raise InvalidVcfError(blob.strip())
    if "error: reference not found" in lower:
        raise CheckBuildExecutionFailed(
            f"Reference file missing: {blob.strip()}",
            stdout=stdout,
            stderr=stderr,
            returncode=returncode,
        )
    if returncode != 0 and not allow_nonzero:
        raise CheckBuildExecutionFailed(
            f"check_build exited with status {returncode}: {stderr.strip() or stdout.strip()}",
            stdout=stdout,
            stderr=stderr,
            returncode=returncode,
        )


def _parse_detect_output(stdout: str, stderr: str, returncode: int) -> BuildResult:
    _expected_failure(stdout, stderr, returncode, allow_nonzero=False)
    m = _DETECT_RE.search(stdout)
    if not m:
        raise CheckBuildExecutionFailed(
            "Could not parse --detect output. Has the check_build CLI format changed?",
            stdout=stdout,
            stderr=stderr,
            returncode=returncode,
        )
    hg19_rate = float(m.group(1))
    hg19_matched = int(m.group(2))
    hg19_lines = int(m.group(3))
    hg38_rate = float(m.group(4))
    hg38_matched = int(m.group(5))
    hg38_lines = int(m.group(6))
    return BuildResult(
        hg19_match_rate=hg19_rate,
        hg38_match_rate=hg38_rate,
        hg19_lines=hg19_lines,
        hg38_lines=hg38_lines,
        hg19_mismatches=hg19_lines - hg19_matched,
        hg38_mismatches=hg38_lines - hg38_matched,
    )


def _parse_verify_output(stdout: str, stderr: str, returncode: int) -> VerificationResult:
    # Exit code 1 just means mismatches found — that's a normal outcome.
    _expected_failure(stdout, stderr, returncode, allow_nonzero=True)
    matches = list(_VERIFY_HG_RE.finditer(stdout))
    if not matches:
        raise CheckBuildExecutionFailed(
            "Could not parse verify_both output.",
            stdout=stdout,
            stderr=stderr,
            returncode=returncode,
        )
    hg19_lines = hg19_mm = hg38_lines = hg38_mm = 0
    for m in matches:
        ref = m.group(1).lower()
        lines = int(m.group(2))
        mm = int(m.group(3))
        if ref == "hg19":
            hg19_lines, hg19_mm = lines, mm
        else:
            hg38_lines, hg38_mm = lines, mm
    # Mismatch details are anything between the "===" banner and the
    # "Verification Summary:" line — handy for debugging but optional.
    details: Tuple[str, ...] = ()
    if "Verification Summary" in stdout:
        head, _, _ = stdout.partition("Verification Summary")
        # Skip the leading "Verifying against..." line and trailing banner.
        lines = [
            ln.strip()
            for ln in head.splitlines()
            if ln.strip()
            and not ln.lstrip().startswith("=")
            and not ln.startswith("Verifying against")
        ]
        details = tuple(lines)
    return VerificationResult(
        hg19_lines=hg19_lines,
        hg19_mismatches=hg19_mm,
        hg38_lines=hg38_lines,
        hg38_mismatches=hg38_mm,
        mismatch_details=details,
    )


def _parse_single_output(
    reference: Reference, stdout: str, stderr: str, returncode: int
) -> Tuple[int, int]:
    _expected_failure(stdout, stderr, returncode, allow_nonzero=True)
    needle = reference.value.lower()
    for m in _VERIFY_HG_RE.finditer(stdout):
        if m.group(1).lower() == needle:
            return int(m.group(2)), int(m.group(3))
    raise CheckBuildExecutionFailed(
        f"Could not find {needle} summary line in single-reference output.",
        stdout=stdout,
        stderr=stderr,
        returncode=returncode,
    )
