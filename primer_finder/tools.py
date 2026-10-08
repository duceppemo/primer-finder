"""Running the external programs, with their command lines and output in the log."""

from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from primer_finder import PrimerFinderError

log = logging.getLogger(__name__)

# The conda package that provides each program, for the "missing program" message.
PACKAGES = {
    "kmc": "kmc",
    "primer3_core": "primer3",
    "kmc_tools": "kmc",
    "skesa": "skesa",
    "spades.py": "spades",
    "minimap2": "minimap2",
    "blastn": "blast",
    "makeblastdb": "blast",
}

# How to ask each program for its version: blast wants one dash, and KMC prints its banner with no
# argument at all.
VERSION_FLAGS = {
    "kmc": (),
    "kmc_tools": (),
    "blastn": ("-version",),
    "makeblastdb": ("-version",),
    "primer3_core": ("--about",),  # --version prints nothing at all
}

# How many lines of a failed program's output to put in the error message.
ERROR_LINES = 15


class ToolError(PrimerFinderError):
    """An external program failed."""


def which(name: str) -> str | None:
    """The path of an executable, or None."""
    return shutil.which(name)


def require(names: Sequence[str]) -> None:
    """Check that every program is on PATH, or raise one error listing all the missing ones."""
    missing = [name for name in names if which(name) is None]
    if missing:
        packages = sorted({PACKAGES.get(name, name) for name in missing})
        raise PrimerFinderError(
            f"Required program(s) not found on PATH: {', '.join(missing)}. "
            f"Install them with: conda install -c conda-forge -c bioconda {' '.join(packages)}"
        )


class ToolTimeout(ToolError):
    """An external program was still running when its time was up."""


def run(cmd: Sequence[str | Path], *, stdout_path: Path | None = None, cwd: Path | None = None,
        timeout: float | None = None) -> str:
    """Run a program, waiting for it to finish.

    Its output is captured and logged; with `stdout_path`, its standard output goes to that file instead.
    With `timeout`, a program still running after that many seconds is killed and ToolTimeout raised:
    some inputs make a program search for an answer that does not exist, for as long as it is allowed to.
    Returns the captured output. Raises ToolError if the program fails or is not found.
    """
    argv = [str(part) for part in cmd]
    log.debug("Running: %s", shlex.join(argv))
    stdout = stdout_path.open("wb") if stdout_path else subprocess.PIPE
    try:
        process = subprocess.run(argv, stdout=stdout, stderr=subprocess.PIPE, cwd=cwd, timeout=timeout)
    except FileNotFoundError as exc:
        raise ToolError(f"Program not found: {argv[0]}") from exc
    except subprocess.TimeoutExpired as exc:
        raise ToolTimeout(
            f"{argv[0]} was still running after {timeout:g} s and was stopped:\n  {shlex.join(argv)}"
        ) from exc
    finally:
        if stdout_path:
            stdout.close()
    output = (process.stderr or b"").decode(errors="replace")
    if not stdout_path:
        output = (process.stdout or b"").decode(errors="replace") + output
    if output.strip():
        log.debug("%s output:\n%s", argv[0], output.strip())
    if process.returncode != 0:
        tail = "\n".join(output.strip().splitlines()[-ERROR_LINES:])
        raise ToolError(
            f"{argv[0]} failed (exit code {process.returncode}):\n  {shlex.join(argv)}"
            + (f"\n{tail}" if tail else "")
        )
    return output


def version(name: str, *flags: str) -> str:
    """The version a program reports, on one line, or "unknown". Never raises: this is for the record only."""
    if which(name) is None:
        return "not found"
    attempts = [[name, flag] for flag in flags or VERSION_FLAGS.get(name, ("--version",))] or [[name]]
    for argv in attempts:
        try:
            process = subprocess.run(argv, capture_output=True, timeout=60)
        except (OSError, subprocess.SubprocessError):  # pragma: no cover - a program that cannot be started
            continue
        text = (process.stdout + process.stderr).decode(errors="replace")
        for line in text.splitlines():
            if line.strip():
                return line.strip()
    return "unknown"


def versions(names: Sequence[str]) -> dict[str, str]:
    """What every program reports as its version, for run_info.json."""
    return {name: version(name) for name in names}
