"""Fixtures: fasta files, genome folders, and the stub programs on PATH."""

from __future__ import annotations

import json
import os
import stat
import sys
from pathlib import Path

import pytest

from primer_finder.pipeline import Settings
from tests import stub_programs


def write_fasta(path: Path, records: dict[str, str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records.items()))
    return path


@pytest.fixture
def fasta(tmp_path):
    return lambda name, records: write_fasta(tmp_path / name, records)


@pytest.fixture
def genomes(tmp_path):
    """Two inclusion and two exclusion genomes, in their folders."""
    for number in (1, 2):
        write_fasta(tmp_path / "inclusion" / f"inclusion_{number}.fasta", {"chr": "ACGT" * 50})
        write_fasta(tmp_path / "exclusion" / f"exclusion_{number}.fasta", {"chr": "TTGA" * 50})
    return tmp_path / "inclusion", tmp_path / "exclusion"


@pytest.fixture
def stubs(tmp_path, monkeypatch):
    """Put a launcher for every stub program first on PATH. Returns a function to set the scenario, and
    the file every call is appended to."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    tests_dir = str(Path(stub_programs.__file__).resolve().parents[1])
    for name in stub_programs.PROGRAMS:
        launcher = bin_dir / name
        launcher.write_text(
            f"#!{sys.executable}\n"
            f"import sys; sys.path.insert(0, {tests_dir!r})\n"
            f"from tests.stub_programs import main\n"
            f"raise SystemExit(main({name!r}))\n"
        )
        launcher.chmod(launcher.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    calls = tmp_path / "calls.log"
    monkeypatch.setenv("PRIMER_FINDER_CALLS", str(calls))

    def scenario(**values) -> Path:
        monkeypatch.setenv("PRIMER_FINDER_STUBS", json.dumps(values))
        return calls

    scenario()
    return scenario


@pytest.fixture
def settings(genomes, tmp_path):
    inclusion, exclusion = genomes
    defaults = {
        "inclusion": inclusion,
        "exclusion": exclusion,
        "output": tmp_path / "out",
        "threads": 1,
        "memory_gb": 2,
    }
    return lambda **changes: Settings(**{**defaults, **changes})
