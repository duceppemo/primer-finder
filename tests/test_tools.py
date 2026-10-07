"""Running external programs."""

from __future__ import annotations

import stat
import sys

import pytest

from primer_finder import PrimerFinderError, tools


def script(tmp_path, name, body):
    path = tmp_path / name
    path.write_text(f"#!{sys.executable}\n{body}\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def test_require_lists_every_missing_program_and_its_package(monkeypatch, tmp_path):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(PrimerFinderError) as error:
        tools.require(["kmc", "kmc_tools", "blastn"])
    message = str(error.value)
    assert "kmc, kmc_tools, blastn" in message
    assert "conda install -c conda-forge -c bioconda blast kmc" in message  # One package per program


def test_require_is_happy_when_everything_is_there(monkeypatch, tmp_path):
    script(tmp_path, "kmc", "pass")
    monkeypatch.setenv("PATH", str(tmp_path))
    assert tools.require(["kmc"]) is None


def test_run_captures_the_output(tmp_path):
    program = script(tmp_path, "hello", "import sys; print('out'); print('err', file=sys.stderr)")
    assert tools.run([program]) == "out\nerr\n"


def test_run_writes_the_standard_output_to_a_file(tmp_path):
    program = script(tmp_path, "hello", "import sys; print('data'); print('log', file=sys.stderr)")
    out = tmp_path / "out.txt"
    assert tools.run([program], stdout_path=out) == "log\n"
    assert out.read_text() == "data\n"


def test_a_failing_program_raises_with_the_last_lines_of_its_output(tmp_path):
    program = script(tmp_path, "bad", "import sys; print('why it failed', file=sys.stderr); sys.exit(3)")
    with pytest.raises(tools.ToolError, match="(?s)bad failed .exit code 3.*why it failed"):
        tools.run([program])


def test_a_program_that_does_not_exist(tmp_path):
    with pytest.raises(tools.ToolError, match="Program not found"):
        tools.run([tmp_path / "nothing"])


def test_version_never_raises(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    assert tools.version("nothing") == "not found"
    script(tmp_path, "quiet", "pass")
    assert tools.version("quiet") == "unknown"
    script(tmp_path, "talks", "print('talks 1.2.3')")
    assert tools.version("talks") == "talks 1.2.3"
    assert tools.versions(["talks"]) == {"talks": "talks 1.2.3"}


def test_the_version_flag_each_program_wants(tmp_path, monkeypatch):
    """blast wants -version, and KMC prints its banner with no argument."""
    monkeypatch.setenv("PATH", str(tmp_path))
    script(tmp_path, "blastn", "import sys; print('USAGE') if '--version' in sys.argv else print('blastn: 2.17.0+')")
    script(tmp_path, "kmc", "import sys; print('K-Mer Counter (KMC) ver. 3.2.4' if len(sys.argv) == 1 else 'Usage:')")
    assert tools.version("blastn") == "blastn: 2.17.0+"
    assert tools.version("kmc").startswith("K-Mer Counter")
    assert tools.version("blastn", "--version") == "USAGE"  # An explicit flag wins
