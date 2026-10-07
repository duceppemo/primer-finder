"""How many CPUs and how much memory the run may use."""

from __future__ import annotations

from primer_finder import system


def test_usable_cpus_is_at_least_one():
    assert system.usable_cpus() >= 1


def test_default_memory_is_at_least_one_gb():
    assert system.default_memory_gb() >= 1


def test_the_control_group_limit_is_read(tmp_path):
    limit = tmp_path / "memory.max"
    limit.write_text("2147483648\n")  # 2 GB
    assert system.cgroup_memory_limit((str(limit),)) == 2 * 1024**3


def test_no_control_group_limit(tmp_path):
    missing = tmp_path / "gone"
    unlimited = tmp_path / "memory.max"
    unlimited.write_text("max\n")
    assert system.cgroup_memory_limit((str(missing), str(unlimited))) is None


def test_the_files_of_a_cgroup_v2_hierarchy_are_listed(tmp_path):
    """A job scheduler sets the limit on the job's own control group, or on a parent of it."""
    proc = tmp_path / "cgroup"
    proc.write_text("0::/user.slice/user-1000.slice/job.scope\n")
    files = system.cgroup_memory_files(proc, tmp_path / "fs")
    assert files[:2] == system.CGROUP_MEMORY_FILES
    assert [name[len(str(tmp_path / "fs")):] for name in files[2:]] == [
        "/user.slice/user-1000.slice/job.scope/memory.max",
        "/user.slice/user-1000.slice/memory.max",
        "/user.slice/memory.max",
        "/memory.max",
    ]


def test_the_files_of_a_cgroup_v1_hierarchy_are_listed(tmp_path):
    proc = tmp_path / "cgroup"
    proc.write_text("12:memory,cpu:/slurm/job_42\n3:cpuset:/elsewhere\n")
    files = system.cgroup_memory_files(proc, tmp_path / "fs")
    assert [name[len(str(tmp_path / "fs")):] for name in files[2:]] == [
        "/memory/slurm/job_42/memory.limit_in_bytes",
        "/memory/slurm/memory.limit_in_bytes",
        "/memory/memory.limit_in_bytes",
    ]
    assert not any("cpuset" in name for name in files)


def test_a_missing_or_odd_proc_cgroup_is_ignored(tmp_path):
    assert system.cgroup_memory_files(tmp_path / "gone", tmp_path) == system.CGROUP_MEMORY_FILES
    odd = tmp_path / "cgroup"
    odd.write_text("nonsense\n\n")
    assert system.cgroup_memory_files(odd, tmp_path) == system.CGROUP_MEMORY_FILES


def test_the_smallest_limit_of_the_hierarchy_wins(tmp_path):
    big, small = tmp_path / "big", tmp_path / "small"
    big.write_text("8000000000")
    small.write_text("2000000000")
    assert system.cgroup_memory_limit((str(big), str(small))) == 2000000000


def test_a_control_group_limit_caps_the_default(monkeypatch):
    monkeypatch.setattr(system, "cgroup_memory_limit", lambda *_: 4 * 10**9)
    assert system.default_memory_gb() == 3  # 85% of 4 GB


def test_a_machine_whose_memory_cannot_be_read(monkeypatch):
    monkeypatch.setattr(system, "total_memory", lambda: None)
    assert system.default_memory_gb() == 8
