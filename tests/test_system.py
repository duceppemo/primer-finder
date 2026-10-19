"""How many CPUs and how much memory the run may use."""

from __future__ import annotations

import os

import pytest

from primer_finder import PrimerFinderError, system


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


def test_a_tight_control_group_still_gives_kmc_enough(monkeypatch):
    """KMC refuses -m below 2, so a small machine or a tight cgroup must not produce 1. This holds the
    floor here to the one KMC actually enforces."""
    from primer_finder.kmers import MIN_MEMORY_GB

    monkeypatch.setattr(system, "cgroup_memory_limit", lambda *_: 10**9)  # 1 GB: 85% of it is 0 GB
    assert system.default_memory_gb() == MIN_MEMORY_GB
    assert MIN_MEMORY_GB == 2


def test_a_machine_whose_memory_cannot_be_read(monkeypatch):
    monkeypatch.setattr(system, "total_memory", lambda: None)
    assert system.default_memory_gb() == 8


def test_a_control_group_cpu_quota_is_read(tmp_path):
    """cgroup v2 writes "<quota> <period>" in cpu.max, where the period is 100,000 microseconds."""
    quota = tmp_path / "cpu.max"
    quota.write_text("400000 100000\n")
    assert system.cgroup_cpu_limit((str(quota),)) == 4
    quota.write_text("max 100000\n")  # The format of a group with no quota, as this machine has
    assert system.cgroup_cpu_limit((str(quota),)) is None
    quota.write_text("50000 100000\n")  # Half a CPU: there is no half a thread
    assert system.cgroup_cpu_limit((str(quota),)) == 1
    quota.write_text("\n")  # And a file with nothing in it is not a limit either
    assert system.cgroup_cpu_limit((str(quota),)) is None


def test_a_cgroup_v1_cpu_quota_is_read_with_its_period(tmp_path):
    """v1 keeps the quota and the period in two files, and -1 means no quota."""
    quota, period = tmp_path / "cpu.cfs_quota_us", tmp_path / "cpu.cfs_period_us"
    period.write_text("50000\n")
    quota.write_text("600000\n")
    assert system.cgroup_cpu_limit((str(quota),)) == 12
    quota.write_text("-1\n")
    assert system.cgroup_cpu_limit((str(quota),)) is None
    quota.write_text("600000\n")
    period.unlink()  # Without the period there is nothing to divide by, so no limit is claimed
    assert system.cgroup_cpu_limit((str(quota),)) is None


def test_the_cpu_quota_caps_the_usable_cpus(monkeypatch):
    """`sched_getaffinity` sees a cpuset but not a quota, so -t used to be the whole machine under a
    scheduler that limits CPU *time*."""
    mask = len(os.sched_getaffinity(0))
    monkeypatch.setattr(system, "cgroup_cpu_limit", lambda *_: 4)
    assert system.usable_cpus() == min(4, mask)
    monkeypatch.setattr(system, "cgroup_cpu_limit", lambda *_: 100000)
    assert system.usable_cpus() == mask  # Never more than the affinity mask
    monkeypatch.setattr(system, "cgroup_cpu_limit", lambda *_: None)
    assert system.usable_cpus() >= 1


def test_the_files_of_a_cpu_quota_hierarchy_are_listed(tmp_path):
    proc = tmp_path / "cgroup"
    proc.write_text("0::/slurm/job_42\n")
    files = system.cgroup_cpu_files(proc, tmp_path / "fs")
    assert files[:2] == system.CGROUP_CPU_FILES
    assert [name[len(str(tmp_path / "fs")):] for name in files[2:]] == [
        "/slurm/job_42/cpu.max", "/slurm/cpu.max", "/cpu.max"]


def test_two_runs_cannot_share_an_output_folder(tmp_path):
    """They would count kmers into the same database and blast against the same one."""
    folder = tmp_path / "out"
    with system.lock_output(folder):
        assert (folder / system.LOCK_NAME).read_text().startswith(f"pid {os.getpid()} ")
        with pytest.raises(PrimerFinderError, match="Another primer-finder design is using"), \
                system.lock_output(folder, "design"):
            pass
    # And the claim is dropped on the way out, so the next run may take it
    with system.lock_output(folder):
        pass
