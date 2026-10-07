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


def test_a_control_group_limit_caps_the_default(monkeypatch):
    monkeypatch.setattr(system, "cgroup_memory_limit", lambda *_: 4 * 10**9)
    assert system.default_memory_gb() == 3  # 85% of 4 GB


def test_a_machine_whose_memory_cannot_be_read(monkeypatch):
    monkeypatch.setattr(system, "total_memory", lambda: None)
    assert system.default_memory_gb() == 8
