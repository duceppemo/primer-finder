"""How many CPUs and how much memory this process may use."""

from __future__ import annotations

import os
from pathlib import Path

# Where a control group's memory limit is, when this process sees its own cgroup as the root (containers).
CGROUP_MEMORY_FILES = ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes")  # v2, v1
CGROUP_ROOT = Path("/sys/fs/cgroup")
PROC_CGROUP = Path("/proc/self/cgroup")
MEMORY_FRACTION = 0.85  # Of the total memory, as the previous versions did


def usable_cpus() -> int:
    """The CPUs primer-finder may use: all of them, or those a job scheduler or `taskset` gave it."""
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0)) or 1
    return os.cpu_count() or 1  # pragma: no cover - macOS


def cgroup_memory_files(proc_cgroup: Path = PROC_CGROUP, root: Path = CGROUP_ROOT) -> tuple[str, ...]:
    """Every file that could hold a memory limit for this process: the roots, then this process's own
    control group and each of its parents.

    A job scheduler (Slurm, systemd) puts the job in a control group below the root, so the limit is not in
    `/sys/fs/cgroup/memory.max` but in `/sys/fs/cgroup<path from /proc/self/cgroup>/memory.max`, and it may
    be set on a parent of it. In a container the process sees its cgroup as the root instead, which is what
    the first two files cover.
    """
    files = list(CGROUP_MEMORY_FILES)
    try:
        lines = proc_cgroup.read_text().splitlines()
    except OSError:  # pragma: no cover - not Linux
        return tuple(files)
    for line in lines:
        parts = line.split(":", 2)
        if len(parts) != 3:
            continue
        _, controllers, path = parts
        if controllers == "":  # cgroup v2: one line, "0::<path>"
            base, limit = root, "memory.max"
        elif "memory" in controllers.split(","):  # cgroup v1: one line per controller
            base, limit = root / "memory", "memory.limit_in_bytes"
        else:
            continue
        group = base / path.strip("/") if path.strip("/") else base
        while True:
            files.append(str(group / limit))
            if group == base or base not in group.parents:
                break
            group = group.parent
    return tuple(dict.fromkeys(files))  # In order, without duplicates


def cgroup_memory_limit(files: tuple[str, ...] | None = None) -> int | None:
    """The smallest memory limit any of this process's control groups sets, in bytes, or None."""
    limits = []
    for path in files if files is not None else cgroup_memory_files():
        try:
            value = Path(path).read_text().strip()
        except OSError:
            continue
        if value.isdigit() and int(value) < 1 << 60:  # "max" or a huge number: no limit
            limits.append(int(value))
    return min(limits) if limits else None


def total_memory() -> int | None:
    """The memory this process may use, in bytes: the physical memory, or a control group's limit."""
    try:
        total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):  # pragma: no cover
        return None
    limit = cgroup_memory_limit()
    return min(total, limit) if limit else total


def default_memory_gb() -> int:
    """85% of the usable memory, in whole GB, as KMC and the assemblers expect it."""
    total = total_memory()
    if total is None:  # pragma: no cover
        return 8
    return max(1, int(total * MEMORY_FRACTION / 1e9))
