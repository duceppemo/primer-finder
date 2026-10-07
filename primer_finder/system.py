"""How many CPUs and how much memory this process may use."""

from __future__ import annotations

import os
from pathlib import Path

CGROUP_MEMORY_FILES = ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes")  # v2, v1
MEMORY_FRACTION = 0.85  # Of the total memory, as the previous versions did


def usable_cpus() -> int:
    """The CPUs primer-finder may use: all of them, or those a job scheduler or `taskset` gave it."""
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0)) or 1
    return os.cpu_count() or 1  # pragma: no cover - macOS


def cgroup_memory_limit(files: tuple[str, ...] = CGROUP_MEMORY_FILES) -> int | None:
    """The memory limit of this process's control group (a job scheduler's or a container's), in bytes."""
    for path in files:
        try:
            value = Path(path).read_text().strip()
        except OSError:
            continue
        if value.isdigit() and int(value) < 1 << 60:  # "max" or a huge number: no limit
            return int(value)
    return None


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
