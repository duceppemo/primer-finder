"""What this process may use: how many CPUs, how much memory, and one run per output folder."""

from __future__ import annotations

import contextlib
import fcntl
import os
import socket
import time
from collections.abc import Iterator
from pathlib import Path

from primer_finder import PrimerFinderError

# Where a control group's memory limit is, when this process sees its own cgroup as the root (containers).
CGROUP_MEMORY_FILES = ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes")  # v2, v1
# And its CPU quota, the same way round.
CGROUP_CPU_FILES = ("/sys/fs/cgroup/cpu.max", "/sys/fs/cgroup/cpu/cpu.cfs_quota_us")  # v2, v1
CGROUP_ROOT = Path("/sys/fs/cgroup")
PROC_CGROUP = Path("/proc/self/cgroup")
MEMORY_FRACTION = 0.85  # Of the total memory, as the previous versions did
LOCK_NAME = ".primer_finder.lock"


def usable_cpus() -> int:
    """The CPUs primer-finder may use: all of them, those a job scheduler or `taskset` gave it, or what a
    control group's CPU quota allows.

    `sched_getaffinity` sees a cpuset -- which CPUs the job may run on -- but not a *quota*, the share of
    them it may use, and a scheduler may set either. A quota of 400,000 microseconds per 100,000 is four
    CPUs' worth of time on however many the affinity mask holds, so handing `-t 64` to KMC and blast there
    oversubscribes the job by sixteen times rather than running it faster. The share is rounded down, and
    never below 1.
    """
    # macOS has no sched_getaffinity, and no cgroups either
    cpus = (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()) or 1
    quota = cgroup_cpu_limit()
    return max(1, min(cpus, quota)) if quota is not None else cpus


def cgroup_files(controller: str, v2_name: str, v1_name: str, roots: tuple[str, ...],
                 proc_cgroup: Path = PROC_CGROUP, root: Path = CGROUP_ROOT) -> tuple[str, ...]:
    """Every file that could hold one kind of limit for this process: the roots, then this process's own
    control group and each of its parents.

    A job scheduler (Slurm, systemd) puts the job in a control group below the root, so the limit is not in
    `/sys/fs/cgroup/memory.max` but in `/sys/fs/cgroup<path from /proc/self/cgroup>/memory.max`, and it may
    be set on a parent of it. In a container the process sees its cgroup as the root instead, which is what
    `roots` covers.
    """
    files = list(roots)
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
            base, limit = root, v2_name
        elif controller in controllers.split(","):  # cgroup v1: one line per controller
            base, limit = root / controller, v1_name
        else:
            continue
        group = base / path.strip("/") if path.strip("/") else base
        while True:
            files.append(str(group / limit))
            if group == base or base not in group.parents:
                break
            group = group.parent
    return tuple(dict.fromkeys(files))  # In order, without duplicates


def cgroup_memory_files(proc_cgroup: Path = PROC_CGROUP, root: Path = CGROUP_ROOT) -> tuple[str, ...]:
    """Every file that could hold a memory limit for this process."""
    return cgroup_files("memory", "memory.max", "memory.limit_in_bytes", CGROUP_MEMORY_FILES,
                        proc_cgroup, root)


def cgroup_cpu_files(proc_cgroup: Path = PROC_CGROUP, root: Path = CGROUP_ROOT) -> tuple[str, ...]:
    """Every file that could hold a CPU quota for this process."""
    return cgroup_files("cpu", "cpu.max", "cpu.cfs_quota_us", CGROUP_CPU_FILES, proc_cgroup, root)


def cgroup_cpu_limit(files: tuple[str, ...] | None = None) -> int | None:
    """The smallest CPU quota any of this process's control groups sets, as whole CPUs, or None.

    cgroup v2 writes the quota and the period on one line of `cpu.max`, as "max 100000" where there is no
    limit; v1 writes the quota in `cpu.cfs_quota_us`, as -1 for no limit, and the period in
    `cpu.cfs_period_us` beside it. A quota smaller than one whole CPU gives 1, since there is no such
    thing as a third of a thread.
    """
    limits = []
    for path in files if files is not None else cgroup_cpu_files():
        file = Path(path)
        try:
            value = file.read_text().split()
            quota = int(value[0])  # "max" is not a number, and a group with no quota says "max" or -1
            period = (int(value[1]) if len(value) > 1  # v2 keeps both on the one line
                      else int((file.parent / "cpu.cfs_period_us").read_text()))  # v1, in the file beside
        except (OSError, ValueError, IndexError):
            continue
        if quota > 0 and period > 0:
            limits.append(max(1, quota // period))
    return min(limits) if limits else None


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
    # KMC refuses anything under 2 GB, so a small machine or a tight cgroup must not produce 1. The floor
    # is kmers.MIN_MEMORY_GB; it is written out here rather than imported to keep this module free of the
    # pipeline, and a test holds the two together.
    return max(2, int(total * MEMORY_FRACTION / 1e9))


@contextlib.contextmanager
def lock_output(folder: Path, what: str = "run") -> Iterator[None]:
    """Claim `folder` for this run, and hold the claim for as long as the block runs.

    Two runs sharing an output folder write over each other: they count kmers into the same
    `1_kmers/inclusion.kmc_suf`, assemble into the same contigs file and blast against the same database,
    so each one reads some of the other's work and the answer belongs to neither. It is an easy mistake to
    make -- one terminal per group, the same `-o` -- and nothing said so. The claim is an `flock` on a file
    in the folder, which the kernel drops when the process ends however it ends, so a crashed run leaves
    no stale lock to clear out by hand.
    """
    folder.mkdir(parents=True, exist_ok=True)
    handle = os.open(folder / LOCK_NAME, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            held = ""
            with contextlib.suppress(OSError):  # pragma: no cover - the file was just opened
                held = os.pread(handle, 200, 0).decode(errors="replace").strip()
            raise PrimerFinderError(
                f"Another primer-finder {what} is using {folder}"
                + (f" ({held})" if held else "")
                + ". Two runs sharing an output folder overwrite each other's kmer databases, contigs and "
                  "blast databases, so this one is refused rather than left to produce an answer that "
                  "belongs to neither. Give it an output folder of its own."
            ) from error
        os.truncate(handle, 0)
        os.pwrite(handle, f"pid {os.getpid()} on {socket.gethostname()}, started "
                          f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n".encode(), 0)
        yield
    finally:
        os.close(handle)  # Which releases the lock, so a second run in this process may take it
