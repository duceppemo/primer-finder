"""primer-finder: find group-specific kmers to design selective qPCR assays."""

__version__ = "1.7.0"
__author__ = "Marc-Olivier Duceppe"


class PrimerFinderError(Exception):
    """A user-facing error: bad input, a missing program, or a failed external command."""


class NothingSelectiveEnough(PrimerFinderError):
    """No region met the rule the run was asking for.

    Its own class because it is the one error worth answering by asking for less rather than by stopping:
    `find` catches it to try the last-resort pass. A failed program (`tools.ToolError`) is also a
    `PrimerFinderError` and must not be mistaken for it.
    """
