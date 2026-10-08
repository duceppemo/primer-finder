"""primer-finder: find group-specific kmers to design selective qPCR assays."""

__version__ = "1.4.0"
__author__ = "Marc-Olivier Duceppe"


class PrimerFinderError(Exception):
    """A user-facing error: bad input, a missing program, or a failed external command."""
