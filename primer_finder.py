#!/usr/bin/env python3
"""Backward-compatible entry point: `python primer_finder.py ...` is the same as `primer-finder ...`."""

import sys

from primer_finder.cli import main

if __name__ == "__main__":
    sys.exit(main())
