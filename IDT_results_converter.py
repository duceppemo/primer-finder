#!/usr/bin/env python3
"""Deprecated: use `primer-finder idt <sheet> <output.fasta> [prefix]` instead.

Kept so that `python IDT_results_converter.py sheet.xlsx assays.fasta prefix` keeps working.
"""

import sys

from primer_finder.cli import main

if __name__ == "__main__":
    print(__doc__.splitlines()[0], file=sys.stderr)
    sys.exit(main(["idt", *sys.argv[1:]]))
