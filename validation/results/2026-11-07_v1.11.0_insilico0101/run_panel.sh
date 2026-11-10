#!/usr/bin/env bash
# Design the Xylella panel with primer-finder 1.11.0 + insilicoPCR 0.10.1, on the SAME find results the
# 0.8.0 arms of ../2026-11-03_v1.10.0_insilico080/ and ../2026-11-03_3prime_dependence/ used (that
# record's run_panel.sh builds the work folder), so the two insilicoPCR releases can be joined per assay.
# The four targets run side by side, 8 threads each.
set -uo pipefail

W=/tmp/claude-1000/-home-bioinfo-prog-primer-finder/1df804a1-6989-4610-bef1-ce42c6758ed4/scratchpad/xylella-v1.10.0
O=/tmp/claude-1000/-home-bioinfo-prog-primer-finder/fe284be5-8be0-4894-9428-65d4739733cf/scratchpad/xylella-0101
ENVBIN=/home/bioinfo/miniconda3/envs/primer-finder_dev/bin
export PATH=$ENVBIN:$PATH
I0101=/home/bioinfo/prog/insilicoPCR-linux-x64   # the v0.10.1 linux-x64 release, sha256-checked

mkdir -p "$O"
for t in multiplex fastidiosa pauca morus; do
    ( $ENVBIN/python /home/bioinfo/prog/primer-finder/primer_finder.py design "$W/$t/results" -o "$O/$t" \
          -t 8 --insilico-pcr "$I0101" > "$O/$t.log" 2>&1
      echo "$t exit $? $(date +%T)" ) &
done
wait
