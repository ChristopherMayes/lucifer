#!/bin/bash
#
# tau_pass_table.sh prices the Metal backend's pass along tau under the spectral
# propagator, per pass, with run.device_timing = true, on three windows at one plane a
# slice: the 24-slice and 96-slice SASE windows of the device table and the 1280-slice
# one, each one Aramis segment at ngrid 256 with comb_ds_save = -1. Each input file runs
# on the device three times and the table takes the run with the fewest busy seconds. A
# row gives the window's planes along tau with its pad, the lines a threadgroup carries
# through the pass, the walk, the busy seconds the command buffers clock, the pass's own
# seconds and encoders, the pass's share of the sum over the passes, and the passes and
# transforms a step. The seconds are the instrument's, since every pass pays an encoder
# boundary that the production path does not (docs/performance.md, perf-device).
#
# Usage:
#   ./tau_pass_table.sh [--exe <path>] [--workers N] [--work-dir <path>]
#
# The production lucifer is preferred, found as in run_fel_validation.sh. The work
# directory is kept. The input files are the *.toml files beside this script, and the
# lattice comes from tests/bmad. A build with no device backend refuses the device runs,
# which this script reports and stops on.
set -o pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../../util/locations.sh
source "$SCRIPT_DIR/../../util/locations.sh"
EXE=""; WORKERS=""; WORK_DIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --exe)       EXE="$2";      shift 2 ;;
    --workers)   WORKERS="$2";  shift 2 ;;
    --work-dir)  WORK_DIR="$2"; shift 2 ;;
    -h|--help)   sed -n '2,20p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
if [[ -z "$EXE" ]]; then
  for candidate in "$PRODUCTION_BIN/lucifer" "$DEBUG_BIN/lucifer"; do
    [[ -x "$candidate" ]] && { EXE="$candidate"; break; }
  done
fi
if [[ -z "$EXE" || ! -x "$EXE" ]]; then
  echo "Error: no lucifer binary found (use --exe)." >&2; exit 1
fi
[[ "$EXE" == /* ]] || EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"
if [[ -z "$WORKERS" ]]; then
  WORKERS="$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || nproc 2>/dev/null)"
fi
if [[ -z "$WORK_DIR" ]]; then
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/tau_pass_table.XXXXXX")"
fi
mkdir -p "$WORK_DIR"
cp "$SCRIPT_DIR/../bmad/aramis_1seg.bmad" "$WORK_DIR/"
cd "$WORK_DIR" || exit 1

CASES="sase_24x4096 sase_96x8192 sase_1280x1024"

# walk_of <log> prints the walk row's seconds from the footer's Timing block.
walk_of () {
  awk '/^ Timing/ {t = 1; next} t && $1 == "walk" {print $3; exit}' "$1"
}
# busy_of <log> prints the device-busy seconds the backend reports.
busy_of () {
  sed -n 's/.*device busy *\([0-9.]*\) s.*/\1/p' "$1" | head -n 1
}
# steps_of <log> prints the steps the device encoded.
steps_of () {
  sed -n 's/.*Device: \([0-9]*\) steps encoded.*/\1/p' "$1" | head -n 1
}
# pass_of <log> <name> prints a per-pass line's seconds and encoders.
pass_of () {
  sed -n "s/^ *$2 *\([0-9.]*\) s over \([0-9]*\) encoders\./\1 \2/p" "$1" | head -n 1
}
# six_of <log> prints the sum over the passes the footer states.
six_of () {
  sed -n 's/.*the six together *\([0-9.]*\) s against.*/\1/p' "$1" | head -n 1
}
# planes_of <log> and lines_of <log> read the header's device line.
planes_of () {
  sed -n 's/.*, \([0-9]*\) window planes and \([0-9]*\) pad planes under.*/\1 \2/p' "$1" | \
    head -n 1 | awk '{print $1 + $2}'
}
lines_of () {
  sed -n 's/.*, \([0-9]*\) lines a threadgroup along tau\..*/\1/p' "$1" | head -n 1
}

echo "=============================================================================="
echo "The pass along tau on the device, per pass, with run.device_timing on"
echo "=============================================================================="
echo "  lucifer:  $EXE"
echo "  workers:  $WORKERS (CPU threads)"
echo "  workdir:  $WORK_DIR"
echo
status=0

# measure <case> runs one window three times on the device and prints its row.
measure () {
  local c="$1" best="" busy walk rep log tau six planes lines steps passes
  for rep in 1 2 3; do
    sed "s/^out_root = \"$c\"/out_root = \"${c}_tau$rep\"/; s/^\[run\]$/[run]\ndevice = \"metal\"\ndevice_timing = true/; s/^\[field\]$/[field]\npropagator = \"spectral\"\nn_wavelength = 3/" \
      "$SCRIPT_DIR/$c.toml" > "${c}_tau$rep.toml"
    log="${c}_tau$rep.log"
    if ! env OMP_NUM_THREADS="$WORKERS" "$EXE" "${c}_tau$rep.toml" > "$log" 2>&1; then
      echo "  $c: the device run failed; see $WORK_DIR/$log" >&2; status=1; continue
    fi
    busy="$(busy_of "$log")"
    if [[ -z "$best" ]] || awk -v a="$busy" -v b="$(busy_of "$best")" 'BEGIN { exit !(a < b) }'; then best="$log"; fi
  done
  [[ -n "$best" ]] || return
  walk="$(walk_of "$best")"; busy="$(busy_of "$best")"; steps="$(steps_of "$best")"
  read -r tau passes <<< "$(pass_of "$best" tau)"
  six="$(six_of "$best")"; planes="$(planes_of "$best")"; lines="$(lines_of "$best")"
  awk -v c="$c" -v p="$planes" -v l="$lines" -v w="$walk" -v b="$busy" -v t="$tau" -v e="$passes" \
      -v s="$six" -v n="$steps" \
      'BEGIN { printf "  %-16s %7d %6d %9.3f s %9.3f s %9.4f s %6d %7.1f%% %6.2f %6.2f\n", c, p, l, w, b, t, e, 100 * t / s, e / n, 2 * e / n }'
}

printf '  %-16s %7s %6s %11s %11s %12s %6s %8s %6s %6s\n' "case" "planes" "lines" "walk" "busy" "tau" "enc" "share" "pass/st" "tf/st"
for c in $CASES; do
  measure "$c"
done
echo
echo "Outputs kept in: $WORK_DIR"
exit $status
