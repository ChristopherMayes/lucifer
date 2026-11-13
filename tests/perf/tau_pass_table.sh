#!/bin/bash
#
# tau_pass_table.sh measures the GPU time of the Metal backend's transform along tau
# under the spectral propagator, pass by pass, with run.device_timing = true, on three
# windows at one field plane a slice: the 24-slice and 96-slice SASE windows of the device
# table and the 1280-slice one, each one Aramis segment at ngrid 256 with comb_ds_save =
# -1. Each input file runs on the device three times and the table takes the run with
# the least GPU time. A row gives the transform's length (the window's planes and the
# pad's), the lines (transverse grid points) one threadgroup transforms together, the
# walk, the GPU time the command buffers record, the transform's own GPU time and its
# encoders, its share of the sum over all passes, and the passes and transforms a step.
# The timing instrument adds a timestamped encoder boundary at every pass, which the
# production path does not have, so the seconds are longer than a production run's
# (docs/performance.md, perf-device).
#
# Usage:
#   ./tau_pass_table.sh [--exe <path>] [--workers N] [--work-dir <path>] [--python <path>]
#
# The production lucifer is preferred, found as in run_fel_validation.sh. The work
# directory is kept. The input files are the *.toml files beside this script, and the
# lattice comes from tests/bmad. A build with no device backend refuses the device runs,
# which this script reports and stops on.
set -o pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../../util/locations.sh
source "$SCRIPT_DIR/../../util/locations.sh"
EXE=""; WORKERS=""; WORK_DIR=""; PYTHON=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --exe)       EXE="$2";      shift 2 ;;
    --workers)   WORKERS="$2";  shift 2 ;;
    --work-dir)  WORK_DIR="$2"; shift 2 ;;
    --python)    PYTHON="$2";   shift 2 ;;
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
if [[ -z "$PYTHON" ]]; then
  for candidate in "${LUCIFER_PYTHON:-}" \
                   "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
                   "$(command -v python3 2>/dev/null)"; do
    if [[ -x "$candidate" ]]; then PYTHON="$candidate"; break; fi
  done
fi
if [[ -z "$PYTHON" || ! -x "$PYTHON" ]]; then
  echo "Error: no python3 with the lucifer package found (use --python)." >&2; exit 1
fi
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

# value_of <log> <group/name> ... prints values of the stats file the run of that log
# wrote, <root>.stats.h5 beside <root>.log.
value_of () {
  local stats="${1%.log}.stats.h5"; shift
  "$PYTHON" -m lucifer.stats "$stats" "$@"
}
# walk_of <log> prints the walk's seconds, and busy_of <log> the device-busy seconds.
walk_of () { value_of "$1" timing/walk; }
busy_of () { value_of "$1" timing/device_busy; }
# steps_of <log> prints the steps the device encoded, pass_of <log> <name> a kind of
# pass's seconds and encoders, and six_of <log> the sum over the kinds of pass.
steps_of () { value_of "$1" timing/device_steps; }
pass_of () { value_of "$1" "timing/device_pass_$2" "timing/device_pass_$2_encoders"; }
six_of () {
  value_of "$1" timing/device_pass_{transverse,push,zero,deposit,filter,solve,quiver,kick,spont,tau} | \
    awk '{s = 0; for (i = 1; i <= NF; i++) s += $i; print s}'
}
# planes_of <log> prints the window's planes with the pad's, and lines_of <log> the
# lines a threadgroup transforms along tau.
planes_of () { value_of "$1" setup/n_field_plane setup/n_pad_free setup/n_pad_sponge | awk '{print $1 + $2 + $3}'; }
lines_of () { value_of "$1" setup/device_tau_lines; }

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
