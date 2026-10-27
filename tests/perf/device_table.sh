#!/bin/bash
#
# device_table.sh measures the Metal backend against the CPU path on the four input files of the
# device table in docs/performance.md, one Aramis segment at ngrid 256 with comb_ds_save
# = -1, so both sides pay stats at the element end only. Each input file runs three times on
# the CPU at the machine's performance-core count and three times on the device, and the
# table takes the best walk of each side and the device-busy seconds the backend reports
# from its command-buffer timestamps.
#
# Usage:
#   ./device_table.sh [--exe <path>] [--workers N] [--work-dir <path>]
#
# The production lucifer is preferred, found as in run_fel_validation.sh. The work
# directory is kept. The input files are the *.toml files beside this script, and the lattice
# comes from tests/bmad. A build with no device backend refuses the device runs, which
# this script reports and stops on.
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
    -h|--help)   sed -n '2,17p' "${BASH_SOURCE[0]}"; exit 0 ;;
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
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/device_table.XXXXXX")"
fi
mkdir -p "$WORK_DIR"
cp "$SCRIPT_DIR/../bmad/aramis_1seg.bmad" "$WORK_DIR/"
cd "$WORK_DIR" || exit 1

CASES="steady_1x8192 sase_24x4096 sase_96x8192 floor_1x1024"

# walk_of <log> prints the walk row's seconds from the footer's Timing block.
walk_of () {
  awk '/^ Timing/ {t = 1; next} t && $1 == "walk" {print $3; exit}' "$1"
}
# busy_of <log> prints the device-busy seconds the backend reports.
busy_of () {
  sed -n 's/.*device busy *\([0-9.]*\) s.*/\1/p' "$1" | head -n 1
}

echo "=============================================================================="
echo "The device against the CPU: best of three walks on each side"
echo "=============================================================================="
echo "  lucifer:  $EXE"
echo "  workers:  $WORKERS (CPU threads)"
echo "  workdir:  $WORK_DIR"
echo
printf '  %-28s %14s %14s %14s %8s\n' "case" "CPU, ${WORKERS} thr" "device wall" "device busy" "ratio"
status=0
for c in $CASES; do
  best_cpu=""; best_dev=""; busy=""
  for rep in 1 2 3; do
    sed "s/^out_root = \"$c\"/out_root = \"${c}_cpu$rep\"/" "$SCRIPT_DIR/$c.toml" > "${c}_cpu$rep.toml"
    if ! env OMP_NUM_THREADS="$WORKERS" "$EXE" "${c}_cpu$rep.toml" > "${c}_cpu$rep.log" 2>&1; then
      echo "  $c: the CPU run failed; see $WORK_DIR/${c}_cpu$rep.log" >&2; status=1; continue
    fi
    w="$(walk_of "${c}_cpu$rep.log")"
    if [[ -z "$best_cpu" ]] || awk -v a="$w" -v b="$best_cpu" 'BEGIN { exit !(a < b) }'; then best_cpu="$w"; fi
  done
  for rep in 1 2 3; do
    # The device switch goes into [run], right after the table's heading.
    sed "s/^out_root = \"$c\"/out_root = \"${c}_dev$rep\"/; s/^\[run\]$/[run]\ndevice = \"metal\"/" \
      "$SCRIPT_DIR/$c.toml" > "${c}_dev$rep.toml"
    if ! env OMP_NUM_THREADS="$WORKERS" "$EXE" "${c}_dev$rep.toml" > "${c}_dev$rep.log" 2>&1; then
      echo "  $c: the device run failed; see $WORK_DIR/${c}_dev$rep.log" >&2; status=1; continue
    fi
    w="$(walk_of "${c}_dev$rep.log")"
    if [[ -z "$best_dev" ]] || awk -v a="$w" -v b="$best_dev" 'BEGIN { exit !(a < b) }'; then
      best_dev="$w"; busy="$(busy_of "${c}_dev$rep.log")"
    fi
  done
  if [[ -n "$best_cpu" && -n "$best_dev" ]]; then
    ratio="$(awk -v a="$best_cpu" -v b="$best_dev" 'BEGIN { printf "%.1f", a / b }')"
    printf '  %-28s %12s s %12s s %12s s %7sx\n' "$c" "$best_cpu" "$best_dev" "${busy:-?}" "$ratio"
  fi
done
echo
echo "Outputs kept in: $WORK_DIR"
exit $status
