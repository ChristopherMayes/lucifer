#!/bin/bash
#
# device_table.sh measures the Metal backend against the CPU path on the four input files of the
# device table in docs/performance.md, one Aramis segment at ngrid 256 with comb_ds_save
# = -1, so both sides pay stats at the element end only. Each input file runs three times on
# the CPU at the machine's performance-core count and three times on the device, and the
# table takes the best walk of each side and the device-busy seconds the backend reports
# from its command-buffer timestamps. The same four run again under the spectral
# propagator at one plane a slice, and the 96-slice window once more with the field at
# its derived sampling, six slices a plane, which is the second table of that page.
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
# The spectral rows: each case under field.propagator = "spectral" with the field's
# plane the window's slice, and the 96-slice window again with the plane left to derive.
SPECTRAL_CASES="steady_1x8192:spectral sase_24x4096:spectral sase_96x8192:spectral floor_1x1024:spectral sase_96x8192:coarse"

# variant_of <toml> <variant> rewrites one case's input file for a variant: spectral adds
# the propagator and holds the field's plane at the window's slice, coarse adds the
# propagator alone, and slices leaves the file as it is.
variant_of () {
  case "$2" in
    spectral)
      if grep -q '^n_wavelength = ' "$1"; then
        nwl="$(sed -n 's/^n_wavelength = //p' "$1" | head -n 1)"
        sed -i.bak "s/^\[field\]$/[field]\npropagator = \"spectral\"\nn_wavelength = $nwl/" "$1"
      else
        sed -i.bak 's/^\[field\]$/[field]\npropagator = "spectral"/' "$1"
      fi
      rm -f "$1.bak" ;;
    coarse)
      sed -i.bak 's/^\[field\]$/[field]\npropagator = "spectral"/' "$1"; rm -f "$1.bak" ;;
  esac
}

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
status=0

# measure <case> <variant> runs one case three times a side and prints its row.
measure () {
  local c="$1" v="$2" tag="$1" best_cpu="" best_dev="" busy="" w ratio rep
  [[ "$v" == slices ]] || tag="${c}_$v"
  for rep in 1 2 3; do
    sed "s/^out_root = \"$c\"/out_root = \"${tag}_cpu$rep\"/" "$SCRIPT_DIR/$c.toml" > "${tag}_cpu$rep.toml"
    variant_of "${tag}_cpu$rep.toml" "$v"
    if ! env OMP_NUM_THREADS="$WORKERS" "$EXE" "${tag}_cpu$rep.toml" > "${tag}_cpu$rep.log" 2>&1; then
      echo "  $tag: the CPU run failed; see $WORK_DIR/${tag}_cpu$rep.log" >&2; status=1; continue
    fi
    w="$(walk_of "${tag}_cpu$rep.log")"
    if [[ -z "$best_cpu" ]] || awk -v a="$w" -v b="$best_cpu" 'BEGIN { exit !(a < b) }'; then best_cpu="$w"; fi
  done
  for rep in 1 2 3; do
    # The device switch goes into [run], right after the table's heading.
    sed "s/^out_root = \"$c\"/out_root = \"${tag}_dev$rep\"/; s/^\[run\]$/[run]\ndevice = \"metal\"/" \
      "$SCRIPT_DIR/$c.toml" > "${tag}_dev$rep.toml"
    variant_of "${tag}_dev$rep.toml" "$v"
    if ! env OMP_NUM_THREADS="$WORKERS" "$EXE" "${tag}_dev$rep.toml" > "${tag}_dev$rep.log" 2>&1; then
      echo "  $tag: the device run failed; see $WORK_DIR/${tag}_dev$rep.log" >&2; status=1; continue
    fi
    w="$(walk_of "${tag}_dev$rep.log")"
    if [[ -z "$best_dev" ]] || awk -v a="$w" -v b="$best_dev" 'BEGIN { exit !(a < b) }'; then
      best_dev="$w"; busy="$(busy_of "${tag}_dev$rep.log")"
    fi
  done
  if [[ -n "$best_cpu" && -n "$best_dev" ]]; then
    ratio="$(awk -v a="$best_cpu" -v b="$best_dev" 'BEGIN { printf "%.1f", a / b }')"
    printf '  %-28s %12s s %12s s %12s s %7sx\n' "$tag" "$best_cpu" "$best_dev" "${busy:-?}" "$ratio"
  fi
}

printf '  %-28s %14s %14s %14s %8s\n' "case" "CPU, ${WORKERS} thr" "device wall" "device busy" "ratio"
for c in $CASES; do
  measure "$c" slices
done
echo
echo "The spectral propagator, the same four at one plane a slice and the 96-slice window"
echo "with the plane derived:"
printf '  %-28s %14s %14s %14s %8s\n' "case" "CPU, ${WORKERS} thr" "device wall" "device busy" "ratio"
for cv in $SPECTRAL_CASES; do
  measure "${cv%%:*}" "${cv##*:}"
done
echo
echo "Outputs kept in: $WORK_DIR"
exit $status
