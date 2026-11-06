#!/bin/bash
#
# unavg_spectral_table.sh measures the unaveraged mode's record step under the spectral
# propagator against the slice propagator, on the CPU and on the Metal backend, on the
# sandwich probe of the validation suite: an averaged segment, a pipe, an unaveraged
# segment, a pipe and an averaged segment (tests/bmad/unavg_sandwich.bmad), over 8 slices
# of one wavelength with 512 macroparticles a slice at ngrid 64, the input file
# unavg_sandwich_8x512.toml beside this script. Each of the four cases runs three times
# and the table takes the best walk of each, with the GPU time the backend reports from
# its command-buffer timestamps beside the device walks.
#
# Usage:
#   ./unavg_spectral_table.sh [--exe <path>] [--workers N] [--work-dir <path>]
#
# The production lucifer is preferred, found as in run_fel_validation.sh. The work
# directory is kept. A build with no device backend refuses the device runs, which this
# script reports and stops on.
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
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/unavg_spectral_table.XXXXXX")"
fi
mkdir -p "$WORK_DIR"
cp "$SCRIPT_DIR/../bmad/unavg_sandwich.bmad" "$WORK_DIR/"
cd "$WORK_DIR" || exit 1

# walk_of <log> prints the walk row's seconds from the footer's Timing block.
walk_of () {
  awk '/^ Timing/ {t = 1; next} t && $1 == "walk" {print $3; exit}' "$1"
}
# busy_of <log> prints the GPU seconds the backend reports.
busy_of () {
  sed -n 's/.*device busy *\([0-9.]*\) s.*/\1/p' "$1" | head -n 1
}

echo "=============================================================================="
echo "The unaveraged sandwich under the two propagators: best of three walks"
echo "=============================================================================="
echo "  lucifer:  $EXE"
echo "  workers:  $WORKERS (CPU threads)"
echo "  workdir:  $WORK_DIR"
echo
status=0

# measure <side> <propagator> runs one case three times and prints its row. The
# propagator and the device go on the command line after the input file.
measure () {
  local side="$1" prop="$2" tag="${1}_${2}" best="" busy="" w rep log
  local -a over=()
  [[ "$prop" == spectral ]] && over+=('field.propagator="spectral"' 'field.n_wavelength=1')
  [[ "$side" == device ]] && over+=('run.device="metal"')
  for rep in 1 2 3; do
    log="${tag}$rep.log"
    if ! env OMP_NUM_THREADS="$WORKERS" "$EXE" "$SCRIPT_DIR/unavg_sandwich_8x512.toml" \
         "run.out_root=${tag}$rep" "${over[@]}" > "$log" 2>&1; then
      echo "  $tag: the run failed; see $WORK_DIR/$log" >&2; status=1; continue
    fi
    w="$(walk_of "$log")"
    if [[ -z "$best" ]] || awk -v a="$w" -v b="$best" 'BEGIN { exit !(a < b) }'; then
      best="$w"; busy="$(busy_of "$log")"
    fi
  done
  printf '  %-10s %-10s %10s s %12s s\n' "$side" "$prop" "$best" "${busy:--}"
}

printf '  %-10s %-10s %12s %14s\n' "side" "propagator" "walk" "GPU time"
for side in CPU device; do
  for prop in slices spectral; do
    measure "$side" "$prop"
  done
done
echo
echo "Outputs kept in: $WORK_DIR"
exit $status
