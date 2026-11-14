#!/bin/bash
#
# flash1_cell_table.sh times Lucifer against Genesis4 on one cell of the FLASH1 example:
# one 165-period planar segment and its 0.8 m intersection, 5.3045 m, in steps of one
# period. Both codes start from the same particles and the same dark field, which
# Genesis4 makes once (flash1_cell/prep.in): 528 slices of 12 wavelengths, 8192 particles
# a slice, a flat 2500 A, on a 256-point grid. No wake, no space charge, the source
# filter off, and both codes write their final beam and field.
#
# The runs, one after another and each with the machine to itself, in rounds:
#   Genesis4 with the ADI field solver and with the FFT one, at one MPI rank a worker.
#   Lucifer on the CPU at one thread a worker, under run.fft_plan = wisdom and estimate,
#   each with Genesis4's transport and interlude models and with its own.
# The table gives each run's best wall time, as /usr/bin/time measures it, and for
# Lucifer the field solve and the particle push of that run (docs/performance.md,
# perf-flash1-cell).
#
# Usage:
#   ./flash1_cell_table.sh [--exe <path>] [--genesis <path>] [--mpirun <path>]
#                          [--workers N] [--rounds N] [--work-dir <path>] [--python <path>]
#
# The production lucifer is preferred, found as in run_fel_validation.sh. genesis4 and
# its MPI launcher default to those beside the Python interpreter of the bmad-fel-validate
# environment. The work directory is kept, about 4 GB.
set -o pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../../util/locations.sh
source "$SCRIPT_DIR/../../util/locations.sh"
EXE=""; GENESIS="${GENESIS4:-}"; MPIRUN=""; WORKERS=""; ROUNDS=3; WORK_DIR=""; PYTHON=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --exe)       EXE="$2";      shift 2 ;;
    --genesis)   GENESIS="$2";  shift 2 ;;
    --mpirun)    MPIRUN="$2";   shift 2 ;;
    --workers)   WORKERS="$2";  shift 2 ;;
    --rounds)    ROUNDS="$2";   shift 2 ;;
    --work-dir)  WORK_DIR="$2"; shift 2 ;;
    --python)    PYTHON="$2";   shift 2 ;;
    -h|--help)   sed -n '2,24p' "${BASH_SOURCE[0]}"; exit 0 ;;
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
[[ -n "$GENESIS" ]] || GENESIS="$(dirname "$PYTHON")/genesis4"
[[ -n "$MPIRUN" ]] || MPIRUN="$(dirname "$GENESIS")/mpirun"
for tool in "$GENESIS" "$MPIRUN"; do
  [[ -x "$tool" ]] || { echo "Error: not found: $tool (use --genesis and --mpirun)." >&2; exit 1; }
done
if [[ -z "$WORKERS" ]]; then
  WORKERS="$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || nproc 2>/dev/null)"
fi
if [[ -z "$WORK_DIR" ]]; then
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/lucifer-flash1-cell.XXXXXX")"
fi
mkdir -p "$WORK_DIR" && cd "$WORK_DIR" || exit 1
cp "$SCRIPT_DIR"/flash1_cell/* .
export OMP_NUM_THREADS="$WORKERS"

echo "lucifer:  $EXE"
echo "genesis4: $GENESIS"
echo "workers:  $WORKERS    rounds: $ROUNDS    work directory: $WORK_DIR"

# The shared start, made once and not timed. The slice count divides over the ranks.
if [[ ! -f flash1-start.beam.h5 || ! -f zero-start.wf.h5 ]]; then
  "$MPIRUN" -n "$WORKERS" "$GENESIS" prep.in > prep.log 2>&1 \
    || { echo "Error: the Genesis4 run that makes the start failed, see $WORK_DIR/prep.log" >&2; exit 1; }
  convert="$SCRIPT_DIR/../scripts/convert_genesis.py"
  "$PYTHON" "$convert" to-openpmd flash1-start.par.h5 flash1-start.beam.h5 > convert.log 2>&1 \
    && "$PYTHON" "$convert" to-openpmd zero-start.fld.h5 zero-start.wf.h5 >> convert.log 2>&1 \
    || { echo "Error: the conversion to openPMD failed, see $WORK_DIR/convert.log" >&2; exit 1; }
fi

timed () {   # <name> <command...>: the wall seconds to <name>.wall, appended
  local name="$1"; shift
  /usr/bin/time -p "$@" > "$name.log" 2> "$name.time" \
    || { echo "Error: $name failed, see $WORK_DIR/$name.log" >&2; exit 1; }
  awk '/^real/ {print $2}' "$name.time" >> "$name.wall"
  echo "  $name  $(tail -n 1 "$name.wall") s"
}

rm -f ./*.wall
for round in $(seq "$ROUNDS"); do
  echo "round $round"
  timed g4_adi "$MPIRUN" -n "$WORKERS" "$GENESIS" genesis_adi.in
  timed g4_fft "$MPIRUN" -n "$WORKERS" "$GENESIS" genesis_fft.in
  for model in genesis bmad; do
    for plan in wisdom estimate; do
      name="luc_${model}_${plan}"
      timed "$name" "$EXE" lucifer.toml run.out_root="$name" run.transport_model="$model" \
            run.interlude_model="$model" run.fft_plan="$plan"
      # The phases of the fastest round so far, from the stats file of that run.
      if [[ "$(sort -n "$name.wall" | head -n 1)" == "$(tail -n 1 "$name.wall")" ]]; then
        "$PYTHON" - "$name" > "$name.phases" <<'PYEOF'
import sys
import h5py
with h5py.File(sys.argv[1] + ".stats.h5") as h:
    get = lambda k: float(h["timing/" + k][()].ravel()[0])
    plan = h["setup/fft_plan"][()]
    print(f"{get('field_solve'):.1f} s | {get('particle_push'):.1f} s | "
          f"{plan.decode() if isinstance(plan, bytes) else plan}")
PYEOF
      fi
    done
  done
done

best () { sort -n "$1.wall" | head -n 1; }
echo
echo "| run | wall | field solve | particle push | plan in use |"
echo "|---|---|---|---|---|"
echo "| Genesis4, ADI solver, $WORKERS ranks | $(best g4_adi) s | | | |"
echo "| Genesis4, FFT solver, $WORKERS ranks | $(best g4_fft) s | | | |"
for model in genesis bmad; do
  [[ "$model" == genesis ]] && label="Genesis4's models" || label="its own models"
  for plan in wisdom estimate; do
    name="luc_${model}_${plan}"
    echo "| Lucifer, $label, fft_plan = $plan, $WORKERS threads | $(best "$name") s | $(cat "$name.phases") |"
  done
done
