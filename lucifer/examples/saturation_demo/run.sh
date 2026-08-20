#!/bin/bash
#
# The saturation demo, thin runner: every input is a real file in THIS directory
# (Genesis decks, Bmad namelists, the unaveraged wrapper lattice) -- read them, edit
# them, rerun. One practical SASE case -- Genesis's own Benchmark1-SASE configuration,
# the full 57 m Aramis line, dark start, 96 slices x 2048 particles -- tracked three
# ways from IDENTICAL initial dumps, each on the machine's full performance-core count:
# Genesis4 (MPI), the Bmad averaged mode (the bmad_standard default), and the Bmad
# unaveraged mode (~32x the averaged cost; the progress lines show it working).
#
# Methodology: an untimed Genesis prep run (same ranks, same seed, no &track) writes
# the shared initial dumps; wall times come from one external clock (/usr/bin/time -p);
# check_agreement.py must pass before any timing or report is produced. The window is
# fixed at 96 slices (sample = 3), divisible by any reasonable core count, so the input
# decks stay static. Measured results: the "The saturation demo" section of
# bsim/fel/README.md.
#
# Usage:  ./run.sh [--workers N] [--genesis <path>] [--exe <path>] [--mpirun <path>]
#                  [--python <path>] [--work-dir <path>]
# Outputs land in ./output (or --work-dir); the summary report: output/sat-demo-report.pdf.

set -o pipefail

DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BMAD_ROOT="$(cd "$DEMO_DIR/../../../.." && pwd)"
SCRIPTS="$BMAD_ROOT/bsim/fel/tests/scripts"

GENESIS="$HOME/Code/GitHub/Genesis-1.3-Version4/build-metal/genesis4"
EXE="";  PYTHON="";  MPIRUN="";  WORKERS=""
WORK_DIR="$DEMO_DIR/output"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genesis)  GENESIS="$2";  shift 2 ;;
    --exe)      EXE="$2";      shift 2 ;;
    --python)   PYTHON="$2";   shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --workers)  WORKERS="$2";  shift 2 ;;
    --mpirun)   MPIRUN="$2";   shift 2 ;;
    -h|--help)  sed -n '2,21p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$WORKERS" ]]; then
  [[ "$(uname)" == "Darwin" ]] && WORKERS="$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null)"
  [[ -z "$WORKERS" ]] && WORKERS="$(nproc 2>/dev/null)"
fi
[[ -z "$WORKERS" || "$WORKERS" -lt 2 ]] && { echo "Error: pass --workers N." >&2; exit 1; }
if (( 96 % WORKERS != 0 )); then
  echo "note: 96 slices not divisible by $WORKERS workers; Genesis will pad its window."
fi

[[ -x "$GENESIS" ]] || { echo "Error: genesis4 not found at $GENESIS (use --genesis)." >&2; exit 1; }

# rpath-matched launcher (a PATH mpirun from the wrong MPI aborts the binary on sight).
if [[ -z "$MPIRUN" && "$(uname)" == "Darwin" ]]; then
  while read -r rpath; do
    [[ -x "$rpath/../bin/mpiexec" ]] && { MPIRUN="$rpath/../bin/mpiexec"; break; }
  done < <(otool -l "$GENESIS" 2>/dev/null | awk '/LC_RPATH/ {found=1} found && /path / {print $2; found=0}')
fi
[[ -z "$MPIRUN" ]] && MPIRUN="$(command -v mpiexec || command -v mpirun)"
[[ -x "$MPIRUN" ]] || { echo "Error: no MPI launcher (use --mpirun)." >&2; exit 1; }

if [[ -z "$EXE" ]]; then
  for c in "$BMAD_ROOT/production/bin/fel_track_test" "$BMAD_ROOT/debug/bin/fel_track_test"; do
    [[ -x "$c" ]] && { EXE="$c"; break; }
  done
fi
[[ -x "$EXE" ]] || { echo "Error: fel_track_test not found; build with util/conda_compile." >&2; exit 1; }
[[ "$EXE" == *"/debug/"* ]] && echo "WARNING: debug fel_track_test; timings not comparable to a release Genesis."

if [[ -z "$PYTHON" ]]; then
  for c in "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
           "$HOME/Code/miniforge3/envs/bmad-fel-validate/bin/python3"; do
    [[ -x "$c" ]] && { PYTHON="$c"; break; }
  done
fi
[[ -x "$PYTHON" ]] || { echo "Error: python (numpy+h5py+matplotlib) not found." >&2; exit 1; }

mkdir -p "$WORK_DIR"
cp "$DEMO_DIR"/Aramis.lat "$DEMO_DIR"/aramis.bmad "$DEMO_DIR"/sat_unavg.bmad \
   "$DEMO_DIR"/sat-prep.in "$DEMO_DIR"/sat-genesis.in "$DEMO_DIR"/sat-avg.nml \
   "$DEMO_DIR"/sat-unavg.nml "$WORK_DIR/"
cd "$WORK_DIR" || exit 1
export FI_PROVIDER=tcp

echo "== saturation demo: 96 slices x 2048 particles, 57 m, $WORKERS cores each =="
echo "-- prep: Genesis writes the shared initial dumps (untimed, no tracking)"
"$MPIRUN" -np "$WORKERS" "$GENESIS" sat-prep.in > sat-prep.log 2>&1 || { tail -10 sat-prep.log >&2; exit 1; }

declare -a TIMES
timed () {  # <label> <log> <command...>
  echo "-- $1"
  /usr/bin/time -p "${@:3}" > "$2" 2> "$2.time" || { tail -10 "$2" >&2; exit 1; }
  local wall;  wall="$(awk '/^real/ {print $2}' "$2.time")"
  echo "   wall: ${wall}s"
  TIMES+=("$wall")
}

timed "Genesis4, $WORKERS MPI ranks" sat-genesis.log "$MPIRUN" -np "$WORKERS" "$GENESIS" sat-genesis.in
timed "Bmad averaged, $WORKERS threads" sat-avg.log env OMP_NUM_THREADS="$WORKERS" "$EXE" sat-avg.nml
timed "Bmad unaveraged, $WORKERS threads (~32x; watch its progress lines in sat-unavg.log)" \
      sat-unavg.log env OMP_NUM_THREADS="$WORKERS" "$EXE" sat-unavg.nml
echo

"$PYTHON" "$DEMO_DIR/check_agreement.py" AramisSat.out.h5 sat-avg.diag.txt sat-unavg.diag.txt || exit 1

"$PYTHON" "$SCRIPTS/report_fel_saturation.py" "$WORK_DIR" || exit 1

echo
echo "Outputs in: $WORK_DIR (report: sat-demo-report.pdf; per-run stats: sat-*.stats.h5)"
