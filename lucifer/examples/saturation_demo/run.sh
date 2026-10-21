#!/bin/bash
#
# The saturation demo, thin runner: every input is a real file in THIS directory
# (Genesis decks, Bmad namelists, the unaveraged wrapper lattice) -- read them, edit
# them, rerun. One practical SASE case -- Genesis's own Benchmark1-SASE configuration,
# the full 57 m Aramis line, dark start, 96 slices x 2048 particles -- tracked three
# ways from IDENTICAL initial dumps, each on the machine's full performance-core count:
# Genesis4 (MPI), the Bmad averaged mode (the bmad_standard default), and the Bmad
# unaveraged mode (~32x the averaged cost, the progress lines show it working).
#
# Methodology: an untimed Genesis prep run (same ranks, same seed, no &track) writes
# the shared initial dumps. Wall times come from one external clock (/usr/bin/time -p).
# check_agreement.py must pass before any timing or report is produced. The window is
# fixed at 96 slices (sample = 3), divisible by any reasonable core count, so the input
# decks stay static. Measured results: the "The saturation demo" section of
# lucifer/examples/saturation_demo/README.md.
#
# Usage:  ./run.sh [--workers N] [--genesis <path>] [--exe <path>] [--mpirun <path>]
#                  [--python <path>] [--work-dir <path>] [--show-locations]
# Outputs land in ./output (or --work-dir). The summary report: output/sat-demo-report.pdf.
# --show-locations prints every location the run would use, absolute, checks from a
# scratch directory that the executable, the interpreter and the helper scripts are
# there, and exits before anything runs.

set -o pipefail

DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# The locations, LUCIFER_DIR, DEBUG_BIN, PRODUCTION_BIN and BMAD_SRC, with their
# defaults and the variables that name them elsewhere (util/locations.sh).
source "$DEMO_DIR/../../util/locations.sh"
SCRIPTS="$LUCIFER_DIR/tests/scripts"

# genesis4 comes from conda-forge, in this project's own environment. PATH is not
# searched first on purpose: an unrelated environment on PATH would supply a different
# build, and the comparison levels are tied to the build that produced them.
GENESIS="${GENESIS4:-}"
if [[ -z "$GENESIS" ]]; then
  for candidate in "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/genesis4" \
                   "$(command -v genesis4 2>/dev/null)"; do
    if [[ -x "$candidate" ]]; then GENESIS="$candidate"; break; fi
  done
fi
EXE="";  PYTHON="";  MPIRUN="";  WORKERS="";  SHOW=0
WORK_DIR="$DEMO_DIR/output"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genesis)  GENESIS="$2";  shift 2 ;;
    --exe)      EXE="$2";      shift 2 ;;
    --python)   PYTHON="$2";   shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --workers)  WORKERS="$2";  shift 2 ;;
    --mpirun)   MPIRUN="$2";   shift 2 ;;
    --show-locations) SHOW=1;  shift 1 ;;
    -h|--help)  sed -n '2,24p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$WORKERS" ]]; then
  [[ "$(uname)" == "Darwin" ]] && WORKERS="$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null)"
  [[ -z "$WORKERS" ]] && WORKERS="$(nproc 2>/dev/null)"
fi
[[ -z "$WORKERS" || "$WORKERS" -lt 2 ]] && { echo "Error: pass --workers N." >&2; exit 1; }
if (( 96 % WORKERS != 0 )); then
  echo "note: 96 slices not divisible by $WORKERS workers. Genesis will pad its window."
fi

# Every path is made absolute here, against the directory the script was started in,
# because the run changes into the work directory before it launches anything.
[[ -n "$GENESIS" ]] && GENESIS="$(lucifer_abs "$GENESIS")"
[[ -n "$MPIRUN" ]] && MPIRUN="$(lucifer_abs "$MPIRUN")"
[[ -n "$EXE" ]] && EXE="$(lucifer_abs "$EXE")"
[[ -n "$PYTHON" ]] && PYTHON="$(lucifer_abs "$PYTHON")"
WORK_DIR="$(lucifer_abs "$WORK_DIR")"

[[ $SHOW -eq 1 || -x "$GENESIS" ]] || { echo "Error: genesis4 not found at $GENESIS (use --genesis)." >&2; exit 1; }

# The launcher of the MPI genesis4 was built against: a PATH mpirun from another MPI
# aborts the binary on sight, and an unrelated environment on PATH supplies exactly that.
# A conda genesis4 has its environment's launcher beside it, as mpirun with or without an
# mpiexec, so that is looked for first. Then the binary's LC_RPATH, with conda's
# @loader_path resolved against the binary's own directory, and PATH last.
if [[ -z "$MPIRUN" && -x "$GENESIS" ]]; then
  gdir="$(dirname "$GENESIS")"
  for c in "$gdir/mpiexec" "$gdir/mpirun"; do
    [[ -x "$c" ]] && { MPIRUN="$c"; break; }
  done
fi
if [[ -z "$MPIRUN" && -x "$GENESIS" && "$(uname)" == "Darwin" ]]; then
  while read -r rpath; do
    rpath="${rpath/@loader_path/$(dirname "$GENESIS")}"
    for c in "$rpath/../bin/mpiexec" "$rpath/../bin/mpirun"; do
      [[ -x "$c" ]] && { MPIRUN="$c"; break 2; }
    done
  done < <(otool -l "$GENESIS" 2>/dev/null | awk '/LC_RPATH/ {found=1} found && /path / {print $2; found=0}')
fi
[[ -z "$MPIRUN" ]] && MPIRUN="$(command -v mpiexec || command -v mpirun)"
[[ $SHOW -eq 1 || -x "$MPIRUN" ]] || { echo "Error: no MPI launcher (use --mpirun)." >&2; exit 1; }

if [[ -z "$EXE" ]]; then
  for c in "$PRODUCTION_BIN/lucifer" "$DEBUG_BIN/lucifer"; do
    [[ -x "$c" ]] && { EXE="$c"; break; }
  done
fi
[[ -x "$EXE" ]] || { echo "Error: lucifer not found. Build with lucifer/util/conda_compile." >&2; exit 1; }
[[ "$EXE" == "$DEBUG_BIN/"* ]] && echo "WARNING: debug lucifer. Timings are not comparable to a release Genesis."

if [[ -z "$PYTHON" ]]; then
  for c in "${LUCIFER_PYTHON:-}" \
           "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
           "$(command -v python3 2>/dev/null)"; do
    [[ -x "$c" ]] && { PYTHON="$c"; break; }
  done
fi
[[ -x "$PYTHON" ]] || { echo "Error: python (numpy+h5py+matplotlib) not found." >&2; exit 1; }

# The locations alone. The check runs from a scratch directory, so a path that held only
# relative to where the script started is caught the way the real run would hit it, and
# the real work directory is not created. genesis4 and the MPI launcher are reported and
# not required, since this answers where the run's files are and not whether it can run.
if [[ $SHOW -eq 1 ]]; then
  echo "saturation demo locations:"
  for pair in "demo:$DEMO_DIR" "source:$LUCIFER_DIR" "scripts:$SCRIPTS" \
              "lucifer:$EXE" "python:$PYTHON" "genesis4:${GENESIS:-<not found>}" \
              "mpi launcher:${MPIRUN:-<not found>}" "work dir:$WORK_DIR"; do
    printf '  %-15s %s\n' "${pair%%:*}" "${pair#*:}"
  done
  scratch="$(mktemp -d "${TMPDIR:-/tmp}/sat_locations.XXXXXX")"
  cd "$scratch" || exit 1
  missing=0
  for f in "$EXE" "$PYTHON"; do
    [[ -x "$f" ]] || { echo "  missing after changing directory: $f"; missing=1; }
  done
  for f in "$SCRIPTS/report_fel_saturation.py" "$DEMO_DIR/check_agreement.py" \
           "$DEMO_DIR/Aramis.lat" "$DEMO_DIR/aramis.bmad" "$DEMO_DIR/sat_unavg.bmad" \
           "$DEMO_DIR/sat-prep.in" "$DEMO_DIR/sat-genesis.in" \
           "$DEMO_DIR/lucifer-avg.in" "$DEMO_DIR/lucifer-unavg.in"; do
    [[ -f "$f" ]] || { echo "  missing after changing directory: $f"; missing=1; }
  done
  cd / && rm -rf "$scratch"
  if [[ $missing -ne 0 ]]; then
    echo "saturation demo locations: FAIL"
    exit 1
  fi
  echo "saturation demo locations: all present from $scratch"
  exit 0
fi

mkdir -p "$WORK_DIR"
cp "$DEMO_DIR"/Aramis.lat "$DEMO_DIR"/aramis.bmad "$DEMO_DIR"/sat_unavg.bmad \
   "$DEMO_DIR"/sat-prep.in "$DEMO_DIR"/sat-genesis.in "$DEMO_DIR"/lucifer-avg.in \
   "$DEMO_DIR"/lucifer-unavg.in "$WORK_DIR/"
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
timed "Bmad averaged, $WORKERS threads" sat-avg.log env OMP_NUM_THREADS="$WORKERS" "$EXE" lucifer-avg.in
timed "Bmad unaveraged, $WORKERS threads (~32x; watch its progress lines in sat-unavg.log)" \
      sat-unavg.log env OMP_NUM_THREADS="$WORKERS" "$EXE" lucifer-unavg.in
echo

"$PYTHON" "$DEMO_DIR/check_agreement.py" AramisSat.out.h5 sat-avg.diag.txt sat-unavg.diag.txt || exit 1

"$PYTHON" "$SCRIPTS/report_fel_saturation.py" "$WORK_DIR" || exit 1

echo
echo "Outputs in: $WORK_DIR (report: sat-demo-report.pdf; per-run stats: sat-*.stats.h5)"
