#!/bin/bash
#
# Head-to-head performance benchmark: the Bmad FEL tracker (OpenMP) against Genesis 1.3
# Version 4 (MPI), serial and at the machine's full performance-core count, on the same
# physics. This reproduces the head-to-head measurement of doc/validation.md as one
# command instead of a by-hand procedure.
#
# Methodology, each point deliberate:
#
#   - Both codes run the full 6-FODO Benchmark1-SASE line, time dependent, seeded plus
#     shot noise, 2048 particles per slice.
#   - The window is 4 slices per worker (sample = 3), so the slice count divides the
#     worker count exactly and Genesis does not pad its time window -- padding would
#     charge Genesis for slices the Bmad side does not track.
#   - Genesis writes the initial dumps. The Bmad tracker imports them, so both track
#     identical particles and field. Before any timing is quoted, the final total powers
#     are required to agree at the documented seam level (the Bmad tracker runs its
#     production configuration, interlude_model = "bmad", which differs from Genesis by
#     the priced transport model difference of ~4e-2 -- see doc/validation.md).
#   - Wall times come from one uniform external clock (/usr/bin/time -p), not from each
#     code's self-report. Serial Genesis runs the plain binary (one rank, no launcher),
#     and parallel Genesis runs mpirun -np <workers>. The Bmad tracker runs the same binary
#     with OMP_NUM_THREADS=1 and =<workers>.
#   - Use production builds. A debug lucifer is accepted with a loud warning,
#     and its numbers are not comparable.
#
# Usage:
#   ./run_perf_benchmark.sh [--workers N] [--genesis <path>] [--exe <path>]
#                           [--mpirun <path>] [--python <path>] [--work-dir <path>]
#
# Defaults: workers = the machine's performance-core count (hw.perflevel0.physicalcpu on
# macOS, nproc on Linux). Binaries and python found as in run_fel_benchmark.sh, except
# the PRODUCTION lucifer is preferred. The work directory is kept.

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BMAD_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

GENESIS="$HOME/Code/GitHub/Genesis-1.3-Version4/build-metal/genesis4"
EXE=""
PYTHON=""
WORK_DIR=""
WORKERS=""
MPIRUN=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genesis)  GENESIS="$2";  shift 2 ;;
    --exe)      EXE="$2";      shift 2 ;;
    --python)   PYTHON="$2";   shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --workers)  WORKERS="$2";  shift 2 ;;
    --mpirun)   MPIRUN="$2";   shift 2 ;;
    -h|--help)  sed -n '2,35p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$WORKERS" ]]; then
  if [[ "$(uname)" == "Darwin" ]]; then
    WORKERS="$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null)"
  fi
  [[ -z "$WORKERS" ]] && WORKERS="$(nproc 2>/dev/null)"
fi
if [[ -z "$WORKERS" || "$WORKERS" -lt 2 ]]; then
  echo "Error: could not detect a performance-core count; pass --workers N." >&2
  exit 1
fi

if [[ ! -x "$GENESIS" ]]; then
  echo "Error: genesis4 binary not found at: $GENESIS (use --genesis)." >&2
  exit 1
fi

# The launcher must match the MPI the Genesis binary links (an OpenMPI mpirun aborts an
# MPICH-linked genesis4 with "unsupported PMI version PMIx", and vice versa). On macOS
# the binary's LC_RPATH names the MPI installation it loads from. Prefer the mpiexec
# that lives beside it, and fall back to PATH only if that turns up nothing.
if [[ -z "$MPIRUN" && "$(uname)" == "Darwin" ]]; then
  while read -r rpath; do
    if [[ -x "$rpath/../bin/mpiexec" ]]; then MPIRUN="$rpath/../bin/mpiexec"; break; fi
  done < <(otool -l "$GENESIS" 2>/dev/null | awk '/LC_RPATH/ {found=1} found && /path / {print $2; found=0}')
fi
if [[ -z "$MPIRUN" ]]; then
  MPIRUN="$(command -v mpiexec || command -v mpirun)"
fi
if [[ -z "$MPIRUN" || ! -x "$MPIRUN" ]]; then
  echo "Error: no MPI launcher found; the parallel Genesis run needs one (use --mpirun)." >&2
  exit 1
fi

DEBUG_EXE_WARNING=0
if [[ -z "$EXE" ]]; then
  for candidate in "$BMAD_ROOT/production/bin/lucifer" "$BMAD_ROOT/debug/bin/lucifer"; do
    if [[ -x "$candidate" ]]; then EXE="$candidate"; break; fi
  done
fi
if [[ -z "$EXE" || ! -x "$EXE" ]]; then
  echo "Error: lucifer not found. Build with: cd $BMAD_ROOT && ./util/conda_compile" >&2
  exit 1
fi
[[ "$EXE" == *"/debug/"* ]] && DEBUG_EXE_WARNING=1

if [[ -z "$PYTHON" ]]; then
  for candidate in "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
                   "$HOME/Code/miniforge3/envs/bmad-fel-validate/bin/python3"; do
    if [[ -x "$candidate" ]]; then PYTHON="$candidate"; break; fi
  done
fi
if [[ -z "$PYTHON" || ! -x "$PYTHON" ]]; then
  echo "Error: Python (numpy + h5py) not found; see run_fel_benchmark.sh." >&2
  exit 1
fi

if [[ -z "$WORK_DIR" ]]; then
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/fel_perf.XXXXXX")"
fi
mkdir -p "$WORK_DIR"

NSLICE=$((4 * WORKERS))
SLEN="$(awk -v n="$NSLICE" 'BEGIN { printf "%.6e", n * 3 * 1e-10 }')"

echo "=============================================================================="
echo "FEL performance benchmark: Bmad tracker (OpenMP) vs Genesis4 (MPI)"
echo "=============================================================================="
echo "  workers:     $WORKERS (performance cores)"
echo "  window:      $NSLICE slices (4 per worker; sample = 3)"
echo "  lucifer: $EXE"
echo "  genesis4:    $GENESIS"
echo "  mpirun:      $MPIRUN"
echo "  workdir:     $WORK_DIR"
if [[ $DEBUG_EXE_WARNING -eq 1 ]]; then
  echo
  echo "  WARNING: using a DEBUG lucifer. Its timings are not comparable to a"
  echo "  release Genesis build; compile production (./util/conda_compile) for real numbers."
fi
echo

cp "$SCRIPT_DIR/genesis4/Aramis.lat" "$SCRIPT_DIR/bmad/aramis.bmad" "$WORK_DIR/"
cd "$WORK_DIR" || exit 1
export FI_PROVIDER=tcp

# The deck is Aramis-td.in's configuration with the window sized to the worker count.
# The serial run carries the &write blocks (the dumps the Bmad side imports). The timed
# parallel run drops them so neither code is charged for dump I/O asymmetrically --
# the Bmad runs write only their final dumps, Genesis's serial run writes initial+final.
# Dump I/O is seconds against minutes of tracking. The asymmetry is noted, not modeled.

make_deck () {   # <file> <with_write: 1|0>
  cat > "$1" <<DECK
&setup
rootname=AramisPerf
lattice=Aramis.lat
beamline=ARAMIS
lambda0=1e-10
gamma0=11357.82
delz=0.045000
shotnoise=1
npart = 2048
nbins = 8
beam_global_stat = true
field_global_stat = true
&end

&lattice
zmatch=9.5
&end

&time
slen = $SLEN
sample = 3
&end

&field
power=5e3
dgrid=2.000000e-04
ngrid=255
waist_size=30e-6
&end

&beam
current=3000
delgam=1.000000
ex=4.000000e-07
ey=4.000000e-07
&end
DECK
  if [[ "$2" == "1" ]]; then
    cat >> "$1" <<'DECK'

&write
field = AramisPerf-initial
beam = AramisPerf-initial
&end
DECK
  fi
  cat >> "$1" <<'DECK'

&track
fft_fieldsolver = true
&end
DECK
}

make_deck perf-serial.in 1
make_deck perf-parallel.in 0

cat > perf.nml <<NML
&fel_params
  lat_file = "aramis.bmad"
  global%out_root = "perf"
  global%interlude_model = "bmad"
  global%write_diag = T
/
&fel_beam_init
  beam_file = "AramisPerf-initial.beam.h5"
  nbins = 8
/
&fel_wavefront_init
  field_file = "AramisPerf-initial.wf.h5"
  wavefront_init%lambda0 = 1e-10
  wavefront_init%window_sample = 3
/
NML

# timed <label> <log> <command...>: uniform external wall clock.
declare -a TIMES
timed () {
  local label="$1" log="$2"; shift 2
  echo "--- $label ---"
  /usr/bin/time -p "$@" > "$log" 2> "$log.time"
  local status=$?
  if [[ $status -ne 0 ]]; then
    echo "FAILED (exit $status); log tail:" >&2
    tail -10 "$log" >&2
    exit 1
  fi
  local wall
  wall="$(awk '/^real/ {print $2}' "$log.time")"
  echo "    wall: ${wall}s"
  TIMES+=("$wall")
}

timed "Genesis4, serial (1 rank; writes the shared dumps)" genesis-serial.log \
  "$GENESIS" perf-serial.in
timed "Genesis4, $WORKERS MPI ranks" genesis-parallel.log \
  "$MPIRUN" -np "$WORKERS" "$GENESIS" perf-parallel.in

# The tracker reads openPMD only, so the Genesis dumps convert once here. The conversion
# is outside every timed section: it is harness work, not either code's.

for pair in "AramisPerf-initial.par.h5 AramisPerf-initial.beam.h5" \
            "AramisPerf-initial.fld.h5 AramisPerf-initial.wf.h5"; do
  if ! "$PYTHON" "$SCRIPT_DIR/scripts/convert_genesis.py" to-openpmd $pair; then
    echo "FAIL: could not convert $pair" >&2
    exit 1
  fi
done
timed "lucifer, 1 thread" fel-serial.log \
  env OMP_NUM_THREADS=1 "$EXE" perf.nml
mv perf.diag.txt perf-serial.diag.txt
timed "lucifer, $WORKERS threads" fel-parallel.log \
  env OMP_NUM_THREADS="$WORKERS" "$EXE" perf.nml
echo

# Sanity before quoting numbers: same physics on both sides. The Bmad seam differs from
# Genesis by the priced ~4e-2 transport model difference. A factor-level disagreement
# means the runs are not comparable and the timings are meaningless.

"$PYTHON" - <<EOF
import numpy as np, h5py, sys
d = np.loadtxt("perf-serial.diag.txt")
nslice = $NSLICE
d = d.reshape(-1, nslice, d.shape[1])
p_bmad = d[-1, :, 2].sum()
with h5py.File("AramisPerf.out.h5") as h5:
    p_gen = h5["Field/power"][-1, :].sum()
rel = abs(p_bmad - p_gen) / p_gen
print(f"sanity: exit total power bmad {p_bmad:.4e} W, genesis {p_gen:.4e} W, rel {rel:.2e}")
if rel > 0.1:
    print("FAIL: disagreement beyond the documented seam level; timings not comparable.")
    sys.exit(1)
EOF
[[ $? -ne 0 ]] && { echo "Outputs kept in: $WORK_DIR" >&2; exit 1; }
echo

echo "=============================================================================="
echo "Results ($NSLICE slices x 2048 particles, full line, $WORKERS workers)"
echo "=============================================================================="
"$PYTHON" - <<EOF
g1, gp, f1, fp = [float(t) for t in "${TIMES[@]}".split()]
print(f"  {'':26s} {'serial':>10s} {'$WORKERS workers':>13s} {'speedup':>9s}")
print(f"  {'Genesis 1.3 v4 (MPI)':26s} {g1:>9.1f}s {gp:>11.1f}s {g1/gp:>8.1f}x")
print(f"  {'lucifer (OpenMP)':26s} {f1:>9.1f}s {fp:>11.1f}s {f1/fp:>8.1f}x")
print()
print(f"  lucifer vs Genesis: {g1/f1:.2f}x faster serial, {gp/fp:.2f}x faster at $WORKERS workers")
EOF

echo
echo "Outputs kept in: $WORK_DIR"
