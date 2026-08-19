#!/bin/bash
#
# The saturation demonstration: one practical SASE case -- the Aramis benchmark line,
# dark start, lasing through saturation -- tracked three ways from IDENTICAL initial
# particles and field, each on the machine's full performance-core count:
#
#   1. Genesis 1.3 Version 4, MPI, all cores;
#   2. the Bmad FEL tracker's averaged mode (the bmad_standard default), OpenMP, all cores;
#   3. the Bmad tracker's unaveraged mode (fel_tracking = fel_unaveraged: RK4 through the
#      real undulator field, fc/JJ nowhere in its inputs), OpenMP, all cores.
#
# It ends with one figure (scripts/plot_fel_saturation.py) overlaying gain curves,
# per-slice exit power, bunching, energy statistics and the relative differences, with
# the three wall clocks annotated. The methodology follows run_perf_benchmark.sh:
#
#   - One untimed Genesis "prep" run (same rank count, same seed, no &track) writes the
#     initial dumps; the Bmad runs import them, and the timed Genesis run regenerates
#     the identical beam and field from the same seed. Nobody is charged for the dump I/O.
#   - The window is slices-per-worker * workers slices (sample = 3), so the slice count
#     divides the rank count exactly and Genesis pads nothing.
#   - Wall times come from one external clock (/usr/bin/time -p), never a self-report.
#   - Before any timing or figure is produced, the exit total powers must agree:
#     averaged vs Genesis within the documented seam level, unaveraged vs Genesis within
#     the priced integrator-structure difference. Disagreement fails the demo.
#
# Expect the unaveraged run to cost ~25-35x the averaged one: it takes fel_steps_per_period
# RK4 substeps (default 20) plus a diffraction per substep where the averaged mode takes
# one step per period. That price is the point: the same answer emerges from raw dynamics.
#
# Usage:
#   ./run_saturation_demo.sh [--workers N] [--slices-per-worker N] [--genesis <path>]
#                            [--exe <path>] [--mpirun <path>] [--python <path>]
#                            [--work-dir <path>]
#
# Defaults: workers = performance-core count; 8 slices per worker (96 slices at 12);
# binaries and python found as in run_perf_benchmark.sh (production fel_track_test).

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BMAD_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

GENESIS="$HOME/Code/GitHub/Genesis-1.3-Version4/build-metal/genesis4"
EXE=""
PYTHON=""
WORK_DIR=""
WORKERS=""
MPIRUN=""
SPW=8

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genesis)  GENESIS="$2";  shift 2 ;;
    --exe)      EXE="$2";      shift 2 ;;
    --python)   PYTHON="$2";   shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --workers)  WORKERS="$2";  shift 2 ;;
    --mpirun)   MPIRUN="$2";   shift 2 ;;
    --slices-per-worker) SPW="$2"; shift 2 ;;
    -h|--help)  sed -n '2,38p' "${BASH_SOURCE[0]}"; exit 0 ;;
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

# rpath-matched launcher, as in run_perf_benchmark.sh (a mismatched mpirun aborts).
if [[ -z "$MPIRUN" && "$(uname)" == "Darwin" ]]; then
  while read -r rpath; do
    if [[ -x "$rpath/../bin/mpiexec" ]]; then MPIRUN="$rpath/../bin/mpiexec"; break; fi
  done < <(otool -l "$GENESIS" 2>/dev/null | awk '/LC_RPATH/ {found=1} found && /path / {print $2; found=0}')
fi
if [[ -z "$MPIRUN" ]]; then
  MPIRUN="$(command -v mpiexec || command -v mpirun)"
fi
if [[ -z "$MPIRUN" || ! -x "$MPIRUN" ]]; then
  echo "Error: no MPI launcher found (use --mpirun)." >&2
  exit 1
fi

DEBUG_EXE_WARNING=0
if [[ -z "$EXE" ]]; then
  for candidate in "$BMAD_ROOT/production/bin/fel_track_test" "$BMAD_ROOT/debug/bin/fel_track_test"; do
    if [[ -x "$candidate" ]]; then EXE="$candidate"; break; fi
  done
fi
if [[ -z "$EXE" || ! -x "$EXE" ]]; then
  echo "Error: fel_track_test not found. Build with: cd $BMAD_ROOT && ./util/conda_compile" >&2
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
  echo "Error: Python (numpy + h5py + matplotlib) not found; see run_fel_benchmark.sh." >&2
  exit 1
fi

if [[ -z "$WORK_DIR" ]]; then
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/fel_sat.XXXXXX")"
fi
mkdir -p "$WORK_DIR"

NSLICE=$((SPW * WORKERS))
SLEN="$(awk -v n="$NSLICE" 'BEGIN { printf "%.6e", n * 3 * 1e-10 }')"

echo "=============================================================================="
echo "FEL saturation demo: one SASE case, three trackers, $WORKERS cores each"
echo "=============================================================================="
echo "  window:      $NSLICE slices (sample = 3), full 6-FODO Aramis line, dark start"
echo "  fel_track_test: $EXE"
echo "  genesis4:    $GENESIS"
echo "  mpirun:      $MPIRUN"
echo "  workdir:     $WORK_DIR"
if [[ $DEBUG_EXE_WARNING -eq 1 ]]; then
  echo
  echo "  WARNING: using a DEBUG fel_track_test; timings are not comparable to a release"
  echo "  Genesis build. Compile production (./util/conda_compile) for real numbers."
fi
echo

cp "$SCRIPT_DIR/genesis4/Aramis.lat" "$SCRIPT_DIR/bmad/aramis.bmad" "$WORK_DIR/"
cd "$WORK_DIR" || exit 1
export FI_PROVIDER=tcp

# The deck is the benchmark's Aramis configuration, dark start (power = 0, growth from
# shot noise alone). The prep deck writes the initial dumps and does not track; the
# timed deck tracks and writes nothing but its .out.h5. Same seed, same rank count,
# so the timed run regenerates exactly the particles and field the dumps hold.

make_deck () {   # <file> <rootname> <with_write: 1|0> <with_track: 1|0>
  cat > "$1" <<DECK
&setup
rootname=$2
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
power=0
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
  if [[ "$3" == "1" ]]; then
    cat >> "$1" <<'DECK'

&write
field = AramisSat-initial
beam = AramisSat-initial
&end
DECK
  fi
  if [[ "$4" == "1" ]]; then
    cat >> "$1" <<'DECK'

&track
fft_fieldsolver = true
&end
DECK
  fi
}

make_deck sat-prep.in    AramisSatPrep 1 0
make_deck sat-genesis.in AramisSat     0 1

cat > sat-avg.nml <<'NML'
&fel_track_params
  lat_file = "aramis.bmad"
  beam_file = "AramisSat-initial.par.h5"
  field_file = "AramisSat-initial.fld.h5"
  out_root = "sat-avg"
  interlude_model = "bmad"
&end
NML

# The unaveraged selection lives in the lattice, where element parameters belong.
# The default 2-period ramps carry two priced effects at the segment ends: their
# slippage deficit (~3.3 rad of optical phase per end at these parameters) is
# compensated exactly by the mode's built-in end-of-ramp phase jump (FINDINGS 7.26,
# manual sec:unaveraged), and their reduced coupling length costs ~2% in ln power
# per segment -- real field physics of adiabatic ends, deliberately not hidden.
cat > sat_unavg.bmad <<'LAT'
call, file = aramis.bmad
fel_unaveraged = 1
wiggler::*[FEL_TRACKING] = fel_unaveraged
LAT
sed -e 's/aramis.bmad/sat_unavg.bmad/' -e 's/sat-avg/sat-unavg/' sat-avg.nml > sat-unavg.nml

echo "--- prep: Genesis writes the shared initial dumps (untimed, no tracking) ---"
if ! "$MPIRUN" -np "$WORKERS" "$GENESIS" sat-prep.in > sat-prep.log 2>&1; then
  echo "Prep run FAILED; log tail:" >&2; tail -10 sat-prep.log >&2; exit 1
fi
echo

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

timed "Genesis4, $WORKERS MPI ranks" sat-genesis.log \
  "$MPIRUN" -np "$WORKERS" "$GENESIS" sat-genesis.in
timed "Bmad averaged (bmad_standard default), $WORKERS threads" sat-avg.log \
  env OMP_NUM_THREADS="$WORKERS" "$EXE" sat-avg.nml
timed "Bmad unaveraged (fel_tracking = fel_unaveraged), $WORKERS threads" sat-unavg.log \
  env OMP_NUM_THREADS="$WORKERS" "$EXE" sat-unavg.nml
echo

# The answers must agree before the timings mean anything. Averaged vs Genesis carries
# the documented seam-transport difference (~4e-2 on the benchmark line; measured 4.9e-4
# at saturation, where the power self-limits). The unaveraged mode additionally carries
# the shot-noise radiation channel it physically resolves and the averaged model does
# not track (README "The saturation demo", FINDINGS 7.27): its SASE curve rides ~2%/m
# above the KMR codes, +0.6 ln at the 57 m exit -- measured n-particle- and
# steps-per-period-independent, so it is the loaded shot noise radiating, not a
# numerical artifact. The check holds |ln| <= 1.0: an order-of-magnitude disagreement
# still fails the demo.

"$PYTHON" - <<EOF
import numpy as np, h5py, sys, math
nslice = $NSLICE
def total_exit(fn):
    d = np.loadtxt(fn); d = d.reshape(-1, nslice, d.shape[1])
    return d[-1, :, 2].sum()
with h5py.File("AramisSat.out.h5") as h5:
    p_gen = h5["Field/power"][-1, :].sum()
p_avg, p_unavg = total_exit("sat-avg.diag.txt"), total_exit("sat-unavg.diag.txt")
rel = abs(p_avg - p_gen) / p_gen
lnr = abs(math.log(p_unavg / p_gen))
print(f"check: exit total power  Genesis {p_gen:.4e} W")
print(f"check: Bmad averaged     {p_avg:.4e} W, rel {rel:.2e} (must be <= 0.15)")
print(f"check: Bmad unaveraged   {p_unavg:.4e} W, |ln ratio| {lnr:.3f} (must be <= 1.0)")
if rel > 0.15 or lnr > 1.0:
    print("FAIL: disagreement beyond the documented levels; timings are meaningless.")
    sys.exit(1)
EOF
[[ $? -ne 0 ]] && { echo "Outputs kept in: $WORK_DIR" >&2; exit 1; }
echo

echo "=============================================================================="
echo "Results ($NSLICE slices x 2048 particles, dark start to z = 57 m, $WORKERS cores)"
echo "=============================================================================="
G_T="${TIMES[0]}"; A_T="${TIMES[1]}"; U_T="${TIMES[2]}"
"$PYTHON" - <<EOF
g, a, u = $G_T, $A_T, $U_T
print(f"  {'Genesis 1.3 v4, MPI':34s} {g:>9.1f} s")
print(f"  {'Bmad averaged (bmad_standard)':34s} {a:>9.1f} s   ({g/a:.2f}x vs Genesis)")
print(f"  {'Bmad unaveraged (RK4, real field)':34s} {u:>9.1f} s   ({u/a:.1f}x the averaged mode)")
EOF
echo

"$PYTHON" "$SCRIPT_DIR/scripts/plot_fel_saturation.py" AramisSat.out.h5 \
  sat-avg.diag.txt sat-unavg.diag.txt --times "$G_T,$A_T,$U_T" --workers "$WORKERS" \
  -o sat-demo.png || exit 1

echo
echo "Outputs kept in: $WORK_DIR (figure: $WORK_DIR/sat-demo.png)"
