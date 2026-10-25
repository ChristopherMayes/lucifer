#!/bin/bash
#
# run_perf_benchmark.sh times the Bmad FEL tracker (OpenMP) head to head against Genesis
# 1.3 Version 4 (MPI), serial and at the machine's full performance-core count, on the
# same physics. This reproduces the head-to-head measurement of docs/validation.md as one
# command instead of a by-hand procedure.
#
# Each point of the method is deliberate:
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
#     the priced transport model difference of ~4e-2 -- see docs/validation.md).
#   - Wall times come from one uniform external clock (/usr/bin/time -p), not from each
#     code's self-report. Serial Genesis runs the plain binary (one rank, no launcher),
#     and parallel Genesis runs mpirun -np <workers>. The Bmad tracker runs the same binary
#     with OMP_NUM_THREADS=1 and =<workers>.
#   - Use production builds. A debug lucifer is accepted with a loud warning,
#     and its numbers are not comparable.
#
# The --phases mode answers a different question, where this code's own time goes. It
# needs no genesis4 and no MPI, since it compares the code against itself. Two runs at
# production-like parameters give the per-phase table the tracker's footer prints, and a
# thread sweep over the same configuration gives the scaling curve and the serial
# fraction Amdahl's law implies. The phases partition the walk in the tracker itself
# (code/fel_timer_mod.f90), so the table needs no instrument here beyond the footer.
#
# --phases cannot see the phases that interleave per slice inside one parallel region,
# which are the deposit against its FFT and the field gather against the RK4. Those come
# from a sampling profiler, and docs/performance.md records that split.
#
# Usage:
#   ./run_perf_benchmark.sh [--workers N] [--genesis <path>] [--exe <path>]
#                           [--mpirun <path>] [--python <path>] [--work-dir <path>]
#   ./run_perf_benchmark.sh --phases [--workers N] [--exe <path>] [--python <path>]
#                           [--work-dir <path>] [--big-slices N] [--npart N]
#                           [--ngrid N] [--no-sweep] [--grid-scan "N1 N2 ..."]
#
# --ngrid sets the transverse grid of every input file the phases mode writes, 256 by default.
# --no-sweep leaves the thread sweep out, for a particle count whose profile alone is
# wanted. --grid-scan runs the small configuration at each listed grid three times at the
# worker count, takes the best walk of each, and prints the transform-length table in
# place of the profile and the sweep.
#
# Defaults: workers = the machine's performance-core count (hw.perflevel0.physicalcpu on
# macOS, nproc on Linux). The binaries and Python are found as in run_fel_benchmark.sh,
# except that the production lucifer is preferred. The work directory is kept.

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# The next line sources the locations, LUCIFER_DIR, DEBUG_BIN, PRODUCTION_BIN and
# BMAD_SRC, with their defaults and the names that set them elsewhere (util/locations.sh).
source "$SCRIPT_DIR/../util/locations.sh"

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
EXE=""
PYTHON=""
WORK_DIR=""
WORKERS=""
MPIRUN=""
PHASES=0
BIG_SLICES=504
NPART=2048
NGRID=256
SWEEP_ON=1
GRID_SCAN=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genesis)     GENESIS="$2";    shift 2 ;;
    --exe)         EXE="$2";        shift 2 ;;
    --python)      PYTHON="$2";     shift 2 ;;
    --work-dir)    WORK_DIR="$2";   shift 2 ;;
    --workers)     WORKERS="$2";    shift 2 ;;
    --mpirun)      MPIRUN="$2";     shift 2 ;;
    --phases)      PHASES=1;        shift 1 ;;
    --big-slices)  BIG_SLICES="$2"; shift 2 ;;
    --npart)       NPART="$2";      shift 2 ;;
    --ngrid)       NGRID="$2";      shift 2 ;;
    --no-sweep)    SWEEP_ON=0;      shift 1 ;;
    --grid-scan)   GRID_SCAN="$2";  shift 2 ;;
    -h|--help)  sed -n '2,54p' "${BASH_SOURCE[0]}"; exit 0 ;;
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

if [[ $PHASES -eq 0 && ! -x "$GENESIS" ]]; then
  echo "Error: genesis4 binary not found at: $GENESIS (use --genesis)." >&2
  exit 1
fi

# The launcher must match the MPI the Genesis binary links (an OpenMPI mpirun aborts an
# MPICH-linked genesis4 with "unsupported PMI version PMIx", and vice versa). On macOS
# the binary's LC_RPATH names the MPI installation it loads from. The mpiexec that lives
# beside it is preferred, and PATH is the fallback only if that turns up nothing.
if [[ $PHASES -eq 1 ]]; then
  MPIRUN="not needed in the phases mode"
fi
if [[ -z "$MPIRUN" && "$(uname)" == "Darwin" ]]; then
  while read -r rpath; do
    if [[ -x "$rpath/../bin/mpiexec" ]]; then MPIRUN="$rpath/../bin/mpiexec"; break; fi
  done < <(otool -l "$GENESIS" 2>/dev/null | awk '/LC_RPATH/ {found=1} found && /path / {print $2; found=0}')
fi
if [[ -z "$MPIRUN" ]]; then
  MPIRUN="$(command -v mpiexec || command -v mpirun)"
fi
if [[ $PHASES -eq 0 && ( -z "$MPIRUN" || ! -x "$MPIRUN" ) ]]; then
  echo "Error: no MPI launcher found; the parallel Genesis run needs one (use --mpirun)." >&2
  exit 1
fi

DEBUG_EXE_WARNING=0
if [[ -z "$EXE" ]]; then
  for candidate in "$PRODUCTION_BIN/lucifer" "$DEBUG_BIN/lucifer"; do
    if [[ -x "$candidate" ]]; then EXE="$candidate"; break; fi
  done
fi
if [[ -z "$EXE" || ! -x "$EXE" ]]; then
  echo "Error: lucifer not found. Build with: $LUCIFER_DIR/util/conda_compile" >&2
  exit 1
fi
# Every path is absolute from here on, since the run changes directory.
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"
[[ "$EXE" == *debug* ]] && DEBUG_EXE_WARNING=1

if [[ -z "$PYTHON" ]]; then
  for candidate in "${LUCIFER_PYTHON:-}" \
                   "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
                   "$(command -v python3 2>/dev/null)"; do
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

# ------------------------------------------------------------------------------------
# The phases mode measures where this code's own time goes, from the footer table the
# tracker prints, plus a thread sweep for the serial fraction. Genesis is not involved, so the
# starting state is the tracker's own shot-noise quiet start rather than an imported
# dump: nothing here compares two codes, so nothing here needs identical particles.

if [[ $PHASES -eq 1 ]]; then
  echo "=============================================================================="
  echo "FEL phase profile: where the walk spends its time"
  echo "=============================================================================="
  echo "  lucifer:     $EXE"
  echo "  workers:     $WORKERS (performance cores)"
  echo "  particles:   $NPART per slice"
  echo "  grid:        $NGRID points"
  echo "  workdir:     $WORK_DIR"
  if [[ $DEBUG_EXE_WARNING -eq 1 ]]; then
    echo
    echo "  WARNING: using a DEBUG lucifer. Its phase shares are not the production"
    echo "  build's: the per-file -O3 on fel_track_mod is a production-only setting,"
    echo "  so the particle path is charged more here than it costs in production."
  fi
  echo

  cp "$SCRIPT_DIR/bmad/aramis.bmad" "$WORK_DIR/"
  cd "$WORK_DIR" || exit 1

  # Each slice count gets one input file, the SASE example's configuration with the window resized.
  # A time window of nslice slices at sample = 3 is nslice * 3 * lambda0 long, which is
  # the same arithmetic the head-to-head above does for Genesis.
  #
  # The input file states its slice spacing. An unstated one is derived from the beam and the
  # gain, which on this line is 38 wavelengths rather than 3, so the same window_length
  # would build 8 slices where this profile wants 96. The slice count is what the two
  # columns compare, so it is stated and then checked below.
  #
  # The bunch charge does not scale with the window, so a longer window is a lower
  # current and a different gain. That is deliberate and it does not matter here: the
  # work per step is the same whatever the charge, which the measurement confirms (the
  # phase shares agree to 0.1% across a factor of 5.25 in slice count). Read the two
  # columns as the same work at two sizes rather than as the same physics.

  phase_input () {   # <file> <nslice> <out_root>
    local slen sigz
    slen="$(awk -v n="$2" 'BEGIN { printf "%.6e", n * 3 * 1e-10 }')"
    sigz="$(awk -v n="$2" 'BEGIN { printf "%.6e", n * 3 * 1e-10 / 2 }')"
    cat > "$1" <<INPUT
lattice = "aramis.bmad"

[run]
out_root = "$3"

[window]
window_length = $slen
n_wavelength = 3

[beam]
n_particle = $NPART
bunch_charge = 2.881993782512e-13
z_distribution = "GRID"
z_min = -$sigz
z_max = $sigz
sig_pz = 8.804506566858e-5
a_norm_emit = 4e-7
b_norm_emit = 4e-7
shot_noise = true

[field]
lambda0 = 1e-10
seed_power = 0
grid_n_pts = $NGRID
grid_half_width = 2e-4
INPUT
  }

  # phase_run <label> <nslice> <threads> [<suffix>] does one run and parses out its
  # footer table. The suffix keeps a scan's repeated runs apart.
  phase_run () {
    local label="$1" ns="$2" nt="$3" root log
    root="phase_${ns}_${nt}${4:-}"
    log="$root.log"
    phase_input "$root.toml" "$ns" "$root"
    echo "--- $label: $ns slices, $nt thread(s) ---"
    if ! env OMP_NUM_THREADS="$nt" "$EXE" "$root.toml" > "$log" 2>&1; then
      echo "FAILED; log tail:" >&2
      tail -20 "$log" >&2
      exit 1
    fi
    if ! grep -q '^ Timing' "$log"; then
      echo "FAILED: the run printed no Timing block. An old binary, or a phase left open." >&2
      exit 1
    fi

    # The profile compares two slice counts, so a run that built a different one is
    # measuring a different configuration. The header says what it built.

    local built
    built="$(awk '/^ Beam / {for (i = 1; i <= NF; i++) if ($i == "slices" || $i == "slice") print $(i-1)}' "$log")"
    if [[ "$built" != "$ns" ]]; then
      echo "FAILED: asked for $ns slices and the run built ${built:-none}." >&2
      echo "The input file states window.n_wavelength; a derived spacing would change the count." >&2
      exit 1
    fi
  }

  # The grid scan: the small configuration at each grid, three runs each, the best walk.
  # The transform's cost per grid point is what the table isolates, since the particle
  # work is the same at every grid.

  if [[ -n "$GRID_SCAN" ]]; then
    for g in $GRID_SCAN; do
      NGRID="$g"
      for rep in 1 2 3; do
        phase_run "grid scan" "$((8 * WORKERS))" "$WORKERS" "_g${g}_r${rep}"
      done
    done
    "$PYTHON" - "$WORK_DIR" "$((8 * WORKERS))" "$WORKERS" "$NPART" $GRID_SCAN <<'SCAN'
import re, sys, os
from collections import Counter

work, small, workers, npart = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
grids = [int(g) for g in sys.argv[5:]]

def walk_of(path):
    for line in open(path):
        m = re.match(r'^\s{2,}walk\s+\d+\s+(\d+\.\d+)', line.rstrip())
        if m:
            return float(m.group(1))
    return None

def factors(n):
    out, d = [], 2
    while d * d <= n:
        while n % d == 0:
            out.append(d); n //= d
        d += 1
    if n > 1:
        out.append(n)
    if len(out) == 1:
        return 'prime'
    c = Counter(out)
    return ' x '.join(f'{p}^{k}' if k > 1 else f'{p}' for p, k in sorted(c.items()))

print('==============================================================================')
print(f'Transform length ({small} slices x {npart} particles, full line, {workers} threads, best of three)')
print('==============================================================================')
print(f"    {'ngrid':>6s} {'factors':>12s} {'points':>8s} {'walk [s]':>9s} {'seconds a point':>16s}")
for g in grids:
    walks = [walk_of(os.path.join(work, f'phase_{small}_{workers}_g{g}_r{r}.log')) for r in (1, 2, 3)]
    walks = [w for w in walks if w is not None]
    if not walks:
        print(f'    {g:>6d}  no walk found')
        continue
    best = min(walks)
    print(f'    {g:>6d} {factors(g):>12s} {g * g:>8d} {best:>9.3f} {best / (g * g):>16.3e}')
SCAN
    exit 0
  fi

  phase_run "profile" "$((8 * WORKERS))" "$WORKERS"
  phase_run "profile" "$BIG_SLICES" "$WORKERS"

  # The thread sweep runs the smaller configuration, whose 1-thread run is affordable.
  # It takes powers of two up to the worker count, then the worker count itself.

  if [[ $SWEEP_ON -eq 1 ]]; then
    SWEEP=""
    n=1
    while [[ $n -lt $WORKERS ]]; do SWEEP="$SWEEP $n"; n=$((n * 2)); done
    SWEEP="$SWEEP $WORKERS"
    for nt in $SWEEP; do
      [[ "$nt" == "$WORKERS" && -f "phase_$((8 * WORKERS))_$WORKERS.log" ]] && continue
      phase_run "sweep" "$((8 * WORKERS))" "$nt"
    done
  fi
  echo

  "$PYTHON" - "$WORK_DIR" "$((8 * WORKERS))" "$BIG_SLICES" "$WORKERS" "$NPART" <<'EOF'
import re, sys, glob, os

work, small, big, workers, npart = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])

# The footer's own table, read back rather than re-derived: name, regions, seconds,
# share of the walk. The name field can hold spaces, so the numbers anchor the parse.
ROW = re.compile(r'^\s{2,}(?P<name>\S.*?)\s+(?P<reg>\d+)?\s+(?P<sec>\d+\.\d+)(?:\s+(?P<frac>[-\d.]+)%)?\s*$')

def read_phases(path):
    rows, inblock = [], False
    for line in open(path):
        if line.startswith(' Timing'):
            inblock = True
            continue
        if not inblock:
            continue
        if line.startswith(' Wrote') or line.startswith('====') or not line.strip():
            break
        m = ROW.match(line.rstrip())
        if m:
            rows.append((m.group('name'), int(m.group('reg')) if m.group('reg') else None,
                         float(m.group('sec')), float(m.group('frac')) if m.group('frac') else None))
    return rows

def table(ns, nt, title):
    path = os.path.join(work, f'phase_{ns}_{nt}.log')
    rows = read_phases(path)
    if not rows:
        print(f'  {title}: no Timing block in {path}')
        return
    print(f'  {title}')
    print(f"    {'phase':22s} {'regions':>10s} {'seconds':>10s} {'of walk':>9s}")
    for name, reg, sec, frac in rows:
        r = '' if reg is None else str(reg)
        f = '' if frac is None else f'{frac:.1f}%'
        print(f'    {name:22s} {r:>10s} {sec:>10.3f} {f:>9s}')
    print()

print('==============================================================================')
print(f'Phase profile ({npart} particles per slice, full 6-FODO line, {workers} threads)')
print('==============================================================================')
table(small, workers, f'{small} slices')
table(big,   workers, f'{big} slices')

# The scaling curve, from the walk row of each sweep run, and the serial fraction
# Amdahl's law implies from the best speedup: S = 1/(f + (1-f)/n) inverts to
# f = (n/S - 1)/(n - 1). It is a bound on this configuration, not a constant of the
# code: a serial region that shrinks with slice count reads differently at another size.

walks = {}
for path in sorted(glob.glob(os.path.join(work, f'phase_{small}_*.log'))):
    nt = int(re.search(r'_(\d+)\.log$', path).group(1))
    for name, reg, sec, frac in read_phases(path):
        if name == 'walk':
            walks[nt] = sec

if len(walks) > 1:
    base = walks[min(walks)]
    print('==============================================================================')
    print(f'Thread scaling ({small} slices x {npart} particles, full line)')
    print('==============================================================================')
    print(f"    {'threads':>8s} {'walk [s]':>10s} {'speedup':>9s} {'efficiency':>11s} {'implied serial':>15s}")
    for nt in sorted(walks):
        sp = base / walks[nt]
        ser = ''
        if nt > 1:
            f = (nt / sp - 1) / (nt - 1)
            ser = f'{100 * f:.1f}%'
        print(f'    {nt:>8d} {walks[nt]:>10.1f} {sp:>8.2f}x {100 * sp / nt:>10.0f}% {ser:>15s}')
    print()
EOF
  status=$?
  echo "Outputs kept in: $WORK_DIR"
  exit $status
fi

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
  echo "  release Genesis build; compile production (util/conda_compile) for real numbers."
fi
echo

cp "$SCRIPT_DIR/genesis4/Aramis.lat" "$SCRIPT_DIR/bmad/aramis.bmad" "$WORK_DIR/"
cd "$WORK_DIR" || exit 1
export FI_PROVIDER=tcp

# The input file is Aramis-td.in's configuration with the window sized to the worker count.
# The serial run carries the &write blocks (the dumps the Bmad side imports). The timed
# parallel run drops them so neither code is charged for dump I/O asymmetrically --
# the Bmad runs write only their final dumps, Genesis's serial run writes initial+final.
# Dump I/O is seconds against minutes of tracking. The asymmetry is noted, not modeled.

make_input () {   # <file> <with_write: 1|0>
  cat > "$1" <<INPUT
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
ngrid=256
waist_size=30e-6
&end

&beam
current=3000
delgam=1.000000
ex=4.000000e-07
ey=4.000000e-07
&end
INPUT
  if [[ "$2" == "1" ]]; then
    cat >> "$1" <<'INPUT'

&write
field = AramisPerf-initial
beam = AramisPerf-initial
&end
INPUT
  fi
  cat >> "$1" <<'INPUT'

&track
fft_fieldsolver = true
&end
INPUT
}

make_input perf-serial.in 1
make_input perf-parallel.in 0

cat > perf.toml <<INPUT
lattice = "aramis.bmad"

[run]
out_root = "perf"
interlude_model = "bmad"
write_diag = true

[beam]
file = "AramisPerf-initial.beam.h5"
beamlet_size = 8

[field]
file = "AramisPerf-initial.wf.h5"
lambda0 = 1e-10
INPUT

# timed <label> <log> <command...> times a command on one uniform external wall clock.
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
# is outside every timed section: it is the benchmark's own work, not either code's.

for pair in "AramisPerf-initial.par.h5 AramisPerf-initial.beam.h5" \
            "AramisPerf-initial.fld.h5 AramisPerf-initial.wf.h5"; do
  if ! "$PYTHON" "$SCRIPT_DIR/scripts/convert_genesis.py" to-openpmd $pair; then
    echo "FAIL: could not convert $pair" >&2
    exit 1
  fi
done
timed "lucifer, 1 thread" fel-serial.log \
  env OMP_NUM_THREADS=1 "$EXE" perf.toml
mv perf.diag.txt perf-serial.diag.txt
timed "lucifer, $WORKERS threads" fel-parallel.log \
  env OMP_NUM_THREADS="$WORKERS" "$EXE" perf.toml
echo

# The physics on both sides must agree before any number is quoted. The Bmad seam differs from
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
