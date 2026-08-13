#!/bin/bash
#
# Measure the coarse-step question of the design brief's section 8.3: how does the gain
# curve converge with the undulator integration step delz, up to SIMPLEX's matched
# configuration of twelve periods per step with slippage exactly one slice per step?
#
# Method: Genesis generates ONE time-dependent initial state (Aramis-td-s12.in: 32 slices
# of spacing 12*lambda0, shot noise on), and the Bmad tracker runs the full 6-FODO line
# from that same dump at delz = 1, 2, 3, 6 and 12 undulator periods. Identical initial
# conditions -- one shot-noise realization -- make the run-to-run differences pure
# integration error, not SASE statistics. The comparison samples total power (summed over
# slices) at the twelve undulator-segment exits, positions every run records regardless
# of delz, against the finest run; a convergence order is fitted from successive step
# halvings/triplings. This is a measurement, not a gate: results are recorded in
# bsim/fel/README.md.
#
#   ./run_delz_sweep.sh [--genesis <path>] [--exe <path>] [--python <path>] [--work-dir <path>]
#
# Same prerequisites as run_fel_benchmark.sh.

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BMAD_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

GENESIS="$HOME/Code/GitHub/Genesis-1.3-Version4/build-metal/genesis4"
EXE=""
PYTHON=""
WORK_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genesis)  GENESIS="$2";  shift 2 ;;
    --exe)      EXE="$2";      shift 2 ;;
    --python)   PYTHON="$2";   shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    -h|--help)  sed -n '2,20p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

if [[ ! -x "$GENESIS" ]]; then
  echo "Error: genesis4 binary not found at: $GENESIS" >&2
  exit 1
fi

if [[ -z "$EXE" ]]; then
  for candidate in "$BMAD_ROOT/debug/bin/fel_track_test" "$BMAD_ROOT/production/bin/fel_track_test"; do
    if [[ -x "$candidate" ]]; then EXE="$candidate"; break; fi
  done
fi
if [[ -z "$EXE" || ! -x "$EXE" ]]; then
  echo "Error: fel_track_test not found. Build with: cd $BMAD_ROOT && BUILD_PRODUCTION=N ./util/conda_compile" >&2
  exit 1
fi

if [[ -z "$PYTHON" ]]; then
  for candidate in "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
                   "$HOME/Code/miniforge3/envs/bmad-fel-validate/bin/python3"; do
    if [[ -x "$candidate" ]]; then PYTHON="$candidate"; break; fi
  done
fi
if [[ -z "$PYTHON" || ! -x "$PYTHON" ]]; then
  echo "Error: Python with numpy/h5py not found (bmad-fel-validate env)." >&2
  exit 1
fi

if [[ -z "$WORK_DIR" ]]; then
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/fel_delz_sweep.XXXXXX")"
fi
mkdir -p "$WORK_DIR"

echo "delz sweep workdir: $WORK_DIR"
cp "$SCRIPT_DIR/genesis4/sweep/Aramis-td-s12.in" "$SCRIPT_DIR/genesis4/Aramis.lat" "$SCRIPT_DIR/bmad/aramis.bmad" "$WORK_DIR/"
cd "$WORK_DIR" || exit 1

export FI_PROVIDER=tcp

echo "--- Genesis: generating the shared initial state (sample = 12) ---------------"
if ! "$GENESIS" Aramis-td-s12.in > genesis-s12.log 2>&1; then
  echo "Genesis dump generation FAILED; log tail:" >&2
  tail -20 genesis-s12.log >&2
  exit 1
fi

# delz in undulator periods (lambdau = 0.015 m). 1 period is the reference.

for np in 1 2 3 6 12; do
  delz=$(echo "$np * 0.015" | bc -l)
  cat > "sweep_p$np.nml" <<NML
&fel_track_params
  lat_file = "aramis.bmad"
  beam_file = "AramisS12-initial.par.h5"
  field_file = "AramisS12-initial.fld.h5"
  out_root = "sweep_p$np"
  gamma0 = 11357.82
  delz = $delz
  interlude_model = "bmad"
&end
NML
  echo "--- fel_track_test: delz = $np period(s) --------------------------------------"
  if ! "$EXE" "sweep_p$np.nml" > "sweep_p$np.log" 2>&1; then
    echo "fel_track_test delz=$delz FAILED; log tail:" >&2
    tail -20 "sweep_p$np.log" >&2
    exit 1
  fi
done

"$PYTHON" - "$WORK_DIR" <<'EOF'
import sys
import numpy as np

w = sys.argv[1]
periods = [1, 2, 3, 6, 12]
lambdau, l_und = 0.015, 3.99

# Total power (summed over slices) per record, from the diag files. Element-end records
# exist in every run; undulator-segment exits are located by z value.

def curve(np_periods):
    d = np.loadtxt(f"{w}/sweep_p{np_periods}.diag.txt")
    nslice = int(d[:, 1].max())
    d = d.reshape(-1, nslice, d.shape[1])
    return d[:, 0, 0], d[:, :, 2].sum(axis=1)

z_ref, p_ref = curve(periods[0])

# Undulator segment exits: 12 segments in the 6-FODO line.
fodo = [3.99, 0.44, 0.08, 0.24, 3.99, 0.44, 0.08, 0.24]
z_exits, z = [], 0.0
for _ in range(6):
    for i, l in enumerate(fodo):
        z += l
        if i in (0, 4):
            z_exits.append(z)
z_exits = np.array(z_exits)

def at_exits(zc, pc):
    idx = [np.argmin(np.abs(zc - ze)) for ze in z_exits]
    assert all(abs(zc[i] - ze) < 1e-9 for i, ze in zip(idx, z_exits))
    return pc[idx]

ref = at_exits(z_ref, p_ref)
print()
print("Total power at the 12 undulator-segment exits, relative to delz = 1 period:")
print(f"{'delz':>10} {'max |dP/P|':>12} {'at saturation':>14}")
errs = {}
for np_p in periods[1:]:
    zc, pc = curve(np_p)
    p = at_exits(zc, pc)
    rel = np.abs(p - ref) / ref
    errs[np_p] = rel
    print(f"{np_p:>7} pd {rel.max():>12.3e} {rel[-1]:>14.3e}")

print()
print("Apparent convergence order p from successive pairs (err ~ delz^p), on the")
print("max-over-exits error (dominated by the exponential-gain region, where the error")
print("is a quasi-systematic gain-length shift; the saturation error is oscillatory and")
print("fits no clean order):")
for b, a in [(2, 3), (3, 6), (6, 12)]:
    ea, eb = errs[a].max(), errs[b].max()
    if eb > 0:
        print(f"  {b} pd -> {a} pd: p = {np.log(ea/eb)/np.log(a/b):.2f}")
EOF
echo
echo "Outputs kept in: $WORK_DIR"
