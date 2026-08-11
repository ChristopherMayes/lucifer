#!/bin/bash
#
# Validate the Bmad FEL tracker against Genesis 1.3 Version 4: steady state and
# time dependent with slippage.
#
# One command runs everything: Genesis over the Benchmark1-SASE lattice (writing the
# initial dumps both codes start from), Genesis over a single undulator segment importing
# the same dumps, the same pair again time dependent (32 slices, Aramis-td.in), the Bmad
# tracker over every configuration (each full line twice, once with the Bmad seam and
# once with the transcribed Genesis interlude model), a thread-count-independence rerun
# of the time-dependent single segment (8 threads must reproduce 1 thread bit for bit),
# and the tiered comparison printing the largest relative difference of each tier.
#
#   ./run_fel_benchmark.sh
#
# Prerequisites:
#
#   1. Bmad built, so debug/bin/fel_track_test (or production/bin) exists:
#        BUILD_PRODUCTION=N ./util/conda_compile        # from the bmad-ecosystem root
#
#   2. A genesis4 binary (CPU is fine). Default location below; override with --genesis.
#      Built per the Genesis repository's instructions; the benchmark needs FFTW support
#      compiled in, since it runs with fft_fieldsolver = true.
#
#   3. Python with numpy and h5py: the bmad-fel-validate environment
#      (conda env create -f ../../wavefront/tests/environment.yml).
#
# Options:
#   --genesis <path>    genesis4 binary. Default: ~/Code/GitHub/Genesis-1.3-Version4/build-metal/genesis4
#   --exe <path>        fel_track_test binary. Default: debug then production.
#   --python <path>     Python interpreter. Default: the bmad-fel-validate conda env.
#   --work-dir <path>   Where to run. Default: a temporary directory (kept on failure).
#
# Exit status is zero only if every tier passes its tolerance.

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
    -h|--help)  sed -n '2,32p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

# The comparison is against Genesis; without the binary there is nothing to compare and
# the only correct behavior is to fail loudly, not to skip.

if [[ ! -x "$GENESIS" ]]; then
  echo "Error: genesis4 binary not found at: $GENESIS" >&2
  echo "Point at one with --genesis <path>." >&2
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
  echo "Error: Python not found. Create the environment:" >&2
  echo "  conda env create -f $BMAD_ROOT/bsim/wavefront/tests/environment.yml" >&2
  exit 1
fi

KEEP_WORK_DIR=1
if [[ -z "$WORK_DIR" ]]; then
  WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/fel_benchmark.XXXXXX")"
  KEEP_WORK_DIR=0
fi
mkdir -p "$WORK_DIR"

echo "=============================================================================="
echo "FEL steady-state benchmark"
echo "=============================================================================="
echo "  fel_track_test: $EXE"
echo "  genesis4:    $GENESIS"
echo "  python:      $PYTHON"
echo "  workdir:     $WORK_DIR"
echo

cp "$SCRIPT_DIR/Aramis-ss.in" "$SCRIPT_DIR/Aramis.lat" \
   "$SCRIPT_DIR/Aramis-1seg.in" "$SCRIPT_DIR/Aramis-1seg.lat" \
   "$SCRIPT_DIR/Aramis-td.in" "$SCRIPT_DIR/Aramis-td-1seg.in" \
   "$SCRIPT_DIR/aramis.bmad" "$SCRIPT_DIR/aramis_1seg.bmad" "$WORK_DIR/"

cd "$WORK_DIR" || exit 1

# Without FI_PROVIDER=tcp the MPI runtime's provider search adds tens of seconds
# (FINDINGS.md / design brief 11.1).
export FI_PROVIDER=tcp

echo "--- Genesis: full line (writes the shared initial dumps) ---------------------"
if ! "$GENESIS" Aramis-ss.in > genesis-full.log 2>&1; then
  echo "Genesis full-line run FAILED; log tail:" >&2
  tail -20 genesis-full.log >&2
  exit 1
fi
tail -3 genesis-full.log
echo

echo "--- Genesis: single segment (imports the same dumps) -------------------------"
if ! "$GENESIS" Aramis-1seg.in > genesis-1seg.log 2>&1; then
  echo "Genesis single-segment run FAILED; log tail:" >&2
  tail -20 genesis-1seg.log >&2
  exit 1
fi
tail -3 genesis-1seg.log
echo

echo "--- Genesis: full line, time dependent (writes the TD initial dumps) ---------"
if ! "$GENESIS" Aramis-td.in > genesis-td.log 2>&1; then
  echo "Genesis time-dependent run FAILED; log tail:" >&2
  tail -20 genesis-td.log >&2
  exit 1
fi
tail -3 genesis-td.log
echo

echo "--- Genesis: single segment, time dependent (imports the TD dumps) -----------"
if ! "$GENESIS" Aramis-td-1seg.in > genesis-td1seg.log 2>&1; then
  echo "Genesis time-dependent single-segment run FAILED; log tail:" >&2
  tail -20 genesis-td1seg.log >&2
  exit 1
fi
tail -3 genesis-td1seg.log
echo

# make_nml <nml> <lattice> <out_root> <interlude_model> <dump_root> [extra]

make_nml () {
  cat > "$1" <<NML
&fel_track_params
  lat_file = "$2"
  beam_file = "$5-initial.par.h5"
  field_file = "$5-initial.fld.h5"
  out_root = "$3"
  gamma0 = 11357.82
  delz = 0.045
  und_aw = 0.84853
  und_lambdau = 0.015
  und_kx = 0.5, und_ky = 0.5
  und_helical = T
  interlude_model = "$4"
${6:+  $6}
&end
NML
}

make_nml tier1.nml  aramis_1seg.bmad tier1  bmad    Aramis
make_nml tier2.nml  aramis.bmad      tier2  bmad    Aramis
make_nml tier2g.nml aramis.bmad      tier2g genesis Aramis
make_nml tier1s.nml aramis_1seg.bmad tier1s bmad    Aramis "split_weights = T"
make_nml td1.nml    aramis_1seg.bmad td1    bmad    AramisTD
make_nml td2.nml    aramis.bmad      td2    bmad    AramisTD
make_nml td2g.nml   aramis.bmad      td2g   genesis AramisTD

# The documented tier numbers are single-thread runs; results must not depend on the
# thread count, and the explicit gate for that follows the loop.

export OMP_NUM_THREADS=1

for tier in tier1 tier2 tier2g tier1s td1 td2 td2g; do
  echo "--- fel_track_test: $tier -------------------------------------------------------"
  if ! "$EXE" $tier.nml > fel-$tier.log 2>&1; then
    echo "fel_track_test $tier FAILED; log tail:" >&2
    tail -20 fel-$tier.log >&2
    exit 1
  fi
  tail -4 fel-$tier.log
  echo
done

# Thread-count independence: the time-dependent single-segment configuration rerun with
# eight threads must reproduce the one-thread run BIT FOR BIT -- each slice's arithmetic
# is independent of which thread runs it, so any difference at all is a race. The diag
# file is compared byte for byte; the dumps dataset by dataset, exactly (HDF5 object
# headers carry timestamps, so whole-file cmp would false-alarm).

echo "--- thread independence: td1 with OMP_NUM_THREADS=8 --------------------------"
sed 's/out_root = "td1"/out_root = "td1t8"/' td1.nml > td1t8.nml
if ! OMP_NUM_THREADS=8 "$EXE" td1t8.nml > fel-td1t8.log 2>&1; then
  echo "fel_track_test td1t8 FAILED; log tail:" >&2
  tail -20 fel-td1t8.log >&2
  exit 1
fi
THREADS_OK=1
if ! cmp -s td1.diag.txt td1t8.diag.txt; then
  echo "FAIL: diag files differ between 1 and 8 threads" >&2
  THREADS_OK=0
fi
if ! "$PYTHON" - <<'PYEOF'
import sys
import h5py
import numpy as np

def identical(fa, fb):
    bad = []
    with h5py.File(fa) as a, h5py.File(fb) as b:
        names_a, names_b = [], []
        a.visit(lambda n: names_a.append(n) if isinstance(a[n], h5py.Dataset) else None)
        b.visit(lambda n: names_b.append(n) if isinstance(b[n], h5py.Dataset) else None)
        if names_a != names_b:
            return [f"dataset lists differ: {len(names_a)} vs {len(names_b)}"]
        for n in names_a:
            if not np.array_equal(a[n][...], b[n][...]):
                bad.append(n)
    return bad

ok = True
for fa, fb in (("td1-final.fld.h5", "td1t8-final.fld.h5"),
               ("td1-final.par.h5", "td1t8-final.par.h5")):
    bad = identical(fa, fb)
    if bad:
        ok = False
        print(f"FAIL: {fa} vs {fb}: {len(bad)} datasets differ (first: {bad[0]})")
sys.exit(0 if ok else 1)
PYEOF
then
  THREADS_OK=0
fi
if [[ $THREADS_OK -ne 1 ]]; then
  echo "FAIL: thread-count independence violated; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
echo "  1-thread and 8-thread runs are bit-identical (diag byte-equal, dumps dataset-equal)"
echo

"$PYTHON" "$SCRIPT_DIR/compare_fel.py" "$WORK_DIR"
STATUS=$?

if [[ $STATUS -eq 0 && $KEEP_WORK_DIR -eq 0 ]]; then
  rm -rf "$WORK_DIR"
elif [[ $STATUS -ne 0 ]]; then
  echo "Outputs kept in: $WORK_DIR"
fi

exit $STATUS
