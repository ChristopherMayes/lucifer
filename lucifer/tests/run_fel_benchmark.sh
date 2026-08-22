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
#   1. Bmad built, so debug/bin/lucifer (or production/bin) exists:
#        BUILD_PRODUCTION=N ./util/conda_compile        # from the bmad-ecosystem root
#
#   2. A genesis4 binary (CPU is fine). Default location below; override with --genesis.
#      Built per the Genesis repository's instructions; the benchmark needs FFTW support
#      compiled in, since it runs with fft_fieldsolver = true.
#
#   3. Python with numpy and h5py: the bmad-fel-validate environment
#      (conda env create -f ../wavefront/tests/environment.yml).
#
# Options:
#   --genesis <path>    genesis4 binary. Default: ~/Code/GitHub/Genesis-1.3-Version4/build-metal/genesis4
#   --exe <path>        lucifer binary. Default: debug then production.
#   --python <path>     Python interpreter. Default: the bmad-fel-validate conda env.
#   --work-dir <path>   Where to run. Default: a temporary directory (kept on failure).
#
# Exit status is zero only if every tier passes its tolerance.

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BMAD_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

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
  for candidate in "$BMAD_ROOT/debug/bin/lucifer" "$BMAD_ROOT/production/bin/lucifer"; do
    if [[ -x "$candidate" ]]; then EXE="$candidate"; break; fi
  done
fi
if [[ -z "$EXE" || ! -x "$EXE" ]]; then
  echo "Error: lucifer not found. Build with: cd $BMAD_ROOT && BUILD_PRODUCTION=N ./util/conda_compile" >&2
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
  echo "  conda env create -f $BMAD_ROOT/lucifer/wavefront/tests/environment.yml" >&2
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
echo "  lucifer: $EXE"
echo "  genesis4:    $GENESIS"
echo "  python:      $PYTHON"
echo "  workdir:     $WORK_DIR"
echo

# Inputs are grouped by consumer and configuration (genesis4/<config>/*.in with the
# shared .lat files one level up; bmad/*.bmad); the run itself is flat in WORK_DIR, so
# the decks' internal lattice= references need no paths.
cp "$SCRIPT_DIR"/genesis4/*.lat "$SCRIPT_DIR"/genesis4/*/*.in \
   "$SCRIPT_DIR"/bmad/*.bmad "$WORK_DIR/"

cd "$WORK_DIR" || exit 1

# Without FI_PROVIDER=tcp the MPI runtime's provider search adds tens of seconds to
# every Genesis launch (measured; the tcp provider is always sufficient on one node).
export FI_PROVIDER=tcp

# Every section prints its wall time, so a regression in test cost is visible from the
# harness output itself.

BENCH_T_LAST=$SECONDS
section_time () {   # <label>
  echo "  [time: $1 $((SECONDS - BENCH_T_LAST)) s]"
  BENCH_T_LAST=$SECONDS
}

# The Genesis reference runs form three independent chains -- [ss -> 1seg],
# [td -> td-1seg, td-sc, td-wake], [td-sase] -- which run concurrently (same
# decks, same single-process runs; only the wall clock changes).

echo "--- Genesis references: three chains, concurrently ---------------------------"
run_genesis () {   # <deck> ...: the decks of one chain, in order
  for deck in "$@"; do
    if ! "$GENESIS" $deck.in > genesis-$deck.log 2>&1; then
      echo "Genesis $deck run FAILED; log tail:" >&2
      tail -20 genesis-$deck.log >&2
      return 1
    fi
  done
}
run_genesis Aramis-ss Aramis-1seg &
GEN_PID_SS=$!
run_genesis Aramis-td Aramis-td-1seg Aramis-td-sc Aramis-td-wake &
GEN_PID_TD=$!
run_genesis Aramis-td-sase &
GEN_PID_SASE=$!
GEN_OK=1
wait $GEN_PID_SS   || GEN_OK=0
wait $GEN_PID_TD   || GEN_OK=0
wait $GEN_PID_SASE || GEN_OK=0
if [[ $GEN_OK -ne 1 ]]; then
  echo "FAIL: a Genesis reference chain failed (see above)" >&2
  exit 1
fi
for log in genesis-Aramis-ss genesis-Aramis-td genesis-Aramis-td-sase; do
  tail -3 $log.log
done
section_time genesis-references
echo

# make_nml <nml> <lattice> <out_root> <interlude_model> <dump_root> [params_extra] [beam_extra]
# The three input groups (manual sec:program): extra &fel_params content (wake/sc) as
# argument 6, extra &fel_beam_init content (check knobs) as argument 7.

make_nml () {
  cat > "$1" <<NML
&fel_params
  lat_file = "$2"
  global%out_root = "$3"
  global%interlude_model = "$4"
  global%write_diag = T
${6:+  $6}
/
&fel_beam_init
  beam_file = "$5-initial.par.h5"
${7:+  $7}
/
&fel_wavefront_init
  field_file = "$5-initial.fld.h5"
/
NML
}

make_nml tier1.nml  aramis_1seg_val.bmad tier1  bmad    Aramis
make_nml tier2.nml  aramis_val.bmad      tier2  bmad    Aramis
make_nml tier2g.nml aramis_val.bmad      tier2g genesis Aramis
make_nml tier1s.nml aramis_1seg_val.bmad tier1s bmad    Aramis "" "split_weights = T"
make_nml tier1u.nml aramis_1seg_unavg.bmad tier1u bmad Aramis
make_nml td1.nml    aramis_1seg_val.bmad td1    bmad    AramisTD
make_nml td2.nml    aramis_val.bmad      td2    bmad    AramisTD
make_nml td2g.nml   aramis_val.bmad      td2g   genesis AramisTD
make_nml tdsase.nml aramis_val.bmad      tdsase genesis AramisTDSASE
make_nml tdsc.nml   aramis_1seg_val.bmad tdsc   genesis AramisTD "sc%rmax = 250e-6
  sc%nz = 2
  sc%nphi = 1
  sc%longrange = T"
make_nml tdwk.nml   aramis_1seg_val.bmad tdwk   genesis AramisTD "wake%on = T
  wake%radius = 2.5e-3
  wake%conductivity = 5.813e7
  wake%relaxation = 8.1e-6
  wake%gap = 0.5e-3
  wake%lgap = 0.015
  wake%hrough = 100e-9
  wake%lrough = 100e-6"

# Assertion checks: a lattice whose FEL element is missing b_max, missing l_period, or
# uses a fieldmap field_calc (deliverable 9, brief 7.5), or that carries Bmad wakes on
# any element (the slice-at-a-time seam cannot apply them meaningfully), must be
# REFUSED BY NAME -- the failure message names the offending attribute and element -- not passed
# through to fail downstream with an unrelated message (missing b_max) or a segfault in
# the parse-time reference tracking (fieldmap). Each check mutates one attribute of the
# real single-segment lattice and requires both a nonzero exit and the by-name message.

echo "--- FEL-element assertion checks (refusal by name) -----------------------------"
grep -v "b_max" aramis_1seg.bmad                                   > refusal_bmax.bmad
sed 's/l_period = 0.015, //' aramis_1seg.bmad                      > refusal_lperiod.bmad
sed 's/field_calc = helical_model/field_calc = fieldmap/' aramis_1seg.bmad > refusal_fieldmap.bmad

GATES_OK=1
run_assert_refusal () {   # <name> <by-name message fragment>
  make_nml refusal_$1.nml refusal_$1.bmad refusal_$1 bmad Aramis
  if "$EXE" refusal_$1.nml > fel-refusal_$1.log 2>&1; then
    echo "FAIL: refusal_$1 lattice was accepted (exit 0); it must be refused" >&2
    GATES_OK=0
  elif ! grep -q "$2" fel-refusal_$1.log; then
    echo "FAIL: refusal_$1 refused, but not by name; log tail:" >&2
    tail -5 fel-refusal_$1.log >&2
    GATES_OK=0
  else
    echo "  refusal_$1: refused by name ($(grep "$2" fel-refusal_$1.log | head -1 | cut -c1-60)...)"
  fi
}
cat > refusal_lrwake.bmad <<'LAT'
call, file = aramis_1seg.bmad
PW: pipe, l = 0.1, lr_wake = {mode = {2e5, 0.1, 1e-5, 0.3, 2, 0.7}}
SEGW: line = (UND, PW)
use, SEGW
LAT

cat > refusal_zmax.bmad <<'LAT'
call, file = aramis_1seg.bmad
PW: pipe, l = 0.1, sr_wake = {amp_scale = 1.0, scale_with_length = T, z_max = 1e-12,
  longitudinal = {-8266.6e2, 26.6, 46089.2, 1.578966/twopi, none}}
SEGW: line = (UND, PW)
use, SEGW
LAT

run_assert_refusal bmax     "zero b_max"
run_assert_refusal lperiod  "zero l_period"
run_assert_refusal fieldmap "field_calc must be planar_model"
run_assert_refusal lrwake   "lr (multi-bunch) wakes are not supported"
run_assert_refusal zmax     "sr wake z_max can handle"
if [ "$GATES_OK" -ne 1 ]; then
  echo "FAIL: FEL-element assertion checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time refusals
echo

# The documented tier numbers are single-thread runs; results must not depend on the
# thread count, and the explicit check for that follows the loop.

export OMP_NUM_THREADS=1

# The tiers are independent single-thread processes over the shared dumps; they
# run concurrently (12 performance cores hold all eleven) and their logs print
# in the usual order afterwards. Same configurations, same single-thread
# arithmetic -- the documented tier numbers are unchanged.

TIERS="tier1 tier1u tier2 tier2g tier1s td1 td2 td2g tdsase tdsc tdwk"
TIER_PIDS=""
for tier in $TIERS; do
  "$EXE" $tier.nml > fel-$tier.log 2>&1 &
  TIER_PIDS="$TIER_PIDS $tier:$!"
done
TIERS_OK=1
for entry in $TIER_PIDS; do
  tier=${entry%%:*}
  if ! wait "${entry##*:}"; then
    echo "lucifer $tier FAILED; log tail:" >&2
    tail -20 fel-$tier.log >&2
    TIERS_OK=0
  fi
done
if [[ $TIERS_OK -ne 1 ]]; then exit 1; fi
for tier in $TIERS; do
  echo "--- lucifer: $tier -------------------------------------------------------"
  tail -4 fel-$tier.log
  echo
done
section_time tiers-all-eleven

# Thread-count independence: the time-dependent single-segment configuration rerun with
# eight threads must reproduce the one-thread run BIT FOR BIT -- each slice's arithmetic
# is independent of which thread runs it, so any difference at all is a race. The diag
# file is compared byte for byte; the dumps dataset by dataset, exactly (HDF5 object
# headers carry timestamps, so whole-file cmp would false-alarm).

echo "--- thread independence: td1 with OMP_NUM_THREADS=8 --------------------------"
sed 's/out_root = "td1"/out_root = "td1t8"/' td1.nml > td1t8.nml
if ! OMP_NUM_THREADS=8 "$EXE" td1t8.nml > fel-td1t8.log 2>&1; then
  echo "lucifer td1t8 FAILED; log tail:" >&2
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
section_time thread-independence
echo

# Shot-noise checks (deliverable 6). The statistical check is self-referenced (FINDINGS
# 6.9: Genesis cannot represent weighted noise); the SASE startup cross-check pits the
# two codes' fully independent loaders and RNGs against each other at the level the
# noise sets, the startup power.

echo "--- shot-noise statistical check (weighted Fawley loading) ---------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_shot_noise.py" --exe "$EXE" --workdir "$WORK_DIR" --seeds 15; then
  echo "FAIL: shot-noise statistics; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time shot-noise
echo

echo "--- SASE startup cross-check against Genesis's loader -------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_sase_startup.py" --exe "$EXE" --genesis "$GENESIS" --workdir "$WORK_DIR" --seeds 4; then
  echo "FAIL: SASE startup level; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time sase-startup
echo

# Seam-wake checks (deliverable 11): element sr wakes across the whole window --
# closed-form pseudomode ramp, exact causality with the d8 direction cross-check, the
# z_long kernel cross-validation against the deliverable-8 wake model (first-principles
# tight, resolved-beam at the derived boundary bound), split-weight invariance and
# thread determinism. Self-referenced; needs no Genesis.

echo "--- seam-wake checks -------------------------------------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_seam_wake.py" --exe "$EXE" --workdir "$WORK_DIR"; then
  echo "FAIL: seam-wake checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time seam-wake
echo

# Distribution-import checks (deliverable 10): the bunch_struct resampler transcribed
# from Genesis's SDDSBeam.cpp -- exact where no RNG enters (the per-slice current
# profile against Genesis importing the SAME distribution file; the match transform
# hitting its Twiss targets; split-weight invariance; thread determinism), statistical
# where the resampling RNG forces it (slice Twiss recovery, startup power cross-code).

echo "--- distribution-import checks --------------------------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_import.py" --exe "$EXE" --genesis "$GENESIS" --workdir "$WORK_DIR" --seeds 4; then
  echo "FAIL: distribution-import checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time import
echo

# Slice-migration checks (deliverable 7): conservation under heavy migration, exact
# phase continuity of the moves, and no-op bit identity (self-referenced, FINDINGS 6.9
# -- Genesis migrates only under one4one).

echo "--- slice-migration checks ------------------------------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_migration.py" --exe "$EXE" --workdir "$WORK_DIR"; then
  echo "FAIL: migration checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time migration
echo

# Collective-effects self-referenced checks (deliverable 8): exact energy bookkeeping of
# the wake's eloss on a cold dark beam, and the stale-wake structural check (the
# convolution must follow the currents under migration).

echo "--- collective-effects checks ---------------------------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_collective.py" --exe "$EXE" --workdir "$WORK_DIR"; then
  echo "FAIL: collective checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time collective
echo

# Unaveraged-mode checks (deliverable 13; fel-physics.tex sec:unaveraged): the energy
# ledger, ballistic conservation and ramp handoff, fc measured against the closed
# forms in both limits and at h = 3, step-size convergence, and the priced gain-curve
# comparison against the averaged mode. The fc/faw leak grep is part of the check: the
# unaveraged path must not touch the averaged coupling quantities it measures.

echo "--- unaveraged-mode checks ------------------------------------------------------"
if grep -n "fel_und_coupling\|faw" "$SCRIPT_DIR/../code/fel_unaveraged_mod.f90" | grep -v "^[0-9]*: *!"; then
  echo "FAIL: averaged coupling quantities (fc/faw) leaked into the unaveraged path" >&2
  exit 1
fi
echo "--- fc/faw leak grep: clean"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_unaveraged.py" --exe "$EXE" --latdir "$SCRIPT_DIR/bmad" --workdir "$WORK_DIR"; then
  echo "FAIL: unaveraged checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time unaveraged

echo
echo "--- spontaneous-emission checks (FEL modes vs Bmad radiation vs analytic) -----"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_spontaneous.py" --exe "$EXE" --latdir "$SCRIPT_DIR/bmad" --workdir "$WORK_DIR"; then
  echo "FAIL: spontaneous checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time spontaneous

echo
echo "--- two-polarization checks (vector radiation, tilt, crossed undulator) -------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_two_polarization.py" --exe "$EXE" --latdir "$SCRIPT_DIR/bmad" --workdir "$WORK_DIR"; then
  echo "FAIL: two-polarization checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time two-polarization

echo
echo "--- harmonic field-set + openPMD wavefront checks ------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_harmonics.py" "$WORK_DIR/harmonics" --exe "$EXE" --genesis "$GENESIS"; then
  echo "FAIL: harmonic/openPMD checks; outputs kept in: $WORK_DIR/harmonics" >&2
  exit 1
fi
section_time harmonics

echo
echo "--- phasing checks (autophase, z_offset knob, absolute mode, chicanes) ---------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_phasing.py" "$WORK_DIR/phasing" --exe "$EXE" --genesis "$GENESIS"; then
  echo "FAIL: phasing checks; outputs kept in: $WORK_DIR/phasing" >&2
  exit 1
fi
section_time phasing

echo
echo "--- coherent-source checks (SIMPLEX hybrid: limit, claim, guards, refusals) ----"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_coherent.py" "$WORK_DIR/coherent" --exe "$EXE"; then
  echo "FAIL: coherent-source checks; outputs kept in: $WORK_DIR/coherent" >&2
  exit 1
fi
section_time coherent-source

echo
echo "--- program-structure checks (library contract, the comb, the window) ---------"
SMOKE="${EXE%lucifer}lucifer_smoke_test"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_program.py" "$WORK_DIR/program" --exe "$EXE" --smoke "$SMOKE"; then
  echo "FAIL: program-structure checks; outputs kept in: $WORK_DIR/program" >&2
  exit 1
fi
section_time program-structure

echo
echo "--- diagnostic-output checks (stats file, dumps, escaped-field bank) ----------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_diagnostics.py" --exe "$EXE" --latdir "$SCRIPT_DIR/bmad" --workdir "$WORK_DIR"; then
  echo "FAIL: diagnostic checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time diagnostics
echo

"$PYTHON" "$SCRIPT_DIR/scripts/compare_fel.py" "$WORK_DIR"
STATUS=$?
section_time tier-comparison

if [[ $STATUS -eq 0 && $KEEP_WORK_DIR -eq 0 ]]; then
  rm -rf "$WORK_DIR"
  echo "Work directory removed (pass --work-dir <path> to keep every run's outputs,"
  echo "e.g. to overlay a tier's curves: scripts/plot_fel_compare.py tdsase.diag.txt AramisTDSASE.out.h5)"
else
  echo "Outputs kept in: $WORK_DIR"
  echo "Overlay any tier's curves from there, e.g.:"
  echo "  $SCRIPT_DIR/scripts/plot_fel_compare.py tdsase.diag.txt AramisTDSASE.out.h5"
fi

exit $STATUS
