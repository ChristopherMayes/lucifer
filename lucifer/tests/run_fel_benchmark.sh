#!/bin/bash
#
# Validate the Bmad FEL tracker against Genesis 1.3 Version 4: steady state and
# time dependent with slippage.
#
# One command runs everything: Genesis over the Benchmark1-SASE lattice (writing the
# initial dumps both codes start from), a conversion of those dumps to openPMD, which is
# the only format the tracker reads, Genesis over a single undulator segment importing the
# same dumps, the same pair again time dependent (32 slices, Aramis-td.in), the Bmad
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
#   2. A genesis4 binary (CPU is fine). Default location below. Override with --genesis.
#      Built per the Genesis repository's instructions. The benchmark needs FFTW support
#      compiled in, since it runs with fft_fieldsolver = true.
#
#   3. Python with numpy and h5py: the bmad-fel-validate environment
#      (conda env create -f ../wavefront/tests/environment.yml).
#
#   4. An openPMD-beamphysics checkout carrying beamphysics/wavefront/openpmd.py, at
#      ../openPMD-beamphysics by default. It converts the Genesis reference dumps at the
#      boundary, and the harmonics section round-trips the wavefront files through that
#      class in both directions. Both refuse by name if the checkout is missing.
#
# Options:
#   --genesis <path>    genesis4 binary. Default: $GENESIS4, else genesis4 on PATH.
#   --exe <path>        lucifer binary. Default: debug then production.
#   --python <path>     Python interpreter. Default: the bmad-fel-validate conda env.
#   --work-dir <path>   Where to run. Default: a temporary directory (kept on failure).
#   --beamphysics <p>   openPMD-beamphysics checkout. Default: sibling of bmad-ecosystem.
#   --results <path>    Write a machine-readable results file (tiers, check sections,
#                       build flavor) for doc generation. Default: none written.
#
# Two of these may run at once, which is how the keystone runs them: one pass per
# build, each with its own --work-dir. They share only a source tree they read.
#
# The Genesis reference dumps are cached, since they are a pure function of the
# reference binary and the decks that make them. The cache sits under
# $LUCIFER_GENESIS_CACHE, or ~/.cache/lucifer/genesis-refs, keyed by the binary's
# bytes, every deck genesis reads and the pinned reference version. Each run prints
# whether it hit or missed. Deleting the cache is always safe and costs one cold run.
#
# Every section runs on every invocation. There is no quicker mode, on purpose: a
# cheap run that checks less is what this harness exists to prevent.
#
# Exit status is zero only if every tier passes its tolerance.

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BMAD_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

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
BEAMPHYSICS=""
RESULTS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --genesis)  GENESIS="$2";  shift 2 ;;
    --exe)      EXE="$2";      shift 2 ;;
    --python)   PYTHON="$2";   shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --beamphysics) BEAMPHYSICS="$2"; shift 2 ;;
    --results)  RESULTS="$2";  shift 2 ;;
    -h|--help)  sed -n '2,41p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

# The comparison is against Genesis. Without the binary there is nothing to compare and
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
  for candidate in "${LUCIFER_PYTHON:-}" \
                   "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
                   "$(command -v python3 2>/dev/null)"; do
    if [[ -x "$candidate" ]]; then PYTHON="$candidate"; break; fi
  done
fi
if [[ -z "$PYTHON" || ! -x "$PYTHON" ]]; then
  echo "Error: Python not found. Create the environment:" >&2
  echo "  conda env create -f $BMAD_ROOT/lucifer/wavefront/tests/environment.yml" >&2
  exit 1
fi

# openPMD-beamphysics. The boundary conversion of the Genesis reference dumps needs it,
# so it is located here rather than left to each script's own default.

if [[ -z "$BEAMPHYSICS" ]]; then
  for candidate in "${OPENPMD_BEAMPHYSICS:-}" "$BMAD_ROOT/../openPMD-beamphysics"; do
    if [[ -d "$candidate/beamphysics/wavefront" ]]; then
      BEAMPHYSICS="$(cd "$candidate" && pwd)"
      break
    fi
  done
fi
if [[ -z "$BEAMPHYSICS" || ! -d "$BEAMPHYSICS/beamphysics/wavefront" ]]; then
  echo "Error: openPMD-beamphysics checkout not found (looked for beamphysics/wavefront)." >&2
  echo "Pass it with --beamphysics <path>." >&2
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
# Which reference this is. The recorded levels belong to one version of it, and another
# moves them while every check still passes, so this is a refusal rather than a note:
# 4.6.15 both carries the CODATA electron rest energy and seeds its noise per global
# slice index, and digits recorded against anything earlier do not compare. The version
# comes from the binary's own banner, which it prints before complaining about the
# missing input file, so a binary named by --genesis identifies itself the same way a
# packaged one does.
GENESIS_VERSION_WANTED=4.6.15
GENESIS_VERSION="$("$GENESIS" 2>&1 | sed -n 's/.*Version \([0-9][0-9.]*\).*/\1/p' | head -1)"
echo "               version $GENESIS_VERSION"
if [[ "$GENESIS_VERSION" != "$GENESIS_VERSION_WANTED" ]]; then
  echo "Error: the recorded levels belong to genesis4 $GENESIS_VERSION_WANTED, and this reference reports" >&2
  echo "       ${GENESIS_VERSION:-no version at all}: $GENESIS" >&2
  echo "       Install the reference with:" >&2
  echo "         conda install -n bmad-fel-validate -c conda-forge \\" >&2
  echo "                 'genesis4=$GENESIS_VERSION_WANTED=mpi_openmpi*'" >&2
  exit 1
fi

# Which electron rest energy it was compiled with. Releases to v4.6.14 carried the
# pre-CODATA value, which loosened the transcription levels without failing anything.
# The version check above already excludes those, so this line is a statement of what
# the reference carries rather than a warning about what it might.
GENESIS_EEV="$("$PYTHON" "$SCRIPT_DIR/scripts/genesis_constants.py" "$GENESIS" 2>/dev/null)"
echo "               $GENESIS_EEV"
echo "  python:      $PYTHON"
echo "  beamphysics: $BEAMPHYSICS"
echo "  workdir:     $WORK_DIR"
echo

# The results file is for doc generation, so it records what is reproducible: which
# build flavor ran, which sections passed, and each tier's level. No timing, no path,
# no host: those would make the generated documentation differ between runs.
if [[ -n "$RESULTS" ]]; then
  case "$EXE" in
    */debug/bin/*)      BUILD_FLAVOR=debug ;;
    */production/bin/*) BUILD_FLAVOR=production ;;
    *)                  BUILD_FLAVOR=unknown ;;
  esac
  : > "$RESULTS"
  echo "build|$BUILD_FLAVOR" >> "$RESULTS"
fi

# Inputs are grouped by consumer and configuration (genesis4/<config>/*.in with the
# shared .lat files one level up, bmad/*.bmad). The run itself is flat in WORK_DIR, so
# the decks' internal lattice= references need no paths.
cp "$SCRIPT_DIR"/genesis4/*.lat "$SCRIPT_DIR"/genesis4/*/*.in \
   "$SCRIPT_DIR"/bmad/*.bmad "$WORK_DIR/"

cd "$WORK_DIR" || exit 1

# Without FI_PROVIDER=tcp the MPI runtime's provider search adds tens of seconds to
# every Genesis launch (measured, the tcp provider is always sufficient on one node).
export FI_PROVIDER=tcp

# Every section prints its wall time, so a regression in test cost is visible from the
# harness output itself.

BENCH_T_LAST=$SECONDS
section_time () {   # <label>
  echo "  [time: $1 $((SECONDS - BENCH_T_LAST)) s]"
  BENCH_T_LAST=$SECONDS
  # Reaching here means the section passed: every failure exits before it. The timing
  # is machine dependent and deliberately not recorded.
  if [[ -n "$RESULTS" ]]; then echo "section|$1|pass" >> "$RESULTS"; fi
}

# The Genesis reference runs form three independent chains -- [ss -> 1seg],
# [td -> td-1seg, td-sc, td-wake], [td-sase] -- which run concurrently (same
# decks, same single-process runs, so only the wall clock changes).

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

# The references are a pure function of the reference binary and the decks that make
# them, and they cost 131 s on both passes of every run, so they are cached. The key
# names all three things that can change the answer: the binary's own bytes, every
# deck genesis reads, and the version the harness pins, which is what a set generated
# by some other reference would fail. The cache holds whatever the chains wrote, found
# by diffing the directory rather than by a list of names that could fall behind the
# decks. It lives outside the repository, is safe to delete at any time, and the run
# says which way it went, because a cache that silently serves stale physics is worse
# than no cache. A reader requires the .complete marker rather than the directory,
# since the two build passes may run at once: the set is staged under a private name
# with its marker written last and moved into place whole, so a half-filled directory
# is never mistaken for a reference set.
GEN_CACHE_ROOT="${LUCIFER_GENESIS_CACHE:-${XDG_CACHE_HOME:-$HOME/.cache}/lucifer/genesis-refs}"
if command -v shasum >/dev/null 2>&1; then SHA="shasum -a 256"; else SHA="sha256sum"; fi
GEN_KEY="$( { echo "genesis $GENESIS_VERSION"
              $SHA "$GENESIS"
              ls -1 *.lat *.in | sort | while read -r f; do $SHA "$f"; done
            } | $SHA | cut -c1-16 )"
GEN_CACHE="$GEN_CACHE_ROOT/$GEN_KEY"

if [[ -f "$GEN_CACHE/.complete" ]]; then
  cp -p "$GEN_CACHE"/* . 2>/dev/null
  echo "  reference cache HIT  $GEN_KEY  ($GEN_CACHE_ROOT)"
else
  echo "  reference cache MISS $GEN_KEY, generating"
  ls -1 > .gen_before
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
  # Populate a private directory and move it into place, because the two build
  # passes may run concurrently: a reader tests for the directory, so it must not
  # exist until it is complete, or a concurrent pass copies half a reference set and
  # compares against it. The rename is atomic within one filesystem, and a loser of
  # the race simply finds the winner's set already there.
  GEN_STAGE="$GEN_CACHE.stage.$$"
  mkdir -p "$GEN_STAGE"
  ls -1 | grep -v '^\.gen_before$' | comm -13 .gen_before - | while read -r f; do
    [[ -f "$f" ]] && cp -p "$f" "$GEN_STAGE"/
  done
  rm -f .gen_before
  touch "$GEN_STAGE/.complete"
  if [[ -f "$GEN_CACHE/.complete" ]]; then
    rm -rf "$GEN_STAGE"
    echo "  reference cache filled by a concurrent run, discarded ours"
  elif mv "$GEN_STAGE" "$GEN_CACHE" 2>/dev/null; then
    echo "  reference cache saved $GEN_KEY"
  else
    rm -rf "$GEN_STAGE"
    echo "  reference cache could not be written, continuing without it"
  fi
fi
for log in genesis-Aramis-ss genesis-Aramis-td genesis-Aramis-td-sase; do
  tail -3 $log.log
done
section_time genesis-references
echo

# The boundary. The tracker reads openPMD only, so each chain's initial dumps are
# converted once, here, and every tier reads the converted file. Genesis keeps reading its
# own dumps, so both codes still start from the same particles and the same field: the
# conversion is exact (validated at 2.2e-15 for particles and 2.1e-16 for fields by
# convert_genesis.py round-trip, and a converted read matches a direct read at 1.7e-15).

echo "--- convert the reference dumps to openPMD (the tracker reads openPMD only) ---"
for root in Aramis AramisTD AramisTDSASE; do
  for pair in "$root-initial.par.h5 $root-initial.beam.h5" \
              "$root-initial.fld.h5 $root-initial.wf.h5"; do
    if ! "$PYTHON" "$SCRIPT_DIR/scripts/convert_genesis.py" to-openpmd $pair \
            --pyrepo "$BEAMPHYSICS"; then
      echo "FAIL: could not convert $pair" >&2
      exit 1
    fi
  done
done
section_time convert-reference-dumps
echo

# The window each chain's dumps were sliced on, read from the Genesis dump itself. An
# openPMD file carries the slice partition and not the radiation it was sliced on, so a
# deck has to state the wavelength and the sample, and taking them from the reference file
# is one truth rather than two.

read_window () {   # <genesis .par.h5>  ->  "<lambda0> <sample>"
  "$PYTHON" - "$1" <<'WINEOF'
import sys
import h5py
with h5py.File(sys.argv[1]) as f:
    lam = float(f["slicelength"][0])
    spacing = float(f["slicespacing"][0])
sample = round(spacing / lam)
assert abs(sample * lam - spacing) < 1e-9 * spacing, "the slice spacing is not a whole sample"
print(f"{lam:.12e} {sample}")
WINEOF
}

# make_nml <nml> <lattice> <out_root> <interlude_model> <dump_root> [params_extra] [beam_extra]
# Every deck reads the converted openPMD dumps and writes openPMD, the only format the
# tracker speaks. The three input groups (doc/user-guide.md): extra &fel_params content
# (wake/sc) as argument 6, extra &fel_beam_init content (check knobs) as argument 7.
# transport_model = genesis on every tier this builds: a transcription comparison needs
# Genesis4's own transverse maps, so what is left over is transcription fidelity rather
# than a transport model difference. It is a no-op for the unaveraged tier, which
# integrates the field and has no map to choose.
# nbins is the beamlet size the Genesis decks load with, which no dump format carries.

make_nml () {
  local lam sample
  read -r lam sample <<<"$(read_window "$5-initial.par.h5")"
  cat > "$1" <<NML
&fel_params
  lat_file = "$2"
  global%out_root = "$3"
  global%interlude_model = "$4"
  global%transport_model = "genesis"
  global%write_diag = T
${6:+  $6}
/
&fel_beam_init
  beam_file = "$5-initial.beam.h5"
  beamlet_size = 8
${7:+  $7}
/
&fel_wavefront_init
  field_file = "$5-initial.wf.h5"
  wavefront_init%lambda0 = $lam
  wavefront_init%window_sample = $sample
/
NML
}

make_nml tier1.nml  aramis_1seg.bmad tier1  bmad    Aramis
make_nml tier2.nml  aramis.bmad      tier2  bmad    Aramis
make_nml tier2g.nml aramis.bmad      tier2g genesis Aramis
make_nml tier1s.nml aramis_1seg.bmad tier1s bmad    Aramis "" "split_weights = T"
make_nml tier1u.nml aramis_1seg_unavg.bmad tier1u bmad Aramis
make_nml td1.nml    aramis_1seg.bmad td1    bmad    AramisTD
make_nml td2.nml    aramis.bmad      td2    bmad    AramisTD
make_nml td2g.nml   aramis.bmad      td2g   genesis AramisTD
make_nml tdsase.nml aramis.bmad      tdsase genesis AramisTDSASE
# The space-charge tier says where space charge acts on the element, as any Bmad program
# does, and turns on Bmad's master switch in the deck. The solver's numbers stay in
# space_charge%.

cat > aramis_1seg_sc.bmad <<'LAT'
call, file = aramis_1seg.bmad
wiggler::*[SPACE_CHARGE_METHOD] = slice
LAT

make_nml tdsc.nml   aramis_1seg_sc.bmad tdsc   genesis AramisTD "space_charge%rmax = 250e-6
  space_charge%nz = 2
  space_charge%nphi = 1
  space_charge%longrange = T
  bmad_com%csr_and_space_charge_on = T"
make_nml tdwk.nml   aramis_1seg.bmad tdwk   genesis AramisTD "chamber_wake%on = T
  chamber_wake%radius = 2.5e-3
  chamber_wake%conductivity = 5.813e7
  chamber_wake%relaxation = 8.1e-6
  chamber_wake%gap = 0.5e-3
  chamber_wake%lgap = 0.015
  chamber_wake%hrough = 100e-9
  chamber_wake%lrough = 100e-6"

# Assertion checks: a lattice whose FEL element is missing b_max, missing l_period, or
# uses a fieldmap field_calc, or that carries Bmad wakes on
# any element (the slice-at-a-time seam cannot apply them meaningfully), must be
# Refused by name -- the failure message names the offending attribute and element -- not passed
# through to fail downstream with an unrelated message (missing b_max) or a segfault in
# the parse-time reference tracking (fieldmap). Each check mutates one attribute of the
# real single-segment lattice and requires both a nonzero exit and the by-name message.
#
# The fieldmap case is refused by Bmad rather than by this program: an FEL method is not
# valid on a wiggler whose field is a map, so the lattice does not parse. That is one step
# earlier than the program's own assertion, and the assertion stays as the second line of
# defense for a lattice built through the API rather than parsed.

echo "--- FEL-element assertion checks (refusal by name) -----------------------------"
grep -v "b_max" aramis_1seg.bmad                                   > refusal_bmax.bmad
sed 's/l_period = 0.015, //' aramis_1seg.bmad                      > refusal_lperiod.bmad
sed 's/field_calc = helical_model/field_calc = fieldmap/' aramis_1seg.bmad > refusal_fieldmap.bmad

CHECKS_OK=1
run_assert_refusal () {   # <name> <by-name message fragment> [extra &fel_params lines]
  make_nml refusal_$1.nml refusal_$1.bmad refusal_$1 bmad Aramis "${3:-}"
  if "$EXE" refusal_$1.nml > fel-refusal_$1.log 2>&1; then
    echo "FAIL: refusal_$1 lattice was accepted (exit 0); it must be refused" >&2
    CHECKS_OK=0
  elif ! grep -q "$2" fel-refusal_$1.log; then
    echo "FAIL: refusal_$1 refused, but not by name; log tail:" >&2
    tail -5 fel-refusal_$1.log >&2
    CHECKS_OK=0
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

# Space charge asked for with neither term configured: the solve would cost its full
# price and return an exact zero, so which term was meant is worth asking. The other
# space-charge refusal, a three-dimensional method on an FEL element, is covered by the
# collective checks.

cat > refusal_scnone.bmad <<'LAT'
call, file = aramis_1seg.bmad
wiggler::*[SPACE_CHARGE_METHOD] = slice
LAT

run_assert_refusal bmax     "ZERO B_MAX"
run_assert_refusal lperiod  "ZERO L_PERIOD"
run_assert_refusal fieldmap "NOT A VALID TRACKING_METHOD"
run_assert_refusal lrwake   "LR (MULTI-BUNCH) WAKES ARE NOT SUPPORTED"
run_assert_refusal zmax     "Z_MAX CAN HANDLE"
run_assert_refusal scnone   "NEITHER SPACE-CHARGE TERM IS CONFIGURED" \
                            "bmad_com%csr_and_space_charge_on = T"

# tracking_method = custom on a wiggler means some other program's tracking, so this
# program must not claim the element. With the only wiggler tracked that way there is no
# FEL segment left, which is refused by name: the alternative, silently tracking it as an
# FEL element, is what the named methods exist to prevent.

sed 's/tracking_method = fel_averaged/tracking_method = custom/' aramis_1seg.bmad > custom_seam.bmad
make_nml custom_seam.nml custom_seam.bmad custom_seam bmad Aramis
if "$EXE" custom_seam.nml > fel-custom_seam.log 2>&1; then
  echo "FAIL: a custom-tracked wiggler was claimed as an FEL segment" >&2
  CHECKS_OK=0
elif ! grep -q "NO FEL ELEMENT" fel-custom_seam.log; then
  echo "FAIL: custom-tracked wiggler refused, but not by name; log tail:" >&2
  tail -5 fel-custom_seam.log >&2
  CHECKS_OK=0
else
  echo "  custom_seam: a custom-tracked wiggler is not an FEL segment, refused by name"
fi

if [ "$CHECKS_OK" -ne 1 ]; then
  echo "FAIL: FEL-element assertion checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time refusals
echo

# The documented tier numbers are single-thread runs. Results must not depend on the
# thread count, and the explicit check for that follows the loop.

export OMP_NUM_THREADS=1

# The tiers are independent single-thread processes over the shared dumps. They
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
# eight threads must reproduce the one-thread run bit for bit -- each slice's arithmetic
# is independent of which thread runs it, so any difference at all is a race. The diag
# file is compared byte for byte. The dumps dataset by dataset, exactly (HDF5 object
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
for fa, fb in (("td1-final.wf.h5", "td1t8-final.wf.h5"),
               ("td1-final.beam.h5", "td1t8-final.beam.h5")):
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

# Shot-noise checks. The statistical check is self-referenced
# (no cross-code reference exists: Genesis cannot represent weighted noise). The SASE startup cross-check pits the
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

# Seam-wake checks: element sr wakes across the whole window --
# closed-form pseudomode ramp, exact causality with the d8 direction cross-check, the
# z_long kernel cross-validation against the wake model (first-principles
# tight, resolved-beam at the derived boundary bound), split-weight invariance and
# thread determinism. Self-referenced. Needs no Genesis.

echo "--- seam-wake checks -------------------------------------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_seam_wake.py" --exe "$EXE" --workdir "$WORK_DIR"; then
  echo "FAIL: seam-wake checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time seam-wake
echo

# Distribution-import checks: the bunch_struct resampler transcribed
# from Genesis's SDDSBeam.cpp -- exact where no RNG enters (the per-slice current
# profile against Genesis importing the SAME distribution file, the match transform
# hitting its Twiss targets, split-weight invariance, thread determinism), statistical
# where the resampling RNG forces it (slice Twiss recovery, startup power cross-code).

echo "--- distribution-import checks --------------------------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_import.py" --exe "$EXE" --genesis "$GENESIS" --workdir "$WORK_DIR" --seeds 4; then
  echo "FAIL: distribution-import checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time import
echo

# Slice-migration checks: conservation under heavy migration, exact
# phase continuity of the moves, and no-op bit identity (self-referenced
# -- Genesis migrates only under one4one).

echo "--- slice-migration checks ------------------------------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_migration.py" --exe "$EXE" --workdir "$WORK_DIR"; then
  echo "FAIL: migration checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time migration
echo

# Collective-effects self-referenced checks: exact energy bookkeeping of
# the wake's eloss on a cold dark beam, and the stale-wake structural check (the
# convolution must follow the currents under migration).

echo "--- collective-effects checks ---------------------------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_collective.py" --exe "$EXE" --workdir "$WORK_DIR"; then
  echo "FAIL: collective checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time collective
echo

# FP32 lockstep checks: the single-precision particle path's divergence from the
# FP64 reference lands inside its recorded ceilings, the instrument is read-only on
# the FP64 physics (diag byte-identical on vs off), its mutation hook moves the
# recorded level so the check can fail, freerun measures the compounding rate, and
# configurations the twin does not cover are refused by name.

echo "--- FP32 lockstep checks --------------------------------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_fp32.py" --exe "$EXE" --workdir "$WORK_DIR"; then
  echo "FAIL: FP32 lockstep checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time fp32-lockstep
echo

# Device backend checks: the Metal backend judged by the lockstep instrument with
# the device in the twin's role, inside its recorded ceilings; the exact-wrap
# assertion on the fixed-point phase; read-only on the FP64 physics; the perturbed
# kernel constant moves a recorded level so the check can fail; the freerun twin
# and the resident production run land in the same end-to-end band; and everything
# the backend does not cover is refused by name.

echo "--- Device backend checks -------------------------------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_device.py" --exe "$EXE" --workdir "$WORK_DIR"; then
  echo "FAIL: device backend checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time device
echo

# Unaveraged-mode checks (fel-physics.md sec-unaveraged): the energy
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
echo "--- particle dump format checks (openPMD and Genesis .par) ---------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_beam_format.py" --exe "$EXE" --workdir "$WORK_DIR"; then
  echo "FAIL: beam-format checks; outputs kept in: $WORK_DIR" >&2
  exit 1
fi
section_time beam-format

echo
echo "--- harmonic field-set + openPMD wavefront checks ------------------------------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_harmonics.py" "$WORK_DIR/harmonics" --exe "$EXE" --genesis "$GENESIS" --pyrepo "$BEAMPHYSICS"; then
  echo "FAIL: harmonic/openPMD checks; outputs kept in: $WORK_DIR/harmonics" >&2
  exit 1
fi
section_time harmonics

echo
echo "--- phasing checks (autophase, z_offset knob, absolute mode, chicanes) ---------"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_phasing.py" "$WORK_DIR/phasing" --exe "$EXE" --genesis "$GENESIS" --pyrepo "$BEAMPHYSICS"; then
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
echo "--- input-reference conformance (the stated defaults against the declarations) ---"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_input_reference.py"; then
  echo "FAIL: the input reference disagrees with the struct declarations" >&2
  exit 1
fi
section_time input-reference

# The examples are the feature list's evidence, so the two must agree. The generated
# pages quote the decks verbatim, and regenerating them here is what keeps a page from
# drifting: a deck edited without a regeneration leaves a diff the keystone refuses.

echo
echo "--- examples conformance (the feature matrix, the directories, the pages) -----"
if ! "$PYTHON" "$SCRIPT_DIR/scripts/report_examples.py" \
       --examples "$BMAD_ROOT/lucifer/examples" \
       --out "$BMAD_ROOT/lucifer/doc/generated/examples" > /dev/null; then
  echo "FAIL: the example pages could not be generated" >&2
  exit 1
fi
if ! "$PYTHON" "$SCRIPT_DIR/scripts/check_examples.py"; then
  echo "FAIL: the examples, the feature matrix and the generated pages disagree" >&2
  exit 1
fi
section_time examples

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

COMPARE_ARGS=()
if [[ -n "$RESULTS" ]]; then COMPARE_ARGS+=(--results "$RESULTS"); fi
"$PYTHON" "$SCRIPT_DIR/scripts/compare_fel.py" "$WORK_DIR" "${COMPARE_ARGS[@]}"
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
