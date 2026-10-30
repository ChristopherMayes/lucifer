#!/bin/bash
# run_examples.sh runs every example and saves the figure the documentation shows.
#
# It does two jobs in one pass. It answers whether all the examples still run, printing
# one line per input file and exiting nonzero if any input file failed. It also writes the figures
# the generated pages reference, one per example, into docs/generated/examples/.
#
# The figures are committed. Building an FEL run at documentation build time would be
# wrong, and a page whose figure is missing is worse than a page with an old one. Their
# bytes are not compared: a plotting library upgrade rewrites every PNG without changing
# any physics, so check_examples.py asserts that each figure exists and the regeneration
# check covers the Markdown only.
#
# saturation_demo is not run here. It needs the genesis4 binary and its own clock, and it
# has run.sh for that. Where run.sh finds its files is checked here, through its
# --show-locations option and from a directory outside both checkouts.
#
# Usage:
#
#   run_examples.sh [--exe <lucifer>] [--python <python3>] [--dpi <n>] [--no-figures]
#                   [--jobs <n>]
#
# Run it from anywhere. Outputs land in each example's own directory, where .gitignore
# covers them.

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# The next line sources the locations, LUCIFER_DIR, DEBUG_BIN, PRODUCTION_BIN and
# BMAD_SRC, with their defaults and the names that set them elsewhere (util/locations.sh).
source "$SCRIPT_DIR/../util/locations.sh"
FIG_DIR="$LUCIFER_DIR/docs/generated/examples"

EXE=""
PYTHON=""
FIGURES=1
JOBS=0

# The committed figures are read on a page rather than zoomed into, and they live in the
# repository's history, so they are written below the default resolution of
# python -m lucifer.plot.
DPI=110

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exe)        EXE="$2";    shift 2 ;;
    --python)     PYTHON="$2"; shift 2 ;;
    --no-figures) FIGURES=0;   shift ;;
    --dpi)        DPI="$2";     shift 2 ;;
    --jobs)       JOBS="$2";   shift 2 ;;
    -h|--help)    sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Error: unknown argument: $1" >&2; exit 1 ;;
  esac
done

if ! [[ "$JOBS" =~ ^[0-9]+$ ]]; then
  echo "Error: --jobs takes a whole number, and 0 means the default. Got: $JOBS" >&2
  exit 1
fi

# The production build is the default, since the wall times the READMEs give are its own.
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

if [[ -z "$PYTHON" ]]; then
  for candidate in "${LUCIFER_PYTHON:-}" \
                   "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
                   "$(command -v python3 2>/dev/null)"; do
    if [[ -x "$candidate" ]]; then PYTHON="$candidate"; break; fi
  done
fi
# The record of each run and the figures both need the Python with h5py, and the figures
# need the package with its extra, lucifer[plot].
if [[ -z "$PYTHON" || ! -x "$PYTHON" ]]; then
  echo "Error: Python not found, and the figures need h5py and matplotlib." >&2
  echo "  Create it with: conda env create -f $LUCIFER_DIR/tests/environment-validate.yml" >&2
  echo "  Or run with --no-figures to only check that the examples run." >&2
  exit 1
fi

# The block below sets how many directories run at once, and says why.
#
# The unit is a directory, never an input file. The import example's second input file reads the
# openPMD file its first input file writes, so a directory's input files stay in order inside one
# worker. Directories share only the lattice files they read, so several at once change
# no result. The rows are printed after every worker has finished, in the order a
# one-worker run prints them, so the table does not depend on the schedule.
#
# The count is small because one example run is already parallel. Each run takes the
# OpenMP runtime's own thread count, which this script does not set and the run header
# prints (12 on the machine below), and one run alone draws 8.4 cores on average and
# 11.7 at its peak. A fourth worker therefore finds little idle capacity left. The table
# was measured on an M3 Max of 12 performance cores, with the examples alone and no other
# keystone job running:
#
#   workers          1     2     3     4     6
#   seconds        390   293   262   255   248
#   peak RSS (GB)  0.5   0.7   0.9   1.3   1.6
#
# The third worker buys 31 s and the fourth buys 7 s, so the default is one worker per
# four performance cores, which is three here.
#
# The count is also a static allocation against the rest of the keystone, which runs
# five jobs at once (docs/validation.md). The two validation passes take eleven concurrent
# single-thread tier processes each while their tiers run, the regression suite takes
# one and the wavefront validation one. The suite spends 2*11 + JOBS + 2 processes at its
# widest, which is above the core count. A larger share here is taken from the passes
# that set the keystone's critical path.
if [[ "$JOBS" -eq 0 ]]; then
  CORES="$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || true)"
  [[ -n "${CORES:-}" ]] || CORES="$(sysctl -n hw.ncpu 2>/dev/null || true)"
  [[ -n "${CORES:-}" ]] || CORES="$(nproc 2>/dev/null || true)"
  [[ -n "${CORES:-}" ]] || CORES=4
  JOBS=$((CORES / 4))
  [[ $JOBS -lt 1 ]] && JOBS=1
  [[ $JOBS -gt 3 ]] && JOBS=3
fi

cd "$LUCIFER_DIR/examples" || exit 1
mkdir -p "$FIG_DIR"

STATE="$(mktemp -d "${TMPDIR:-/tmp}/run_examples.XXXXXX")"
RECORD="$SCRIPT_DIR/../tests/scripts/example_record.py"
BANK="$SCRIPT_DIR/../tests/scripts/bank_bytes.py"
trap 'rm -rf "$STATE"' EXIT

echo "=============================================================================="
echo " Examples"
echo "   lucifer: $EXE"
echo "   workers: $JOBS directories at once"
[[ $FIGURES -eq 1 ]] && echo "   figures: $FIG_DIR"
echo "------------------------------------------------------------------------------"
printf '%-20s %-26s %6s %6s  %s\n' DIR INPUT EXIT SEC EXIT_LINE

FAILED=0

# The directory list and each directory's input files are settled before anything launches, so
# a directory with no input file is a failure found without running one.
DIRS=()
for dir in */; do
  dir="${dir%/}"
  [[ "$dir" == "saturation_demo" ]] && continue
  # The input files run lucifer.toml first, then the rest in byte order. A directory's input files can
  # depend on each other: import's lucifer_from_file.toml reads the bunch its lucifer.toml
  # writes. The glob's own order follows the locale's collation, under which
  # lucifer_from_file.toml sorts ahead of lucifer.toml, and a bunch file an earlier run left
  # in the tree hides that.
  shopt -s nullglob
  inputs=()
  [[ -f "$dir/lucifer.toml" ]] && inputs=("$dir/lucifer.toml")
  while IFS= read -r d; do
    inputs+=("$d")
  done < <(LC_ALL=C ls -1 "$dir"/lucifer*.toml 2>/dev/null | LC_ALL=C sort | grep -vx "$dir/lucifer.toml")
  shopt -u nullglob
  if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "Error: $dir has no lucifer*.toml input file." >&2
    FAILED=1
    continue
  fi
  printf '%s\n' "${inputs[@]}" > "$STATE/inputs.${#DIRS[@]}"
  DIRS+=("$dir")
done

# run_dir does one directory's work, its input files in order and then its figure. Every input file
# leaves a record carrying the worker's own verdict, so the table and the exit status
# cannot disagree, and an input file that left no record at all is a failure the parent reports
# by name.
run_dir () {   # <dir> <index>
  local dir="$1" idx="$2"
  local k=0 status=0 input_path input log rc dt t0 line verdict primary root stats bank
  local inputs=()
  # The input file list is read into an array first. A while-read loop would hold the list on
  # the loop's stdin, and a tracker that read stdin would eat the input files behind it.
  while IFS= read -r input_path; do inputs+=("$input_path"); done < "$STATE/inputs.$idx"
  for input_path in "${inputs[@]}"; do
    input="$(basename "$input_path")"
    log="${input%.toml}.log"
    t0=$SECONDS
    ( cd "$dir" && "$EXE" "$input" > "$log" 2>&1 < /dev/null )
    rc=$?
    dt=$((SECONDS - t0))
    line="$(grep -m1 '^ Exit' "$dir/$log" 2>/dev/null | sed 's/^ *Exit *//')"
    # A completed run always prints its exit line, so its absence is a failure even when
    # the status says otherwise. A stale binary that rejects a lattice it does not
    # understand would otherwise look like success.
    verdict=ok
    if [[ $rc -ne 0 || -z "$line" ]]; then
      line="$(grep -m1 -E 'ERROR|FATAL' "$dir/$log" 2>/dev/null | head -c 90)"
      verdict=fail
      status=1
    fi
    printf '%s\n%s\n%s\n%s\n' "$verdict" "$rc" "$dt" "$line" > "$STATE/out.$idx.$k"
    # The record holds what the run measured, from its stats file and its log, one row
    # per quantity. The rows of every input file become examples/measured.tsv below, and the
    # pages carry them.
    if [[ "$verdict" == "ok" ]]; then
      if ! "$PYTHON" "$RECORD" "$dir/$input" "$dir/$log" >> "$STATE/rec.$idx" 2> "$STATE/recerr.$idx"; then
        echo "Error: $dir/$input: the record could not be written; see $STATE/recerr.$idx" > "$STATE/note.$idx"
        status=1
      fi
      # The bytes the header states for the escaped-field bank, against the files the
      # run kept, and no pulse or escaped file where the header states no bank.
      if ! bank="$("$PYTHON" "$BANK" "$dir/$log" "$dir" 2>&1)"; then
        printf 'Error: %s/%s: the bank the header states is not the files on disk:\n%s\n' \
               "$dir" "$input" "$bank" > "$STATE/note.$idx"
        status=1
      fi
    fi
    k=$((k + 1))
  done

  # Each example gets one figure, from the primary input file's stats file. The out_root is read from
  # the input file rather than assumed equal to the directory name, because several examples
  # deliberately differ (chicane writes "chicane", coherent_source writes "coherent").
  # It runs in this worker, after this directory's input files, because it reads what they wrote.
  if [[ $FIGURES -eq 1 ]]; then
    primary="$dir/lucifer.toml"
    [[ -f "$primary" ]] || primary="${inputs[0]}"
    root="$(sed -n 's/.*out_root[[:space:]]*=[[:space:]]*"\([^"]*\)".*/\1/p' "$primary" | head -1)"
    stats="$dir/$root.stats.h5"
    if [[ -z "$root" || ! -f "$stats" ]]; then
      echo "Error: $dir: no stats file at $stats to plot." > "$STATE/note.$idx"
      status=1
    elif ! "$PYTHON" -m lucifer.plot "$stats" --dpi "$DPI" -o "$FIG_DIR/$dir.png" \
           > "$dir/plot.log" 2>&1 < /dev/null; then
      echo "Error: $dir: python -m lucifer.plot failed; see $dir/plot.log" > "$STATE/note.$idx"
      status=1
    fi
  fi
  return $status
}

# Each worker holds a slot in a fifo, because this runs under bash 3.2, where "wait -n"
# does not exist. A worker takes a slot before it launches and returns it as it ends, whatever
# it returns, so a failing directory frees its slot like any other. A worker killed by a
# signal returns nothing and the launch loop then waits, which is the tier loop's
# behaviour too and what the keystone's timeout covers (tests/test_keystone.py).
SLOTS="$STATE/slots"
mkfifo "$SLOTS"
exec 9<>"$SLOTS"
islot=0
while [[ $islot -lt $JOBS ]]; do printf '.\n' >&9; islot=$((islot + 1)); done

PIDS=()
idir=0
while [[ $idir -lt ${#DIRS[@]} ]]; do
  IFS= read -r -u 9 slot || { echo "Error: a worker slot could not be taken." >&2; exit 1; }
  ( run_dir "${DIRS[$idir]}" "$idir"; rc=$?; printf '.\n' >&9; exit $rc ) &
  PIDS+=("$!")
  idir=$((idir + 1))
done

# Every worker is waited for by pid, so a worker that ended without leaving records is
# distinguishable from one that left failing records, and neither is mistaken for a pass.
if [[ ${#PIDS[@]} -gt 0 ]]; then
  for pid in "${PIDS[@]}"; do
    if ! wait "$pid"; then FAILED=1; fi
  done
fi
exec 9>&-

idir=0
while [[ $idir -lt ${#DIRS[@]} ]]; do
  dir="${DIRS[$idir]}"
  k=0
  while IFS= read -r input_path; do
    input="$(basename "$input_path")"
    if [[ -f "$STATE/out.$idir.$k" ]]; then
      verdict=""; rc=""; dt=""; line=""
      { read -r verdict; read -r rc; read -r dt; IFS= read -r line || line=""; } \
          < "$STATE/out.$idir.$k"
      [[ "$verdict" == "ok" ]] || FAILED=1
      printf '%-20s %-26s %6s %6s  %s\n' "$dir" "$input" "$rc" "$dt" "$line"
    else
      printf '%-20s %-26s %6s %6s  %s\n' "$dir" "$input" "-" "-" "no outcome recorded"
      FAILED=1
    fi
    k=$((k + 1))
  done < "$STATE/inputs.$idir"
  if [[ -f "$STATE/note.$idir" ]]; then cat "$STATE/note.$idir" >&2; FAILED=1; fi
  idir=$((idir + 1))
done

# The record of every run is assembled in one order, whatever order the workers finished
# in, so two runs of the same binary write the same bytes. The keystone regenerates
# this file and diffs it with the generated pages. The sample output the user guide
# includes comes from the steady-state example's log the same way.
if [[ $FAILED -eq 0 ]]; then
  {
    printf 'directory\tinput\tquantity\tvalue\tunit\n'
    cat "$STATE"/rec.* 2>/dev/null | LC_ALL=C sort
  } > measured.tsv
  if ! "$PYTHON" "$RECORD" --sample steady_state/lucifer.log --out "$LUCIFER_DIR/docs/generated/sample-output.md"; then
    echo "Error: the sample output could not be written from steady_state/lucifer.log" >&2
    FAILED=1
  fi
fi

# The saturation demo's locations are checked from a directory outside both checkouts. The binary
# directories are named under names the defaults never produce, through the two
# variables and relative to that directory, and the executable is given once as a
# relative path and once left to the variables. run.sh makes each absolute and checks it
# from a scratch directory of its own, so a path that held only relative to where it
# started fails here the way the real run would.

echo "------------------------------------------------------------------------------"
echo " saturation_demo: locations from outside both checkouts (no physics runs)"
loc="$(mktemp -d "${TMPDIR:-/tmp}/sat_launch.XXXXXX")"
mkdir -p "$loc/build"
ln -s "$(dirname "$EXE")" "$loc/build/opt-bin"
ln -s "$(dirname "$EXE")" "$loc/build/dbg-bin"
for exe_arg in "--exe build/opt-bin/lucifer" ""; do
  # shellcheck disable=SC2086
  out="$(cd "$loc" && LUCIFER_DEBUG_BIN=build/dbg-bin LUCIFER_PRODUCTION_BIN=build/opt-bin \
         "$LUCIFER_DIR/examples/saturation_demo/run.sh" --show-locations --workers 2 \
         --python "$PYTHON" --work-dir out $exe_arg 2>&1)"
  rc=$?
  got="$(printf '%s\n' "$out" | awk '$1 == "lucifer" {print $2}')"
  wdir="$(printf '%s\n' "$out" | awk '$1 == "work" && $2 == "dir" {print $3}')"
  how="${exe_arg:-from LUCIFER_PRODUCTION_BIN}"
  if [[ $rc -eq 0 && "$got" == "$loc/build/opt-bin/lucifer" && "$wdir" == "$loc/out" ]]; then
    echo "   ok    executable $how resolves to $got"
  else
    echo "   FAIL  executable $how: exit $rc, lucifer '$got', work dir '$wdir'"
    printf '%s\n' "$out" | sed 's/^/         /'
    FAILED=1
  fi
done
rm -rf "$loc"

echo "------------------------------------------------------------------------------"
if [[ $FAILED -ne 0 ]]; then
  echo " FAIL: an example did not run, its figure was not written, or the saturation demo's"
  echo "       locations did not resolve."
  exit 1
fi

# The Python page's plotting examples run on the sase example's run, drawing the page's
# figures into docs/generated/python with figures on and into the scratch directory
# under --no-figures, so the page's code runs on every examples job while its figures,
# like the examples', stay out of the regeneration diff. They are drawn before the link
# check below, since the page links to them.
if [[ $FIGURES -eq 1 ]]; then PY_FIGS="$LUCIFER_DIR/docs/generated/python"; else PY_FIGS="$STATE/python"; fi
if ! "$PYTHON" "$SCRIPT_DIR/../tests/scripts/python_page_figures.py" --stats-dir "$SCRIPT_DIR/sase" \
       --out "$PY_FIGS" > "$STATE/python-figures.log" 2>&1 < /dev/null; then
  echo " FAIL: the Python page's figures could not be drawn:"
  tail -20 "$STATE/python-figures.log"
  exit 1
fi
tail -1 "$STATE/python-figures.log"

# The notebooks under docs/notebooks are executed on the steady-state example, into
# docs/notebooks with figures on and into the scratch directory under --no-figures, so
# every examples job runs them while their stored outputs stay out of the regeneration
# diff, as the figures do. Then their cells are checked: short, shown rather than
# printed, and executed in order.
if [[ $FIGURES -eq 1 ]]; then NB_OUT="$LUCIFER_DIR/docs/notebooks"; else NB_OUT="$STATE/notebooks"; fi
if ! "$PYTHON" "$SCRIPT_DIR/../tests/scripts/run_notebooks.py" --out "$NB_OUT" --exe "$EXE" \
       > "$STATE/notebooks.log" 2>&1 < /dev/null; then
  echo " FAIL: a notebook did not execute:"
  tail -40 "$STATE/notebooks.log"
  exit 1
fi
grep '^run_notebooks:' "$STATE/notebooks.log"
if ! "$PYTHON" "$SCRIPT_DIR/../tests/scripts/check_notebooks.py" --root "$LUCIFER_DIR"; then
  echo " FAIL: a notebook's cells are not in the form the check holds them to."
  exit 1
fi

# Every relative link in the committed Markdown files resolves. The READMEs here are
# what the example pages are generated from, with their links retargeted, so a dead
# link in one reads correctly on the site and is caught nowhere else.
if ! "$PYTHON" "$SCRIPT_DIR/../tests/scripts/check_links.py" --root "$LUCIFER_DIR"; then
  echo " FAIL: a relative link in a committed Markdown file does not resolve."
  exit 1
fi

# And no page, script or comment carries a word the repository has retired.
if ! "$PYTHON" "$SCRIPT_DIR/../tests/scripts/check_retired_words.py" --root "$LUCIFER_DIR"; then
  echo " FAIL: a committed file carries a retired word."
  exit 1
fi

# And every key with a set of names in the key table is compared against those names and no
# other, so the table and the code that acts on it cannot drift apart.
if ! "$PYTHON" "$SCRIPT_DIR/../tests/scripts/check_allowed_values.py" --root "$LUCIFER_DIR"; then
  echo " FAIL: a name-valued key's set and the names the code compares it against disagree."
  exit 1
fi
# The plotting module runs even without figures, on the steady-state example's stats file
# into a scratch image, so a plotter that fails to import or to draw fails this run.
if ! "$PYTHON" -m lucifer.plot "$SCRIPT_DIR/steady_state/steady_state.stats.h5" --dpi 40 \
       -o "$STATE/plot-check.png" > "$STATE/plot-check.log" 2>&1 < /dev/null; then
  echo " FAIL: python -m lucifer.plot failed on the steady-state example:"
  tail -20 "$STATE/plot-check.log"
  exit 1
fi

if [[ $FIGURES -eq 1 ]]; then
  echo " All examples ran. Figures in docs/generated/examples/."
  echo " Regenerate the pages with tests/scripts/report_examples.py after changing an input file."
else
  echo " All examples ran. No figures written (--no-figures)."
fi
echo "=============================================================================="
