#!/bin/bash
# Run every example and save the figure the documentation shows.
#
# Two jobs in one pass. It is the answer to "do all the examples still run", printing one
# line per deck and exiting nonzero if any deck failed. And it produces the figures the
# generated pages reference, one per example, into doc/generated/examples/.
#
# The figures are committed. Building an FEL run at documentation build time would be
# wrong, and a page whose figure is missing is worse than a page with an old one. Their
# bytes are not compared: a plotting library upgrade rewrites every PNG without changing
# any physics, so check_examples.py asserts that each figure exists and the regeneration
# check covers the Markdown only.
#
# saturation_demo is not run here. It needs the genesis4 binary and its own clock, and it
# has run.sh for that.
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
BMAD_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
FIG_DIR="$BMAD_ROOT/lucifer/doc/generated/examples"

EXE=""
PYTHON=""
FIGURES=1
JOBS=0

# The committed figures are read on a page rather than zoomed into, and they live in the
# repository's history, so they are written below plot_fel.py's own default resolution.
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

# Production by default: these timings are the ones the READMEs quote.
if [[ -z "$EXE" ]]; then
  for candidate in "$BMAD_ROOT/production/bin/lucifer" "$BMAD_ROOT/debug/bin/lucifer"; do
    if [[ -x "$candidate" ]]; then EXE="$candidate"; break; fi
  done
fi
if [[ -z "$EXE" || ! -x "$EXE" ]]; then
  echo "Error: lucifer not found. Build with: cd $BMAD_ROOT && ./util/conda_compile" >&2
  exit 1
fi

if [[ -z "$PYTHON" ]]; then
  for candidate in "${LUCIFER_PYTHON:-}" \
                   "$(conda info --base 2>/dev/null)/envs/bmad-fel-validate/bin/python3" \
                   "$(command -v python3 2>/dev/null)"; do
    if [[ -x "$candidate" ]]; then PYTHON="$candidate"; break; fi
  done
fi
if [[ $FIGURES -eq 1 && ( -z "$PYTHON" || ! -x "$PYTHON" ) ]]; then
  echo "Error: Python not found, and the figures need h5py and matplotlib." >&2
  echo "  Create it with: conda env create -f $BMAD_ROOT/lucifer/wavefront/tests/environment.yml" >&2
  echo "  Or run with --no-figures to only check that the examples run." >&2
  exit 1
fi

# How many directories run at once, and why the number is what it is.
#
# The unit is a directory, never a deck. The import example's second deck reads the
# openPMD file its first deck writes, so a directory's decks stay in order inside one
# worker. Directories share only the lattice files they read, so several at once change
# no result. The rows are printed after every worker has finished, in the order a
# one-worker run prints them, so the table does not depend on the schedule.
#
# The count is small because one example run is already parallel. Each run takes the
# OpenMP runtime's own thread count, which this script does not set and the run header
# prints (12 on the machine below), and one run alone draws 8.4 cores on average and
# 11.7 at its peak. A fourth worker therefore finds little idle capacity left. Measured
# on an M3 Max of 12 performance cores, the examples alone with no other keystone job
# running, 2026-09-13:
#
#   workers          1     2     3     4     6
#   seconds        390   293   262   255   248
#   peak RSS (GB)  0.5   0.7   0.9   1.3   1.6
#
# The third worker buys 31 s and the fourth buys 7 s, so the default is one worker per
# four performance cores, which is three here.
#
# The count is also a static allocation against the rest of the keystone, which runs
# five jobs at once (doc/validation.md). The two benchmark passes take eleven concurrent
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

cd "$BMAD_ROOT/lucifer/examples" || exit 1
mkdir -p "$FIG_DIR"

STATE="$(mktemp -d "${TMPDIR:-/tmp}/run_examples.XXXXXX")"
trap 'rm -rf "$STATE"' EXIT

echo "=============================================================================="
echo " Examples"
echo "   lucifer: $EXE"
echo "   workers: $JOBS directories at once"
[[ $FIGURES -eq 1 ]] && echo "   figures: $FIG_DIR"
echo "------------------------------------------------------------------------------"
printf '%-20s %-26s %6s %6s  %s\n' DIR DECK EXIT SEC EXIT_LINE

FAILED=0

# The directory list and each directory's decks are settled before anything launches, so
# a directory with no deck is a failure found without running one.
DIRS=()
for dir in */; do
  dir="${dir%/}"
  [[ "$dir" == "saturation_demo" ]] && continue
  shopt -s nullglob
  decks=("$dir"/lucifer*.in)
  shopt -u nullglob
  if [[ ${#decks[@]} -eq 0 ]]; then
    echo "Error: $dir has no lucifer*.in deck." >&2
    FAILED=1
    continue
  fi
  printf '%s\n' "${decks[@]}" > "$STATE/decks.${#DIRS[@]}"
  DIRS+=("$dir")
done

# One directory's work: its decks in order, then its figure. Every deck leaves a record
# carrying the worker's own verdict, so the table and the exit status cannot disagree,
# and a deck that left no record at all is a failure the parent reports by name.
run_dir () {   # <dir> <index>
  local dir="$1" idx="$2"
  local k=0 status=0 deck_path deck log rc dt t0 line verdict primary root stats
  local decks=()
  # The deck list is read into an array first. A while-read loop would hold the list on
  # the loop's stdin, and a tracker that read stdin would eat the decks behind it.
  while IFS= read -r deck_path; do decks+=("$deck_path"); done < "$STATE/decks.$idx"
  for deck_path in "${decks[@]}"; do
    deck="$(basename "$deck_path")"
    log="${deck%.in}.log"
    t0=$SECONDS
    ( cd "$dir" && "$EXE" "$deck" > "$log" 2>&1 < /dev/null )
    rc=$?
    dt=$((SECONDS - t0))
    line="$(grep -m1 '^ Exit' "$dir/$log" 2>/dev/null | sed 's/^ *Exit *//')"
    # A completed run always prints its exit line, so its absence is a failure even when
    # the status says otherwise. A stale binary rejecting a lattice it does not understand
    # looked like success here until this check existed.
    verdict=ok
    if [[ $rc -ne 0 || -z "$line" ]]; then
      line="$(grep -m1 -E 'ERROR|FATAL' "$dir/$log" 2>/dev/null | head -c 90)"
      verdict=fail
      status=1
    fi
    printf '%s\n%s\n%s\n%s\n' "$verdict" "$rc" "$dt" "$line" > "$STATE/out.$idx.$k"
    k=$((k + 1))
  done

  # One figure per example, from the primary deck's stats file. The out_root is read from
  # the deck rather than assumed equal to the directory name, because several examples
  # deliberately differ (chicane writes "chicane", coherent_source writes "coherent").
  # It runs in this worker, after this directory's decks, because it reads what they wrote.
  if [[ $FIGURES -eq 1 ]]; then
    primary="$dir/lucifer.in"
    [[ -f "$primary" ]] || primary="${decks[0]}"
    root="$(sed -n 's/.*out_root[[:space:]]*=[[:space:]]*"\([^"]*\)".*/\1/p' "$primary" | head -1)"
    stats="$dir/$root.stats.h5"
    if [[ -z "$root" || ! -f "$stats" ]]; then
      echo "Error: $dir: no stats file at $stats to plot." > "$STATE/note.$idx"
      status=1
    elif ! "$PYTHON" plot_fel.py "$stats" --dpi "$DPI" -o "$FIG_DIR/$dir.png" \
           > "$dir/plot.log" 2>&1 < /dev/null; then
      echo "Error: $dir: plot_fel.py failed; see $dir/plot.log" > "$STATE/note.$idx"
      status=1
    fi
  fi
  return $status
}

# A slot per worker, held in a fifo, because this runs under bash 3.2 where "wait -n" does
# not exist. A worker takes a slot before it launches and returns it as it ends, whatever
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
  while IFS= read -r deck_path; do
    deck="$(basename "$deck_path")"
    if [[ -f "$STATE/out.$idir.$k" ]]; then
      verdict=""; rc=""; dt=""; line=""
      { read -r verdict; read -r rc; read -r dt; IFS= read -r line || line=""; } \
          < "$STATE/out.$idir.$k"
      [[ "$verdict" == "ok" ]] || FAILED=1
      printf '%-20s %-26s %6s %6s  %s\n' "$dir" "$deck" "$rc" "$dt" "$line"
    else
      printf '%-20s %-26s %6s %6s  %s\n' "$dir" "$deck" "-" "-" "no outcome recorded"
      FAILED=1
    fi
    k=$((k + 1))
  done < "$STATE/decks.$idir"
  if [[ -f "$STATE/note.$idir" ]]; then cat "$STATE/note.$idir" >&2; FAILED=1; fi
  idir=$((idir + 1))
done

echo "------------------------------------------------------------------------------"
if [[ $FAILED -ne 0 ]]; then
  echo " FAIL: an example did not run, or its figure was not written."
  exit 1
fi
if [[ $FIGURES -eq 1 ]]; then
  echo " All examples ran. Figures in doc/generated/examples/."
  echo " Regenerate the pages with tests/scripts/report_examples.py after changing a deck."
else
  echo " All examples ran. No figures written (--no-figures)."
fi
echo "=============================================================================="
