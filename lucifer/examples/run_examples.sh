#!/bin/bash
# Run every example and save the figure the documentation shows.
#
# Two jobs in one pass. It is the answer to "do all the examples still run", printing one
# line per deck and failing on the first nonzero exit. And it produces the figures the
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

# The committed figures are read on a page rather than zoomed into, and they live in the
# repository's history, so they are written below plot_fel.py's own default resolution.
DPI=110

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exe)        EXE="$2";    shift 2 ;;
    --python)     PYTHON="$2"; shift 2 ;;
    --no-figures) FIGURES=0;   shift ;;
    --dpi)        DPI="$2";     shift 2 ;;
    -h|--help)    sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Error: unknown argument: $1" >&2; exit 1 ;;
  esac
done

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

cd "$BMAD_ROOT/lucifer/examples" || exit 1
mkdir -p "$FIG_DIR"

echo "=============================================================================="
echo " Examples"
echo "   lucifer: $EXE"
[[ $FIGURES -eq 1 ]] && echo "   figures: $FIG_DIR"
echo "------------------------------------------------------------------------------"
printf '%-20s %-26s %6s %6s  %s\n' DIR DECK EXIT SEC EXIT_LINE

FAILED=0

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

  # The example's own figure comes from its primary deck, whose out_root names the file.
  for deck_path in "${decks[@]}"; do
    deck="$(basename "$deck_path")"
    log="${deck%.in}.log"
    t0=$SECONDS
    ( cd "$dir" && "$EXE" "$deck" > "$log" 2>&1 )
    rc=$?
    dt=$((SECONDS - t0))
    line="$(grep -m1 '^ Exit' "$dir/$log" 2>/dev/null | sed 's/^ *Exit *//')"
    if [[ $rc -ne 0 ]]; then
      line="$(grep -m1 'ERROR' "$dir/$log" 2>/dev/null | head -c 90)"
      FAILED=1
    fi
    printf '%-20s %-26s %6s %6s  %s\n' "$dir" "$deck" "$rc" "$dt" "$line"
  done

  # One figure per example, from the primary deck's stats file. The out_root is read from
  # the deck rather than assumed equal to the directory name, because several examples
  # deliberately differ (chicane writes "chicane", coherent_source writes "coherent").
  if [[ $FIGURES -eq 1 ]]; then
    primary="$dir/lucifer.in"
    [[ -f "$primary" ]] || primary="${decks[0]}"
    root="$(sed -n 's/.*out_root[[:space:]]*=[[:space:]]*"\([^"]*\)".*/\1/p' "$primary" | head -1)"
    stats="$dir/$root.stats.h5"
    if [[ -z "$root" || ! -f "$stats" ]]; then
      echo "Error: $dir: no stats file at $stats to plot." >&2
      FAILED=1
    else
      if ! "$PYTHON" plot_fel.py "$stats" --dpi "$DPI" -o "$FIG_DIR/$dir.png" \
             > "$dir/plot.log" 2>&1; then
        echo "Error: $dir: plot_fel.py failed; see $dir/plot.log" >&2
        FAILED=1
      fi
    fi
  fi
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
