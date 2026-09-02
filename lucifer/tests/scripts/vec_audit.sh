#!/bin/bash
#
# What the vectorizer did to the hot files, and what stopped it where it stopped.
#
# The compile command comes from the production build's own build.make rather than from
# a copy of the flags kept here. A copy drifts, and a vectorization report taken under
# flags the build does not use describes a binary nobody runs. That also means the
# production tree must be built first: this script reads its makefile.
#
# Nothing is written into the build tree. The report, the object and the module files
# all go to the output directory, so a build in progress is unaffected and the audit
# leaves the tree as it found it.
#
# The report flag is -fopt-info-vec-all rather than -fopt-info-vec-missed, because the
# question is which loops vectorize as well as what stopped the ones that did not. A
# missed-only report cannot say whether a hot loop is in the first list.
#
# Usage:
#   ./vec_audit.sh [--out <dir>] [--file <name.f90>]...
#
# Defaults: the three files that hold the per-particle loops, into a temporary
# directory whose path is printed. doc/performance.md records the standing summary.

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LUCIFER="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_MAKE="$LUCIFER/production/CMakeFiles/lucifer.dir/build.make"
FLAGS_MAKE="$LUCIFER/production/CMakeFiles/lucifer.dir/flags.make"

OUT=""
FILES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out)  OUT="$2";        shift 2 ;;
    --file) FILES+=("$2");   shift 2 ;;
    -h|--help) sed -n '2,22p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

if [[ ${#FILES[@]} -eq 0 ]]; then
  FILES=(fel_track_mod.f90 fel_unaveraged_mod.f90 fel_collective_mod.f90)
fi

for f in "$BUILD_MAKE" "$FLAGS_MAKE"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: $f is missing. Build production first: ./util/conda_compile" >&2
    exit 1
  fi
done

[[ -z "$OUT" ]] && OUT="$(mktemp -d "${TMPDIR:-/tmp}/fel_vec.XXXXXX")"
mkdir -p "$OUT" || exit 1

# The flags file defines Fortran_FLAGS and Fortran_INCLUDES as make variables, and the
# per-file recipe in build.make interpolates them. Read both, then expand by hand: the
# alternative is running make, which would rebuild into the tree.

read_var () {   # <file> <name>
  sed -n "s/^$2 = //p" "$1" | head -1
}

F_FLAGS="$(read_var "$FLAGS_MAKE" Fortran_FLAGS)"
F_INCLUDES="$(read_var "$FLAGS_MAKE" Fortran_INCLUDES)"
F_DEFINES="$(read_var "$FLAGS_MAKE" Fortran_DEFINES)"

# -J in the build's flags points at the production module directory. The audit must not
# write there, so the module output is redirected and the original -J becomes an -I: the
# module files it holds are exactly what these files compile against.

MOD_DIR="$(echo "$F_FLAGS" | sed -n 's/.*-J\([^ ]*\).*/\1/p')"
F_FLAGS="$(echo "$F_FLAGS" | sed 's|-J[^ ]*||')"
[[ -n "$MOD_DIR" ]] && F_INCLUDES="$F_INCLUDES -I$MOD_DIR"

echo "=============================================================================="
echo "Vectorization audit of the FEL hot files"
echo "=============================================================================="
echo "  flags from:  $FLAGS_MAKE"
echo "  output:      $OUT"
echo

status=0
for src in "${FILES[@]}"; do

  # The per-file recipe carries this file's own options, which is how the production -O3
  # on fel_track_mod reaches the audit. Take the compiler and everything after it from
  # the recipe line, drop its make variables (already expanded above) and its -o.

  recipe="$(grep -m1 -A 4 "^CMakeFiles/lucifer.dir/code/$src.o: CMakeFiles" "$BUILD_MAKE" | \
            grep -m1 'gfortran')"
  if [[ -z "$recipe" ]]; then
    echo "  $src: no compile recipe in build.make. Not a source file of this library?" >&2
    status=1
    continue
  fi
  FC="$(echo "$recipe" | awk '{print $1}')"
  # Everything between the make variables and the -c, which is this file's own options.
  own="$(echo "$recipe" | sed 's|.*\$(Fortran_FLAGS)||; s| -c .*||')"

  rpt="$OUT/${src%.f90}.vec.txt"
  echo "--- $src ---"
  # shellcheck disable=SC2086
  "$FC" $F_DEFINES $F_INCLUDES $F_FLAGS $own \
        -J"$OUT" -fopt-info-vec-all="$rpt" \
        -c "$LUCIFER/code/$src" -o "$OUT/${src%.f90}.o" 2> "$OUT/${src%.f90}.err"
  if [[ $? -ne 0 ]]; then
    echo "  compile FAILED; stderr tail:" >&2
    tail -5 "$OUT/${src%.f90}.err" >&2
    status=1
    continue
  fi

  n_ok="$(grep -c 'loop vectorized' "$rpt" 2>/dev/null)"
  n_no="$(grep -c 'not vectorized' "$rpt" 2>/dev/null)"
  echo "  vectorized: $n_ok    not vectorized: $n_no    report: $rpt"

  # The reasons, most common first. A count is what makes the report a shopping list:
  # one reason behind twenty loops is one change.

  echo "  reasons the vectorizer gave, by count:"
  grep 'not vectorized' "$rpt" 2>/dev/null | sed 's/.*not vectorized: *//' | \
        sort | uniq -c | sort -rn | head -8 | sed 's/^/    /'
  echo
done

echo "Reports kept in: $OUT"
exit $status
