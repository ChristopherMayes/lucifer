# The locations the Lucifer scripts share, for sourcing from bash. Each has a default
# that is the layout inside a Bmad checkout, so a Lucifer source tree there runs with
# nothing set.
#
#   LUCIFER_DIR     The Lucifer source: this file's own tree.
#   DEBUG_BIN       Lucifer's debug binaries. Master.cmake places a project's binaries
#   PRODUCTION_BIN  in debug/ and production/ beside the project's source, in the Bmad
#                   checkout and out of it alike, so both default to the tree beside
#                   this one, and LUCIFER_DEBUG_BIN and LUCIFER_PRODUCTION_BIN name them
#                   elsewhere.
#   BMAD_SRC        The Bmad checkout, which carries the regression suite and Tao:
#                   BMAD_DIR, else DIST_BASE_DIR as util/dist_source_me exports it, else
#                   the parent of the Lucifer source. A variable naming a directory with
#                   no util/dist_source_me is passed over with a note, since
#                   dist_source_me uses BMAD_DIR as a loop variable and leaves it at
#                   the checkout's util/ once sourced. Empty when no candidate is a
#                   checkout, and the script that needs it refuses then.
#
# A relative value is taken against the directory the script was started in, before
# any script changes directory, and every value is absolute from here on.

lucifer_abs () {   # <path>
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *)  printf '%s\n' "$PWD/$1" ;;
  esac
}

LUCIFER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
lucifer_parent="$(cd "$LUCIFER_DIR/.." && pwd)"
DEBUG_BIN="$(lucifer_abs "${LUCIFER_DEBUG_BIN:-$lucifer_parent/debug/bin}")"
PRODUCTION_BIN="$(lucifer_abs "${LUCIFER_PRODUCTION_BIN:-$lucifer_parent/production/bin}")"

BMAD_SRC=""
for lucifer_var in BMAD_DIR DIST_BASE_DIR; do
  lucifer_candidate="${!lucifer_var:-}"
  [[ -n "$lucifer_candidate" ]] || continue
  lucifer_candidate="$(lucifer_abs "$lucifer_candidate")"
  if [[ -f "$lucifer_candidate/util/dist_source_me" ]]; then
    BMAD_SRC="$lucifer_candidate"
    break
  fi
  echo "note: $lucifer_var=$lucifer_candidate holds no util/dist_source_me and is passed over" >&2
done
if [[ -z "$BMAD_SRC" && -f "$lucifer_parent/util/dist_source_me" ]]; then
  BMAD_SRC="$lucifer_parent"
fi
unset lucifer_var lucifer_candidate lucifer_parent
