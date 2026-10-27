# Sets Lucifer's version and writes the module that carries it.
#
# CMakeLists.txt includes this file at configure time, and the target lucifer_version
# runs it as a script at every build, so a commit or a tag made after the configure shows
# in the next build without a reconfigure. The version is the first of:
#
#   1. LUCIFER_VERSION, when it is given. A packaged build passes its own version here.
#   2. git describe --tags --dirty --match "v[0-9]*", in a checkout whose top level is
#      this source tree: v0.1.0 on the tag, v0.1.0-3-g<hash> three commits on, and
#      0.0.0-g<hash> in a checkout with no such tag. A tree with uncommitted changes ends
#      in -dirty. A source tree unpacked inside some other repository is not described.
#   3. The describe-name line of .git_archival.txt. .gitattributes marks the file
#      export-subst, so git archive, and a release archive on GitHub with it, fills the
#      line from the tagged commit. An unexpanded $Format placeholder counts as absent.
#   4. 0.0.0-unknown, with a note.
#
# The module is rewritten only when the version changes, so a build at an unchanged
# version recompiles nothing.
#
# Takes:
#   LUCIFER_SOURCE_DIR     The Lucifer source root.
#   LUCIFER_VERSION_FILE   The module to write, from src/fel_version_mod.f90.in.
#   LUCIFER_VERSION        The version to use, or empty.
#   GIT_EXECUTABLE         git, or empty. Found here when unset.
#
# Sets:
#   LUCIFER_VERSION_STRING   The version.
#   LUCIFER_VERSION_ORIGIN   Which of the four gave it.

if (NOT GIT_EXECUTABLE)
  find_program (GIT_EXECUTABLE git)
endif ()

set (LUCIFER_VERSION_STRING "")

if (LUCIFER_VERSION)
  set (LUCIFER_VERSION_STRING "${LUCIFER_VERSION}")
  set (LUCIFER_VERSION_ORIGIN "LUCIFER_VERSION")
endif ()

if (NOT LUCIFER_VERSION_STRING AND GIT_EXECUTABLE)
  execute_process (COMMAND "${GIT_EXECUTABLE}" rev-parse --show-toplevel
                   WORKING_DIRECTORY "${LUCIFER_SOURCE_DIR}"
                   OUTPUT_VARIABLE _top RESULT_VARIABLE _err
                   OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if (_err EQUAL 0)
    get_filename_component (_top "${_top}" REALPATH)
    get_filename_component (_src "${LUCIFER_SOURCE_DIR}" REALPATH)
  endif ()
  if (_err EQUAL 0 AND _top STREQUAL _src)
    execute_process (COMMAND "${GIT_EXECUTABLE}" describe --tags --dirty --match "v[0-9]*"
                     WORKING_DIRECTORY "${LUCIFER_SOURCE_DIR}"
                     OUTPUT_VARIABLE _describe RESULT_VARIABLE _err
                     OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if (NOT _err EQUAL 0)
      execute_process (COMMAND "${GIT_EXECUTABLE}" describe --always --dirty
                       WORKING_DIRECTORY "${LUCIFER_SOURCE_DIR}"
                       OUTPUT_VARIABLE _describe RESULT_VARIABLE _err
                       OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      set (_describe "0.0.0-g${_describe}")
    endif ()
    if (_err EQUAL 0)
      set (LUCIFER_VERSION_STRING "${_describe}")
      set (LUCIFER_VERSION_ORIGIN "git describe")
    endif ()
  endif ()
endif ()

if (NOT LUCIFER_VERSION_STRING AND EXISTS "${LUCIFER_SOURCE_DIR}/.git_archival.txt")
  file (STRINGS "${LUCIFER_SOURCE_DIR}/.git_archival.txt" _line REGEX "^describe-name:")
  string (REGEX REPLACE "^describe-name:[ ]*" "" _line "${_line}")
  string (STRIP "${_line}" _line)
  if (_line AND NOT _line MATCHES "\\$Format")
    set (LUCIFER_VERSION_STRING "${_line}")
    set (LUCIFER_VERSION_ORIGIN ".git_archival.txt")
  endif ()
endif ()

if (NOT LUCIFER_VERSION_STRING)
  set (LUCIFER_VERSION_STRING "0.0.0-unknown")
  set (LUCIFER_VERSION_ORIGIN "fallback")
  message (STATUS "No version from LUCIFER_VERSION, git or .git_archival.txt, so the "
           "version is ${LUCIFER_VERSION_STRING}")
endif ()

set (_old "")
if (EXISTS "${LUCIFER_VERSION_FILE}")
  file (READ "${LUCIFER_VERSION_FILE}" _old)
endif ()
configure_file ("${LUCIFER_SOURCE_DIR}/src/fel_version_mod.f90.in" "${LUCIFER_VERSION_FILE}" @ONLY)
file (READ "${LUCIFER_VERSION_FILE}" _new)
if (NOT CMAKE_SCRIPT_MODE_FILE OR NOT _new STREQUAL _old)
  message (STATUS "Lucifer version ${LUCIFER_VERSION_STRING}, from ${LUCIFER_VERSION_ORIGIN}")
endif ()
