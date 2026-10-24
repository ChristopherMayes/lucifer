# Finds an installed Bmad for Lucifer to build against.
#
# LUCIFER_BMAD_ROOT names the tree, a distribution build's production/ or debug/ root, which holds
# lib/, modules/ and include/, or a prefix of the same shape with the modules under
# include/. The variable comes from the cache, else from the environment, where the build
# script and the keystone set it from BMAD_DIR and the flavor. Nothing is searched for
# outside that tree and the libraries it needs beside it: a Bmad found by accident on a
# path would build a binary whose digits belong to no declared release.
#
# Sets:
#   Bmad_FOUND
#   Bmad_MODULE_DIR      The directory holding bmad.mod and sim_utils.mod.
#   Bmad_INCLUDE_DIR     Bmad's include/ when the tree has one.
#   Bmad_LIBRARIES       The archives and libraries to link, in link order.
#
# The link list is the one Bmad's own build hands a user program through its
# ACC_BMAD_LINK_LIBS variable, read from the production build of 2026-10-23, with the three
# Bmad archives first. Each entry carries the reason it is there, so a library that Bmad
# stops needing can be dropped with its line.

include(FindPackageHandleStandardArgs)

if (NOT LUCIFER_BMAD_ROOT AND DEFINED ENV{LUCIFER_BMAD_ROOT})
  set (LUCIFER_BMAD_ROOT "$ENV{LUCIFER_BMAD_ROOT}" CACHE PATH "The Bmad tree to build against")
endif ()
if (NOT LUCIFER_BMAD_ROOT)
  message (FATAL_ERROR "LUCIFER_BMAD_ROOT is not set. Name the Bmad tree to build against, "
           "a distribution build's production or debug root or a prefix of that shape.")
endif ()
get_filename_component (LUCIFER_BMAD_ROOT "${LUCIFER_BMAD_ROOT}" ABSOLUTE)

find_path (Bmad_MODULE_DIR bmad.mod
  PATHS "${LUCIFER_BMAD_ROOT}/modules" "${LUCIFER_BMAD_ROOT}/include" NO_DEFAULT_PATH)
find_path (Bmad_INCLUDE_DIR NAMES bmad_std_typedefs.h sim_utils.h
  PATHS "${LUCIFER_BMAD_ROOT}/include" NO_DEFAULT_PATH)

# The three archives Lucifer calls into: the Bmad library, its utilities and PTC's
# forest, which bmad's element tracking reaches. The archive is taken ahead of the shared
# library of the same name, as the distribution's link line did, so the binary carries
# Bmad and runs without the tree it was built against.
set (_bmad_own bmad sim_utils forest)
set (Bmad_LIBRARIES)
foreach (_lib ${_bmad_own})
  find_library (Bmad_${_lib}_LIBRARY NAMES lib${_lib}.a ${_lib} PATHS "${LUCIFER_BMAD_ROOT}/lib" NO_DEFAULT_PATH)
  list (APPEND Bmad_LIBRARIES "${Bmad_${_lib}_LIBRARY}")
endforeach ()

# What the archives call, found where the toolchain's libraries are: beside the tree for a
# distribution build that built its own packages, and in the prefix the compilers come
# from otherwise, and nowhere else, since the system's own lapack and blas are not the
# ones the Bmad tree was built against. Each with its reason.
#   xrlf03, xrl        xraylib, the photon cross sections bmad's x-ray elements use.
#   fgsl, gsl, gslcblas  the GNU Scientific Library and its Fortran binding, sim_utils'
#                      special functions and integrators.
#   lapack95, lapack, blas  the linear algebra of the mat6 and Twiss code.
#   fftw3, fftw3_omp   the transforms sim_utils and Lucifer's field solver take.
set (_bmad_deps xrlf03 xrl fgsl gsl gslcblas lapack95 lapack blas fftw3 fftw3_omp)
set (_dep_paths "${LUCIFER_BMAD_ROOT}/lib" "${LUCIFER_BMAD_ROOT}/../packages/production/lib"
     "${LUCIFER_BMAD_ROOT}/../packages/debug/lib")
if (DEFINED ENV{CONDA_PREFIX})
  list (APPEND _dep_paths "$ENV{CONDA_PREFIX}/lib")
endif ()
foreach (_lib ${_bmad_deps})
  find_library (Bmad_${_lib}_LIBRARY NAMES ${_lib} PATHS ${_dep_paths} NO_DEFAULT_PATH)
  if (Bmad_${_lib}_LIBRARY)
    list (APPEND Bmad_LIBRARIES "${Bmad_${_lib}_LIBRARY}")
  else ()
    message (FATAL_ERROR "Bmad needs ${_lib} and it was not found beside ${LUCIFER_BMAD_ROOT} "
             "or in the environment's lib.")
  endif ()
endforeach ()

find_package_handle_standard_args (Bmad REQUIRED_VARS Bmad_MODULE_DIR Bmad_bmad_LIBRARY
  Bmad_sim_utils_LIBRARY Bmad_forest_LIBRARY)
