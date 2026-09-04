cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set (LIBNAME lucifer)

IF ($ENV{ACC_CONDA_BUILD})
  SET (INC_DIRS
    $ENV{ACC_CONDA_PATH}/include
  )
ELSE ()
  SET (INC_DIRS
  )
ENDIF ()

set (SRC_DIRS
  code
)

# The FEL step's particle path is transcendental-bound (perf goal 2026-08-21):
# -O3 + the paired-sincos shim together take the per-particle-step cost from 152
# to ~108 ns (neither pays alone). This is a NAMED VALUE CHANGE at the 1-ulp level
# (gfortran's sin/cos intrinsics differ from libm's by 1 ulp on ~2e-6 of arguments,
# and -O3 re-contracts FP): every benchmark tier re-measured and re-recorded.
set_source_files_properties(code/fel_track_mod.f90 PROPERTIES COMPILE_OPTIONS "-O3")

# The FP32 field twin (code/fel_fp32_mod.f90) transforms with FFTW's single-precision
# interface. Not every toolchain carries it: the conda environment does, and the
# off-site distribution's from-source FFTW is double precision only, where an
# unconditional -lfftw3f fails to link. Detected here rather than assumed, and this
# runs before Master.cmake, so the search path is composed the way that file composes
# its own link directories. Without the library nothing links against it and
# fp32_check refuses it; every other path is untouched.

find_library (LUCIFER_FFTW3F_LIBRARY NAMES fftw3f
  PATHS
    $ENV{ACC_CONDA_PATH}/lib
    $ENV{DIST_BASE_DIR}/production/lib
    $ENV{DIST_BASE_DIR}/debug/lib
    $ENV{DIST_BASE_DIR}/packages/production/lib
    $ENV{DIST_BASE_DIR}/packages/debug/lib
    ${CMAKE_SOURCE_DIR}/../production/lib
    ${CMAKE_SOURCE_DIR}/../debug/lib
    $ENV{BMAD_USER_LIB_DIRS}
)

if (LUCIFER_FFTW3F_LIBRARY)
  set (LUCIFER_FFTW3F_LIBS fftw3f)
  set_source_files_properties(code/fel_fp32_mod.f90 PROPERTIES
                              COMPILE_DEFINITIONS LUCIFER_HAVE_FFTW3F)
  message (STATUS "Single-precision FFTW found, so fp32_check is available")
else ()
  set (LUCIFER_FFTW3F_LIBS)
  message (STATUS "No single-precision FFTW (libfftw3f), so fp32_check is refused")
endif ()

# The device seam (device/lucifer_device.h). The Metal backend is ARC-managed
# Objective-C++ against the Metal framework, so it needs a Clang-family toolchain and
# not merely macOS: this build system may hand OBJCXX to a GNU compiler, which
# rejects -fobjc-arc outright. Decided here, used after Master.cmake has made the
# target. Forcing a different compiler for the one file is deliberately not attempted:
# two C++ runtimes in one binary is the failure the reference backend's own manual
# records against a hijacked HDF5 wrapper.

set (LUCIFER_DEVICE_BACKEND stub)
if (APPLE)
  enable_language(OBJCXX)
  find_library (METAL_LIBRARY Metal)
  find_library (FOUNDATION_LIBRARY Foundation)
  if (CMAKE_OBJCXX_COMPILER_ID MATCHES "Clang" AND METAL_LIBRARY AND FOUNDATION_LIBRARY)
    set (LUCIFER_DEVICE_BACKEND metal)
  endif ()
endif ()

set (EXE_SPECS
  cmake.lucifer
  cmake.lucifer_smoke_test
  cmake.lucifer_advance_bench
  cmake.wavefront_drift_test
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)

# Exactly one implementation of the device seam joins the library, per the decision
# made above. A build that cannot carry the Metal backend takes the stub, so the same
# tree builds everywhere and a deck asking for a device then gets an error rather than
# a quiet CPU fallback. device/ is deliberately outside SRC_DIRS: Master.cmake globs
# whole directories, and these two files must never both compile.

if (LUCIFER_DEVICE_BACKEND STREQUAL metal)
  set (CMAKE_OBJCXX_STANDARD 17)
  set (CMAKE_OBJCXX_STANDARD_REQUIRED ON)
  target_sources (lucifer PRIVATE device/lucifer_metal.mm)
  set_source_files_properties (device/lucifer_metal.mm PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
  target_link_libraries (lucifer ${METAL_LIBRARY} ${FOUNDATION_LIBRARY})
  if (TARGET lucifer-shared)
    target_sources (lucifer-shared PRIVATE device/lucifer_metal.mm)
    target_link_libraries (lucifer-shared ${METAL_LIBRARY} ${FOUNDATION_LIBRARY})
  endif ()
  message (STATUS "Device backend: Metal")
else ()
  target_sources (lucifer PRIVATE device/lucifer_device_stub.c)
  if (TARGET lucifer-shared)
    target_sources (lucifer-shared PRIVATE device/lucifer_device_stub.c)
  endif ()
  message (STATUS "Device backend: none, so global%device is refused")
endif ()
