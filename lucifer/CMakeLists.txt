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

# -O3 with the paired-sincos shim takes the FEL step's particle path from 152 to about
# 108 ns a particle-step. Neither pays alone (validation.md, "The particle-path cost").
# The change is at the 1-ulp level. gfortran's sin/cos intrinsics differ from libm's by
# 1 ulp on about 2e-6 of arguments, and -O3 re-contracts floating point.
set_source_files_properties(code/fel_track_mod.f90 PROPERTIES COMPILE_OPTIONS "-O3")

# The FP32 field twin (code/fel_fp32_mod.f90) transforms with FFTW's single-precision
# interface. Not every toolchain carries it. The conda environment does. A from-source
# FFTW in the distribution is double precision only, and there an unconditional -lfftw3f
# fails to link. So the library is searched for here. This runs before Master.cmake, so
# the search path is composed the way that file composes its link directories.
# Without the library nothing links against it and fp32_check is refused.

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
# Objective-C++ against the Metal framework. So it needs a Clang-family toolchain, and
# macOS alone is not enough: this build system may hand OBJCXX to a GNU compiler, which
# rejects -fobjc-arc. The backend is chosen here and used after Master.cmake makes the
# target. A different compiler is *not* forced for the one file, since that would put
# two C++ runtimes in one binary.

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

# Exactly one implementation of the device seam joins the library. A build that cannot
# carry the Metal backend takes the stub. So the same tree builds everywhere, and a deck
# that asks for a device gets an error instead of a quiet CPU fallback. device/ is kept
# out of SRC_DIRS: Master.cmake globs whole directories, and these two files must never
# both compile.

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
