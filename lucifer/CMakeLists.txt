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

set (EXE_SPECS
  cmake.lucifer
  cmake.lucifer_smoke_test
  cmake.lucifer_advance_bench
  cmake.wavefront_drift_test
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
