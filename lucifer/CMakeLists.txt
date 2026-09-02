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

# The device seam (device/lucifer_device.h): exactly one implementation file joins
# the library, chosen here and nowhere else. Darwin compiles the Metal backend
# (Objective-C++, so its own language and the two frameworks) and every other
# platform compiles a stub that refuses by name, so the same tree builds everywhere
# and a deck asking for a device on a stub build gets an error, not a fallback.
# device/ is deliberately outside SRC_DIRS: Master.cmake globs whole directories,
# and these two files must never both compile.

if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  enable_language(OBJCXX)
  set (CMAKE_OBJCXX_STANDARD 17)
  set (CMAKE_OBJCXX_STANDARD_REQUIRED ON)
  find_library (METAL_LIBRARY Metal REQUIRED)
  find_library (FOUNDATION_LIBRARY Foundation REQUIRED)
  target_sources (lucifer PRIVATE device/lucifer_metal.mm)
  set_source_files_properties (device/lucifer_metal.mm PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
  target_link_libraries (lucifer ${METAL_LIBRARY} ${FOUNDATION_LIBRARY})
  if (TARGET lucifer-shared)
    target_sources (lucifer-shared PRIVATE device/lucifer_metal.mm)
    target_link_libraries (lucifer-shared ${METAL_LIBRARY} ${FOUNDATION_LIBRARY})
  endif ()
else ()
  target_sources (lucifer PRIVATE device/lucifer_device_stub.c)
  if (TARGET lucifer-shared)
    target_sources (lucifer-shared PRIVATE device/lucifer_device_stub.c)
  endif ()
endif ()
