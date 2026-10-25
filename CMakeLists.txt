# Lucifer, an FEL tracker inside Bmad, built against an installed Bmad.
#
#   cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DLUCIFER_BMAD_ROOT=<bmad>/production
#   cmake --build build-release --parallel
#
# LUCIFER_BMAD_ROOT names the Bmad tree to build against, a distribution build's production/ or
# debug/ root or a prefix of the same shape (cmake/FindBmad.cmake). The binaries land in
# the build directory's bin/, the Fortran modules in its modules/ and the library in its
# lib/. util/conda_compile runs these lines for both configurations.
#
# The flags are the ones Bmad's distribution build gave this program as a user program,
# since the recorded levels of docs/validation.md were measured under them: gfortran at
# -O2 for Release, with the FEL step at -O3, and -O0 with bounds checks for Debug.

cmake_minimum_required (VERSION 3.20)
project (lucifer LANGUAGES Fortran C)

list (APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Release or Debug" FORCE)
endif ()
if (NOT CMAKE_BUILD_TYPE MATCHES "^(Release|Debug)$")
  message (FATAL_ERROR "CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}; the two this program "
           "records levels for are Release and Debug.")
endif ()

# The compiler flags. The common set is what every Fortran file of the distribution's
# user programs is compiled with. Release holds the level the production levels were
# recorded at, and Debug the bounds-checked tree the regression suite and the debug pass
# run on. CMake's own -O3 for Release is replaced.

set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fno-range-check -fdollar-ok -fbacktrace -ffree-line-length-none -fPIC")
set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
set (CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fno-range-check -fbounds-check -Wuninitialized")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -fPIC")
set (CMAKE_C_FLAGS_RELEASE "-O2")
set (CMAKE_C_FLAGS_DEBUG "-g")
set (CMAKE_OBJCXX_FLAGS_RELEASE "")
set (CMAKE_OBJCXX_FLAGS_DEBUG "-g")

set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# The conda environment the build runs in holds toml-f's CMake configuration and HDF5's,
# so its prefix joins the search path when it is active.
if (DEFINED ENV{CONDA_PREFIX})
  list (APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif ()

find_package (Bmad REQUIRED)
find_package (OpenMP REQUIRED COMPONENTS Fortran)
find_package (HDF5 REQUIRED COMPONENTS Fortran HL)
# The deck is TOML, read with toml-f (src/fel_input_mod.f90, src/fel_schema_mod.f90).
find_package (toml-f REQUIRED)

# The FP32 field twin (src/fel_fp32_mod.f90) transforms with FFTW's single-precision
# interface, which not every toolchain carries. Without the library nothing links against
# it and fp32_check is refused.

find_library (LUCIFER_FFTW3F_LIBRARY NAMES fftw3f
  PATHS "${LUCIFER_BMAD_ROOT}/lib" "${LUCIFER_BMAD_ROOT}/../packages/production/lib" "$ENV{CONDA_PREFIX}/lib"
  NO_DEFAULT_PATH)
if (LUCIFER_FFTW3F_LIBRARY)
  message (STATUS "Single-precision FFTW found, so fp32_check is available")
else ()
  message (STATUS "No single-precision FFTW (libfftw3f), so fp32_check is refused")
endif ()

# The Metal backend is Objective-C++ against the Metal framework, so it needs a
# Clang-family toolchain, and macOS alone is not enough. A build that cannot carry it
# takes the stub, so the same tree builds everywhere and a deck that asks for a device
# gets a refusal instead of a quiet CPU fallback. Exactly one of the two joins the
# library.

set (LUCIFER_DEVICE_BACKEND stub)
if (APPLE)
  include (CheckLanguage)
  check_language (OBJCXX)
  if (CMAKE_OBJCXX_COMPILER)
    enable_language (OBJCXX)
    find_library (METAL_LIBRARY Metal)
    find_library (FOUNDATION_LIBRARY Foundation)
    if (CMAKE_OBJCXX_COMPILER_ID MATCHES "Clang" AND METAL_LIBRARY AND FOUNDATION_LIBRARY)
      set (LUCIFER_DEVICE_BACKEND metal)
    endif ()
  endif ()
endif ()

# The library, liblucifer.a: every module under src/, the paired-sincos shim, and the
# device seam. The target is lucifer_lib since the tracker's executable takes the name.

file (GLOB LUCIFER_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.f90")
add_library (lucifer_lib STATIC ${LUCIFER_SOURCES} src/fel_sincos.c)
set_target_properties (lucifer_lib PROPERTIES OUTPUT_NAME lucifer)

# -O3 with the paired-sincos shim takes the FEL step's particle path from 152 to about
# 108 ns a particle-step, and neither pays alone (docs/validation.md). The change is at
# the 1-ulp level, since gfortran's sin and cos differ from libm's by 1 ulp on about 2e-6
# of arguments and -O3 re-contracts floating point. It applies in both configurations.
set_source_files_properties (src/fel_track_mod.f90 PROPERTIES COMPILE_OPTIONS "-O3")

if (LUCIFER_FFTW3F_LIBRARY)
  set_source_files_properties (src/fel_fp32_mod.f90 PROPERTIES COMPILE_DEFINITIONS LUCIFER_HAVE_FFTW3F)
  target_link_libraries (lucifer_lib PUBLIC "${LUCIFER_FFTW3F_LIBRARY}")
endif ()

if (LUCIFER_DEVICE_BACKEND STREQUAL metal)
  set (CMAKE_OBJCXX_STANDARD 17)
  set (CMAKE_OBJCXX_STANDARD_REQUIRED ON)
  target_sources (lucifer_lib PRIVATE src/device/lucifer_metal.mm)
  set_source_files_properties (src/device/lucifer_metal.mm PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
  target_link_libraries (lucifer_lib PUBLIC ${METAL_LIBRARY} ${FOUNDATION_LIBRARY})
  message (STATUS "Device backend: Metal")
else ()
  target_sources (lucifer_lib PRIVATE src/device/lucifer_device_stub.c)
  message (STATUS "Device backend: none, so run.device is refused")
endif ()

target_include_directories (lucifer_lib PUBLIC "${Bmad_MODULE_DIR}" "${CMAKE_SOURCE_DIR}/src/device")
if (Bmad_INCLUDE_DIR)
  target_include_directories (lucifer_lib PUBLIC "${Bmad_INCLUDE_DIR}")
endif ()
# HDF5 is named through the two lists FindHDF5 fills rather than its hdf5:: targets.
# With an HDF5 that ships a CMake configuration, conda-forge's among them, the lists hold
# that configuration's targets, which carry the Fortran stub libraries hdf5_hl_f90cstub
# and hdf5_f90cstub that bmad's hdf5_interface calls into. The hdf5:: targets name the
# two front libraries alone and the link of a static libbmad.a then fails.
target_link_libraries (lucifer_lib PUBLIC ${Bmad_LIBRARIES}
                       ${HDF5_Fortran_HL_LIBRARIES} ${HDF5_Fortran_LIBRARIES}
                       toml-f::toml-f OpenMP::OpenMP_Fortran ${CMAKE_DL_LIBS})
if (APPLE)
  # The distribution linked its binaries with this, and the recorded levels carry it.
  target_link_options (lucifer_lib PUBLIC "-Wl,-no_compact_unwind")
endif ()

# The four executables: the tracker, the smoke test and its namelist twin, the advance
# benchmark and the wavefront propagator's test driver.

foreach (exe lucifer lucifer_smoke_test lucifer_advance_bench wavefront_drift_test)
  add_executable (${exe} "programs/${exe}.f90")
  target_link_libraries (${exe} PRIVATE lucifer_lib)
endforeach ()
