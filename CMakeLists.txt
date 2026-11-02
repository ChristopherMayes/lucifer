# Lucifer, an FEL tracker inside Bmad, built against an installed Bmad.
#
#   cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DLUCIFER_BMAD_ROOT=<bmad>/production
#   cmake --build build-release --parallel
#
# LUCIFER_BMAD_ROOT names the Bmad tree to build against, a distribution build's production/ or
# debug/ root or a prefix of the same shape (cmake/FindBmad.cmake). The binaries land in
# the build directory's bin/, the Fortran modules in its modules/ and the library in its
# lib/. The library is shared, liblucifer.so or .dylib, linked against Bmad's shared
# libraries, and -DBUILD_SHARED_LIBS=OFF builds the archive liblucifer.a against Bmad's
# archives instead. util/conda_compile runs these lines for both configurations.
# cmake --install with --prefix <dir> then places the tracker in <dir>/bin, the library in
# <dir>/lib, its modules in <dir>/include/lucifer and LuciferConfig.cmake in
# <dir>/lib/cmake/lucifer, so another program builds against it with find_package (Lucifer).
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

# The library's kind. Shared is the default, and what the conda-forge package carries.
# FindBmad reads the same switch to take Bmad's shared libraries or its archives.
option (BUILD_SHARED_LIBS "Build liblucifer as a shared library, else as an archive" ON)

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
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# The conda environment the build runs in holds toml-f's CMake configuration and HDF5's,
# so its prefix joins the search path when it is active.
if (DEFINED ENV{CONDA_PREFIX})
  list (APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif ()

include (GNUInstallDirs)
find_package (Bmad REQUIRED)
find_package (OpenMP REQUIRED COMPONENTS Fortran)
find_package (HDF5 REQUIRED COMPONENTS Fortran HL)
# The input file is TOML, read with toml-f (src/fel_input_mod.f90, src/fel_schema_mod.f90).
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
# takes the stub, so the same tree builds everywhere and an input file that asks for a device
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

# The version (cmake/LuciferVersion.cmake): LUCIFER_VERSION when given, else git describe,
# else the release archive's .git_archival.txt, else 0.0.0-unknown. It is written into
# fel_version_mod.f90 in the build directory, here and again at every build by the target
# lucifer_version, so each build tree carries the version it was built at.

set (LUCIFER_VERSION "" CACHE STRING "The version to build as, in place of git's")
find_package (Git QUIET)
set (LUCIFER_SOURCE_DIR "${CMAKE_SOURCE_DIR}")
set (LUCIFER_VERSION_FILE "${CMAKE_BINARY_DIR}/generated/fel_version_mod.f90")
include (LuciferVersion)
add_custom_target (lucifer_version
  COMMAND "${CMAKE_COMMAND}" "-DLUCIFER_SOURCE_DIR=${CMAKE_SOURCE_DIR}"
          "-DLUCIFER_VERSION_FILE=${LUCIFER_VERSION_FILE}" "-DLUCIFER_VERSION=${LUCIFER_VERSION}"
          "-DGIT_EXECUTABLE=${GIT_EXECUTABLE}" -P "${CMAKE_SOURCE_DIR}/cmake/LuciferVersion.cmake"
  BYPRODUCTS "${LUCIFER_VERSION_FILE}")

# The library, liblucifer: every module under src/, the version module, the
# paired-sincos shim, and the device seam. The target is lucifer_lib since the tracker's
# executable takes the name. A shared library carries the numeric part of the version,
# 0.1.0 from v0.1.0-3-g<hash>, as its version and the major part as its soversion, and
# 0.0.0 when the version has no numeric part.

file (GLOB LUCIFER_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.f90")
add_library (lucifer_lib ${LUCIFER_SOURCES} "${LUCIFER_VERSION_FILE}" src/fel_sincos.c)
if (LUCIFER_VERSION_STRING MATCHES "^v?([0-9]+)\\.([0-9]+)\\.([0-9]+)")
  set (LUCIFER_VERSION_NUMERIC "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
  set (LUCIFER_VERSION_MAJOR "${CMAKE_MATCH_1}")
else ()
  set (LUCIFER_VERSION_NUMERIC "0.0.0")
  set (LUCIFER_VERSION_MAJOR "0")
endif ()
set_target_properties (lucifer_lib PROPERTIES OUTPUT_NAME lucifer
                       VERSION "${LUCIFER_VERSION_NUMERIC}" SOVERSION "${LUCIFER_VERSION_MAJOR}")
add_dependencies (lucifer_lib lucifer_version)

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

# The device header is the library's own. Bmad's module directory is public, since a
# program that uses Lucifer's modules uses Bmad's, and it is carried by its path into the
# installed configuration. The installed modules are under include/lucifer.
target_include_directories (lucifer_lib PRIVATE "${CMAKE_SOURCE_DIR}/src/device")
target_include_directories (lucifer_lib PUBLIC "${Bmad_MODULE_DIR}"
                            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/lucifer>")
if (Bmad_INCLUDE_DIR)
  target_include_directories (lucifer_lib PUBLIC "${Bmad_INCLUDE_DIR}")
endif ()
# HDF5 is named through the two lists FindHDF5 fills rather than its hdf5:: targets.
# With an HDF5 that ships a CMake configuration, conda-forge's among them, the lists hold
# that configuration's targets, which carry the Fortran stub libraries hdf5_hl_f90cstub
# and hdf5_f90cstub that bmad's hdf5_interface calls into. The hdf5:: targets name the
# two front libraries alone and the link of a static libbmad.a then fails.
#
# Bmad's libraries and OpenMP are public, since a program that uses Lucifer's modules
# uses Bmad's and the modules carry threadprivate variables. HDF5 and toml-f are the
# shared library's own, and public only for the archive, whose user links them itself.
if (BUILD_SHARED_LIBS)
  set (LUCIFER_OWN_SCOPE PRIVATE)
else ()
  set (LUCIFER_OWN_SCOPE PUBLIC)
endif ()
target_link_libraries (lucifer_lib PUBLIC ${Bmad_LIBRARIES} OpenMP::OpenMP_Fortran ${CMAKE_DL_LIBS})
# Bmad's names carry a dollar sign, lucifer_version$ among Lucifer's own, so a program
# that uses the modules compiles with the extension, as every Bmad program does.
target_compile_options (lucifer_lib PUBLIC "$<$<COMPILE_LANGUAGE:Fortran>:-fdollar-ok>")
target_link_libraries (lucifer_lib ${LUCIFER_OWN_SCOPE}
                       ${HDF5_Fortran_HL_LIBRARIES} ${HDF5_Fortran_LIBRARIES} toml-f::toml-f)
if (APPLE)
  # The distribution linked its binaries with this, and the recorded levels carry it.
  target_link_options (lucifer_lib PUBLIC "-Wl,-no_compact_unwind")
endif ()

# The four executables: the tracker, the smoke test and its namelist twin, the
# advance benchmark and the wavefront propagator's test driver.

foreach (exe lucifer lucifer_smoke_test lucifer_advance_bench wavefront_drift_test)
  add_executable (${exe} "programs/${exe}.f90")
  target_link_libraries (${exe} PRIVATE lucifer_lib)
endforeach ()

# What cmake --install places: the tracker in the prefix's bin/, the library in its
# lib/, every module in its include/lucifer/, and LuciferConfig.cmake with its version
# file and the exported target in lib/cmake/lucifer. The three test programs stay in the
# build directory. The installed tracker finds the library beside its own bin/ through
# its rpath, and both keep the directories of the libraries they were linked against,
# the conda environment's lib/ among them, so they run from the prefix as they run from
# the build directory. The configuration finds OpenMP again for the program that uses
# it, HDF5 and toml-f too for the archive, and carries Bmad's libraries and module
# directory by path.

if (APPLE)
  set (LUCIFER_ORIGIN "@loader_path")
else ()
  set (LUCIFER_ORIGIN "$ORIGIN")
endif ()
set_target_properties (lucifer PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE
                       INSTALL_RPATH "${LUCIFER_ORIGIN}/../${CMAKE_INSTALL_LIBDIR}")
set_target_properties (lucifer_lib PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install (TARGETS lucifer RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install (TARGETS lucifer_lib EXPORT LuciferTargets
         LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
         ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install (DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}/"
         DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/lucifer" FILES_MATCHING PATTERN "*.mod")
set (LUCIFER_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/lucifer")
install (EXPORT LuciferTargets NAMESPACE Lucifer:: DESTINATION "${LUCIFER_CMAKE_DIR}")
include (CMakePackageConfigHelpers)
configure_package_config_file (cmake/LuciferConfig.cmake.in
                               "${CMAKE_BINARY_DIR}/generated/LuciferConfig.cmake"
                               INSTALL_DESTINATION "${LUCIFER_CMAKE_DIR}")
write_basic_package_version_file ("${CMAKE_BINARY_DIR}/generated/LuciferConfigVersion.cmake"
                                  VERSION "${LUCIFER_VERSION_NUMERIC}" COMPATIBILITY SameMajorVersion)
install (FILES "${CMAKE_BINARY_DIR}/generated/LuciferConfig.cmake"
               "${CMAKE_BINARY_DIR}/generated/LuciferConfigVersion.cmake"
         DESTINATION "${LUCIFER_CMAKE_DIR}")
