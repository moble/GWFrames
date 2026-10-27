cmake_minimum_required(VERSION 3.20...3.31)
project(GWFrames LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# External libraries.  GSL has a module that ships with CMake; FFTW does
# not, so it is found through pkg-config when that is available, and by
# searching the usual paths otherwise.
find_package(GSL 2.0 REQUIRED)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(FFTW3 QUIET IMPORTED_TARGET fftw3)
endif()
if(TARGET PkgConfig::FFTW3)
  add_library(FFTW3::fftw3 ALIAS PkgConfig::FFTW3)
else()
  find_path(FFTW3_INCLUDE_DIR fftw3.h REQUIRED)
  find_library(FFTW3_LIBRARY NAMES fftw3 REQUIRED)
  add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
  set_target_properties(FFTW3::fftw3 PROPERTIES
    IMPORTED_LOCATION "${FFTW3_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}")
endif()

# The vendored copy of spinsfast is plain C.
file(GLOB SPINSFAST_SOURCES CONFIGURE_DEPENDS extern/spinsfast/code/*.c)
add_library(spinsfast STATIC ${SPINSFAST_SOURCES})
target_include_directories(spinsfast PUBLIC extern/spinsfast/include)
target_link_libraries(spinsfast PUBLIC FFTW3::fftw3)

# GWFrames itself, together with the C++ code it uses from the submodules.
add_library(gwframes STATIC
  extern/Quaternions/Quaternions.cpp
  extern/Quaternions/IntegrateAngularVelocity.cpp
  extern/Quaternions/QuaternionUtilities.cpp
  extern/PostNewtonian/C++/PNEvolution.cpp
  extern/PostNewtonian/C++/PNEvolution_Q.cpp
  extern/PostNewtonian/C++/PNWaveformModes.cpp
  extern/SphericalFunctions/Combinatorics.cpp
  extern/SphericalFunctions/WignerDMatrices.cpp
  extern/SphericalFunctions/SWSHs.cpp
  extern/SpacetimeAlgebra/SpacetimeAlgebra.cpp
  cpp/Utilities.cpp
  cpp/Waveforms.cpp
  cpp/PNWaveforms.cpp
  cpp/WaveformsAtAPointFT.cpp
  cpp/fft.cpp
  cpp/NoiseCurves.cpp
  cpp/Interpolate.cpp
  cpp/Scri.cpp
)
target_include_directories(gwframes PUBLIC
  cpp
  extern
  extern/Quaternions
  extern/SphericalFunctions
  extern/SpacetimeAlgebra
)
target_compile_definitions(gwframes PUBLIC USE_GSL)
target_link_libraries(gwframes PUBLIC spinsfast GSL::gsl FFTW3::fftw3)

# Waveform histories record the code revision: the git commit when building
# from a clone, and the package version otherwise (as when building an sdist).
set(GWFRAMES_CODE_REVISION "${SKBUILD_PROJECT_VERSION}")
find_package(Git QUIET)
if(Git_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    OUTPUT_VARIABLE GWFRAMES_GIT_REVISION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if(GWFRAMES_GIT_REVISION)
    set(GWFRAMES_CODE_REVISION "${GWFRAMES_GIT_REVISION}")
  endif()
endif()
if(NOT GWFRAMES_CODE_REVISION)
  set(GWFRAMES_CODE_REVISION "unknown")
endif()
target_compile_definitions(gwframes PRIVATE "CodeRevision=\"${GWFRAMES_CODE_REVISION}\"")

# The Python extension module, installed inside the `GWFrames` package
find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(_GWFrames
  cpp/bindings/module.cpp
  cpp/bindings/quaternions.cpp
  cpp/bindings/utilities.cpp
  cpp/bindings/waveforms.cpp
  cpp/bindings/scri.cpp
)
target_link_libraries(_GWFrames PRIVATE gwframes)
install(TARGETS _GWFrames LIBRARY DESTINATION GWFrames)
