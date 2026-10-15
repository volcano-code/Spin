cmake_minimum_required(VERSION 3.21)

project(occt_step_roundtrip LANGUAGES CXX)

# This target intentionally consumes the locally pinned OCCT 7.9.3 *build
# tree*.  Its generated package config is install-tree-only, so importing it
# here would accidentally require an unpinned system installation.
# It is a synthetic regression fixture, not a Spin-It Faster paper-case runner.
set(OCCT_BUILD_ROOT "" CACHE PATH "Pinned OCCT 7.9.3 build tree")
if(NOT EXISTS "${OCCT_BUILD_ROOT}/include/opencascade/BRepAlgoAPI_Cut.hxx")
  message(FATAL_ERROR "OCCT_BUILD_ROOT must point to the pinned OCCT 7.9.3 build tree")
endif()
if(NOT EXISTS "${OCCT_BUILD_ROOT}/OpenCASCADETargets.cmake")
  message(FATAL_ERROR "Pinned OCCT build tree has no exported target dependency closure")
endif()
set(OCCT_INCLUDE_DIR "${OCCT_BUILD_ROOT}/include/opencascade")
# This export is generated in the build tree and preserves OCCT's full shared
# library dependency closure.  In contrast, OpenCASCADEConfig.cmake is only
# complete after installation and is intentionally not used here.
include("${OCCT_BUILD_ROOT}/OpenCASCADETargets.cmake")
find_path(GMPXX_INCLUDE_DIR gmpxx.h REQUIRED)
find_library(GMPXX_LIBRARY gmpxx REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)

add_executable(occt_step_roundtrip main.cpp)
target_compile_features(occt_step_roundtrip PRIVATE cxx_std_17)
target_include_directories(occt_step_roundtrip PRIVATE "${OCCT_INCLUDE_DIR}" "${GMPXX_INCLUDE_DIR}")
target_link_libraries(occt_step_roundtrip PRIVATE
  TKDESTEP TKBO TKPrim TKBRep TKTopAlgo "${GMPXX_LIBRARY}" "${GMP_LIBRARY}"
)
