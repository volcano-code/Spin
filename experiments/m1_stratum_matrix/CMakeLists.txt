cmake_minimum_required(VERSION 3.16)
project(M1StratumMatrix LANGUAGES CXX)

include(CTest)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

get_filename_component(MATRIX_ROOT_DEFAULT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(MATRIX_ROOT "${MATRIX_ROOT_DEFAULT}" CACHE PATH "Test-only source root")
# This isolated replay links frozen archives from a Conda-backed dependency
# prefix.  Stage the selected compiler's C++ runtime ahead of that prefix,
# matching the main project's build-tree policy; otherwise the executable can
# load Conda's older libstdc++ before it reaches any stratum assertion.
list(APPEND CMAKE_MODULE_PATH "${MATRIX_ROOT}/cmake")
include(CompilerRuntime)
set(MIS_NATIVE_ROOT "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-plane-family-20260912/source" CACHE PATH "Frozen native source root")
set(MOTION_ROOT "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-general-differential-20260912/source" CACHE PATH "Frozen general differential source root")
set(FULL_RANK_SOLVER_ROOT "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-full-rank-newton-20260912/source" CACHE PATH "Frozen solver source root")
set(CANONICAL_ROOT "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-replay-determinism-20260912/source" CACHE PATH "Frozen canonical source root")
set(FROZEN_GENERAL_BUILD "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-general-differential-20260912/build" CACHE PATH "Hash-verified general build")
set(FROZEN_CANONICAL_BUILD "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-replay-determinism-20260912/acceptance-build" CACHE PATH "Hash-verified canonical build")
set(IGL_ROOT "${MATRIX_ROOT}/libigl" CACHE PATH "libigl source root")

foreach(required MIS_NATIVE_ROOT MOTION_ROOT FULL_RANK_SOLVER_ROOT CANONICAL_ROOT FROZEN_GENERAL_BUILD FROZEN_CANONICAL_BUILD IGL_ROOT)
  if(NOT IS_DIRECTORY "${${required}}")
    message(FATAL_ERROR "${required} must bind an existing validated directory")
  endif()
endforeach()

set(FROZEN_CANONICAL_ARCHIVE "${FROZEN_CANONICAL_BUILD}/libspin_canonical_replay.a")
set(FROZEN_GENERAL_ARCHIVE "${FROZEN_GENERAL_BUILD}/libspin_general_differential.a")
set(FROZEN_NATIVE_ARCHIVE "${FROZEN_GENERAL_BUILD}/native/libmake_it_stand_core.a")
foreach(archive FROZEN_CANONICAL_ARCHIVE FROZEN_GENERAL_ARCHIVE FROZEN_NATIVE_ARCHIVE)
  if(NOT EXISTS "${${archive}}")
    message(FATAL_ERROR "${archive} must bind an existing validated archive")
  endif()
endforeach()

find_package(CGAL CONFIG REQUIRED)
find_library(MIS_BOOST_THREAD_LIBRARY boost_thread REQUIRED)
find_library(MIS_BOOST_SYSTEM_LIBRARY boost_system REQUIRED)
find_library(MIS_GMP_LIBRARY gmp REQUIRED)
find_library(MIS_MPFR_LIBRARY mpfr REQUIRED)

add_library(stratum_dependencies INTERFACE)
target_include_directories(stratum_dependencies INTERFACE
  "${MIS_NATIVE_ROOT}/include"
  "${MOTION_ROOT}/include"
  "${FULL_RANK_SOLVER_ROOT}/include"
  "${CANONICAL_ROOT}/include")
target_include_directories(stratum_dependencies SYSTEM INTERFACE
  "${IGL_ROOT}/include"
  "${IGL_ROOT}/external/eigen")
target_compile_features(stratum_dependencies INTERFACE cxx_std_17)
target_link_libraries(stratum_dependencies INTERFACE
  "${FROZEN_CANONICAL_ARCHIVE}"
  "${FROZEN_GENERAL_ARCHIVE}"
  "${FROZEN_NATIVE_ARCHIVE}"
  "${MIS_BOOST_THREAD_LIBRARY}"
  "${MIS_BOOST_SYSTEM_LIBRARY}"
  CGAL::CGAL
  "${MIS_GMP_LIBRARY}"
  "${MIS_MPFR_LIBRARY}")

add_executable(quadric_stratum_matrix "${MATRIX_ROOT}/tests/quadric_stratum_matrix.cpp")
target_link_libraries(quadric_stratum_matrix PRIVATE stratum_dependencies)

add_test(NAME stratum-matrix
  COMMAND quadric_stratum_matrix "${CMAKE_CURRENT_BINARY_DIR}/stratum-matrix.json")
add_test(NAME stratum-demand-unsupported-j
  COMMAND quadric_stratum_matrix --mode demand-unsupported-j)
set_tests_properties(stratum-demand-unsupported-j PROPERTIES WILL_FAIL TRUE)
