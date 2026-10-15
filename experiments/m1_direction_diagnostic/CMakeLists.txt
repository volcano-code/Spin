cmake_minimum_required(VERSION 3.16)
project(M1DirectionDiagnostic LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
get_filename_component(CAMPAIGN_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
foreach(required MIS_NATIVE_ROOT MOTION_ROOT FULL_RANK_SOLVER_ROOT CANONICAL_ROOT FROZEN_GENERAL_BUILD FROZEN_CANONICAL_BUILD IGL_ROOT)
  if(NOT IS_DIRECTORY "${${required}}")
    message(FATAL_ERROR "${required} must bind a validated directory")
  endif()
endforeach()
find_package(CGAL CONFIG REQUIRED)
find_library(MIS_BOOST_THREAD_LIBRARY boost_thread REQUIRED)
find_library(MIS_BOOST_SYSTEM_LIBRARY boost_system REQUIRED)
find_path(MIS_JSON_INCLUDE nlohmann/json.hpp REQUIRED)
add_executable(direction_probe "${CMAKE_CURRENT_LIST_DIR}/direction_probe.cpp")
target_compile_features(direction_probe PRIVATE cxx_std_17)
target_include_directories(direction_probe PRIVATE
  "${MIS_NATIVE_ROOT}/include" "${MOTION_ROOT}/include"
  "${FULL_RANK_SOLVER_ROOT}/include" "${CANONICAL_ROOT}/include")
target_include_directories(direction_probe SYSTEM PRIVATE
  "${IGL_ROOT}/external/eigen" "${MIS_JSON_INCLUDE}")
target_link_libraries(direction_probe PRIVATE
  "${FROZEN_CANONICAL_BUILD}/libspin_canonical_replay.a"
  "${FROZEN_GENERAL_BUILD}/libspin_general_differential.a"
  "${FROZEN_GENERAL_BUILD}/native/libmake_it_stand_core.a"
  "${MIS_BOOST_THREAD_LIBRARY}" "${MIS_BOOST_SYSTEM_LIBRARY}" CGAL::CGAL)

target_compile_definitions(direction_probe PRIVATE MIS_FROZEN_SOLVER_CPP="${MIS_NATIVE_ROOT}/src/spin/QuadricHollowingSolver.cpp")
