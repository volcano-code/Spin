cmake_minimum_required(VERSION 3.16)
project(SpinGeneralDifferentialContract LANGUAGES CXX)
include(CTest)
get_filename_component(MOTION_ROOT_DEFAULT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(MOTION_ROOT "${MOTION_ROOT_DEFAULT}" CACHE PATH "New-module source tree (may be a frozen snapshot)")
set(MIS_NATIVE_ROOT "${MOTION_ROOT}" CACHE PATH "Native source tree (may be a frozen evidence snapshot)")
add_subdirectory("${MIS_NATIVE_ROOT}" native EXCLUDE_FROM_ALL)
if(NOT MAKE_IT_STAND_BUILD_SPIN_IT_FASTER OR NOT MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
  message(FATAL_ERROR "General differential requires native Faster and exact geometry")
endif()
add_library(spin_general_differential STATIC
  "${MOTION_ROOT}/src/spin/PolyhedralMomentDifferential.cpp"
  "${MOTION_ROOT}/src/spin/TruncationRootCertificate.cpp"
  "${MOTION_ROOT}/src/spin/QuadricSourceMotion.cpp"
  "${MOTION_ROOT}/src/spin/FullRankMomentEvaluator.cpp")
# In a snapshot build existing native headers must match the compiled core.
target_include_directories(spin_general_differential BEFORE PUBLIC
  "${MIS_NATIVE_ROOT}/include" "${MOTION_ROOT}/include")
target_link_libraries(spin_general_differential PUBLIC make_it_stand_core)
target_compile_features(spin_general_differential PUBLIC cxx_std_17)
foreach(name IN ITEMS polyhedral_moment_differential truncation_root_certificate quadric_source_motion full_rank_moment)
  add_executable(${name}_regression "${MOTION_ROOT}/tests/${name}_regression.cpp")
  target_link_libraries(${name}_regression PRIVATE spin_general_differential)
  add_test(NAME general-${name} COMMAND ${name}_regression)
endforeach()
