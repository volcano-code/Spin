cmake_minimum_required(VERSION 3.16)
project(M1FinalAuditTests LANGUAGES CXX)
include(CTest)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
get_filename_component(AUDIT_ROOT_DEFAULT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(AUDIT_ROOT "${AUDIT_ROOT_DEFAULT}" CACHE PATH "Test-only source root")
foreach(required MIS_NATIVE_ROOT MOTION_ROOT FULL_RANK_SOLVER_ROOT CANONICAL_ROOT FROZEN_GENERAL_BUILD FROZEN_CANONICAL_BUILD IGL_ROOT)
  if(NOT IS_DIRECTORY "${${required}}")
    message(FATAL_ERROR "${required} must bind an existing validated directory")
  endif()
endforeach()
find_package(CGAL CONFIG REQUIRED)
find_library(MIS_BOOST_THREAD_LIBRARY boost_thread REQUIRED)
find_library(MIS_BOOST_SYSTEM_LIBRARY boost_system REQUIRED)
add_library(audit_dependencies INTERFACE)
target_include_directories(audit_dependencies INTERFACE
  "${MIS_NATIVE_ROOT}/include" "${MOTION_ROOT}/include"
  "${FULL_RANK_SOLVER_ROOT}/include" "${CANONICAL_ROOT}/include")
target_include_directories(audit_dependencies SYSTEM INTERFACE
  "${IGL_ROOT}/include" "${IGL_ROOT}/external/eigen")
target_compile_features(audit_dependencies INTERFACE cxx_std_17)
target_link_libraries(audit_dependencies INTERFACE
  "${FROZEN_CANONICAL_BUILD}/libspin_canonical_replay.a"
  "${FROZEN_GENERAL_BUILD}/libspin_general_differential.a"
  "${FROZEN_GENERAL_BUILD}/native/libmake_it_stand_core.a"
  "${MIS_BOOST_THREAD_LIBRARY}" "${MIS_BOOST_SYSTEM_LIBRARY}" CGAL::CGAL)
# This compilation is test-only. It retains the exact frozen audit body and
# changes dependency symbols, not the checks, tolerances or numerical solver.
set(AUDIT_SUBJECT_SOURCE "${FULL_RANK_SOLVER_ROOT}/src/spin/FullRankHollowingSolver.cpp"
  CACHE FILEPATH "Frozen audit source; override only for named mutation controls")
add_library(audit_subject OBJECT "${AUDIT_SUBJECT_SOURCE}")
target_link_libraries(audit_subject PRIVATE audit_dependencies)
target_compile_definitions(audit_subject PRIVATE
  evaluateFullRankMaterialMoments=auditTestMaterialEvaluation
  solveSpinMomentSystem=auditTestNumericalSolve
  makeFullRankMomentEvaluator=auditTestMakeEvaluator
  solveFullRankQuadricHollowing=auditTestHollowingSolve
  integrateFullRankBoundaryMoments=auditTestBoundaryIntegration)
add_executable(final_audit_regression "${AUDIT_ROOT}/tests/final_audit_regression.cpp"
  $<TARGET_OBJECTS:audit_subject>)
target_link_libraries(final_audit_regression PRIVATE audit_dependencies)
foreach(mode control frame missing-j nonfinite-j moments inward stale-lambda invalid-mu missing-quadric residual)
  add_test(NAME final-audit-${mode} COMMAND final_audit_regression ${mode})
endforeach()
