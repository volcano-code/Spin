cmake_minimum_required(VERSION 3.16)
project(SpinCanonicalReplay LANGUAGES CXX)
include(CTest)
get_filename_component(CANONICAL_ROOT_DEFAULT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(CANONICAL_ROOT "${CANONICAL_ROOT_DEFAULT}" CACHE PATH "New canonical module root")
set(MIS_NATIVE_ROOT "" CACHE PATH "Frozen native source root")
set(MOTION_ROOT "" CACHE PATH "Frozen differential source root")
set(FULL_RANK_SOLVER_ROOT "" CACHE PATH "Frozen solver companion source root")
set(FROZEN_GENERAL_BUILD "" CACHE PATH "Validated read-only differential build")
set(FROZEN_FULL_RANK_BUILD "" CACHE PATH "Validated read-only original solver companion build")
set(IGL_ROOT "" CACHE PATH "Native build libigl root")
foreach(required MIS_NATIVE_ROOT MOTION_ROOT FULL_RANK_SOLVER_ROOT FROZEN_GENERAL_BUILD IGL_ROOT)
  if(NOT IS_DIRECTORY "${${required}}")
    message(FATAL_ERROR "${required} must name an existing frozen/validated dependency")
  endif()
endforeach()
find_package(CGAL CONFIG REQUIRED)
find_library(MIS_BOOST_THREAD_LIBRARY boost_thread REQUIRED)
find_library(MIS_BOOST_SYSTEM_LIBRARY boost_system REQUIRED)
add_library(canonical_dependencies INTERFACE)
target_include_directories(canonical_dependencies INTERFACE
  "${MIS_NATIVE_ROOT}/include" "${MOTION_ROOT}/include"
  "${FULL_RANK_SOLVER_ROOT}/include" "${CANONICAL_ROOT}/include")
target_include_directories(canonical_dependencies SYSTEM INTERFACE
  "${IGL_ROOT}/include" "${IGL_ROOT}/external/eigen")
target_link_libraries(canonical_dependencies INTERFACE CGAL::CGAL
  "${FROZEN_GENERAL_BUILD}/libspin_general_differential.a"
  "${FROZEN_GENERAL_BUILD}/native/libmake_it_stand_core.a"
  "${MIS_BOOST_THREAD_LIBRARY}" "${MIS_BOOST_SYSTEM_LIBRARY}")
target_compile_features(canonical_dependencies INTERFACE cxx_std_17)

add_library(canonical_material_impl OBJECT "${MOTION_ROOT}/src/spin/FullRankMomentEvaluator.cpp")
target_link_libraries(canonical_material_impl PRIVATE canonical_dependencies)
target_compile_definitions(canonical_material_impl PRIVATE
  evaluateFullRankMaterialMoments=canonicalFullRankMaterialImplementation
  differentiatePolyhedralMaterialMoments=differentiateCanonicalPolyhedralMaterialMoments)
add_library(canonical_solver_impl OBJECT "${FULL_RANK_SOLVER_ROOT}/src/spin/FullRankHollowingSolver.cpp")
target_link_libraries(canonical_solver_impl PRIVATE canonical_dependencies)
target_compile_definitions(canonical_solver_impl PRIVATE
  integrateFullRankBoundaryMoments=integrateCanonicalFullRankBoundaryMoments
  evaluateFullRankMaterialMoments=evaluateCanonicalFullRankMaterialMoments
  makeFullRankMomentEvaluator=canonicalFullRankAdapterImplementation
  solveFullRankQuadricHollowing=canonicalFullRankSolverImplementation)
add_library(spin_canonical_replay STATIC
  "${CANONICAL_ROOT}/src/spin/FasterTheoremConditions.cpp"
  "${CANONICAL_ROOT}/src/spin/CanonicalPolyhedralMomentDifferential.cpp"
  "${CANONICAL_ROOT}/src/spin/CanonicalFullRankHollowing.cpp"
  $<TARGET_OBJECTS:canonical_material_impl> $<TARGET_OBJECTS:canonical_solver_impl>)
target_link_libraries(spin_canonical_replay PUBLIC canonical_dependencies)

add_executable(canonical_polyhedral_replay_regression "${CANONICAL_ROOT}/tests/canonical_polyhedral_replay_regression.cpp")
target_compile_definitions(canonical_polyhedral_replay_regression PRIVATE
  differentiatePolyhedralMaterialMoments=differentiateCanonicalPolyhedralMaterialMoments)
target_link_libraries(canonical_polyhedral_replay_regression PRIVATE spin_canonical_replay)
add_test(NAME canonical-permutation COMMAND canonical_polyhedral_replay_regression)
add_executable(canonical_full_rank_regression "${FULL_RANK_SOLVER_ROOT}/tests/full_rank_hollowing_regression.cpp")
target_compile_definitions(canonical_full_rank_regression PRIVATE
  integrateFullRankBoundaryMoments=integrateCanonicalFullRankBoundaryMoments
  evaluateFullRankMaterialMoments=evaluateCanonicalFullRankMaterialMoments
  makeFullRankMomentEvaluator=makeCanonicalFullRankMomentEvaluator
  solveFullRankQuadricHollowing=solveCanonicalFullRankQuadricHollowing)
target_link_libraries(canonical_full_rank_regression PRIVATE spin_canonical_replay)
foreach(mode adapter chain negative solver integration modes offtarget-yoyo audit replay)
  add_test(NAME full-rank-${mode} COMMAND canonical_full_rank_regression ${mode})
endforeach()
foreach(name polyhedral_moment_differential full_rank_moment)
  add_executable(canonical_${name} "${MOTION_ROOT}/tests/${name}_regression.cpp")
  target_compile_definitions(canonical_${name} PRIVATE
    evaluateFullRankMaterialMoments=evaluateCanonicalFullRankMaterialMoments
    differentiatePolyhedralMaterialMoments=differentiateCanonicalPolyhedralMaterialMoments)
  target_link_libraries(canonical_${name} PRIVATE spin_canonical_replay)
  add_test(NAME general-${name} COMMAND canonical_${name})
endforeach()
foreach(name truncation_root_certificate quadric_source_motion)
  add_test(NAME general-${name} COMMAND "${FROZEN_GENERAL_BUILD}/${name}_regression")
endforeach()
add_executable(canonical_solver_process_replay "${CANONICAL_ROOT}/tests/canonical_solver_process_replay.cpp")
target_compile_definitions(canonical_solver_process_replay PRIVATE
  integrateFullRankBoundaryMoments=integrateCanonicalFullRankBoundaryMoments
  evaluateFullRankMaterialMoments=evaluateCanonicalFullRankMaterialMoments
  makeFullRankMomentEvaluator=makeCanonicalFullRankMomentEvaluator
  solveFullRankQuadricHollowing=solveCanonicalFullRankQuadricHollowing
  M1_FROZEN_REGRESSION="${FULL_RANK_SOLVER_ROOT}/tests/full_rank_hollowing_regression.cpp")
target_link_libraries(canonical_solver_process_replay PRIVATE spin_canonical_replay)
if(NOT EXISTS "${FROZEN_FULL_RANK_BUILD}/libspin_full_rank_solver.a")
  message(FATAL_ERROR "FROZEN_FULL_RANK_BUILD must bind the existing solver for coexistence testing")
endif()
add_executable(canonical_legacy_coexistence "${CANONICAL_ROOT}/tests/canonical_legacy_coexistence.cpp")
target_compile_definitions(canonical_legacy_coexistence PRIVATE
  M1_FROZEN_REGRESSION="${FULL_RANK_SOLVER_ROOT}/tests/full_rank_hollowing_regression.cpp")
target_link_libraries(canonical_legacy_coexistence PRIVATE spin_canonical_replay
  "${FROZEN_FULL_RANK_BUILD}/libspin_full_rank_solver.a")
add_test(NAME canonical-legacy-coexistence COMMAND canonical_legacy_coexistence)
