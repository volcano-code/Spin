cmake_minimum_required(VERSION 3.16)

get_filename_component(
    SOURCE_DIR
    "${CMAKE_CURRENT_LIST_DIR}/../.."
    ABSOLUTE
)
get_filename_component(SOURCE_DIR_REAL "${SOURCE_DIR}" REALPATH)

if(NOT DEFINED MAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT
   OR MAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT STREQUAL "")
    message(FATAL_ERROR "MAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT is required")
endif()
if(NOT IS_ABSOLUTE "${MAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT}"
   OR NOT IS_DIRECTORY "${MAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT}")
    message(FATAL_ERROR
        "MAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT must be an existing absolute directory"
    )
endif()
get_filename_component(
    TEST_BINARY_ROOT_REAL
    "${MAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT}"
    REALPATH
)
if(TEST_BINARY_ROOT_REAL STREQUAL SOURCE_DIR_REAL)
    message(FATAL_ERROR "dependency test binary root must not be the source root")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/allocate_retained_test_directory.cmake")
make_it_stand_allocate_retained_test_directory(
    "${TEST_BINARY_ROOT_REAL}"
    RETAINED_RUN_DIR
)
message("SPIN_IT_2014_DEPENDENCY_RETAINED_DIR=${RETAINED_RUN_DIR}")
file(WRITE "${RETAINED_RUN_DIR}/contract.status" "RUNNING\n")

function(run_retained_child name)
    string(JOIN "\n" retained_command ${ARGN})
    file(WRITE "${RETAINED_RUN_DIR}/${name}.command" "${retained_command}\n")
    file(WRITE "${RETAINED_RUN_DIR}/${name}.status" "RUNNING\n")
    execute_process(
        COMMAND ${ARGN}
        RESULT_VARIABLE child_status
        OUTPUT_FILE "${RETAINED_RUN_DIR}/${name}.stdout.log"
        ERROR_FILE "${RETAINED_RUN_DIR}/${name}.stderr.log"
    )
    file(WRITE "${RETAINED_RUN_DIR}/${name}.status" "${child_status}\n")
    set(${name}_status "${child_status}" PARENT_SCOPE)
endfunction()

function(fail_retained_contract failure_message)
    file(WRITE "${RETAINED_RUN_DIR}/contract.status" "FAILED\n")
    message(FATAL_ERROR "${failure_message}\nRetained evidence: ${RETAINED_RUN_DIR}")
endfunction()

set(DISABLED_BINARY_DIR "${RETAINED_RUN_DIR}/disabled")
run_retained_child(
    configure-disabled
        "${CMAKE_COMMAND}"
        -S "${SOURCE_DIR}"
        -B "${DISABLED_BINARY_DIR}"
        -G Ninja
        -DMAKE_IT_STAND_BUILD_VIEWER=OFF
        -DBUILD_TESTING=OFF
        -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF
        -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF
        -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=OFF
        -DLIBIGL_SKIP_DOWNLOAD=ON
)
if(NOT "${configure-disabled_status}" STREQUAL "0")
    fail_retained_contract("disabled Spin-It 2014 configuration failed")
endif()
file(READ "${DISABLED_BINARY_DIR}/build.ninja" disabled_graph)
if(disabled_graph MATCHES "src/spin2014|src/api/SpinIt2014")
    fail_retained_contract("disabled build contains Spin-It 2014 implementation")
endif()

set(ENABLED_BINARY_DIR "${RETAINED_RUN_DIR}/enabled")
run_retained_child(
    configure-enabled
        "${CMAKE_COMMAND}"
        -S "${SOURCE_DIR}"
        -B "${ENABLED_BINARY_DIR}"
        -G Ninja
        -DMAKE_IT_STAND_BUILD_VIEWER=OFF
        -DBUILD_TESTING=ON
        -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF
        -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF
        -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=ON
        -DLIBIGL_SKIP_DOWNLOAD=ON
)
if(NOT "${configure-enabled_status}" STREQUAL "0")
    fail_retained_contract("enabled Spin-It 2014 configuration failed")
endif()
file(READ "${ENABLED_BINARY_DIR}/build.ninja" enabled_graph)
if(NOT enabled_graph MATCHES "src/spin2014/AdaptiveVoxelSpinSolver.cpp")
    fail_retained_contract("enabled build omits Spin-It 2014 implementation")
endif()
if(enabled_graph MATCHES "igl_cgal|external/cgal|libgmp|libmpfr")
    fail_retained_contract("Spin-It 2014-only build unexpectedly requires CGAL")
endif()
run_retained_child(
    build-enabled
    "${CMAKE_COMMAND}" --build "${ENABLED_BINARY_DIR}" --target spin-it-2014-audited -j3
)
if(NOT "${build-enabled_status}" STREQUAL "0")
    fail_retained_contract("Spin-It 2014-only core did not build")
endif()
run_retained_child(audited "${ENABLED_BINARY_DIR}/spin-it-2014-audited")
if(NOT "${audited_status}" STREQUAL "0")
    fail_retained_contract("2014 ON / reachability OFF linked audited solver failed")
endif()
file(READ "${RETAINED_RUN_DIR}/audited.stdout.log" audited_stdout)
file(WRITE "${RETAINED_RUN_DIR}/contract.status" "PASSED\n")
message(STATUS "2014 ON / reachability OFF actual linked audit fixture: ${audited_stdout}")
