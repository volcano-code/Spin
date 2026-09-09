cmake_minimum_required(VERSION 3.16)

get_filename_component(
    SOURCE_DIR
    "${CMAKE_CURRENT_LIST_DIR}/../.."
    ABSOLUTE
)

set(DISABLED_BINARY_DIR "${SOURCE_DIR}/build/spin-2014-disabled-contract")
file(REMOVE_RECURSE "${DISABLED_BINARY_DIR}")
execute_process(
    COMMAND
        "${CMAKE_COMMAND}"
        -S "${SOURCE_DIR}"
        -B "${DISABLED_BINARY_DIR}"
        -G Ninja
        -DMAKE_IT_STAND_BUILD_VIEWER=OFF
        -DBUILD_TESTING=OFF
        -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF
        -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF
        -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=OFF
    RESULT_VARIABLE disabled_status
    OUTPUT_VARIABLE disabled_stdout
    ERROR_VARIABLE disabled_stderr
)
if(NOT disabled_status EQUAL 0)
    message(FATAL_ERROR
        "disabled Spin-It 2014 configuration failed\n"
        "stdout:\n${disabled_stdout}\nstderr:\n${disabled_stderr}"
    )
endif()
file(READ "${DISABLED_BINARY_DIR}/build.ninja" disabled_graph)
if(disabled_graph MATCHES "src/spin2014|src/api/SpinIt2014")
    message(FATAL_ERROR "disabled build contains Spin-It 2014 implementation")
endif()

set(ENABLED_BINARY_DIR "${SOURCE_DIR}/build/spin-2014-enabled-contract")
file(REMOVE_RECURSE "${ENABLED_BINARY_DIR}")
execute_process(
    COMMAND
        "${CMAKE_COMMAND}"
        -S "${SOURCE_DIR}"
        -B "${ENABLED_BINARY_DIR}"
        -G Ninja
        -DMAKE_IT_STAND_BUILD_VIEWER=OFF
        -DBUILD_TESTING=ON
        -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF
        -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF
        -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=ON
    RESULT_VARIABLE enabled_status
    OUTPUT_VARIABLE enabled_stdout
    ERROR_VARIABLE enabled_stderr
)
if(NOT enabled_status EQUAL 0)
    message(FATAL_ERROR
        "enabled Spin-It 2014 configuration failed\n"
        "stdout:\n${enabled_stdout}\nstderr:\n${enabled_stderr}"
    )
endif()
file(READ "${ENABLED_BINARY_DIR}/build.ninja" enabled_graph)
if(NOT enabled_graph MATCHES "src/spin2014/AdaptiveVoxelSpinSolver.cpp")
    message(FATAL_ERROR "enabled build omits Spin-It 2014 implementation")
endif()
if(enabled_graph MATCHES "igl_cgal|external/cgal|libgmp|libmpfr")
    message(FATAL_ERROR "Spin-It 2014-only build unexpectedly requires CGAL")
endif()
execute_process(
    COMMAND "${CMAKE_COMMAND}" --build "${ENABLED_BINARY_DIR}" --target spin-it-2014-audited -j3
    RESULT_VARIABLE build_status
    OUTPUT_VARIABLE build_stdout
    ERROR_VARIABLE build_stderr
)
if(NOT build_status EQUAL 0)
    message(FATAL_ERROR
        "Spin-It 2014-only core did not build\n"
        "stdout:\n${build_stdout}\nstderr:\n${build_stderr}"
    )
endif()
execute_process(COMMAND "${ENABLED_BINARY_DIR}/spin-it-2014-audited"
    RESULT_VARIABLE audited_status OUTPUT_VARIABLE audited_stdout ERROR_VARIABLE audited_stderr)
if(NOT audited_status EQUAL 0)
    message(FATAL_ERROR "2014 ON / reachability OFF linked audited solver failed\n${audited_stdout}\n${audited_stderr}")
endif()
message(STATUS "2014 ON / reachability OFF actual linked audit fixture: ${audited_stdout}")
