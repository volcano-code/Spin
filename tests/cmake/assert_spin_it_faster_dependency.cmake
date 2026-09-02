cmake_minimum_required(VERSION 3.16)

get_filename_component(
    SOURCE_DIR
    "${CMAKE_CURRENT_LIST_DIR}/../.."
    ABSOLUTE
)

set(DISABLED_BINARY_DIR "${SOURCE_DIR}/build/spin-faster-disabled-contract")
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
        -DMAKE_IT_STAND_CGAL_GIT_TAG=definitely-not-the-pinned-revision
    RESULT_VARIABLE disabled_status
    OUTPUT_VARIABLE disabled_stdout
    ERROR_VARIABLE disabled_stderr
)
if(NOT disabled_status EQUAL 0)
    message(FATAL_ERROR
        "disabled Faster configuration consulted the CGAL pin\n"
        "stdout:\n${disabled_stdout}\n"
        "stderr:\n${disabled_stderr}"
    )
endif()
file(READ "${DISABLED_BINARY_DIR}/CMakeCache.txt" disabled_cache)
if(NOT disabled_cache MATCHES "LIBIGL_WITH_CGAL:BOOL=OFF")
    message(FATAL_ERROR "disabled Faster configuration enabled CGAL")
endif()
file(READ "${DISABLED_BINARY_DIR}/build.ninja" disabled_build_graph)
if(disabled_build_graph MATCHES
   "igl_cgal|external/cgal|libgmp|libmpfr")
    message(FATAL_ERROR
        "disabled Faster build graph contains a CGAL/GMP/MPFR dependency"
    )
endif()

set(INVALID_BINARY_DIR "${SOURCE_DIR}/build/spin-faster-invalid-revision")
file(REMOVE_RECURSE "${INVALID_BINARY_DIR}")
execute_process(
    COMMAND
        "${CMAKE_COMMAND}"
        -S "${SOURCE_DIR}"
        -B "${INVALID_BINARY_DIR}"
        -G Ninja
        -DMAKE_IT_STAND_BUILD_VIEWER=OFF
        -DBUILD_TESTING=OFF
        -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF
        -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=ON
        -DMAKE_IT_STAND_CGAL_GIT_TAG=definitely-not-the-pinned-revision
    RESULT_VARIABLE invalid_status
    OUTPUT_VARIABLE invalid_stdout
    ERROR_VARIABLE invalid_stderr
)
if(invalid_status EQUAL 0)
    message(FATAL_ERROR
        "Faster configuration accepted an unpinned CGAL revision\n"
        "stdout:\n${invalid_stdout}\n"
        "stderr:\n${invalid_stderr}"
    )
endif()
