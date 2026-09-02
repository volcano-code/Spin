cmake_minimum_required(VERSION 3.16)

get_filename_component(
    SOURCE_DIR
    "${CMAKE_CURRENT_LIST_DIR}/../.."
    ABSOLUTE
)
set(BINARY_DIR "${SOURCE_DIR}/build/reachability-invalid-revision")
file(REMOVE_RECURSE "${BINARY_DIR}")

execute_process(
    COMMAND
        "${CMAKE_COMMAND}"
        -S "${SOURCE_DIR}"
        -B "${BINARY_DIR}"
        -DMAKE_IT_STAND_BUILD_VIEWER=OFF
        -DBUILD_TESTING=OFF
        -DMAKE_IT_STAND_BUILD_REACHABILITY=ON
        -DMAKE_IT_STAND_HIGHS_GIT_TAG=definitely-not-the-pinned-revision
    RESULT_VARIABLE configure_status
    OUTPUT_VARIABLE configure_stdout
    ERROR_VARIABLE configure_stderr
)

if(configure_status EQUAL 0)
    message(FATAL_ERROR
        "reachability configuration accepted an unpinned HiGHS revision\n"
        "stdout:\n${configure_stdout}\n"
        "stderr:\n${configure_stderr}"
    )
endif()
