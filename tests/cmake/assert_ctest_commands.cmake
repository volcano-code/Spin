if(NOT DEFINED BINARY_DIR)
    message(FATAL_ERROR "BINARY_DIR is required")
endif()
if(NOT DEFINED CTEST_COMMAND)
    message(FATAL_ERROR "CTEST_COMMAND is required")
endif()

execute_process(
    COMMAND "${CTEST_COMMAND}" --test-dir "${BINARY_DIR}"
        --show-only=json-v1
    RESULT_VARIABLE status
    OUTPUT_VARIABLE document
    ERROR_VARIABLE error_output
)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "CTest registration query failed: ${error_output}")
endif()

string(JSON test_count ERROR_VARIABLE json_error LENGTH "${document}" tests)
if(json_error)
    message(FATAL_ERROR "CTest registration JSON is invalid: ${json_error}")
endif()
if(test_count EQUAL 0)
    message(FATAL_ERROR "CTest registered no tests")
endif()

math(EXPR last_test "${test_count} - 1")
foreach(index RANGE 0 ${last_test})
    string(JSON test_name GET "${document}" tests ${index} name)
    string(
        JSON command_count ERROR_VARIABLE command_error
        LENGTH "${document}" tests ${index} command
    )
    if(command_error OR command_count EQUAL 0)
        message(FATAL_ERROR "CTest test '${test_name}' has an empty command")
    endif()
    string(JSON executable GET "${document}" tests ${index} command 0)
    if(executable STREQUAL "")
        message(FATAL_ERROR "CTest test '${test_name}' has a blank executable")
    endif()
endforeach()
