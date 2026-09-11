cmake_minimum_required(VERSION 3.16)

function(make_it_stand_allocate_retained_test_directory root output_variable)
    if(root STREQUAL "" OR output_variable STREQUAL "")
        message(FATAL_ERROR "retained test allocation requires a root and output variable")
    endif()
    if(NOT IS_ABSOLUTE "${root}" OR NOT IS_DIRECTORY "${root}")
        message(FATAL_ERROR "retained test root must be an existing absolute directory: ${root}")
    endif()

    get_filename_component(root_real "${root}" REALPATH)
    set(runs_directory "${root_real}/spin-2014-dependency-contract-runs")
    if(EXISTS "${runs_directory}" OR IS_SYMLINK "${runs_directory}")
        if(NOT IS_DIRECTORY "${runs_directory}")
            message(FATAL_ERROR "retained test run namespace is not a directory: ${runs_directory}")
        endif()
        get_filename_component(runs_real "${runs_directory}" REALPATH)
        if(NOT runs_real STREQUAL runs_directory)
            message(FATAL_ERROR
                "retained test run namespace resolves outside its reserved path: ${runs_directory}"
            )
        endif()
    else()
        file(MAKE_DIRECTORY "${runs_directory}")
        get_filename_component(runs_real "${runs_directory}" REALPATH)
    endif()

    file(
        LOCK "${root_real}/.spin-2014-dependency-contract-allocator.lock"
        GUARD FUNCTION
        TIMEOUT 30
    )
    string(RANDOM LENGTH 16 ALPHABET 0123456789abcdef run_token)
    set(run_directory "${runs_real}/run-${run_token}")
    if(EXISTS "${run_directory}" OR IS_SYMLINK "${run_directory}")
        message(FATAL_ERROR "retained test allocation collision: ${run_directory}")
    endif()
    file(MAKE_DIRECTORY "${run_directory}")
    get_filename_component(run_real "${run_directory}" REALPATH)
    file(RELATIVE_PATH run_relative "${root_real}" "${run_real}")
    if(run_relative STREQUAL ".." OR run_relative MATCHES "^\\.\\./" OR IS_ABSOLUTE "${run_relative}")
        message(FATAL_ERROR "retained test directory escaped its root: ${run_real}")
    endif()
    set(${output_variable} "${run_real}" PARENT_SCOPE)
endfunction()
