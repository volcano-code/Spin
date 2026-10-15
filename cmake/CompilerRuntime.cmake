# Build-tree runtime selection must agree with the selected GNU compiler.
# A dependency prefix (e.g. Conda GMP) may otherwise add an older libstdc++
# ahead of the system search path. Python may already have loaded that runtime;
# extension tests additionally retain their explicit preload policy.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    execute_process(COMMAND "${CMAKE_CXX_COMPILER}" -print-file-name=libstdc++.so.6
        RESULT_VARIABLE _spin_runtime_result OUTPUT_VARIABLE _spin_runtime
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT _spin_runtime_result STREQUAL "0" OR NOT IS_ABSOLUTE "${_spin_runtime}" OR NOT EXISTS "${_spin_runtime}")
        message(FATAL_ERROR "Selected GNU compiler's libstdc++ runtime cannot be resolved")
    endif()
    # Do not prepend the entire system library directory: that could also
    # replace explicitly selected GMP/MPFR dependencies with other versions.
    # Stage only this compiler runtime in the private build tree.
    set(_spin_runtime_dir "${CMAKE_BINARY_DIR}/compiler-runtime")
    file(MAKE_DIRECTORY "${_spin_runtime_dir}")
    configure_file("${_spin_runtime}" "${_spin_runtime_dir}/libstdc++.so.6" COPYONLY)
    list(PREPEND CMAKE_BUILD_RPATH "${_spin_runtime_dir}")
    list(REMOVE_DUPLICATES CMAKE_BUILD_RPATH)
    message(STATUS "Build runtime from selected GNU compiler: ${_spin_runtime}")
endif()
