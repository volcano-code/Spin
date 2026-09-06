cmake_minimum_required(VERSION 3.16)
get_filename_component(SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
# Every invocation preserves earlier evidence. Do not remove or reuse old builds.
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef suffix)
set(root "${SOURCE_DIR}/build/exact-geometry-config-${suffix}")
foreach(mode IN ITEMS disabled exact-only bad-pin)
    set(exact ON)
    set(pin f7c3c8212b56c0d6dae63787efc99093f4383415)
    if(mode STREQUAL "disabled")
        set(exact OFF)
        set(pin intentionally-invalid)
    elseif(mode STREQUAL "bad-pin")
        set(pin intentionally-invalid)
    endif()
    execute_process(COMMAND "${CMAKE_COMMAND}" -S "${SOURCE_DIR}" -B "${root}/${mode}" -G Ninja
        -DMAKE_IT_STAND_BUILD_VIEWER=OFF -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF
        -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=OFF
        -DMAKE_IT_STAND_BUILD_PYTHON=OFF -DBUILD_TESTING=OFF -DLIBIGL_SKIP_DOWNLOAD=ON
        "-DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=${exact}" "-DMAKE_IT_STAND_CGAL_GIT_TAG=${pin}"
        RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(mode STREQUAL "bad-pin")
        if(status EQUAL 0 OR NOT "${error}" MATCHES "requires pinned CGAL revision")
            message(FATAL_ERROR "exact-only invalid pin was not rejected: ${output}\n${error}")
        endif()
    else()
        if(NOT status EQUAL 0)
            message(FATAL_ERROR "${mode} configuration failed: ${output}\n${error}")
        endif()
        file(READ "${root}/${mode}/build.ninja" graph)
        if(mode STREQUAL "disabled")
            if(graph MATCHES "igl_cgal|external/cgal|libgmp|libmpfr|ExactBoundaryAudit.cpp|ExactMaterialRegion.cpp")
                message(FATAL_ERROR "disabled exact feature acquired geometry dependencies")
            endif()
        else()
            # A static-library-only graph has no final link step; CGAL header
            # paths establish its backend here. Runtime linkage is tested in the
            # separate 2014=ON/Faster=OFF/Python=ON/exact=ON build.
            if(NOT graph MATCHES "ExactBoundaryAudit.cpp" OR NOT graph MATCHES "ExactMaterialRegion.cpp" OR NOT graph MATCHES "external/cgal")
                message(FATAL_ERROR "exact-only feature did not acquire the common exact backend")
            endif()
            if(graph MATCHES "src/spin/Quadric|src/spin2014/")
                message(FATAL_ERROR "exact-only feature acquired a solver")
            endif()
        endif()
    endif()
    message(STATUS "exact geometry configuration ${mode}: expected outcome; preserved ${root}/${mode}")
endforeach()

# Missing locally pinned CGAL must fail before any legacy download path. A
# minimal libigl fixture reaches precisely the optional-backend preflight.
set(fixture "${root}/missing-cgal-libigl")
file(MAKE_DIRECTORY "${fixture}/include/igl" "${fixture}/cmake")
file(COPY "${SOURCE_DIR}/libigl/include/igl/read_triangle_mesh.h" DESTINATION "${fixture}/include/igl")
file(COPY "${SOURCE_DIR}/libigl/cmake/LibiglDownloadExternal.cmake" DESTINATION "${fixture}/cmake")
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${SOURCE_DIR}" -B "${root}/missing-cgal" -G Ninja
    "-DMAKE_IT_STAND_LIBIGL_DIR=${fixture}" -DMAKE_IT_STAND_BUILD_VIEWER=OFF
    -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF
    -DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON -DBUILD_TESTING=OFF -DLIBIGL_SKIP_DOWNLOAD=ON
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status EQUAL 0 OR NOT "${error}" MATCHES "Exact geometry checks require the existing pinned CGAL checkout")
    message(FATAL_ERROR "missing exact dependency was not rejected at the offline preflight: ${output}\n${error}")
endif()
message(STATUS "missing exact CGAL rejected before dependency loading; preserved ${root}/missing-cgal")
