include_guard(GLOBAL)

set(
    MAKE_IT_STAND_CGAL_PINNED_REVISION
    "f7c3c8212b56c0d6dae63787efc99093f4383415"
)
set(
    MAKE_IT_STAND_CGAL_GIT_TAG
    "${MAKE_IT_STAND_CGAL_PINNED_REVISION}"
    CACHE STRING
    "Exact CGAL revision used by Spin-It Faster mesh Booleans"
)
if(NOT MAKE_IT_STAND_CGAL_GIT_TAG STREQUAL
       MAKE_IT_STAND_CGAL_PINNED_REVISION)
    message(FATAL_ERROR
        "Spin-It Faster requires pinned CGAL revision "
        "${MAKE_IT_STAND_CGAL_PINNED_REVISION}; received "
        "${MAKE_IT_STAND_CGAL_GIT_TAG}"
    )
endif()

# The parent repository pins libigl, whose download declaration is the source
# of truth for the CGAL checkout.  Refuse to proceed if that declaration ever
# drifts instead of silently fetching another geometry kernel.
set(
    MAKE_IT_STAND_LIBIGL_DOWNLOAD_DECLARATION
    "${MAKE_IT_STAND_LIBIGL_DIR}/cmake/LibiglDownloadExternal.cmake"
)
file(READ
    "${MAKE_IT_STAND_LIBIGL_DOWNLOAD_DECLARATION}"
    MAKE_IT_STAND_LIBIGL_DOWNLOAD_BYTES
)
string(FIND
    "${MAKE_IT_STAND_LIBIGL_DOWNLOAD_BYTES}"
    "GIT_TAG        ${MAKE_IT_STAND_CGAL_PINNED_REVISION}"
    MAKE_IT_STAND_CGAL_PIN_OFFSET
)
if(MAKE_IT_STAND_CGAL_PIN_OFFSET EQUAL -1)
    message(FATAL_ERROR
        "Pinned libigl no longer declares the required CGAL revision "
        "${MAKE_IT_STAND_CGAL_PINNED_REVISION}"
    )
endif()

set(LIBIGL_WITH_CGAL ON CACHE BOOL "Use pinned CGAL" FORCE)


function(make_it_stand_verify_spin_it_faster_dependencies)
    if(NOT TARGET igl::cgal)
        message(FATAL_ERROR "Spin-It Faster requires target igl::cgal")
    endif()
    if(NOT TARGET CGAL::CGAL OR NOT TARGET CGAL::CGAL_Core)
        message(FATAL_ERROR
            "Spin-It Faster requires CGAL::CGAL and CGAL::CGAL_Core"
        )
    endif()

    find_package(Git REQUIRED QUIET)
    set(
        MAKE_IT_STAND_CGAL_SOURCE_DIR
        "${MAKE_IT_STAND_LIBIGL_DIR}/external/cgal"
    )
    execute_process(
        COMMAND
            "${GIT_EXECUTABLE}" -C "${MAKE_IT_STAND_CGAL_SOURCE_DIR}"
            rev-parse HEAD
        RESULT_VARIABLE MAKE_IT_STAND_CGAL_HEAD_STATUS
        OUTPUT_VARIABLE MAKE_IT_STAND_CGAL_HEAD
        ERROR_VARIABLE MAKE_IT_STAND_CGAL_HEAD_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT MAKE_IT_STAND_CGAL_HEAD_STATUS EQUAL 0)
        message(FATAL_ERROR
            "Unable to audit downloaded CGAL checkout: "
            "${MAKE_IT_STAND_CGAL_HEAD_ERROR}"
        )
    endif()
    if(NOT MAKE_IT_STAND_CGAL_HEAD STREQUAL
           MAKE_IT_STAND_CGAL_PINNED_REVISION)
        message(FATAL_ERROR
            "Downloaded CGAL revision ${MAKE_IT_STAND_CGAL_HEAD} does not "
            "match ${MAKE_IT_STAND_CGAL_PINNED_REVISION}"
        )
    endif()

    find_program(MAKE_IT_STAND_PKG_CONFIG_EXECUTABLE pkg-config REQUIRED)
    execute_process(
        COMMAND "${MAKE_IT_STAND_PKG_CONFIG_EXECUTABLE}" --modversion gmp
        RESULT_VARIABLE MAKE_IT_STAND_GMP_STATUS
        OUTPUT_VARIABLE MAKE_IT_STAND_GMP_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND "${MAKE_IT_STAND_PKG_CONFIG_EXECUTABLE}" --modversion mpfr
        RESULT_VARIABLE MAKE_IT_STAND_MPFR_STATUS
        OUTPUT_VARIABLE MAKE_IT_STAND_MPFR_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(NOT MAKE_IT_STAND_GMP_STATUS EQUAL 0 OR
       NOT MAKE_IT_STAND_MPFR_STATUS EQUAL 0)
        message(FATAL_ERROR "Unable to audit GMP/MPFR package versions")
    endif()
    if(NOT MAKE_IT_STAND_GMP_VERSION STREQUAL "6.3.0" OR
       NOT MAKE_IT_STAND_MPFR_VERSION STREQUAL "4.2.1")
        message(FATAL_ERROR
            "Spin-It Faster requires audited GMP 6.3.0 and MPFR 4.2.1; "
            "found GMP ${MAKE_IT_STAND_GMP_VERSION}, "
            "MPFR ${MAKE_IT_STAND_MPFR_VERSION}"
        )
    endif()

    set(
        MAKE_IT_STAND_CGAL_RESOLVED_REVISION
        "${MAKE_IT_STAND_CGAL_HEAD}"
        CACHE INTERNAL "Resolved Spin-It Faster CGAL revision" FORCE
    )
    set(
        MAKE_IT_STAND_GMP_RESOLVED_VERSION
        "${MAKE_IT_STAND_GMP_VERSION}"
        CACHE INTERNAL "Resolved GMP version" FORCE
    )
    set(
        MAKE_IT_STAND_MPFR_RESOLVED_VERSION
        "${MAKE_IT_STAND_MPFR_VERSION}"
        CACHE INTERNAL "Resolved MPFR version" FORCE
    )
endfunction()
