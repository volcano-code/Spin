include_guard(GLOBAL)

include(FetchContent)

set(
    MAKE_IT_STAND_HIGHS_PINNED_REVISION
    "7df0786de3088c832297e5ed821db236d8fab281"
)
set(
    MAKE_IT_STAND_HIGHS_GIT_TAG
    "${MAKE_IT_STAND_HIGHS_PINNED_REVISION}"
    CACHE STRING
    "Exact HiGHS revision used by carving reachability"
)
if(NOT MAKE_IT_STAND_HIGHS_GIT_TAG STREQUAL
       MAKE_IT_STAND_HIGHS_PINNED_REVISION)
    message(FATAL_ERROR
        "Carving reachability requires pinned HiGHS revision "
        "${MAKE_IT_STAND_HIGHS_PINNED_REVISION}; received "
        "${MAKE_IT_STAND_HIGHS_GIT_TAG}"
    )
endif()

set(MAKE_IT_STAND_SAVED_BUILD_TESTING "${BUILD_TESTING}")
set(BUILD_TESTING OFF CACHE BOOL "Disable dependency tests" FORCE)
set(FAST_BUILD ON CACHE BOOL "Build only the HiGHS library" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Use static dependencies" FORCE)
set(BUILD_CXX_EXE OFF CACHE BOOL "Do not build the HiGHS CLI" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build HiGHS examples" FORCE)

FetchContent_Declare(
    highs
    GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
    GIT_TAG "${MAKE_IT_STAND_HIGHS_GIT_TAG}"
    GIT_SHALLOW FALSE
    GIT_PROGRESS FALSE
)
FetchContent_MakeAvailable(highs)

set(
    BUILD_TESTING
    "${MAKE_IT_STAND_SAVED_BUILD_TESTING}"
    CACHE BOOL
    "Build Make-It-Stand tests"
    FORCE
)

if(NOT TARGET highs::highs)
    message(FATAL_ERROR "Pinned HiGHS did not provide target highs::highs")
endif()
