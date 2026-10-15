include_guard(GLOBAL)

# Explicit opt-in: historical builds keep their original arithmetic route.
# No download, install, search-path fallback or global runtime modification.
set(MAKE_IT_STAND_RANK_TWO_FLINT_ROOT "" CACHE PATH
    "Audited isolated FLINT 3.3.1 static build for authenticated root factors")
if(MAKE_IT_STAND_RANK_TWO_FLINT_ROOT)
    set(_factor_library "${MAKE_IT_STAND_RANK_TWO_FLINT_ROOT}/libflint.a")
    set(_factor_header "${MAKE_IT_STAND_RANK_TWO_FLINT_ROOT}/src/flint.h")
    if(NOT EXISTS "${_factor_library}" OR NOT EXISTS "${_factor_header}")
        message(FATAL_ERROR "Native factor selection requires the explicit audited static build")
    endif()
    file(SHA256 "${_factor_library}" _factor_hash)
    if(NOT _factor_hash STREQUAL "461f465567572bad6a6083d8f7e7a428c99c3e2105f770c2ab4be6f49cf72486")
        message(FATAL_ERROR "Native factor static library differs from the audited FLINT build")
    endif()
    file(READ "${_factor_header}" _factor_header_bytes)
    string(FIND "${_factor_header_bytes}" "#define FLINT_VERSION \"3.3.1\"" _factor_version)
    if(_factor_version EQUAL -1)
        message(FATAL_ERROR "Native factor headers must belong to FLINT 3.3.1")
    endif()
    add_library(rank_two_native_factor INTERFACE)
    target_include_directories(rank_two_native_factor INTERFACE
        "${MAKE_IT_STAND_RANK_TWO_FLINT_ROOT}/src")
    target_compile_definitions(rank_two_native_factor INTERFACE
        MAKE_IT_STAND_RANK_TWO_NATIVE_FACTOR=1)
    # Static FLINT references MPFR/GMP: declare those AFTER the archive rather
    # than relying on CGAL's earlier shared-library entries and --as-needed.
    target_link_libraries(rank_two_native_factor INTERFACE "${_factor_library}"
        "${MPFR_LIBRARIES}" "${GMP_LIBRARIES}")
    # Every adapter compilation already consumes this exact-arithmetic target,
    # including the isolated kernel tests; no parallel adapter implementation.
    get_target_property(_factor_core_target CGAL::CGAL_Core ALIASED_TARGET)
    if(NOT _factor_core_target)
        set(_factor_core_target CGAL::CGAL_Core)
    endif()
    target_link_libraries(${_factor_core_target} INTERFACE rank_two_native_factor)
endif()
