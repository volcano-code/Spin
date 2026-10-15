# Shared registration for the main build and the dependency-light test harness.
function(make_it_stand_register_geometry_diagnostics source_root python)
    if(NOT EXISTS "${python}")
        message(FATAL_ERROR "Geometry diagnostics require an existing Python interpreter")
    endif()
    execute_process(COMMAND "${python}" -c "import numpy"
        RESULT_VARIABLE numpy_result OUTPUT_QUIET ERROR_VARIABLE numpy_error)
    if(NOT "${numpy_result}" STREQUAL "0")
        message(FATAL_ERROR "Geometry diagnostics require NumPy: ${numpy_error}")
    endif()
    set(modules
        acute_wedge_query material_cut_certificate boundary_cap_query boundary_cap_invariance
        closest_set_adversarial closest_set_query source_star_query boundary_link_query
        segment_material_query segment_material_broadphase near_feature_witness
        parallel_plane_clearance teapot_contact_protocol)
    foreach(module IN LISTS modules)
        set(test_file "${source_root}/tests/python/test_${module}.py")
        if(NOT EXISTS "${test_file}")
            message(FATAL_ERROR "Missing required geometry diagnostic test: ${test_file}")
        endif()
        set(name "spin-geometry-${module}-python")
        add_test(NAME "${name}" COMMAND "${python}" -B
            "${source_root}/tests/python/run_required_geometry_diagnostic.py"
            "${source_root}" "test_${module}")
        set_tests_properties("${name}" PROPERTIES
            WORKING_DIRECTORY "${source_root}"
            ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;MKL_NUM_THREADS=1"
            LABELS "verification;spin2014;geometry-diagnostic;python" TIMEOUT 120)
    endforeach()
endfunction()
