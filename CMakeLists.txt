cmake_minimum_required(VERSION 3.16...3.28)

project(
    make_it_stand_reproduction
    VERSION 0.1.0
    LANGUAGES CXX
)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0")
    message(FATAL_ERROR
        "The pinned libigl revision is not compatible with CMake 4.x. "
        "Use the verified CMake 3.28.3 toolchain documented in "
        "docs/reproduction_environment.md."
    )
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(MAKE_IT_STAND_BUILD_VIEWER "Build the interactive GLFW viewer" ON)

# Keep dependency resolution anchored to the git submodule. A different
# checkout must be selected explicitly through MAKE_IT_STAND_LIBIGL_DIR.
set(
    MAKE_IT_STAND_LIBIGL_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/libigl"
    CACHE PATH
    "Path to the pinned libigl checkout"
)

if(NOT EXISTS "${MAKE_IT_STAND_LIBIGL_DIR}/include/igl/read_triangle_mesh.h")
    message(FATAL_ERROR
        "Pinned libigl checkout not found at ${MAKE_IT_STAND_LIBIGL_DIR}. "
        "Run: git submodule update --init --recursive"
    )
endif()

set(LIBIGL_DIR "${MAKE_IT_STAND_LIBIGL_DIR}" CACHE PATH "libigl root" FORCE)
option(LIBIGL_WITH_OPENGL "Use OpenGL" ${MAKE_IT_STAND_BUILD_VIEWER})
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW" ${MAKE_IT_STAND_BUILD_VIEWER})
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui" ${MAKE_IT_STAND_BUILD_VIEWER})
set(LIBIGL_WITH_MOSEK OFF CACHE BOOL "Do not use proprietary MOSEK" FORCE)
set(
    LIBIGL_USE_STATIC_LIBRARY
    OFF
    CACHE BOOL
    "Use the verified header-only libigl core configuration"
    FORCE
)

find_package(LIBIGL REQUIRED QUIET)

set(MAKE_IT_STAND_CORE_SOURCES
    src/carving/CarvingDiagnostics.cpp
    src/carving/PlaneCarvingPolicy.cpp
    src/carving/CarvingSolver.cpp
    src/center_of_mass.cpp
    src/deformation/BoundedBiharmonicWeights.cpp
    src/deformation/DeformationObjective.cpp
    src/deformation/DeformationOptimizer.cpp
    src/deformation/DeformationSystem.cpp
    src/deformation/GridDeformationDomain.cpp
    src/deformation/HandleSet.cpp
    src/deformation/LaplacianEditingEnergy.cpp
    src/deformation/LinearBlendSkinning.cpp
    src/extraction/VoidMeshExtractor.cpp
    src/grid_util.cpp
    src/inner_void_mesh.cpp
    src/make_it_stand.cpp
    src/physics/BalanceEnergy.cpp
    src/physics/BalanceObjective.cpp
    src/physics/MassProperties.cpp
    src/physics/SurfaceMassDifferential.cpp
    src/stability/SupportRegion.cpp
    src/voxel/Voxel.cpp
    src/voxel/VoxelGrid.cpp
    src/voxel/VoxelSurface.cpp
    src/voxel/VoxelTopology.cpp
)

set(MAKE_IT_STAND_CORE_HEADERS
    include/carving/CarvingDiagnostics.h
    include/carving/PlaneCarvingPolicy.h
    include/carving/CarvingSolver.h
    include/center_of_mass.h
    include/deformation/BoundedBiharmonicWeights.h
    include/deformation/DeformationObjective.h
    include/deformation/DeformationOptimizer.h
    include/deformation/DeformationSystem.h
    include/deformation/GridDeformationDomain.h
    include/deformation/HandleSet.h
    include/deformation/LaplacianEditingEnergy.h
    include/deformation/LinearBlendSkinning.h
    include/extraction/VoidMeshExtractor.h
    include/grid_util.h
    include/inner_void_mesh.h
    include/list3d.h
    include/make_it_stand.h
    include/physics/BalanceEnergy.h
    include/physics/BalanceObjective.h
    include/physics/MassProperties.h
    include/physics/SurfaceMassDifferential.h
    include/stability/SupportRegion.h
    include/voxel/Voxel.h
    include/voxel/VoxelGrid.h
    include/voxel/VoxelSurface.h
    include/voxel/VoxelTopology.h
)

add_library(
    make_it_stand_core
    STATIC
    ${MAKE_IT_STAND_CORE_SOURCES}
    ${MAKE_IT_STAND_CORE_HEADERS}
)
add_library(make_it_stand::core ALIAS make_it_stand_core)

target_compile_features(make_it_stand_core PUBLIC cxx_std_17)
set_target_properties(make_it_stand_core PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(
    make_it_stand_core
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(make_it_stand_core PUBLIC igl::core)

if(MAKE_IT_STAND_BUILD_VIEWER)
    add_executable(make_it_stand_viewer main.cpp)
    set_target_properties(
        make_it_stand_viewer
        PROPERTIES OUTPUT_NAME "make-it-stand"
    )
    target_compile_features(make_it_stand_viewer PRIVATE cxx_std_17)
    set_target_properties(make_it_stand_viewer PROPERTIES CXX_EXTENSIONS OFF)
    target_include_directories(
        make_it_stand_viewer
        PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    target_link_libraries(
        make_it_stand_viewer
        PRIVATE
            make_it_stand::core
            igl::opengl
            igl::opengl_glfw
            igl::opengl_glfw_imgui
    )
endif()

include(CTest)

if(BUILD_TESTING)
    add_executable(
        mass-properties-regression
        tests/mass_properties_regression.cpp
    )
    target_compile_features(
        mass-properties-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        mass-properties-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        mass-properties-regression
        PRIVATE make_it_stand::core
    )

    function(add_mass_properties_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:mass-properties-regression>"
                "${scenario}"
        )

        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r2;physics"
                TIMEOUT 30
        )
    endfunction()

    add_mass_properties_test(mass-surface-analytic surface)
    add_mass_properties_test(mass-voxel-analytic voxel)
    add_mass_properties_test(mass-deletion-exhaustive deletion)
    add_mass_properties_test(balance-energy-analytic energy)

    add_executable(
        stability-regression
        tests/stability_regression.cpp
    )
    target_compile_features(
        stability-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        stability-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        stability-regression
        PRIVATE make_it_stand::core
    )

    function(add_stability_regression_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:stability-regression>"
                "${scenario}"
        )

        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r3;physics;stability"
                TIMEOUT 30
        )
    endfunction()

    add_stability_regression_test(
        support-polygon-analytic
        support
    )
    add_stability_regression_test(
        safe-region-analytic
        safe
    )
    add_stability_regression_test(
        standing-target-analytic
        standing-target
    )
    add_stability_regression_test(
        suspension-target-analytic
        suspension
    )
    add_stability_regression_test(
        multi-target-energy-analytic
        multi
    )

    add_executable(
        objective-integration-regression
        tests/objective_integration_regression.cpp
    )
    target_compile_features(
        objective-integration-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        objective-integration-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        objective-integration-regression
        PRIVATE make_it_stand::core
    )

    add_test(
        NAME multi-target-solver-boundary
        COMMAND "$<TARGET_FILE:objective-integration-regression>"
            "${CMAKE_CURRENT_SOURCE_DIR}/data/bunny.off"
    )
    set_tests_properties(
        multi-target-solver-boundary
        PROPERTIES
            LABELS "r3;r4;integration;carving"
            TIMEOUT 30
    )

    add_executable(
        plane-carving-policy-regression
        tests/plane_carving_policy_regression.cpp
    )
    target_compile_features(
        plane-carving-policy-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        plane-carving-policy-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        plane-carving-policy-regression
        PRIVATE make_it_stand::core
    )

    function(add_plane_carving_policy_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:plane-carving-policy-regression>"
                "${scenario}"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r4;carving;physics"
                TIMEOUT 30
        )
    endfunction()

    add_plane_carving_policy_test(cutting-plane-analytic planes)
    add_plane_carving_policy_test(multi-plane-score-analytic multi-score)
    add_plane_carving_policy_test(carving-order-deterministic ordering)
    add_plane_carving_policy_test(global-prefix-transient transient-prefix)
    add_plane_carving_policy_test(global-prefix-middle middle-prefix)
    add_plane_carving_policy_test(global-prefix-zero no-improvement)
    add_plane_carving_policy_test(global-prefix-shortest-tie tie-prefix)
    add_plane_carving_policy_test(global-prefix-scale-invariant scale-prefix)
    add_plane_carving_policy_test(global-prefix-over-1000 over-thousand)

    add_executable(
        deformation-regression
        tests/deformation_regression.cpp
    )
    target_compile_features(deformation-regression PRIVATE cxx_std_17)
    set_target_properties(
        deformation-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        deformation-regression
        PRIVATE make_it_stand::core
    )

    function(add_deformation_regression_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:deformation-regression>" "${scenario}"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r5;deformation"
                TIMEOUT 60
        )
    endfunction()

    add_deformation_regression_test(bbw-path-analytic bbw-path)
    add_deformation_regression_test(bbw-active-bound-analytic bbw-active-bound)
    add_deformation_regression_test(bbw-contract-analytic bbw-contract)
    add_deformation_regression_test(handle-parameter-contract handles)
    add_deformation_regression_test(lbs-equation-9-10 lbs)
    add_deformation_regression_test(grid-weight-trilinear domain-trilinear)
    add_deformation_regression_test(deformation-geometry-sync geometry-sync)
    add_deformation_regression_test(deformation-invalid-inputs invalid)

    add_executable(
        deformation-optimization-regression
        tests/deformation_optimization_regression.cpp
    )
    target_compile_features(
        deformation-optimization-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        deformation-optimization-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        deformation-optimization-regression
        PRIVATE make_it_stand::core
    )

    function(add_deformation_optimization_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:deformation-optimization-regression>"
                "${scenario}"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r6;deformation;physics;optimization"
                TIMEOUT 60
        )
    endfunction()

    add_deformation_optimization_test(
        laplacian-editing-analytic
        laplacian
    )
    add_deformation_optimization_test(
        surface-mass-gradient-analytic
        mass-gradient
    )
    add_deformation_optimization_test(
        deformation-objective-gradient-analytic
        objective-gradient
    )
    add_deformation_optimization_test(
        paper-line-search-state-machine
        line-search
    )
    add_deformation_optimization_test(
        deformation-optimizer-end-to-end
        optimizer-e2e
    )

    add_executable(
        voxel-constraints-regression
        tests/voxel_constraints_regression.cpp
    )
    target_compile_features(
        voxel-constraints-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        voxel-constraints-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        voxel-constraints-regression
        PRIVATE make_it_stand::core
    )

    add_test(
        NAME minimum-wall-thickness-analytic
        COMMAND "$<TARGET_FILE:voxel-constraints-regression>" wall
    )
    add_test(
        NAME voxel-topology-report
        COMMAND "$<TARGET_FILE:voxel-constraints-regression>" topology
    )
    set_tests_properties(
        minimum-wall-thickness-analytic
        voxel-topology-report
        PROPERTIES
            LABELS "r4;carving;manufacturing"
            TIMEOUT 30
    )

    add_executable(carving-regression tests/carving_regression.cpp)
    target_compile_features(carving-regression PRIVATE cxx_std_17)
    set_target_properties(carving-regression PROPERTIES CXX_EXTENSIONS OFF)
    target_link_libraries(carving-regression PRIVATE make_it_stand::core)

    function(add_carving_regression_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:carving-regression>"
                "${scenario}"
                "${CMAKE_CURRENT_SOURCE_DIR}/data/bunny.off"
                "${CMAKE_CURRENT_BINARY_DIR}/${test_name}.csv"
        )

        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r1;r2;carving"
                TIMEOUT 30
        )
    endfunction()

    # These retain R1/R2 diagnostic and mass-property invariants while the
    # one-step/full cases also exercise the integrated R4 best-prefix policy.
    add_carving_regression_test(
        carving-initial-state
        initial
    )
    add_carving_regression_test(
        carving-one-step
        one-step
    )
    add_carving_regression_test(
        carving-full-mass-consistency
        full-mass-consistency
    )
    add_carving_regression_test(
        carving-no-candidates
        no-candidates
    )

    set_tests_properties(
        carving-one-step
        carving-full-mass-consistency
        PROPERTIES LABELS "r1;r2;r4;carving"
    )

    # Debug builds intentionally retain expensive libigl assertions and the
    # full scenario executes two complete deterministic runs.
    set_tests_properties(
        carving-full-mass-consistency
        PROPERTIES TIMEOUT 120
    )
endif()
