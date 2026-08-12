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

find_package(LIBIGL REQUIRED QUIET)

set(MAKE_IT_STAND_CORE_SOURCES
    src/carving/CarvingDiagnostics.cpp
    src/carving/CarvingSolver.cpp
    src/center_of_mass.cpp
    src/extraction/VoidMeshExtractor.cpp
    src/grid_util.cpp
    src/inner_void_mesh.cpp
    src/make_it_stand.cpp
    src/physics/BalanceEnergy.cpp
    src/physics/MassProperties.cpp
    src/voxel/Voxel.cpp
    src/voxel/VoxelGrid.cpp
    src/voxel/VoxelSurface.cpp
)

set(MAKE_IT_STAND_CORE_HEADERS
    include/carving/CarvingDiagnostics.h
    include/carving/CarvingSolver.h
    include/center_of_mass.h
    include/extraction/VoidMeshExtractor.h
    include/grid_util.h
    include/inner_void_mesh.h
    include/list3d.h
    include/make_it_stand.h
    include/physics/BalanceEnergy.h
    include/physics/MassProperties.h
    include/voxel/Voxel.h
    include/voxel/VoxelGrid.h
    include/voxel/VoxelSurface.h
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

    # These verify diagnostic/control-flow contracts and R2 mass-property
    # invariants. R4 still owns paper-faithful best-prefix carving.
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

    # Debug builds intentionally retain expensive libigl assertions and the
    # full scenario executes two complete deterministic runs.
    set_tests_properties(
        carving-full-mass-consistency
        PROPERTIES TIMEOUT 120
    )
endif()
