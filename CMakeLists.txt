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
include(CompilerObservation)

option(MAKE_IT_STAND_BUILD_VIEWER "Build the interactive GLFW viewer" ON)
option(MAKE_IT_STAND_BUILD_PYTHON "Build the pybind11 Python extension" OFF)
option(
    MAKE_IT_STAND_BUILD_REACHABILITY
    "Build the HiGHS-backed carving reachability analyzer"
    ON
)
option(
    MAKE_IT_STAND_BUILD_SPIN_IT_FASTER
    "Build the Spin-It Faster formula-reproduction solver"
    OFF
)
option(
    MAKE_IT_STAND_BUILD_SPIN_IT_2014
    "Build the adaptive-voxel Spin-It 2014 formula-reproduction solver"
    OFF
)
option(
    MAKE_IT_STAND_ENABLE_PAPER_PYTHON_TESTS
    "Run pure-Python paper/beta compatibility contract tests"
    OFF
)
option(
    MAKE_IT_STAND_ENABLE_PYBULLET_TESTS
    "Run optional PyBullet/trimesh simulation contract tests"
    OFF
)
set(
    MAKE_IT_STAND_RELEASED_BETA_ROOT
    ""
    CACHE PATH
    "Optional root of the read-only released Make-It-Stand beta archive"
)

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
if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
    include(SpinItFasterDependencies)
else()
    set(LIBIGL_WITH_CGAL OFF CACHE BOOL "Use CGAL" FORCE)
endif()

set(LIBIGL_WITH_MOSEK OFF CACHE BOOL "Do not use proprietary MOSEK" FORCE)
set(
    LIBIGL_USE_STATIC_LIBRARY
    OFF
    CACHE BOOL
    "Use the verified header-only libigl core configuration"
    FORCE
)

find_package(LIBIGL REQUIRED QUIET)

if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
    make_it_stand_verify_spin_it_faster_dependencies()
endif()

if(MAKE_IT_STAND_BUILD_REACHABILITY)
    include(HighsDependency)
endif()

set(MAKE_IT_STAND_CORE_SOURCES
    src/api/BalanceModel.cpp
    src/carving/CarvingDiagnostics.cpp
    src/carving/PlaneCarvingPolicy.cpp
    src/carving/CarvingSolver.cpp
    src/center_of_mass.cpp
    src/deformation/BoundedBiharmonicWeights.cpp
    src/deformation/DeformedVoxelMassField.cpp
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
    src/manufacturing/ManufacturingGeometry.cpp
    src/manufacturing/PartitionManifest.cpp
    src/manufacturing/SlicerValidation.cpp
    src/optimization/AlternatingBalanceOptimizer.cpp
    src/physics/BalanceEnergy.cpp
    src/physics/BalanceObjective.cpp
    src/physics/MassProperties.cpp
    src/physics/RawMassMoments10.cpp
    src/physics/SurfaceMassDifferential.cpp
    src/stability/SupportRegion.cpp
    src/stability/BalanceTargetSchedule.cpp
    src/voxel/Voxel.cpp
    src/voxel/VoxelDomain.cpp
    src/voxel/VoxelGrid.cpp
    src/voxel/MaterialOccupancy.cpp
    src/voxel/VoxelSurface.cpp
    src/voxel/VoxelTopology.cpp
    src/viewer/ResultView.cpp
)

set(MAKE_IT_STAND_CORE_HEADERS
    include/api/BalanceModel.h
    include/carving/CarvingDiagnostics.h
    include/carving/PlaneCarvingPolicy.h
    include/carving/CarvingSolver.h
    include/center_of_mass.h
    include/deformation/BoundedBiharmonicWeights.h
    include/deformation/DeformedVoxelMassField.h
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
    include/viewer/ResultView.h
    include/manufacturing/ManufacturingGeometry.h
    include/manufacturing/PartitionManifest.h
    include/manufacturing/SlicerValidation.h
    include/optimization/AlternatingBalanceOptimizer.h
    include/physics/BalanceEnergy.h
    include/physics/BalanceObjective.h
    include/physics/MassProperties.h
    include/physics/RawMassMoments10.h
    include/physics/SurfaceMassDifferential.h
    include/stability/SupportRegion.h
    include/stability/BalanceTargetSchedule.h
    include/voxel/Voxel.h
    include/voxel/VoxelDomain.h
    include/voxel/VoxelGrid.h
    include/voxel/MaterialOccupancy.h
    include/voxel/VoxelSurface.h
    include/voxel/VoxelTopology.h
)

if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER OR MAKE_IT_STAND_BUILD_SPIN_IT_2014)
    list(APPEND MAKE_IT_STAND_CORE_SOURCES
        src/spin/SpinDesignRequest.cpp
        src/spin/SpinMomentProblem.cpp
    )
    list(APPEND MAKE_IT_STAND_CORE_HEADERS
        include/spin/SpinDesignRequest.h
        include/spin/SpinMomentProblem.h
    )
endif()

if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
    list(APPEND MAKE_IT_STAND_CORE_SOURCES
        src/spin/QuadricCoefficients.cpp
        src/spin/QuadricMomentEvaluator.cpp
        src/spin/QuadricHollowingSolver.cpp
        src/spin/QuadricTriangulator.cpp
    )
    list(APPEND MAKE_IT_STAND_CORE_HEADERS
        include/spin/QuadricCoefficients.h
        include/spin/QuadricMomentEvaluator.h
        include/spin/QuadricHollowingSolver.h
        include/spin/QuadricTriangulator.h
    )
endif()

if(MAKE_IT_STAND_BUILD_SPIN_IT_2014)
    list(APPEND MAKE_IT_STAND_CORE_SOURCES
        src/spin2014/AdaptiveVoxelDomain.cpp
        src/spin2014/AdaptiveFillFunctional.cpp
        src/spin2014/SpinIt2014MomentProblem.cpp
        src/spin2014/SlqpLbfgsSolver.cpp
        src/spin2014/AdaptiveVoxelSpinSolver.cpp
        src/spin2014/DualMaterialSpinSolver.cpp
        src/api/SpinIt2014.cpp
        src/spin2014/SurfaceRawMomentDifferential.cpp
        src/spin2014/CageSpinDeformation.cpp
        src/spin2014/CageSpinDeformationSolver.cpp
        src/spin2014/AutomaticCageBinding.cpp
    )
    list(APPEND MAKE_IT_STAND_CORE_HEADERS
        include/spin2014/AdaptiveVoxelDomain.h
        include/spin2014/AdaptiveFillFunctional.h
        include/spin2014/SpinIt2014MomentProblem.h
        include/spin2014/SlqpLbfgsSolver.h
        include/spin2014/AdaptiveVoxelSpinSolver.h
        include/spin2014/DualMaterialSpinSolver.h
        include/api/SpinIt2014.h
        include/spin2014/SurfaceRawMomentDifferential.h
        include/spin2014/CageSpinDeformation.h
        include/spin2014/CageSpinDeformationSolver.h
        include/spin2014/AutomaticCageBinding.h
    )
endif()

if(MAKE_IT_STAND_BUILD_REACHABILITY)
    list(APPEND MAKE_IT_STAND_CORE_SOURCES
        src/reachability/CarvingReachability.cpp
    )
    list(APPEND MAKE_IT_STAND_CORE_HEADERS
        include/reachability/CarvingReachability.h
    )
endif()

add_library(
    make_it_stand_core
    STATIC
    ${MAKE_IT_STAND_CORE_SOURCES}
    ${MAKE_IT_STAND_CORE_HEADERS}
)
add_library(make_it_stand::core ALIAS make_it_stand_core)

target_compile_features(make_it_stand_core PUBLIC cxx_std_17)
set_target_properties(
    make_it_stand_core
    PROPERTIES
        CXX_EXTENSIONS OFF
        POSITION_INDEPENDENT_CODE ON
)
target_include_directories(
    make_it_stand_core
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_compile_definitions(
    make_it_stand_core
    PUBLIC
        MAKE_IT_STAND_HAS_REACHABILITY=$<BOOL:${MAKE_IT_STAND_BUILD_REACHABILITY}>
        MAKE_IT_STAND_HAS_SPIN_IT_FASTER=$<BOOL:${MAKE_IT_STAND_BUILD_SPIN_IT_FASTER}>
        MAKE_IT_STAND_HAS_SPIN_IT_2014=$<BOOL:${MAKE_IT_STAND_BUILD_SPIN_IT_2014}>
)
target_link_libraries(make_it_stand_core PUBLIC igl::core)
if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
    target_link_libraries(make_it_stand_core PUBLIC igl::cgal)
endif()
if(MAKE_IT_STAND_BUILD_REACHABILITY)
    target_link_libraries(make_it_stand_core PRIVATE highs::highs)
endif()

if(MAKE_IT_STAND_BUILD_PYTHON)
    set(PYBIND11_FINDPYTHON ON)
    find_package(
        Python 3.10
        COMPONENTS Interpreter Development.Module
        REQUIRED
    )
    find_package(Python3 3.10 COMPONENTS Interpreter REQUIRED)
    get_filename_component(
        MAKE_IT_STAND_PYTHON_EXECUTABLE_REAL
        "${Python_EXECUTABLE}"
        REALPATH
    )
    get_filename_component(
        MAKE_IT_STAND_PYTHON3_EXECUTABLE_REAL
        "${Python3_EXECUTABLE}"
        REALPATH
    )
    if(NOT MAKE_IT_STAND_PYTHON_EXECUTABLE_REAL STREQUAL
       MAKE_IT_STAND_PYTHON3_EXECUTABLE_REAL)
        message(FATAL_ERROR
            "Python and Python3 resolved to different interpreters: "
            "${Python_EXECUTABLE} versus ${Python3_EXECUTABLE}"
        )
    endif()

    # Resolve pybind11 from the exact interpreter selected above. This never
    # downloads dependencies and prevents mixing two CPython environments.
    if(NOT pybind11_DIR)
        execute_process(
            COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
            RESULT_VARIABLE MAKE_IT_STAND_PYBIND11_QUERY_STATUS
            OUTPUT_VARIABLE MAKE_IT_STAND_PYBIND11_CMAKE_DIR
            ERROR_VARIABLE MAKE_IT_STAND_PYBIND11_QUERY_ERROR
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(NOT MAKE_IT_STAND_PYBIND11_QUERY_STATUS EQUAL 0)
            message(FATAL_ERROR
                "The selected Python interpreter has no usable pybind11: "
                "${MAKE_IT_STAND_PYBIND11_QUERY_ERROR}"
            )
        endif()
        set(pybind11_DIR "${MAKE_IT_STAND_PYBIND11_CMAKE_DIR}")
    endif()

    find_package(pybind11 3.0.1 EXACT CONFIG REQUIRED)
    pybind11_add_module(
        _make_it_stand
        MODULE
        bindings/pybind_module.cpp
    )
    target_compile_features(_make_it_stand PRIVATE cxx_std_17)
    set_target_properties(
        _make_it_stand
        PROPERTIES
            CXX_EXTENSIONS OFF
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/python"
    )
    target_link_libraries(_make_it_stand PRIVATE make_it_stand::core)

    # Conda may preload an older libstdc++ before importing this GCC-built
    # extension.  Ask the selected compiler for its matching runtime and use
    # it only for extension-backed CTest processes.
    set(MAKE_IT_STAND_PYTHON_TEST_ENV
        "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_BINARY_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
    )
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        execute_process(
            COMMAND "${CMAKE_CXX_COMPILER}" -print-file-name=libstdc++.so.6
            OUTPUT_VARIABLE MAKE_IT_STAND_COMPILER_LIBSTDCXX
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(IS_ABSOLUTE "${MAKE_IT_STAND_COMPILER_LIBSTDCXX}"
           AND EXISTS "${MAKE_IT_STAND_COMPILER_LIBSTDCXX}")
            list(APPEND MAKE_IT_STAND_PYTHON_TEST_ENV
                "LD_PRELOAD=${MAKE_IT_STAND_COMPILER_LIBSTDCXX}"
            )
        endif()
    endif()

    install(
        TARGETS _make_it_stand
        LIBRARY DESTINATION make_it_stand
    )
    install(
        DIRECTORY python/make_it_stand/
        DESTINATION make_it_stand
        FILES_MATCHING PATTERN "*.py"
    )
endif()

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
    if(MAKE_IT_STAND_BUILD_SPIN_IT_2014)
        add_executable(
            spin-it-2014-adaptive-domain-regression
            tests/spin_it_2014_adaptive_domain_regression.cpp
        )
        target_compile_features(
            spin-it-2014-adaptive-domain-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-adaptive-domain-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-adaptive-domain-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-adaptive-domain
            COMMAND "$<TARGET_FILE:spin-it-2014-adaptive-domain-regression>"
        )
        set_tests_properties(
            spin-it-2014-adaptive-domain
            PROPERTIES LABELS "spin-it;spin-it-2014;regression"
        )

        add_executable(
            spin-it-2014-fill-functional-regression
            tests/spin_it_2014_fill_functional_regression.cpp
        )
        target_compile_features(
            spin-it-2014-fill-functional-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-fill-functional-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-fill-functional-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-fill-functional
            COMMAND "$<TARGET_FILE:spin-it-2014-fill-functional-regression>"
        )
        set_tests_properties(
            spin-it-2014-fill-functional
            PROPERTIES LABELS "spin-it;spin-it-2014;regression"
        )

        add_executable(
            spin-it-2014-slqp-regression
            tests/spin_it_2014_slqp_regression.cpp
        )
        target_compile_features(
            spin-it-2014-slqp-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-slqp-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-slqp-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-slqp
            COMMAND "$<TARGET_FILE:spin-it-2014-slqp-regression>"
        )
        set_tests_properties(
            spin-it-2014-slqp
            PROPERTIES LABELS "spin-it;spin-it-2014;solver;regression"
        )

        add_executable(
            spin-it-2014-end-to-end
            tests/spin_it_2014_end_to_end.cpp
        )
        target_compile_features(spin-it-2014-end-to-end PRIVATE cxx_std_17)
        set_target_properties(
            spin-it-2014-end-to-end
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-end-to-end
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-end-to-end
            COMMAND "$<TARGET_FILE:spin-it-2014-end-to-end>"
        )
        set_tests_properties(
            spin-it-2014-end-to-end
            PROPERTIES LABELS "spin-it;spin-it-2014;solver;geometry;regression"
        )

        add_executable(
            spin-it-2014-modes-regression
            tests/spin_it_2014_modes_regression.cpp
        )
        target_compile_features(
            spin-it-2014-modes-regression PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-modes-regression PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-modes-regression PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-modes
            COMMAND "$<TARGET_FILE:spin-it-2014-modes-regression>"
        )
        set_tests_properties(
            spin-it-2014-modes
            PROPERTIES LABELS "spin-it;spin-it-2014;solver;regression"
        )

        add_executable(
            spin-it-2014-objectives-regression
            tests/spin_it_2014_objectives_regression.cpp
        )
        target_compile_features(
            spin-it-2014-objectives-regression PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-objectives-regression PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-objectives-regression PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-objectives
            COMMAND "$<TARGET_FILE:spin-it-2014-objectives-regression>"
        )
        set_tests_properties(
            spin-it-2014-objectives
            PROPERTIES LABELS "spin-it;spin-it-2014;objective;regression"
        )

        add_executable(
            spin-it-2014-resolution-ladder
            tests/spin_it_2014_resolution_ladder.cpp
        )
        target_compile_features(
            spin-it-2014-resolution-ladder PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-resolution-ladder PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-resolution-ladder PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-resolution-ladder
            COMMAND "$<TARGET_FILE:spin-it-2014-resolution-ladder>"
        )
        set_tests_properties(
            spin-it-2014-resolution-ladder
            PROPERTIES
                LABELS "spin-it;spin-it-2014;solver;resolution;regression"
                TIMEOUT 240
        )

        add_executable(
            spin-it-2014-dual-material-regression
            tests/spin_it_2014_dual_material_regression.cpp
        )
        target_compile_features(
            spin-it-2014-dual-material-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-dual-material-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-dual-material-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-dual-material
            COMMAND "$<TARGET_FILE:spin-it-2014-dual-material-regression>"
        )
        set_tests_properties(
            spin-it-2014-dual-material
            PROPERTIES LABELS "spin-it;spin-it-2014;dual-material;regression"
        )

        add_executable(
            spin-it-2014-api-regression
            tests/spin_it_2014_api_regression.cpp
        )
        target_compile_features(spin-it-2014-api-regression PRIVATE cxx_std_17)
        set_target_properties(
            spin-it-2014-api-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-api-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-api
            COMMAND "$<TARGET_FILE:spin-it-2014-api-regression>"
        )
        set_tests_properties(
            spin-it-2014-api
            PROPERTIES LABELS "spin-it;spin-it-2014;api;regression"
        )

        add_executable(
            spin-it-2014-surface-differential-regression
            tests/spin_it_2014_surface_differential_regression.cpp
        )
        target_compile_features(
            spin-it-2014-surface-differential-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-surface-differential-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-surface-differential-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-surface-differential
            COMMAND "$<TARGET_FILE:spin-it-2014-surface-differential-regression>"
        )
        set_tests_properties(
            spin-it-2014-surface-differential
            PROPERTIES LABELS "spin-it;spin-it-2014;deformation;regression"
        )

        add_executable(
            spin-it-2014-cage-functional-regression
            tests/spin_it_2014_cage_functional_regression.cpp
        )
        target_compile_features(
            spin-it-2014-cage-functional-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-cage-functional-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-cage-functional-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-cage-functional
            COMMAND "$<TARGET_FILE:spin-it-2014-cage-functional-regression>"
        )
        set_tests_properties(
            spin-it-2014-cage-functional
            PROPERTIES LABELS "spin-it;spin-it-2014;deformation;regression"
        )

        add_executable(
            spin-it-2014-cage-solver-regression
            tests/spin_it_2014_cage_solver_regression.cpp
        )
        target_compile_features(
            spin-it-2014-cage-solver-regression PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-cage-solver-regression PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-cage-solver-regression PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-cage-solver
            COMMAND "$<TARGET_FILE:spin-it-2014-cage-solver-regression>"
        )
        set_tests_properties(
            spin-it-2014-cage-solver
            PROPERTIES LABELS "spin-it;spin-it-2014;deformation;solver;regression"
        )

        add_executable(
            spin-it-2014-automatic-cage-regression
            tests/spin_it_2014_automatic_cage_regression.cpp
        )
        target_compile_features(
            spin-it-2014-automatic-cage-regression PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-automatic-cage-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-2014-automatic-cage-regression PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-it-2014-automatic-cage
            COMMAND "$<TARGET_FILE:spin-it-2014-automatic-cage-regression>"
        )
        set_tests_properties(
            spin-it-2014-automatic-cage
            PROPERTIES LABELS "spin-it;spin-it-2014;deformation;bbw;regression"
        )

        if(MAKE_IT_STAND_BUILD_PYTHON)
            add_test(
                NAME spin-it-2014-python
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                    test_spin_it_2014.SpinIt2014PythonTests -v
            )
            set_tests_properties(
                spin-it-2014-python
                PROPERTIES
                    LABELS "spin-it;spin-it-2014;python;regression"
                    TIMEOUT 60
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
            add_test(
                NAME spin-it-2014-paper-miniature
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                    test_spin_it_2014_paper_runner.SpinIt2014PaperRunnerNativeTests
                    -v
            )
            set_tests_properties(
                spin-it-2014-paper-miniature
                PROPERTIES
                    LABELS "spin-it;spin-it-2014;python;paper;regression"
                    TIMEOUT 180
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
        endif()
    endif()

    if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
        add_executable(
            spin-design-request-regression
            tests/spin_design_request_regression.cpp
        )
        target_compile_features(
            spin-design-request-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-design-request-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-design-request-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-request
            COMMAND "$<TARGET_FILE:spin-design-request-regression>"
        )
        set_tests_properties(
            spin-faster-request
            PROPERTIES LABELS "spin-it;spin-it-faster;regression"
        )

        add_executable(
            spin-moment-problem-regression
            tests/spin_moment_problem_regression.cpp
        )
        target_compile_features(
            spin-moment-problem-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-moment-problem-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-moment-problem-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-moment-problems
            COMMAND "$<TARGET_FILE:spin-moment-problem-regression>"
        )
        set_tests_properties(
            spin-faster-moment-problems
            PROPERTIES LABELS "spin-it;spin-it-faster;regression"
        )

        add_executable(
            quadric-coefficients-regression
            tests/quadric_coefficients_regression.cpp
        )
        target_compile_features(
            quadric-coefficients-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            quadric-coefficients-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            quadric-coefficients-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-quadric-coefficients
            COMMAND "$<TARGET_FILE:quadric-coefficients-regression>"
        )
        set_tests_properties(
            spin-faster-quadric-coefficients
            PROPERTIES LABELS "spin-it;spin-it-faster;regression"
        )

        add_executable(
            spin-boolean-regression
            tests/spin_boolean_regression.cpp
        )
        target_compile_features(spin-boolean-regression PRIVATE cxx_std_17)
        set_target_properties(
            spin-boolean-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-boolean-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-boolean
            COMMAND "$<TARGET_FILE:spin-boolean-regression>"
        )
        set_tests_properties(
            spin-faster-boolean
            PROPERTIES LABELS "spin-it;spin-it-faster;geometry;regression"
        )

        add_executable(
            quadric-triangulator-regression
            tests/quadric_triangulator_regression.cpp
        )
        target_compile_features(
            quadric-triangulator-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            quadric-triangulator-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            quadric-triangulator-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-quadric-triangulator
            COMMAND "$<TARGET_FILE:quadric-triangulator-regression>"
        )
        set_tests_properties(
            spin-faster-quadric-triangulator
            PROPERTIES LABELS "spin-it;spin-it-faster;geometry;regression"
        )

        add_executable(
            quadric-moment-evaluator-regression
            tests/quadric_moment_evaluator_regression.cpp
        )
        target_compile_features(
            quadric-moment-evaluator-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            quadric-moment-evaluator-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            quadric-moment-evaluator-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-quadric-moments
            COMMAND "$<TARGET_FILE:quadric-moment-evaluator-regression>"
        )
        set_tests_properties(
            spin-faster-quadric-moments
            PROPERTIES LABELS "spin-it;spin-it-faster;geometry;regression"
        )

        add_executable(
            quadric-hollowing-solver-regression
            tests/quadric_hollowing_solver_regression.cpp
        )
        target_compile_features(
            quadric-hollowing-solver-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            quadric-hollowing-solver-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            quadric-hollowing-solver-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-newton
            COMMAND "$<TARGET_FILE:quadric-hollowing-solver-regression>"
        )
        set_tests_properties(
            spin-faster-newton
            PROPERTIES LABELS "spin-it;spin-it-faster;solver;regression"
        )

        add_executable(
            spin-it-faster-end-to-end
            tests/spin_it_faster_end_to_end.cpp
        )
        target_compile_features(spin-it-faster-end-to-end PRIVATE cxx_std_17)
        set_target_properties(
            spin-it-faster-end-to-end
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-faster-end-to-end
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-end-to-end
            COMMAND "$<TARGET_FILE:spin-it-faster-end-to-end>"
        )
        set_tests_properties(
            spin-faster-end-to-end
            PROPERTIES LABELS "spin-it;spin-it-faster;solver;geometry;regression"
        )

        add_executable(
            spin-it-faster-interface-types
            tests/spin_it_faster_interface_types.cpp
        )
        target_compile_features(
            spin-it-faster-interface-types PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-faster-interface-types PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            spin-it-faster-interface-types PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-interface-types
            COMMAND "$<TARGET_FILE:spin-it-faster-interface-types>"
        )
        set_tests_properties(
            spin-faster-interface-types
            PROPERTIES LABELS "spin-it;spin-it-faster;geometry;regression"
        )
    endif()

    add_test(
        NAME reachability-dependency-contract
        COMMAND "${CMAKE_COMMAND}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/assert_reachability_dependency.cmake"
    )

    add_test(
        NAME spin-faster-dependency-contract
        COMMAND "${CMAKE_COMMAND}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/assert_spin_it_faster_dependency.cmake"
    )
    add_test(
        NAME spin-2014-dependency-contract
        COMMAND "${CMAKE_COMMAND}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/assert_spin_it_2014_dependency.cmake"
    )
    set_tests_properties(
        spin-faster-dependency-contract
        PROPERTIES
            LABELS "spin-it;spin-it-faster;build-system"
            TIMEOUT 180
    )
    set_tests_properties(
        spin-2014-dependency-contract
        PROPERTIES
            LABELS "spin-it;spin-it-2014;build-system"
            TIMEOUT 240
    )
    set_tests_properties(
        reachability-dependency-contract
        PROPERTIES
            LABELS "spin-it;reachability;build-system"
            TIMEOUT 60
    )

    if(MAKE_IT_STAND_BUILD_REACHABILITY)
        add_executable(
            carving-reachability-regression
            tests/carving_reachability_regression.cpp
        )
        target_compile_features(
            carving-reachability-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            carving-reachability-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            carving-reachability-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME reachability-input
            COMMAND "$<TARGET_FILE:carving-reachability-regression>" input
        )
        add_test(
            NAME reachability-already-safe
            COMMAND "$<TARGET_FILE:carving-reachability-regression>"
                already-safe
        )
        add_test(
            NAME reachability-feasible
            COMMAND "$<TARGET_FILE:carving-reachability-regression>" feasible
        )
        add_test(
            NAME reachability-infeasible
            COMMAND "$<TARGET_FILE:carving-reachability-regression>"
                infeasible
        )
        add_test(
            NAME reachability-certificate
            COMMAND "$<TARGET_FILE:carving-reachability-regression>"
                certificate
        )
        add_test(
            NAME reachability-exhaustive
            COMMAND "$<TARGET_FILE:carving-reachability-regression>"
                exhaustive
        )
        add_test(
            NAME reachability-sparse-certificate
            COMMAND "$<TARGET_FILE:carving-reachability-regression>"
                sparse-certificate
        )
        set_tests_properties(
            reachability-input
            reachability-already-safe
            reachability-feasible
            reachability-infeasible
            reachability-certificate
            reachability-exhaustive
            reachability-sparse-certificate
            PROPERTIES
                LABELS "spin-it;reachability;regression"
                TIMEOUT 15
        )

        add_executable(
            balance-model-reachability-regression
            tests/balance_model_reachability_regression.cpp
        )
        target_compile_features(
            balance-model-reachability-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            balance-model-reachability-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            balance-model-reachability-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME balance-model-reachability
            COMMAND "$<TARGET_FILE:balance-model-reachability-regression>"
        )
        set_tests_properties(
            balance-model-reachability
            PROPERTIES
                LABELS "spin-it;reachability;api;regression"
                TIMEOUT 120
        )
    endif()

    add_executable(
        result-view-regression
        tests/result_view_regression.cpp
    )
    target_compile_features(result-view-regression PRIVATE cxx_std_17)
    set_target_properties(
        result-view-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        result-view-regression
        PRIVATE make_it_stand::core
    )
    add_test(
        NAME result-view-preserves-sealed-posture
        COMMAND "$<TARGET_FILE:result-view-regression>"
    )
    set_tests_properties(
        result-view-preserves-sealed-posture
        PROPERTIES
            LABELS "viewer;result;regression"
            TIMEOUT 15
    )

    add_executable(
        released-beta-archive-regression
        tests/released_beta_archive_regression.cpp
    )
    target_compile_features(
        released-beta-archive-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        released-beta-archive-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        released-beta-archive-regression
        PRIVATE make_it_stand::core
    )
    add_test(
        NAME released-beta-archive-contract
        COMMAND "$<TARGET_FILE:released-beta-archive-regression>"
    )
    set_tests_properties(
        released-beta-archive-contract
        PROPERTIES
            LABELS "paper;released-beta;archive;regression"
            TIMEOUT 60
    )

    add_test(
        NAME ctest-registration-contract
        COMMAND "${CMAKE_COMMAND}"
            "-DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
            "-DCTEST_COMMAND=${CMAKE_CTEST_COMMAND}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/assert_ctest_commands.cmake"
    )
    set_tests_properties(
        ctest-registration-contract
        PROPERTIES
            LABELS "build-system;registration"
            TIMEOUT 30
    )

    if(MAKE_IT_STAND_ENABLE_PAPER_PYTHON_TESTS)
        find_package(Python3 3.10 COMPONENTS Interpreter REQUIRED)
        execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c "import numpy"
            RESULT_VARIABLE MAKE_IT_STAND_NUMPY_IMPORT_STATUS
            ERROR_VARIABLE MAKE_IT_STAND_NUMPY_IMPORT_ERROR
        )
        if(NOT MAKE_IT_STAND_NUMPY_IMPORT_STATUS EQUAL 0)
            message(FATAL_ERROR
                "Paper Python tests require numpy in ${Python3_EXECUTABLE}: "
                "${MAKE_IT_STAND_NUMPY_IMPORT_ERROR}"
            )
        endif()
        add_test(
            NAME paper-python-contracts
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_official_beta.ReleasedBetaParserTests
                test_paper_compat.PaperTransformTests
                test_paper_evidence.PaperEvidenceTests
                test_paper_trace.PaperTraceTests
                test_experiment_drivers.ExperimentDriverTests
                test_stability_oracle.StabilityOracleTests
                test_spin_it_paper_inventory.SpinItPaperInventoryTests
                test_spin_it_faster_paper_assets.SpinItFasterPaperAssetTests
                test_spin_it_faster_table1_runner.SpinItFasterTable1RunnerTests
                test_spin_it_faster_paper_figures.SpinItFasterPaperFigureTests
                test_spin_it_2014_paper_assets.SpinIt2014PaperAssetTests
                test_spin_it_2014_paper_runner.SpinIt2014PaperRunnerTests
                test_spin_it_complete_release_gate.SpinItCompleteReleaseGateTests
                -v
        )
        set_tests_properties(
            paper-python-contracts
            PROPERTIES
                LABELS "paper;python;released-beta;compatibility"
                TIMEOUT 60
                ENVIRONMENT
                    "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
        )

        add_test(
            NAME spin-it-strict-audit-contract
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_spin_it_strict_audit
                -v
        )
        set_tests_properties(
            spin-it-strict-audit-contract
            PROPERTIES
                LABELS "paper;python;spin-it;strict"
                TIMEOUT 30
                ENVIRONMENT
                    "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
        )

        add_test(
            NAME independent-mass-moments-contract
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_independent_mass_moments
                -v
        )
        set_tests_properties(
            independent-mass-moments-contract
            PROPERTIES
                LABELS "paper;python;spin-it;verification"
                TIMEOUT 30
                ENVIRONMENT
                    "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
        )

        add_test(
            NAME run-evidence-v2-contract
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_run_evidence_v2
                -v
        )
        set_tests_properties(
            run-evidence-v2-contract
            PROPERTIES
                LABELS "paper;python;spin-it;verification"
                TIMEOUT 30
                ENVIRONMENT
                    "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
        )

        add_test(
            NAME compiler-observer-contract
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_compiler_observer.CompilerObserverTests -v
        )
        add_test(
            NAME compiler-observer-cmake-fixture
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_compiler_observer.CompilerObserverCMakeTests -v
        )
        set_tests_properties(
            compiler-observer-contract compiler-observer-cmake-fixture
            PROPERTIES
                LABELS "paper;python;spin-it;verification;compiler-observation"
                RESOURCE_LOCK "compiler-observation-fixtures"
                TIMEOUT 60
                ENVIRONMENT
                    "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
        )

        if(MAKE_IT_STAND_BUILD_PYTHON)
            if(MAKE_IT_STAND_BUILD_SPIN_IT_2014)
                add_test(
                    NAME independent-surface-derivatives-contract
                    COMMAND "${Python3_EXECUTABLE}" -m unittest test_surface_derivatives -v
                )
                set_tests_properties(
                    independent-surface-derivatives-contract
                    PROPERTIES
                        LABELS "paper;python;native;spin-it-2014;independent-verification"
                        TIMEOUT 60
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                )
            endif()
            add_test(
                NAME paper-native-pipeline-contract
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_paper_compat.PaperCompatibilityPipelineTests.test_normalized_optimization_restores_and_flattens_output
                test_paper_compat.PaperCompatibilityPipelineTests.test_sealed_artifact_recomputes_moments_after_support_flattening
                test_paper_compat.PaperCompatibilityPipelineTests.test_invalid_flattened_boundary_cannot_feed_pybullet
                test_paper_compat.PaperCompatibilityPipelineTests.test_sealed_artifact_compacts_unreferenced_cavity_grid_vertices
                -v
            )
            set_tests_properties(
                paper-native-pipeline-contract
                PROPERTIES
                    LABELS "paper;python;native;compatibility"
                    TIMEOUT 180
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )

            if(NOT MAKE_IT_STAND_RELEASED_BETA_ROOT)
                set(MAKE_IT_STAND_RELEASED_BETA_ROOT_CANDIDATE
                    "${CMAKE_CURRENT_SOURCE_DIR}/../reference-only/make-it-stand/beta-extracted/make-it-stand_src_beta"
                )
                if(IS_DIRECTORY "${MAKE_IT_STAND_RELEASED_BETA_ROOT_CANDIDATE}")
                    set(MAKE_IT_STAND_RELEASED_BETA_ROOT
                        "${MAKE_IT_STAND_RELEASED_BETA_ROOT_CANDIDATE}"
                    )
                endif()
            endif()
            if(MAKE_IT_STAND_RELEASED_BETA_ROOT)
                if(NOT IS_DIRECTORY "${MAKE_IT_STAND_RELEASED_BETA_ROOT}")
                    message(FATAL_ERROR
                        "MAKE_IT_STAND_RELEASED_BETA_ROOT is not a directory: "
                        "${MAKE_IT_STAND_RELEASED_BETA_ROOT}"
                    )
                endif()
                add_test(
                    NAME released-beta-native-archive-contracts
                    COMMAND "${Python3_EXECUTABLE}" -m unittest
                        test_official_beta.OfficialArchiveGateTests.test_spheres_cache_contract_and_explicit_gate
                        test_official_beta.OfficialArchiveGateTests.test_horse_true_crossing_remains_fail_closed
                        test_official_beta.OfficialArchiveGateTests.test_horse_occupancy_reconstruction_is_explicit_and_buildable
                        test_official_beta.OfficialArchiveGateTests.test_official_complete_archive_build_matrix
                        -v
                )
                set_tests_properties(
                    released-beta-native-archive-contracts
                    PROPERTIES
                        LABELS "paper;python;native;released-beta;official-assets"
                        TIMEOUT 600
                        ENVIRONMENT
                            "${MAKE_IT_STAND_PYTHON_TEST_ENV};MAKE_IT_STAND_BETA_ROOT=${MAKE_IT_STAND_RELEASED_BETA_ROOT}"
                )
            else()
                message(STATUS
                    "Released-beta archive not found; official asset tests are not registered"
                )
            endif()
        endif()
    endif()

    if(MAKE_IT_STAND_ENABLE_PYBULLET_TESTS)
        find_package(Python3 3.10 COMPONENTS Interpreter REQUIRED)
        execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c "import pybullet, trimesh"
            RESULT_VARIABLE MAKE_IT_STAND_PYBULLET_IMPORT_STATUS
            ERROR_VARIABLE MAKE_IT_STAND_PYBULLET_IMPORT_ERROR
        )
        if(NOT MAKE_IT_STAND_PYBULLET_IMPORT_STATUS EQUAL 0)
            message(FATAL_ERROR
                "PyBullet tests require pybullet and trimesh: "
                "${MAKE_IT_STAND_PYBULLET_IMPORT_ERROR}"
            )
        endif()
        add_test(
            NAME pybullet-material-contracts
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_pybullet_runner.PyBulletMaterialContractTests
                -v
        )
        set_tests_properties(
            pybullet-material-contracts
            PROPERTIES
                LABELS "simulation;python;pybullet;material-contract"
                TIMEOUT 60
                ENVIRONMENT
                    "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
        )
    endif()

    add_executable(
        voxel-domain-regression
        tests/voxel_domain_regression.cpp
    )
    target_compile_features(voxel-domain-regression PRIVATE cxx_std_17)
    set_target_properties(
        voxel-domain-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        voxel-domain-regression
        PRIVATE make_it_stand::core
    )
    foreach(voxel_domain_scenario IN ITEMS domain boundary invalid model)
        add_test(
            NAME "voxel-domain-${voxel_domain_scenario}"
            COMMAND "$<TARGET_FILE:voxel-domain-regression>"
                "${voxel_domain_scenario}"
        )
        set_tests_properties(
            "voxel-domain-${voxel_domain_scenario}"
            PROPERTIES
                LABELS "voxel-domain;occupancy;scaling"
                TIMEOUT 60
        )
    endforeach()

    add_executable(
        balance-model-regression
        tests/balance_model_regression.cpp
    )
    target_compile_features(balance-model-regression PRIVATE cxx_std_17)
    set_target_properties(
        balance-model-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        balance-model-regression
        PRIVATE make_it_stand::core
    )
    add_test(
        NAME r9-balance-model-facade
        COMMAND "$<TARGET_FILE:balance-model-regression>"
    )
    set_tests_properties(
        r9-balance-model-facade
        PROPERTIES
            LABELS "r9;api;regression"
            TIMEOUT 120
    )

    if(MAKE_IT_STAND_BUILD_PYTHON)
        if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
            add_test(
                NAME spin-faster-python-api
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                    discover
                    -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                    -p "test_spin_it_faster.py"
                    -v
            )
            set_tests_properties(
                spin-faster-python-api
                PROPERTIES
                    LABELS "spin-it;spin-it-faster;python;bindings;regression"
                    TIMEOUT 180
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
            add_test(
                NAME spin-faster-table1-miniature
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                    test_spin_it_faster_table1_runner.SpinItFasterTable1NativeRunnerTests
                    -v
            )
            set_tests_properties(
                spin-faster-table1-miniature
                PROPERTIES
                    LABELS "spin-it;spin-it-faster;python;paper;regression"
                    TIMEOUT 180
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
            add_test(
                NAME spin-faster-paper-figures-miniature
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                    test_spin_it_faster_paper_figures.SpinItFasterPaperFigureNativeTests
                    -v
            )
            set_tests_properties(
                spin-faster-paper-figures-miniature
                PROPERTIES
                    LABELS "spin-it;spin-it-faster;python;paper;regression"
                    TIMEOUT 240
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
        endif()
        if(MAKE_IT_STAND_BUILD_REACHABILITY)
            add_test(
                NAME reachability-python-bindings
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                    discover
                    -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                    -p "test_reachability_bindings.py"
                    -v
            )
            set_tests_properties(
                reachability-python-bindings
                PROPERTIES
                    LABELS "spin-it;reachability;python;bindings"
                    TIMEOUT 120
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
        endif()
        add_test(
            NAME r9-python-bindings
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                discover
                -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                -p "test_bindings.py"
                -v
        )
        add_test(
            NAME r9-python-tools
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                discover
                -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                -p "test_python_tools.py"
                -v
        )
        add_test(
            NAME r10-python-evidence
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                discover
                -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                -p "test_r10*.py"
                -v
        )
        set_tests_properties(
            r9-python-bindings
            r9-python-tools
            PROPERTIES
                LABELS "r9;python"
                TIMEOUT 180
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
        )
        set_tests_properties(
            r10-python-evidence
            PROPERTIES
                LABELS "r10;python;evidence"
                TIMEOUT 180
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
        )
    endif()

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
        raw-mass-moments-regression
        tests/raw_mass_moments_regression.cpp
    )
    target_compile_features(
        raw-mass-moments-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        raw-mass-moments-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        raw-mass-moments-regression
        PRIVATE make_it_stand::core
    )

    function(add_raw_mass_moments_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:raw-mass-moments-regression>"
                "${scenario}"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "spin-it;physics;moments"
                TIMEOUT 30
        )
    endfunction()

    add_raw_mass_moments_test(raw-moments-ordering ordering)
    add_raw_mass_moments_test(raw-moments-conversion conversion)
    add_raw_mass_moments_test(raw-moments-transform transform)

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
        surface-mass-triangle-chain-analytic
        triangle-chain-gradient
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
        deformed-voxel-mass-field-regression
        tests/deformed_voxel_mass_field_regression.cpp
    )
    target_compile_features(
        deformed-voxel-mass-field-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        deformed-voxel-mass-field-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        deformed-voxel-mass-field-regression
        PRIVATE make_it_stand::core
    )

    function(add_deformed_voxel_mass_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:deformed-voxel-mass-field-regression>"
                "${scenario}"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r7;deformation;physics;carving"
                TIMEOUT 30
        )
    endfunction()

    add_deformed_voxel_mass_test(
        deformed-voxel-mass-identity-similarity
        identity-similarity
    )
    add_deformed_voxel_mass_test(
        deformed-voxel-mass-cavity-sum
        distorted-cavity-sum
    )
    add_deformed_voxel_mass_test(
        deformed-voxel-mass-invalid
        invalid
    )

    add_executable(
        deformed-carving-solver-regression
        tests/deformed_carving_solver_regression.cpp
    )
    target_compile_features(
        deformed-carving-solver-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        deformed-carving-solver-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        deformed-carving-solver-regression
        PRIVATE make_it_stand::core
    )

    function(add_deformed_carving_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:deformed-carving-solver-regression>"
                "${scenario}"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r7;alternating;carving;transaction"
                TIMEOUT 30
        )
    endfunction()

    add_deformed_carving_test(deformed-carving-identity identity)
    add_deformed_carving_test(deformed-carving-mass-field deformed)
    add_deformed_carving_test(deformed-carving-transaction transaction)
    add_deformed_carving_test(
        deformed-carving-objective-update
        objective-update
    )

    add_executable(
        balance-target-schedule-regression
        tests/balance_target_schedule_regression.cpp
    )
    target_compile_features(
        balance-target-schedule-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        balance-target-schedule-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        balance-target-schedule-regression
        PRIVATE make_it_stand::core
    )

    function(add_balance_target_schedule_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:balance-target-schedule-regression>"
                "${scenario}"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r7;alternating;stability;multi-target"
                TIMEOUT 30
        )
    endfunction()

    add_balance_target_schedule_test(
        alternating-target-refresh
        refresh
    )
    add_balance_target_schedule_test(
        alternating-stability-terminal
        stability
    )
    add_balance_target_schedule_test(
        alternating-multi-target-terminal
        multi
    )
    add_balance_target_schedule_test(
        released-beta-frozen-target
        frozen
    )

    add_executable(
        alternating-block-semantics-regression
        tests/alternating_block_semantics_regression.cpp
    )
    target_compile_features(
        alternating-block-semantics-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        alternating-block-semantics-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        alternating-block-semantics-regression
        PRIVATE make_it_stand::core
    )

    function(add_alternating_block_semantics_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:alternating-block-semantics-regression>"
                "${scenario}"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r7;alternating;transaction"
                TIMEOUT 30
        )
    endfunction()

    add_alternating_block_semantics_test(
        alternating-refill-reset
        refill
    )
    add_alternating_block_semantics_test(
        alternating-fixed-alpha-transaction
        paper-transaction
    )
    add_alternating_block_semantics_test(
        alternating-paper-beta-profile
        profile-discriminator
    )

    add_executable(
        alternating-balance-optimizer-regression
        tests/alternating_balance_optimizer_regression.cpp
    )
    target_compile_features(
        alternating-balance-optimizer-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        alternating-balance-optimizer-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        alternating-balance-optimizer-regression
        PRIVATE make_it_stand::core
    )

    function(add_alternating_balance_optimizer_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:alternating-balance-optimizer-regression>"
                "${scenario}"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r7;alternating;integration;physics"
                TIMEOUT 60
        )
    endfunction()

    add_alternating_balance_optimizer_test(
        alternating-controller-refill-block
        refill-block
    )
    add_alternating_balance_optimizer_test(
        alternating-controller-transaction
        transaction
    )
    add_alternating_balance_optimizer_test(
        alternating-controller-end-to-end
        e2e
    )
    add_alternating_balance_optimizer_test(
        released-beta-outer-loop
        beta-outer
    )
    add_alternating_balance_optimizer_test(
        alternating-invalid-outer-boundary-gate
        invalid-outer
    )

    add_executable(
        r10-ablation-regression
        tests/r10_ablation_regression.cpp
    )
    target_compile_features(r10-ablation-regression PRIVATE cxx_std_17)
    set_target_properties(
        r10-ablation-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        r10-ablation-regression
        PRIVATE make_it_stand::core
    )
    add_test(
        NAME r10-ablation-switches
        COMMAND "$<TARGET_FILE:r10-ablation-regression>"
    )
    set_tests_properties(
        r10-ablation-switches
        PROPERTIES
            LABELS "r10;ablation;alternating;integration;physics"
            TIMEOUT 120
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
    add_test(
        NAME released-beta-unit-cube-voxel-grid
        COMMAND "$<TARGET_FILE:voxel-constraints-regression>" unit-cube
    )
    add_test(
        NAME released-beta-archive-voxel-cache
        COMMAND "$<TARGET_FILE:voxel-constraints-regression>" beta-cache
    )
    set_tests_properties(
        released-beta-unit-cube-voxel-grid
        released-beta-archive-voxel-cache
        PROPERTIES
            LABELS "r4;carving;manufacturing;paper-compat"
            TIMEOUT 30
    )

    set_tests_properties(
        minimum-wall-thickness-analytic
        voxel-topology-report
        PROPERTIES
            LABELS "r4;carving;manufacturing"
            TIMEOUT 30
    )

    add_executable(
        manufacturing-geometry-regression
        tests/manufacturing_geometry_regression.cpp
    )
    target_compile_features(
        manufacturing-geometry-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        manufacturing-geometry-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        manufacturing-geometry-regression
        PRIVATE make_it_stand::core
    )

    function(add_manufacturing_geometry_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:manufacturing-geometry-regression>"
                "${scenario}"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r8;manufacturing;geometry"
                TIMEOUT 30
        )
    endfunction()

    add_manufacturing_geometry_test(r8-boundary-cube cube)
    add_manufacturing_geometry_test(r8-boundary-missing-face missing)
    add_manufacturing_geometry_test(r8-boundary-reversed reversed)
    add_manufacturing_geometry_test(r8-boundary-duplicate duplicate)
    add_manufacturing_geometry_test(r8-boundary-degenerate degenerate)
    add_manufacturing_geometry_test(r8-boundary-double-outer double-outer)
    add_manufacturing_geometry_test(
        r8-boundary-self-intersection
        self-intersection
    )
    add_manufacturing_geometry_test(
        r8-boundary-rounded-sphere-adjacency
        rounded-sphere-adjacency
    )
    add_manufacturing_geometry_test(
        r8-boundary-coplanar-shared-vertex
        coplanar-shared-vertex
    )
    add_manufacturing_geometry_test(r8-prepare-empty prepare-empty)
    add_manufacturing_geometry_test(r8-prepare-cavity prepare-cavity)
    add_manufacturing_geometry_test(r8-prepare-stale-inner stale-inner)
    add_manufacturing_geometry_test(r8-prepare-cell-collapse cell-collapse)
    add_manufacturing_geometry_test(r8-prepare-cell-swap cell-swap)
    add_manufacturing_geometry_test(
        r8-prepare-material-components
        material-components
    )
    add_manufacturing_geometry_test(r8-prepare-containment containment)
    add_manufacturing_geometry_test(r8-prepare-wall-threshold wall-threshold)
    add_manufacturing_geometry_test(r8-source-bound source-bound)
    add_manufacturing_geometry_test(
        r8-source-bound-cavity
        source-bound-cavity
    )
    add_manufacturing_geometry_test(
        r8-source-bound-deformation-terminal
        source-bound-deformation-terminal
    )
    add_manufacturing_geometry_test(
        r10-source-bound-no-carving-terminal
        source-bound-no-carving-terminal
    )
    set_tests_properties(
        r10-source-bound-no-carving-terminal
        PROPERTIES LABELS "r8;r10;manufacturing;geometry;ablation"
    )
    add_manufacturing_geometry_test(
        r8-source-bound-carving-after-deformation
        source-bound-carving-after-deformation
    )
    add_manufacturing_geometry_test(r8-export-roundtrip export-roundtrip)
    add_manufacturing_geometry_test(r8-export-transaction export-transaction)
    add_manufacturing_geometry_test(r8-export-seal export-seal)
    add_manufacturing_geometry_test(r8-source-rejection source-rejection)
    add_manufacturing_geometry_test(r8-raw-substitution raw-substitution)
    add_manufacturing_geometry_test(r8-bound-substitution bound-substitution)
    add_manufacturing_geometry_test(r8-scale-translation scale-translation)

    add_executable(r8-fake-slicer tests/fake_slicer.cpp)
    target_compile_features(r8-fake-slicer PRIVATE cxx_std_17)
    set_target_properties(
        r8-fake-slicer
        PROPERTIES CXX_EXTENSIONS OFF
    )

    add_executable(
        slicer-validation-regression
        tests/slicer_validation_regression.cpp
    )
    target_compile_features(
        slicer-validation-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        slicer-validation-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        slicer-validation-regression
        PRIVATE make_it_stand::core
    )

    function(add_slicer_validation_test test_name scenario)
        add_test(
            NAME "${test_name}"
            COMMAND "$<TARGET_FILE:slicer-validation-regression>"
                "${scenario}"
                "$<TARGET_FILE:r8-fake-slicer>"
        )
        set_tests_properties(
            "${test_name}"
            PROPERTIES
                LABELS "r8;manufacturing;slicer-unit"
                TIMEOUT 15
        )
    endfunction()

    add_test(
        NAME r8-slicer-not-configured
        COMMAND "$<TARGET_FILE:slicer-validation-regression>"
            not-configured
    )
    set_tests_properties(
        r8-slicer-not-configured
        PROPERTIES
            LABELS "r8;manufacturing;slicer-unit"
            TIMEOUT 15
    )
    add_slicer_validation_test(
        r8-slicer-executable-hash-mismatch
        executable-hash-mismatch
    )
    add_slicer_validation_test(
        r8-slicer-profile-hash-mismatch
        profile-hash-mismatch
    )
    add_slicer_validation_test(r8-slicer-version-mismatch version-mismatch)
    add_slicer_validation_test(
        r8-slicer-version-truncation
        version-truncation
    )
    add_slicer_validation_test(r8-slicer-success success)
    add_slicer_validation_test(
        r8-slicer-stable-invocation
        stable-invocation
    )
    add_slicer_validation_test(r8-slicer-repair-rejected repair)
    add_slicer_validation_test(r8-slicer-invalid-toolpath invalid-toolpath)
    add_slicer_validation_test(r8-slicer-nonzero nonzero)
    add_slicer_validation_test(r8-slicer-signal signal)
    add_slicer_validation_test(r8-slicer-timeout timeout)
    add_slicer_validation_test(r8-slicer-shell-safe shell-safe)
    add_slicer_validation_test(r8-slicer-report-json report-json)
    add_slicer_validation_test(r8-slicer-paths paths)
    add_slicer_validation_test(
        r8-slicer-sigchld-ownership
        sigchld-ownership
    )
    add_slicer_validation_test(
        r8-slicer-closed-stdio-launch-failure
        closed-stdio-launch-failure
    )
    add_slicer_validation_test(r8-slicer-capture-limit capture-limit)

    add_executable(
        partition-manifest-regression
        tests/partition_manifest_regression.cpp
    )
    target_compile_features(
        partition-manifest-regression
        PRIVATE cxx_std_17
    )
    set_target_properties(
        partition-manifest-regression
        PROPERTIES CXX_EXTENSIONS OFF
    )
    target_link_libraries(
        partition-manifest-regression
        PRIVATE make_it_stand::core
    )
    add_test(
        NAME r8-partition-manifest
        COMMAND "$<TARGET_FILE:partition-manifest-regression>"
    )
    set_tests_properties(
        r8-partition-manifest
        PROPERTIES
            LABELS "r8;manufacturing;partition-interface"
            TIMEOUT 15
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
