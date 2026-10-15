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
include(CompilerRuntime)

option(MAKE_IT_STAND_BUILD_VIEWER "Build the interactive GLFW viewer" ON)
option(MAKE_IT_STAND_BUILD_PYTHON "Build the pybind11 Python extension" OFF)
set(
    MAKE_IT_STAND_SEALED_SPIN_EXTENSION
    ""
    CACHE FILEPATH
    "Optional exact prebuilt _make_it_stand extension observed by a sealed artifact"
)
option(MAKE_IT_STAND_TEST_USER_SYSTEMD_JOBS
    "Run real Linux user-systemd/cgroup enforcement tests (no backend skip/fallback)" OFF)
option(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS
    "Build the common pinned-CGAL exact represented-boundary auditor" OFF)
option(MAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL
    "Use certified interval predicates with exact Gmpq fallback in the boundary auditor" OFF)
option(MAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE
    "Use exact positive-cone tests for indexed shared-vertex triangle contacts" OFF)
option(MAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS
    "Use conservative exact-projection interval slabs before exact triangle predicates" OFF)
if(MAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS AND NOT MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
    message(FATAL_ERROR "MAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS requires exact geometry checks")
endif()
if(MAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL AND NOT MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
    message(FATAL_ERROR
        "MAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL requires MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS")
endif()
if(MAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE AND NOT MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
    message(FATAL_ERROR
        "MAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE requires MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS")
endif()
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

# Bind Faster evaluator evidence to the actual libigl checkout used by this
# build.  The tree digest includes tracked diffs and the status text includes
# untracked paths; a dirty checkout is therefore never mislabeled as its
# upstream commit alone.
find_package(Git QUIET)
set(MAKE_IT_STAND_LIBIGL_SOURCE_ID "unversioned")
if(GIT_FOUND)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" -C "${MAKE_IT_STAND_LIBIGL_DIR}" rev-parse HEAD
        RESULT_VARIABLE _libigl_revision_result
        OUTPUT_VARIABLE _libigl_revision
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(_libigl_revision_result EQUAL 0 AND _libigl_revision MATCHES "^[0-9a-fA-F]+$")
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" -C "${MAKE_IT_STAND_LIBIGL_DIR}" diff HEAD --binary
            RESULT_VARIABLE _libigl_diff_result
            OUTPUT_VARIABLE _libigl_diff
            ERROR_QUIET
        )
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" -C "${MAKE_IT_STAND_LIBIGL_DIR}" status --porcelain --untracked-files=normal
            RESULT_VARIABLE _libigl_status_result
            OUTPUT_VARIABLE _libigl_status
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(_libigl_diff_result EQUAL 0 AND _libigl_status_result EQUAL 0)
            string(SHA256 _libigl_tree_digest "${_libigl_diff}\n${_libigl_status}")
            set(MAKE_IT_STAND_LIBIGL_SOURCE_ID "revision=${_libigl_revision},tree-sha256=${_libigl_tree_digest}")
        endif()
    endif()
endif()

set(LIBIGL_DIR "${MAKE_IT_STAND_LIBIGL_DIR}" CACHE PATH "libigl root" FORCE)
option(LIBIGL_WITH_OPENGL "Use OpenGL" ${MAKE_IT_STAND_BUILD_VIEWER})
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW" ${MAKE_IT_STAND_BUILD_VIEWER})
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui" ${MAKE_IT_STAND_BUILD_VIEWER})
if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER OR MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
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

if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER OR MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
    make_it_stand_verify_pinned_geometry_dependencies()
    include(RankTwoNativeFactor)
endif()

if(MAKE_IT_STAND_BUILD_REACHABILITY OR MAKE_IT_STAND_BUILD_SPIN_IT_2014)
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
        include/spin/MomentInertiaAlgebra.h
    )
endif()

if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
    list(APPEND MAKE_IT_STAND_CORE_SOURCES
        src/spin/HydrostaticGeometry.cpp
        src/spin/FasterTheoremConditions.cpp
        src/spin/PlanarDiskIntersection.cpp
        src/spin/PlanarFidgetCircleProblem.cpp
        src/spin/PlanarMassMoments.cpp
        src/spin/PlanarSpinMomentProblem.cpp
        src/spin/PlaneMomentEvaluator.cpp
        src/spin/PlaneStaticBalanceSolver.cpp
        src/spin/BuoyantPlaneBalanceSolver.cpp
        src/spin/QuadricCoefficients.cpp
        src/spin/QuadricCoefficientAuthority.cpp
        src/spin/QuadricAuthorityLinearSolve.cpp
        src/spin/QuadricInverseGaussMap.cpp
        src/spin/QuadricMomentEvaluator.cpp
        src/spin/QuadricDiscreteMomentDifferential.cpp
        src/spin/QuadricStratumTangent.cpp
        src/spin/QuadricHollowingSolver.cpp
        src/spin/LegacyPredictionTransaction.cpp
        src/spin/DualMonomialCoordinates.cpp
        src/spin/QuadricTriangulator.cpp
        src/spin/PolyhedralMomentDifferential.cpp
        src/spin/TruncationRootCertificate.cpp
        src/spin/RankTwoSource.cpp
        src/spin/RankTwoAlgebraic.cpp
        src/spin/RankTwoArrangement.cpp
        src/spin/RankTwoNormalMesh.cpp
        src/spin/RankTwoSettings.cpp
        src/spin/RankTwoCarrier.cpp
        src/spin/RankTwoArcSampling.cpp
        src/spin/RankTwoArcStrip.cpp
        src/spin/RankTwoEndpointArc.cpp
        src/spin/RankTwoEndpointFan.cpp
        src/spin/RankTwoEndpointNeighborhood.cpp
        src/spin/RankTwoArcRegistry.cpp
        src/spin/RankTwoFinitePartition.cpp
        src/spin/RankTwoFiniteArcCoverage.cpp
        src/spin/RankTwoCoveredFiniteMesh.cpp
        src/spin/RankTwoFiniteMeshScope.cpp
        src/spin/RankTwoMovingInterior.cpp
        src/spin/RankTwoInteriorInventory.cpp
        src/spin/RankTwoWholeInteriorInventory.cpp
        src/spin/RankTwoAutomaticMaterialMesh.cpp
        src/spin/RankTwoInteriorMesh.cpp
        src/spin/RankTwoMovingMeshEdge.cpp
        src/spin/RankTwoMaterialLifts.cpp
        src/spin/RankTwoStrictMaterialMesh.cpp
        src/spin/RankTwoFiniteTriangulation.cpp
        src/spin/RankTwoRegistryEdge.cpp
        src/spin/RankTwoArcRefinement.cpp
        src/spin/RankTwoChordIncidence.cpp
        src/spin/RankTwoRegistryRotation.cpp
        src/spin/QuadricSourceMotion.cpp
    )
    list(APPEND MAKE_IT_STAND_CORE_HEADERS
        include/spin/HydrostaticGeometry.h
        include/spin/FasterTheoremConditions.h
        include/spin/PlanarDiskIntersection.h
        include/spin/PlanarFidgetCircleProblem.h
        include/spin/PlanarMassMoments.h
        include/spin/PlanarSpinMomentProblem.h
        include/spin/PlaneMomentEvaluator.h
        include/spin/PlaneStaticBalanceSolver.h
        include/spin/BuoyantPlaneBalanceSolver.h
        include/spin/QuadricCoefficients.h
        include/spin/QuadricCoefficientAuthority.h
        include/spin/QuadricAuthorityLinearSolve.h
        include/spin/QuadricInverseGaussMap.h
        include/spin/QuadricMomentEvaluator.h
        include/spin/QuadricDiscreteMomentDifferential.h
        include/spin/QuadricStratumTangent.h
        include/spin/QuadricHollowingSolver.h
        include/spin/QuadricTriangulator.h
        include/spin/PolyhedralMomentDifferential.h
        include/spin/TruncationRootCertificate.h
        include/spin/RankTwoSource.h
        include/spin/QuadricSourceMotion.h
    )
endif()

if(MAKE_IT_STAND_BUILD_SPIN_IT_2014)
    list(APPEND MAKE_IT_STAND_CORE_SOURCES
        src/spin2014/AdaptiveVoxelDomain.cpp
        src/spin2014/AdaptiveFillFunctional.cpp
        src/spin2014/SpinIt2014MomentProblem.cpp
        src/spin2014/QuadricMomentProblem.cpp
        src/spin2014/Spin2014GivensProblem.cpp
        src/spin2014/SlqpLbfgsSolver.cpp
        src/spin2014/SlqpKkt.cpp
        src/spin2014/SlqpAuditLp.cpp
        src/spin2014/SlqpAuditTrace.cpp
        src/spin2014/SlqpIterationStore.cpp
        src/spin2014/SlqpAuditJson.cpp
        src/spin2014/SlqpHistory.cpp
        src/spin2014/SlqpTangentQp.cpp
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

if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER AND MAKE_IT_STAND_BUILD_SPIN_IT_2014)
    list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/QuadricAuthorityMomentProblem.cpp)
    list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/QuadricAuthorityMomentProblem.h)
    list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/QuadricAuthorityNewtonState.h)
endif()

if(MAKE_IT_STAND_BUILD_REACHABILITY)
    list(APPEND MAKE_IT_STAND_CORE_SOURCES
        src/reachability/CarvingReachability.cpp
    )
    list(APPEND MAKE_IT_STAND_CORE_HEADERS
        include/reachability/CarvingReachability.h
    )
endif()

if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
    list(APPEND MAKE_IT_STAND_CORE_SOURCES src/verification/ExactBoundaryAudit.cpp src/verification/ExactMaterialRegion.cpp src/verification/ExactBoundaryDistance.cpp src/verification/ExactTriangleSelfIntersection.cpp)
    if(MAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS)
        set_property(SOURCE src/verification/ExactBoundaryAudit.cpp APPEND PROPERTY
            COMPILE_DEFINITIONS MAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS=1)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/verification/ConservativeDiagonalSlabs.h)
    endif()
    if(MAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL)
        set_property(SOURCE src/verification/ExactBoundaryAudit.cpp APPEND PROPERTY
            COMPILE_DEFINITIONS MAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL=1)
    endif()
    if(MAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE)
        set_property(SOURCE src/verification/ExactBoundaryAudit.cpp APPEND PROPERTY
            COMPILE_DEFINITIONS MAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE=1)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/verification/ExactSharedVertexIntersection.h)
    endif()
    list(APPEND MAKE_IT_STAND_CORE_HEADERS include/verification/ExactBoundaryAudit.h include/verification/ExactMaterialRegion.h include/verification/ExactCuboidUnion.h include/verification/ExactNefBudget.h include/verification/ExactTriangleSelfIntersection.h)
    if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/verification/ExactLibiglBoolean.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/verification/ExactLibiglBoolean.h)
    endif()
    if(MAKE_IT_STAND_BUILD_SPIN_IT_2014)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceSurfaceOctree.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceSurfaceOctree.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceSurfaceCage.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceSurfaceCage.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageGrid.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageGrid.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageFem.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageFem.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageBbw.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageBbw.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageBinding.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageBinding.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceConformingBoundary.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceConformingBoundary.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageMomentFunctional.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageMomentFunctional.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCagePaperFunctional.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCagePaperFunctional.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageSlqp.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageSlqp.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageGeometry.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageGeometry.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageAdaptiveRebind.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageAdaptiveRebind.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageAdaptiveSolver.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageAdaptiveSolver.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageBinaryMaterial.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageBinaryMaterial.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageMaterialRegions.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageMaterialRegions.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceCageComponentClearance.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceCageComponentClearance.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/LinearizedRotationLaplacian.cpp src/spin2014/SourceCageShapeFunctional.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/LinearizedRotationLaplacian.h include/spin2014/SourceCageShapeFunctional.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceAdaptiveDomain.cpp src/spin2014/SourceBinaryMoments.cpp src/spin2014/SourceAdaptiveSpinSolver.cpp src/spin2014/SourceCheckpointCodec.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceAdaptiveDomain.h include/spin2014/SourceBinaryMoments.h include/spin2014/SourceAdaptiveSpinSolver.h)
        list(APPEND MAKE_IT_STAND_CORE_SOURCES src/spin2014/SourceMaterialBoundary.cpp)
        list(APPEND MAKE_IT_STAND_CORE_HEADERS include/spin2014/SourceMaterialBoundary.h)
    endif()
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
        MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS=$<BOOL:${MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS}>
)
target_compile_definitions(
    make_it_stand_core
    PRIVATE
        MAKE_IT_STAND_LIBIGL_SOURCE_ID="${MAKE_IT_STAND_LIBIGL_SOURCE_ID}"
)
target_link_libraries(make_it_stand_core PUBLIC igl::core)
if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER OR MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
    target_link_libraries(make_it_stand_core PUBLIC igl::cgal)
endif()
if(MAKE_IT_STAND_BUILD_REACHABILITY OR MAKE_IT_STAND_BUILD_SPIN_IT_2014)
    target_link_libraries(make_it_stand_core PRIVATE highs::highs)
endif()

# Keep the canonicalized discrete Boolean/moment derivative as an explicit
# Faster solver profile.  Its implementation units are compiled with private
# symbol remapping so the established legacy evaluator and public entrypoints
# remain available unchanged in make_it_stand_core.
if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
    add_library(spin_faster_canonical_material_impl OBJECT
        src/spin/FullRankMomentEvaluator.cpp
    )
    target_link_libraries(spin_faster_canonical_material_impl
        PRIVATE make_it_stand::core
    )
    target_compile_definitions(spin_faster_canonical_material_impl PRIVATE
        evaluateFullRankMaterialMoments=canonicalFullRankMaterialImplementation
        differentiatePolyhedralMaterialMoments=differentiateCanonicalPolyhedralMaterialMoments
    )

    add_library(spin_faster_canonical_solver_impl OBJECT
        src/spin/FullRankHollowingSolver.cpp
    )
    target_link_libraries(spin_faster_canonical_solver_impl
        PRIVATE make_it_stand::core
    )
    target_compile_definitions(spin_faster_canonical_solver_impl PRIVATE
        integrateFullRankBoundaryMoments=integrateCanonicalFullRankBoundaryMoments
        evaluateFullRankMaterialMoments=evaluateCanonicalFullRankMaterialMoments
        makeFullRankMomentEvaluator=canonicalFullRankAdapterImplementation
        solveFullRankQuadricHollowing=canonicalFullRankSolverImplementation
    )

    add_library(make_it_stand_faster_fullrank STATIC
        src/spin/CanonicalPolyhedralMomentDifferential.cpp
        src/spin/CanonicalFullRankHollowing.cpp
        $<TARGET_OBJECTS:spin_faster_canonical_material_impl>
        $<TARGET_OBJECTS:spin_faster_canonical_solver_impl>
    )
    add_library(make_it_stand::faster_fullrank ALIAS make_it_stand_faster_fullrank)
    target_compile_features(make_it_stand_faster_fullrank PUBLIC cxx_std_17)
    set_target_properties(make_it_stand_faster_fullrank PROPERTIES
        CXX_EXTENSIONS OFF
        POSITION_INDEPENDENT_CODE ON
    )
    target_include_directories(make_it_stand_faster_fullrank PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    )
    target_link_libraries(make_it_stand_faster_fullrank
        PUBLIC make_it_stand::core
    )
endif()

if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER AND MAKE_IT_STAND_BUILD_SPIN_IT_2014)
    add_library(make_it_stand_spin2014_quadric STATIC
        src/spin2014/QuadricHollowing.cpp src/spin2014/QuadricStaticBalance.cpp)
    add_library(make_it_stand::spin2014_quadric ALIAS make_it_stand_spin2014_quadric)
    target_compile_features(make_it_stand_spin2014_quadric PUBLIC cxx_std_17)
    set_target_properties(make_it_stand_spin2014_quadric PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(make_it_stand_spin2014_quadric PUBLIC make_it_stand::faster_fullrank)
endif()

if(MAKE_IT_STAND_BUILD_PYTHON)
    set(PYBIND11_FINDPYTHON ON)
    set(MAKE_IT_STAND_PYTHON_COMPONENTS Interpreter Development.Module)
    if(MAKE_IT_STAND_BUILD_VIEWER AND MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
        list(APPEND MAKE_IT_STAND_PYTHON_COMPONENTS Development.Embed)
    endif()
    find_package(
        Python 3.10
        COMPONENTS ${MAKE_IT_STAND_PYTHON_COMPONENTS}
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
    if(MAKE_IT_STAND_BUILD_SPIN_IT_2014 AND MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
        target_sources(_make_it_stand PRIVATE bindings/SourceSpinMaterialBindings.cpp)
        if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
            target_sources(_make_it_stand PRIVATE bindings/Spin2014AuthorityBindings.cpp)
        endif()
    endif()
    set_target_properties(
        _make_it_stand
        PROPERTIES
            CXX_EXTENSIONS OFF
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/python"
    )
    target_link_libraries(_make_it_stand PRIVATE make_it_stand::core)
    if(MAKE_IT_STAND_BUILD_SPIN_IT_FASTER)
        target_link_libraries(_make_it_stand PRIVATE make_it_stand::faster_fullrank)
        if(MAKE_IT_STAND_BUILD_SPIN_IT_2014)
            target_link_libraries(_make_it_stand PRIVATE make_it_stand::spin2014_quadric)
        endif()
    endif()

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
    if(MAKE_IT_STAND_BUILD_PYTHON AND MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
        target_sources(make_it_stand_viewer PRIVATE src/viewer/SealedSpinView.cpp)
        target_link_libraries(make_it_stand_viewer PRIVATE pybind11::embed)
        set(MAKE_IT_STAND_SEALED_SPIN_EXTENSION_FOR_VIEWER "$<TARGET_FILE:_make_it_stand>")
        if(MAKE_IT_STAND_SEALED_SPIN_EXTENSION)
            get_filename_component(
                MAKE_IT_STAND_SEALED_SPIN_EXTENSION_RESOLVED
                "${MAKE_IT_STAND_SEALED_SPIN_EXTENSION}"
                REALPATH
            )
            if(NOT EXISTS "${MAKE_IT_STAND_SEALED_SPIN_EXTENSION_RESOLVED}")
                message(FATAL_ERROR
                    "MAKE_IT_STAND_SEALED_SPIN_EXTENSION does not name an existing native module: "
                    "${MAKE_IT_STAND_SEALED_SPIN_EXTENSION}"
                )
            endif()
            set(MAKE_IT_STAND_SEALED_SPIN_EXTENSION_FOR_VIEWER
                "${MAKE_IT_STAND_SEALED_SPIN_EXTENSION_RESOLVED}")
        else()
            add_dependencies(make_it_stand_viewer _make_it_stand)
        endif()
        target_compile_definitions(make_it_stand_viewer PRIVATE
            MAKE_IT_STAND_SEALED_SPIN_VIEW=1
            IGL_VIEWER_VIEWER_QUIET=1
            SEALED_SPIN_PYTHON="${MAKE_IT_STAND_PYTHON_EXECUTABLE_REAL}"
            SEALED_SPIN_SOURCE_PYTHON="${CMAKE_CURRENT_SOURCE_DIR}/python"
            SEALED_SPIN_EXTENSION="${MAKE_IT_STAND_SEALED_SPIN_EXTENSION_FOR_VIEWER}")
        if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND EXISTS "${MAKE_IT_STAND_COMPILER_LIBSTDCXX}")
            get_filename_component(MAKE_IT_STAND_VIEWER_CXX_RUNTIME "${MAKE_IT_STAND_COMPILER_LIBSTDCXX}" REALPATH)
            get_filename_component(MAKE_IT_STAND_VIEWER_RUNTIME_DIR "${MAKE_IT_STAND_VIEWER_CXX_RUNTIME}" DIRECTORY)
            # Embed's Python library directory may also contain an older C++
            # runtime. The freshly compiled native executable needs its compiler's.
            set_property(TARGET make_it_stand_viewer PROPERTY BUILD_RPATH "${MAKE_IT_STAND_VIEWER_RUNTIME_DIR}")
        endif()
    endif()
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

if(MAKE_IT_STAND_TEST_USER_SYSTEMD_JOBS AND
   NOT (BUILD_TESTING AND MAKE_IT_STAND_BUILD_SPIN_IT_2014 AND
        MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS AND MAKE_IT_STAND_BUILD_PYTHON))
    message(FATAL_ERROR "User-systemd job tests require BUILD_TESTING, SPIN_IT_2014, EXACT_GEOMETRY_CHECKS and PYTHON; refusing an unregistered requested test group")
endif()

if(MAKE_IT_STAND_SEALED_SPIN_EXTENSION AND CMAKE_VERSION VERSION_LESS "3.22")
    message(FATAL_ERROR
        "Using an externally sealed _make_it_stand extension with CTest requires CMake 3.22 or newer "
        "for platform-correct PYTHONPATH modification of the native viewer test."
    )
endif()

if(BUILD_TESTING)
    if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
        add_executable(exact-triangle-self-intersection-regression
            tests/exact_triangle_self_intersection_regression.cpp)
        target_link_libraries(exact-triangle-self-intersection-regression
            PRIVATE make_it_stand::core)
        add_test(NAME exact-triangle-self-intersection-regression
            COMMAND "$<TARGET_FILE:exact-triangle-self-intersection-regression>")
        set_tests_properties(exact-triangle-self-intersection-regression PROPERTIES
            LABELS "verification;exact-geometry;intersection;native" TIMEOUT 120)
        add_executable(exact-partition-clearance-regression tests/exact_partition_clearance_regression.cpp)
        target_link_libraries(exact-partition-clearance-regression PRIVATE make_it_stand::core)
        add_test(NAME exact-partition-clearance-native COMMAND exact-partition-clearance-regression)
        set_tests_properties(exact-partition-clearance-native PROPERTIES LABELS "verification;exact-geometry;native" TIMEOUT 120)
        add_executable(exact-material-region-regression tests/exact_material_region_regression.cpp)
        target_link_libraries(exact-material-region-regression PRIVATE make_it_stand::core)
        add_test(NAME exact-material-region-native COMMAND exact-material-region-regression)
        set_tests_properties(exact-material-region-native PROPERTIES LABELS "verification;exact-geometry;native" TIMEOUT 120)
        add_executable(exact-boundary-audit-regression tests/exact_boundary_audit_regression.cpp)
        target_link_libraries(exact-boundary-audit-regression PRIVATE make_it_stand::core)
        add_test(NAME exact-boundary-audit-native COMMAND exact-boundary-audit-regression)
        set_tests_properties(exact-boundary-audit-native PROPERTIES LABELS "verification;exact-geometry;native" TIMEOUT 60)
        add_test(NAME exact-geometry-dependency-contract
            COMMAND "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/assert_exact_geometry_dependency.cmake")
        set_tests_properties(exact-geometry-dependency-contract PROPERTIES
            LABELS "verification;exact-geometry;cmake" TIMEOUT 120)
        if(MAKE_IT_STAND_BUILD_PYTHON)
            include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/SpinGeometryDiagnostics.cmake")
            make_it_stand_register_geometry_diagnostics(
                "${CMAKE_CURRENT_SOURCE_DIR}" "${Python3_EXECUTABLE}")
            add_test(NAME spin-fixed-boundary-shell-python
                COMMAND "${Python3_EXECUTABLE}" -B
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_required_geometry_diagnostic.py"
                    "${CMAKE_CURRENT_SOURCE_DIR}" test_fixed_boundary_shell)
            set_tests_properties(spin-fixed-boundary-shell-python PROPERTIES
                LABELS "verification;spin-it-2014;python;native;synthetic-regression"
                TIMEOUT 300
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
            add_test(NAME spin-fixed-pivot-critical-oracle-python
                COMMAND "${Python3_EXECUTABLE}" -B
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_required_geometry_diagnostic.py"
                    "${CMAKE_CURRENT_SOURCE_DIR}" test_fixed_pivot_critical_oracle)
            set_tests_properties(spin-fixed-pivot-critical-oracle-python PROPERTIES
                LABELS "verification;spin-it-2014;python;independent-dynamics;synthetic-regression"
                TIMEOUT 120
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
            add_test(NAME spin-torque-free-analytic-python
                COMMAND "${Python3_EXECUTABLE}" -B
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_required_geometry_diagnostic.py"
                    "${CMAKE_CURRENT_SOURCE_DIR}" test_torque_free_analytic)
            set_tests_properties(spin-torque-free-analytic-python PROPERTIES
                LABELS "verification;spin-it-2014;python;independent-dynamics;synthetic-regression"
                TIMEOUT 120
                ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
            add_test(NAME spin-rotational-reference-frames-python
                COMMAND "${Python3_EXECUTABLE}" -B
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_required_geometry_diagnostic.py"
                    "${CMAKE_CURRENT_SOURCE_DIR}" test_rotational_reference_frames)
            set_tests_properties(spin-rotational-reference-frames-python PROPERTIES
                LABELS "verification;spin-it-2014;python;independent-dynamics;synthetic-regression"
                TIMEOUT 120
                ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
            add_test(NAME reconstructed-cuboid-target-frame-python
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_reconstructed_cuboid_target_frame -v)
            set_tests_properties(reconstructed-cuboid-target-frame-python PROPERTIES
                LABELS "verification;spin-it-2014;python;target-frame;synthetic-regression"
                TIMEOUT 60
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
            add_test(NAME spin-dynamics-material-admission-python
                COMMAND "${Python3_EXECUTABLE}" -B
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_required_geometry_diagnostic.py"
                    "${CMAKE_CURRENT_SOURCE_DIR}" test_g7_dynamics_admission)
            set_tests_properties(spin-dynamics-material-admission-python PROPERTIES
                LABELS "verification;spin-it-2014;python;native;independent-dynamics;synthetic-regression"
                TIMEOUT 300
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
            add_test(NAME sealed-spin-view-python
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_sealed_spin_view.SealedSpinViewTests -v)
            set_tests_properties(sealed-spin-view-python PROPERTIES
                LABELS "verification;exact-geometry;python;native;viewer" TIMEOUT 300
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            if(MAKE_IT_STAND_BUILD_VIEWER)
                add_test(NAME sealed-spin-view-native
                    COMMAND "${Python3_EXECUTABLE}" -m unittest test_sealed_spin_view.NativeSealedSpinViewTests -v)
                set_tests_properties(sealed-spin-view-native PROPERTIES
                    LABELS "verification;exact-geometry;native;viewer" TIMEOUT 300
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};MAKE_IT_STAND_SEALED_VIEWER=$<TARGET_FILE:make_it_stand_viewer>;MAKE_IT_STAND_SEALED_VIEWER_EXTENSION=${MAKE_IT_STAND_SEALED_SPIN_EXTENSION_FOR_VIEWER}")
                if(MAKE_IT_STAND_SEALED_SPIN_EXTENSION)
                    get_filename_component(
                        MAKE_IT_STAND_SEALED_SPIN_EXTENSION_TEST_DIR
                        "${MAKE_IT_STAND_SEALED_SPIN_EXTENSION_RESOLVED}"
                        DIRECTORY
                    )
                    set_tests_properties(sealed-spin-view-native PROPERTIES
                        ENVIRONMENT_MODIFICATION
                        "PYTHONPATH=path_list_prepend:${MAKE_IT_STAND_SEALED_SPIN_EXTENSION_TEST_DIR}")
                endif()
            endif()
            add_test(NAME sealed-spin-artifact-python
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_sealed_spin_artifact test_spin_bundle -v)
            set_tests_properties(sealed-spin-artifact-python PROPERTIES
                LABELS "verification;exact-geometry;python;native" TIMEOUT 300
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME exact-partition-clearance-python
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_exact_partition_clearance -v)
            add_test(NAME derived-spin-material-python
                COMMAND "${Python3_EXECUTABLE}" -B
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_required_geometry_diagnostic.py"
                    "${CMAKE_CURRENT_SOURCE_DIR}" "test_derived_spin_material")
            set_tests_properties(derived-spin-material-python PROPERTIES
                LABELS "verification;exact-geometry;python;native" TIMEOUT 180
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME coplanar-boundary-python
                COMMAND "${Python3_EXECUTABLE}" -B
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_required_geometry_diagnostic.py"
                    "${CMAKE_CURRENT_SOURCE_DIR}" "test_coplanar_boundary")
            set_tests_properties(coplanar-boundary-python PROPERTIES
                LABELS "verification;exact-geometry;python" TIMEOUT 180
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME derived-spin-material-runner-python
                COMMAND "${Python3_EXECUTABLE}" -B
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_required_geometry_diagnostic.py"
                    "${CMAKE_CURRENT_SOURCE_DIR}" "test_derived_spin_material_runner")
            set_tests_properties(derived-spin-material-runner-python PROPERTIES
                LABELS "verification;exact-geometry;python;native" TIMEOUT 180
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            set_tests_properties(exact-partition-clearance-python PROPERTIES
                LABELS "verification;exact-geometry;python;native" TIMEOUT 180
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME exact-inter-cavity-clearance-python
                COMMAND "${Python3_EXECUTABLE}" -B
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_required_geometry_diagnostic.py"
                    "${CMAKE_CURRENT_SOURCE_DIR}" "test_inter_cavity_clearance")
            set_tests_properties(exact-inter-cavity-clearance-python PROPERTIES
                LABELS "verification;exact-geometry;python;native;spin2014" TIMEOUT 180
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;MKL_NUM_THREADS=1")
            add_test(NAME spin-it-2014-discrete-cut-plan-python
                COMMAND "${Python3_EXECUTABLE}" -B
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/run_required_geometry_diagnostic.py"
                    "${CMAKE_CURRENT_SOURCE_DIR}" "test_digital_cut_plan")
            set_tests_properties(spin-it-2014-discrete-cut-plan-python PROPERTIES
                LABELS "verification;python;spin2014;digital-cut;synthetic-regression" TIMEOUT 60
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1")
            add_test(NAME exact-material-region-python
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_exact_material_region -v)
            set_tests_properties(exact-material-region-python PROPERTIES
                LABELS "verification;exact-geometry;python;native" TIMEOUT 180
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME exact-boundary-audit-python
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_exact_boundary_audit -v)
            set_tests_properties(exact-boundary-audit-python PROPERTIES
                LABELS "verification;exact-geometry;python;native" TIMEOUT 60
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
        endif()
    endif()
    if(MAKE_IT_STAND_BUILD_SPIN_IT_2014)
        add_executable(spin-it-2014-online-trace-regression tests/spin_it_2014_online_trace_regression.cpp)
        target_link_libraries(spin-it-2014-online-trace-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-online-trace COMMAND spin-it-2014-online-trace-regression)
        if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
            add_executable(spin-it-2014-source-bbw-joint-regression tests/spin_it_2014_source_bbw_joint_regression.cpp)
            add_executable(spin-it-2014-source-conforming-boundary-regression tests/spin_it_2014_source_conforming_boundary_regression.cpp)
            add_executable(spin-it-2014-source-cage-moments-regression tests/spin_it_2014_source_cage_moments_regression.cpp)
            add_executable(spin-it-2014-source-cage-paper-regression tests/spin_it_2014_source_cage_paper_regression.cpp)
            add_executable(spin-it-2014-source-cage-shape-regression tests/spin_it_2014_source_cage_shape_regression.cpp)
            add_executable(spin-it-2014-source-cage-slqp-regression tests/spin_it_2014_source_cage_slqp_regression.cpp)
            add_executable(spin-it-2014-source-cage-geometry-regression tests/spin_it_2014_source_cage_geometry_regression.cpp)
            add_executable(spin-it-2014-source-cage-guarded-slqp-regression tests/spin_it_2014_source_cage_guarded_slqp_regression.cpp)
            add_executable(spin-it-2014-source-cage-asymmetric-joint-regression tests/spin_it_2014_source_cage_asymmetric_joint_regression.cpp)
            add_executable(spin-it-2014-source-cage-adaptive-rebind-regression tests/spin_it_2014_source_cage_adaptive_rebind_regression.cpp)
            add_executable(spin-it-2014-source-cage-adaptive-solver-regression tests/spin_it_2014_source_cage_adaptive_solver_regression.cpp)
            add_executable(spin-it-2014-source-cage-binary-material-regression tests/spin_it_2014_source_cage_binary_material_regression.cpp)
            add_executable(spin-it-2014-source-cage-material-regions-regression tests/spin_it_2014_source_cage_material_regions_regression.cpp)
            add_executable(spin-it-2014-source-cage-component-clearance-regression tests/spin_it_2014_source_cage_component_clearance_regression.cpp)
            target_link_libraries(spin-it-2014-source-cage-component-clearance-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-component-clearance COMMAND spin-it-2014-source-cage-component-clearance-regression)
            set_tests_properties(spin-it-2014-source-cage-component-clearance PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            target_link_libraries(spin-it-2014-source-cage-material-regions-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-material-regions COMMAND spin-it-2014-source-cage-material-regions-regression)
            set_tests_properties(spin-it-2014-source-cage-material-regions PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            target_link_libraries(spin-it-2014-source-cage-binary-material-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-binary-material COMMAND spin-it-2014-source-cage-binary-material-regression)
            set_tests_properties(spin-it-2014-source-cage-binary-material PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            add_test(NAME spin-it-2014-source-cage-joint-index-guard
                COMMAND spin-it-2014-source-cage-binary-material-regression --joint-index-guard)
            set_tests_properties(spin-it-2014-source-cage-joint-index-guard PROPERTIES
                LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 30)
            target_link_libraries(spin-it-2014-source-cage-adaptive-solver-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-adaptive-solver COMMAND spin-it-2014-source-cage-adaptive-solver-regression)
            set_tests_properties(spin-it-2014-source-cage-adaptive-solver PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            target_link_libraries(spin-it-2014-source-cage-adaptive-rebind-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-adaptive-rebind COMMAND spin-it-2014-source-cage-adaptive-rebind-regression)
            set_tests_properties(spin-it-2014-source-cage-adaptive-rebind PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            target_link_libraries(spin-it-2014-source-cage-guarded-slqp-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-guarded-slqp COMMAND spin-it-2014-source-cage-guarded-slqp-regression)
            set_tests_properties(spin-it-2014-source-cage-guarded-slqp PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            add_test(NAME spin-it-2014-source-cage-guarded-slqp-allocation COMMAND spin-it-2014-source-cage-guarded-slqp-regression --allocation)
            set_tests_properties(spin-it-2014-source-cage-guarded-slqp-allocation PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            add_test(NAME spin-it-2014-source-cage-guarded-slqp-move COMMAND spin-it-2014-source-cage-guarded-slqp-regression --move)
            set_tests_properties(spin-it-2014-source-cage-guarded-slqp-move PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            target_link_libraries(spin-it-2014-source-cage-asymmetric-joint-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-asymmetric-joint COMMAND spin-it-2014-source-cage-asymmetric-joint-regression)
            set_tests_properties(spin-it-2014-source-cage-asymmetric-joint PROPERTIES
                LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            if(MAKE_IT_STAND_BUILD_PYTHON AND Python3_EXECUTABLE)
                add_test(NAME spin-it-2014-source-cage-joint-batch-protocol
                    COMMAND "${Python3_EXECUTABLE}" -m unittest discover
                        -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p test_source_cage_joint_batch.py -v)
                set_tests_properties(spin-it-2014-source-cage-joint-batch-protocol PROPERTIES
                    LABELS "spin-it-2014;exact-geometry;python;synthetic-regression"
                    TIMEOUT 180 RUN_SERIAL TRUE
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};SPIN_JOINT_BATCH_BINARY=$<TARGET_FILE:spin-it-2014-source-cage-asymmetric-joint-regression>")
                add_test(NAME spin-it-2014-dual-material-geometry-export
                    COMMAND "${Python3_EXECUTABLE}" -m unittest discover
                        -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p test_dual_material_geometry_export.py -v)
                set_tests_properties(spin-it-2014-dual-material-geometry-export PROPERTIES
                    LABELS "spin-it-2014;exact-geometry;python;synthetic-regression"
                    TIMEOUT 180 RUN_SERIAL TRUE
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};SPIN_DUAL_EXPORT_BINARY=$<TARGET_FILE:spin-it-2014-dual-material-regression>")
                add_test(NAME spin-it-2014-rounding-moment-breakdown
                    COMMAND "${Python3_EXECUTABLE}" -m unittest discover
                        -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p test_rounding_moment_breakdown.py -v)
                set_tests_properties(spin-it-2014-rounding-moment-breakdown PROPERTIES
                    LABELS "spin-it-2014;python;verification" TIMEOUT 30
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            endif()
            target_link_libraries(spin-it-2014-source-cage-geometry-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-geometry COMMAND spin-it-2014-source-cage-geometry-regression)
            set_tests_properties(spin-it-2014-source-cage-geometry PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            target_link_libraries(spin-it-2014-source-cage-slqp-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-slqp COMMAND spin-it-2014-source-cage-slqp-regression)
            set_tests_properties(spin-it-2014-source-cage-slqp PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            target_link_libraries(spin-it-2014-source-cage-shape-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-shape COMMAND spin-it-2014-source-cage-shape-regression)
            set_tests_properties(spin-it-2014-source-cage-shape PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            target_link_libraries(spin-it-2014-source-cage-paper-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-paper COMMAND spin-it-2014-source-cage-paper-regression)
            set_tests_properties(spin-it-2014-source-cage-paper PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600 RUN_SERIAL TRUE)
            add_executable(spin-it-2014-source-cage-paper-allocation-probe tests/spin_it_2014_source_cage_paper_regression.cpp)
            target_compile_definitions(spin-it-2014-source-cage-paper-allocation-probe PRIVATE SPIN2014_PAPER_ALLOCATION_PROBE=1)
            target_link_libraries(spin-it-2014-source-cage-paper-allocation-probe PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-paper-allocation COMMAND spin-it-2014-source-cage-paper-allocation-probe)
            set_tests_properties(spin-it-2014-source-cage-paper-allocation PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 60)
            target_link_libraries(spin-it-2014-source-cage-moments-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-moments COMMAND spin-it-2014-source-cage-moments-regression)
            set_tests_properties(spin-it-2014-source-cage-moments PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600)
            target_link_libraries(spin-it-2014-source-conforming-boundary-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-conforming-boundary COMMAND spin-it-2014-source-conforming-boundary-regression)
            set_tests_properties(spin-it-2014-source-conforming-boundary PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600)
            add_executable(spin-it-2014-source-cage-binding-regression tests/spin_it_2014_source_cage_binding_regression.cpp)
            target_link_libraries(spin-it-2014-source-cage-binding-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-binding COMMAND spin-it-2014-source-cage-binding-regression)
            set_tests_properties(spin-it-2014-source-cage-binding PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600)
            target_link_libraries(spin-it-2014-source-bbw-joint-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-bbw-joint COMMAND spin-it-2014-source-bbw-joint-regression)
            set_tests_properties(spin-it-2014-source-bbw-joint PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600)
            add_executable(spin-it-2014-source-cage-fem-regression tests/spin_it_2014_source_cage_fem_regression.cpp)
            add_executable(spin-it-2014-source-cage-bbw-regression tests/spin_it_2014_source_cage_bbw_regression.cpp)
            target_link_libraries(spin-it-2014-source-cage-bbw-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-bbw COMMAND spin-it-2014-source-cage-bbw-regression)
            set_tests_properties(spin-it-2014-source-cage-bbw PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600)
            target_link_libraries(spin-it-2014-source-cage-fem-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-fem COMMAND spin-it-2014-source-cage-fem-regression)
            set_tests_properties(spin-it-2014-source-cage-fem PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600)
            add_executable(spin-it-2014-source-cage-grid-regression tests/spin_it_2014_source_cage_grid_regression.cpp)
            target_link_libraries(spin-it-2014-source-cage-grid-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-cage-grid COMMAND spin-it-2014-source-cage-grid-regression)
            set_tests_properties(spin-it-2014-source-cage-grid PROPERTIES LABELS "spin-it-2014;exact-geometry;native;synthetic-regression" TIMEOUT 3600)
            add_executable(spin-it-2014-source-surface-cage-regression tests/spin_it_2014_source_surface_cage_regression.cpp)
            target_link_libraries(spin-it-2014-source-surface-cage-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-surface-cage COMMAND spin-it-2014-source-surface-cage-regression)
            set_tests_properties(spin-it-2014-source-surface-cage PROPERTIES LABELS "spin-it-2014;exact-geometry;native" TIMEOUT 3600)
            if(MAKE_IT_STAND_BUILD_PYTHON)
                add_test(NAME source-spin-checkpoint-store-python
                    COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_source_checkpoint_store.py")
                set_tests_properties(source-spin-checkpoint-store-python PROPERTIES
                    LABELS "spin-it-2014;python;checkpoint-foundation" TIMEOUT 120
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTEST_DISABLE_PLUGIN_AUTOLOAD=1")
                add_test(NAME source-spin-run-supervisor-python
                    COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_source_run_supervisor.py")
                set_tests_properties(source-spin-run-supervisor-python PROPERTIES
                    LABELS "spin-it-2014;python;checkpoint-foundation;synthetic-regression" TIMEOUT 120
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTEST_DISABLE_PLUGIN_AUTOLOAD=1")
                if(MAKE_IT_STAND_TEST_USER_SYSTEMD_JOBS)
                    add_test(NAME source-spin-whole-job-python
                        COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                            "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_source_job.py"
                            "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_supervised_teapot_command.py")
                    set_tests_properties(source-spin-whole-job-python PROPERTIES
                        LABELS "spin-it-2014;python;linux-user-systemd;job-resource-gate;synthetic-regression" TIMEOUT 120
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTEST_DISABLE_PLUGIN_AUTOLOAD=1")
                endif()
                add_test(NAME source-spin-online-trace-python
                    COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_source_online_trace.py")
                set_tests_properties(source-spin-online-trace-python PROPERTIES
                    LABELS "spin-it-2014;exact-geometry;python;native;trace-evidence" TIMEOUT 900
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTEST_DISABLE_PLUGIN_AUTOLOAD=1")
                add_test(NAME source-spin-material-python
                    COMMAND "${Python3_EXECUTABLE}" -m unittest test_source_spin_material_bridge test_source_paper_termination -v)
                set_tests_properties(source-spin-material-python PROPERTIES
                    LABELS "spin-it-2014;exact-geometry;python;native" TIMEOUT 900
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
                add_test(NAME source-spin-resumed-material-python
                    COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_source_resumed_material.py")
                set_tests_properties(source-spin-resumed-material-python PROPERTIES
                    LABELS "spin-it-2014;exact-geometry;python;native;checkpoint-foundation"
                    TIMEOUT 900 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTEST_DISABLE_PLUGIN_AUTOLOAD=1")
                add_test(NAME source-spin-trace-stream-python
                    COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_source_trace_stream.py"
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_source_trace_archive.py")
                set_tests_properties(source-spin-trace-stream-python PROPERTIES
                    LABELS "spin-it-2014;exact-geometry;python;native;trace-evidence" TIMEOUT 900
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
                add_test(NAME spin-it-2014-tangent-policy-python
                    COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_native_tangent_policy.py"
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_slqp_direction_replay.py")
                set_tests_properties(spin-it-2014-tangent-policy-python PROPERTIES
                    LABELS "spin-it-2014;python;native;trace-evidence" TIMEOUT 900
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
                add_test(NAME spin-it-2014-finite-box-policy-python
                    COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_native_finite_box_policy.py"
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_finite_box_direction_replay.py")
                set_tests_properties(spin-it-2014-finite-box-policy-python PROPERTIES
                    LABELS "spin-it-2014;python;native;trace-evidence" TIMEOUT 900
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
                add_test(NAME spin-it-2014-primal-box-policy-python
                    COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_native_primal_box_policy.py")
                set_tests_properties(spin-it-2014-primal-box-policy-python PROPERTIES
                    LABELS "spin-it-2014;python;native;trace-evidence" TIMEOUT 900
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            endif()
            add_executable(spin-it-2014-source-material-regression tests/spin_it_2014_source_material_regression.cpp)
            target_link_libraries(spin-it-2014-source-material-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-material COMMAND spin-it-2014-source-material-regression)
            set_tests_properties(spin-it-2014-source-material PROPERTIES LABELS "spin-it-2014;exact-geometry;native" TIMEOUT 3600)
            add_executable(spin-it-2014-source-adaptive-regression tests/spin_it_2014_source_adaptive_regression.cpp)
            add_executable(spin-it-2014-partition-restore tests/spin_it_2014_partition_restore.cpp)
            target_link_libraries(spin-it-2014-partition-restore PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-partition-restore COMMAND spin-it-2014-partition-restore)
            set_tests_properties(spin-it-2014-partition-restore PROPERTIES LABELS "spin-it-2014;exact-geometry;native;checkpoint-foundation" TIMEOUT 120)
            if(MAKE_IT_STAND_BUILD_PYTHON)
                add_executable(spin-it-2014-checkpoint-codec tests/spin_it_2014_checkpoint_codec.cpp)
                target_link_libraries(spin-it-2014-checkpoint-codec PRIVATE make_it_stand::core)
                add_test(NAME spin-it-2014-checkpoint-codec
                    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_source_checkpoint_codec.py"
                        $<TARGET_FILE:spin-it-2014-checkpoint-codec>)
                set_tests_properties(spin-it-2014-checkpoint-codec PROPERTIES
                    LABELS "spin-it-2014;native;checkpoint-foundation;synthetic-regression" TIMEOUT 300
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
                add_executable(spin-it-2014-postsolve-replay tests/spin_it_2014_postsolve_replay.cpp)
                target_link_libraries(spin-it-2014-postsolve-replay PRIVATE make_it_stand::core)
                add_test(NAME spin-it-2014-postsolve-replay
                    COMMAND spin-it-2014-postsolve-replay "${CMAKE_CURRENT_BINARY_DIR}/postsolve-replay")
                set_tests_properties(spin-it-2014-postsolve-replay PROPERTIES
                    LABELS "spin-it-2014;native;postsolve-replay;synthetic-regression" TIMEOUT 180)
                add_test(NAME spin-it-2014-postsolve-replay-python
                    COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                        "${CMAKE_CURRENT_SOURCE_DIR}/experiments/teapot_spin2014/test_postsolve_replay.py")
                set_tests_properties(spin-it-2014-postsolve-replay-python PROPERTIES
                    LABELS "spin-it-2014;python;postsolve-replay;synthetic-regression" TIMEOUT 120
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTEST_DISABLE_PLUGIN_AUTOLOAD=1")
                add_test(NAME source-spin-supervised-native-recovery
                    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_supervised_source_checkpoint.py"
                        $<TARGET_FILE:spin-it-2014-checkpoint-codec>)
                set_tests_properties(source-spin-supervised-native-recovery PROPERTIES
                    LABELS "spin-it-2014;native;checkpoint-foundation;synthetic-regression" TIMEOUT 300
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
                add_executable(spin-it-2014-stage-resume tests/spin_it_2014_stage_resume.cpp)
                target_link_libraries(spin-it-2014-stage-resume PRIVATE make_it_stand::core)
                add_test(NAME spin-it-2014-stage-resume
                    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_stage_resume.py"
                        $<TARGET_FILE:spin-it-2014-stage-resume>)
                set_tests_properties(spin-it-2014-stage-resume PROPERTIES
                    LABELS "spin-it-2014;native;checkpoint-foundation;synthetic-regression" TIMEOUT 240)
                add_executable(spin-it-2014-partition-restore-process tests/spin_it_2014_partition_restore_process.cpp)
                target_link_libraries(spin-it-2014-partition-restore-process PRIVATE make_it_stand::core)
                add_test(NAME spin-it-2014-partition-restore-process COMMAND "${Python3_EXECUTABLE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_partition_restore_process.py"
                    $<TARGET_FILE:spin-it-2014-partition-restore-process>)
                set_tests_properties(spin-it-2014-partition-restore-process PROPERTIES LABELS "spin-it-2014;native;checkpoint-foundation" TIMEOUT 240)
                add_test(NAME source-spin-checkpoint-store-native
                    COMMAND "${Python3_EXECUTABLE}"
                        "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_checkpoint_store_native.py"
                        $<TARGET_FILE:spin-it-2014-partition-restore-process>)
                set_tests_properties(source-spin-checkpoint-store-native PROPERTIES
                    LABELS "spin-it-2014;python;native;checkpoint-foundation;synthetic-regression" TIMEOUT 240
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            endif()
            add_executable(sha256-stream-regression tests/sha256_stream_regression.cpp)
            target_include_directories(sha256-stream-regression PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(sha256-stream-regression PRIVATE cxx_std_17)
            add_test(NAME sha256-stream COMMAND sha256-stream-regression)
            set_tests_properties(sha256-stream PROPERTIES LABELS "spin-it-2014;native;trace-evidence" TIMEOUT 120)
            add_executable(slqp-trace-hash-memory tests/slqp_trace_hash_memory.cpp)
            target_link_libraries(slqp-trace-hash-memory PRIVATE make_it_stand::core)
            add_test(NAME slqp-trace-hash-memory COMMAND slqp-trace-hash-memory)
            set_tests_properties(slqp-trace-hash-memory PROPERTIES LABELS "spin-it-2014;native;trace-evidence;synthetic-regression" TIMEOUT 120)
            if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
                add_test(NAME slqp-trace-hash-memory-capped COMMAND slqp-trace-hash-memory --limited)
                set_tests_properties(slqp-trace-hash-memory-capped PROPERTIES LABELS "spin-it-2014;native;trace-evidence;synthetic-regression" TIMEOUT 120)
            endif()
            target_link_libraries(spin-it-2014-source-adaptive-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-adaptive COMMAND spin-it-2014-source-adaptive-regression)
            add_test(NAME spin-it-2014-source-adaptive-geometry9 COMMAND spin-it-2014-source-adaptive-regression --geometry9)
            set_tests_properties(spin-it-2014-source-adaptive spin-it-2014-source-adaptive-geometry9 PROPERTIES LABELS "spin-it-2014;exact-geometry;native" TIMEOUT 3600)
            add_executable(spin-it-2014-source-surface-octree-regression tests/spin_it_2014_source_surface_octree_regression.cpp)
            target_link_libraries(spin-it-2014-source-surface-octree-regression PRIVATE make_it_stand::core)
            add_test(NAME spin-it-2014-source-surface-octree COMMAND spin-it-2014-source-surface-octree-regression)
            add_test(NAME spin-it-2014-source-surface-octree-level9 COMMAND spin-it-2014-source-surface-octree-regression --level9)
            set_tests_properties(spin-it-2014-source-surface-octree spin-it-2014-source-surface-octree-level9 PROPERTIES LABELS "spin-it-2014;exact-geometry;native" TIMEOUT 3600)
        endif()
        add_executable(spin-it-2014-surface-trace-regression tests/surface_triangle_trace_regression.cpp)
        target_include_directories(spin-it-2014-surface-trace-regression PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
        target_link_libraries(spin-it-2014-surface-trace-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-surface-trace COMMAND spin-it-2014-surface-trace-regression)
        set_tests_properties(spin-it-2014-surface-trace PROPERTIES LABELS "spin-it-2014;native" TIMEOUT 120)
        add_executable(spin-it-2014-surface-value-regression tests/spin_it_2014_surface_value_regression.cpp)
        target_link_libraries(spin-it-2014-surface-value-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-surface-value COMMAND spin-it-2014-surface-value-regression)
        add_test(NAME spin-it-2014-surface-value-allocation COMMAND spin-it-2014-surface-value-regression --allocation)
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
        add_executable(spin-it-2014-audited tests/spin_it_2014_audited.cpp)
        target_link_libraries(spin-it-2014-audited PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-audited COMMAND spin-it-2014-audited)
        add_executable(spin-it-2014-tangent-qp tests/spin_it_2014_tangent_qp.cpp)
        target_link_libraries(spin-it-2014-tangent-qp PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-tangent-qp COMMAND spin-it-2014-tangent-qp)
        add_executable(spin-it-2014-finite-box-qp tests/spin_it_2014_finite_box_qp.cpp)
        add_executable(spin-it-2014-primal-box tests/spin_it_2014_primal_box.cpp)
        target_link_libraries(spin-it-2014-primal-box PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-primal-box COMMAND spin-it-2014-primal-box "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/spin2014_failed_box_108.txt")
        add_executable(spin-it-2014-primal-box-analytic tests/spin_it_2014_finite_box_qp.cpp)
        target_compile_definitions(spin-it-2014-primal-box-analytic PRIVATE TEST_PRIMAL)
        target_link_libraries(spin-it-2014-primal-box-analytic PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-primal-box-analytic COMMAND spin-it-2014-primal-box-analytic)
        add_executable(spin-it-2014-primal-box-integration tests/spin_it_2014_finite_box_integration.cpp)
        target_compile_definitions(spin-it-2014-primal-box-integration PRIVATE TEST_PRIMAL)
        target_link_libraries(spin-it-2014-primal-box-integration PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-primal-box-integration COMMAND spin-it-2014-primal-box-integration)
        add_executable(spin-it-2014-bounded-capacity tests/spin_it_2014_bounded_capacity.cpp)
        add_executable(spin-it-2014-schur-cache tests/spin_it_2014_schur_cache.cpp)
        add_executable(spin-it-2014-schur-fallback tests/spin_it_2014_schur_fallback.cpp)
        target_link_libraries(spin-it-2014-schur-fallback PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-schur-fallback COMMAND spin-it-2014-schur-fallback)
        set_tests_properties(spin-it-2014-schur-fallback PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
        target_link_libraries(spin-it-2014-schur-cache PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-schur-cache COMMAND spin-it-2014-schur-cache)
        set_tests_properties(spin-it-2014-schur-cache PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
        target_link_libraries(spin-it-2014-bounded-capacity PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-bounded-capacity COMMAND spin-it-2014-bounded-capacity)
        set_tests_properties(spin-it-2014-bounded-capacity PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
        # GNU linker wrapping counts real out-of-line metric calls without
        # adding instrumentation to the production solver. Itanium ABI/Linux.
        if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            add_executable(spin-it-2014-column-cache tests/spin_it_2014_column_cache.cpp)
            target_link_libraries(spin-it-2014-column-cache PRIVATE make_it_stand::core)
            target_link_options(spin-it-2014-column-cache PRIVATE
                "LINKER:--wrap=_ZNK8spin201412LbfgsHistory12applyInverseERKN5Eigen6MatrixIdLin1ELi1ELi0ELin1ELi1EEE")
            add_test(NAME spin-it-2014-column-cache COMMAND spin-it-2014-column-cache)
            set_tests_properties(spin-it-2014-column-cache PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
            add_executable(spin-it-2014-schur-allocation tests/spin_it_2014_schur_allocation.cpp)
            target_link_libraries(spin-it-2014-schur-allocation PRIVATE make_it_stand::core)
            target_link_options(spin-it-2014-schur-allocation PRIVATE "LINKER:--wrap=malloc")
            add_test(NAME spin-it-2014-schur-allocation COMMAND spin-it-2014-schur-allocation)
            set_tests_properties(spin-it-2014-schur-allocation PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
            add_executable(spin-it-2014-metric-materialization tests/spin_it_2014_metric_materialization.cpp)
            target_link_libraries(spin-it-2014-metric-materialization PRIVATE make_it_stand::core)
            target_link_options(spin-it-2014-metric-materialization PRIVATE "LINKER:--wrap=malloc")
            add_test(NAME spin-it-2014-metric-materialization COMMAND spin-it-2014-metric-materialization)
            set_tests_properties(spin-it-2014-metric-materialization PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
        endif()
        set_tests_properties(spin-it-2014-primal-box spin-it-2014-primal-box-analytic spin-it-2014-primal-box-integration PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
        add_executable(spin-it-2014-finite-box-integration tests/spin_it_2014_finite_box_integration.cpp)
        target_link_libraries(spin-it-2014-finite-box-integration PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-finite-box-integration COMMAND spin-it-2014-finite-box-integration)
        set_tests_properties(spin-it-2014-finite-box-integration PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
        target_link_libraries(spin-it-2014-finite-box-qp PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-finite-box-qp COMMAND spin-it-2014-finite-box-qp)
        set_tests_properties(spin-it-2014-finite-box-qp PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
        set_tests_properties(spin-it-2014-tangent-qp PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
        add_executable(spin-it-2014-tangent-integration tests/spin_it_2014_tangent_integration.cpp)
        target_link_libraries(spin-it-2014-tangent-integration PRIVATE make_it_stand::core)
        add_test(NAME spin-it-2014-tangent-integration COMMAND spin-it-2014-tangent-integration)
        set_tests_properties(spin-it-2014-tangent-integration PROPERTIES LABELS "spin-it-2014;native;regression" TIMEOUT 120)
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

        add_executable(spin-it-2014-givens-regression tests/spin_it_2014_givens_regression.cpp)
        target_link_libraries(spin-it-2014-givens-regression PRIVATE make_it_stand::core)
        target_compile_features(spin-it-2014-givens-regression PRIVATE cxx_std_17)
        add_test(NAME spin-it-2014-givens COMMAND "$<TARGET_FILE:spin-it-2014-givens-regression>")
        set_tests_properties(spin-it-2014-givens PROPERTIES LABELS "spin-it;spin-it-2014;regression")

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
            spin-it-2014-target-frame-regression
            tests/spin_it_2014_target_frame_regression.cpp
        )
        target_link_libraries(
            spin-it-2014-target-frame-regression PRIVATE make_it_stand::core
        )
        target_compile_features(
            spin-it-2014-target-frame-regression PRIVATE cxx_std_17
        )
        set_target_properties(
            spin-it-2014-target-frame-regression PROPERTIES CXX_EXTENSIONS OFF
        )
        add_test(
            NAME spin-it-2014-target-frame
            COMMAND "$<TARGET_FILE:spin-it-2014-target-frame-regression>"
        )
        set_tests_properties(
            spin-it-2014-target-frame
            PROPERTIES LABELS "spin-it;spin-it-2014;target-frame;regression"
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
            PROPERTIES
                LABELS "spin-it;spin-it-2014;dual-material;regression"
                TIMEOUT 240
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
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
        add_executable(spin-faster-polyhedral-derivative-regression
            tests/polyhedral_moment_differential_regression.cpp
        )
        target_compile_features(spin-faster-polyhedral-derivative-regression
            PRIVATE cxx_std_17
        )
        target_link_libraries(spin-faster-polyhedral-derivative-regression
            PRIVATE make_it_stand::faster_fullrank
        )
        add_test(NAME spin-faster-polyhedral-derivative
            COMMAND "$<TARGET_FILE:spin-faster-polyhedral-derivative-regression>"
        )
        set_tests_properties(spin-faster-polyhedral-derivative PROPERTIES
            LABELS "spin-it;spin-it-faster;derivative;geometry;regression"
        )
        if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
            add_test(NAME spin-faster-rational-material-authority-reuse
                COMMAND "$<TARGET_FILE:spin-faster-polyhedral-derivative-regression>" --rational-authority-reuse)
            set_tests_properties(spin-faster-rational-material-authority-reuse PROPERTIES
                LABELS "spin-it-faster;exact;derivative;regression" TIMEOUT 120)
        endif()

        add_executable(spin-faster-canonical-polyhedral-regression
            tests/canonical_polyhedral_replay_regression.cpp
        )
        target_compile_features(spin-faster-canonical-polyhedral-regression
            PRIVATE cxx_std_17
        )
        target_compile_definitions(spin-faster-canonical-polyhedral-regression
            PRIVATE differentiatePolyhedralMaterialMoments=differentiateCanonicalPolyhedralMaterialMoments
        )
        target_link_libraries(spin-faster-canonical-polyhedral-regression
            PRIVATE make_it_stand::faster_fullrank
        )
        add_test(NAME spin-faster-canonical-polyhedral-ordering
            COMMAND "$<TARGET_FILE:spin-faster-canonical-polyhedral-regression>"
        )
        set_tests_properties(spin-faster-canonical-polyhedral-ordering PROPERTIES
            LABELS "spin-it;spin-it-faster;canonical;geometry;regression"
        )

        add_executable(spin-faster-canonical-fullrank-regression
            tests/full_rank_hollowing_regression.cpp
        )
        target_compile_features(spin-faster-canonical-fullrank-regression
            PRIVATE cxx_std_17
        )
        target_compile_definitions(spin-faster-canonical-fullrank-regression
            PRIVATE
                integrateFullRankBoundaryMoments=integrateCanonicalFullRankBoundaryMoments
                evaluateFullRankMaterialMoments=evaluateCanonicalFullRankMaterialMoments
                makeFullRankMomentEvaluator=makeCanonicalFullRankMomentEvaluator
                solveFullRankQuadricHollowing=solveCanonicalFullRankQuadricHollowing
        )
        target_link_libraries(spin-faster-canonical-fullrank-regression
            PRIVATE make_it_stand::faster_fullrank
        )
        foreach(mode IN ITEMS adapter chain negative solver integration modes offtarget-yoyo audit replay boundary-identity selector-fallback)
            add_test(NAME "spin-faster-canonical-${mode}"
                COMMAND "$<TARGET_FILE:spin-faster-canonical-fullrank-regression>" "${mode}"
            )
            set_tests_properties("spin-faster-canonical-${mode}" PROPERTIES
                LABELS "spin-it;spin-it-faster;canonical;regression"
            )
        endforeach()

        if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
            add_test(NAME spin-faster-canonical-rational-integration
                COMMAND "$<TARGET_FILE:spin-faster-canonical-fullrank-regression>" "rational-integration")
            set_tests_properties(spin-faster-canonical-rational-integration PROPERTIES
                LABELS "spin-it;spin-it-faster;canonical;exact;independent-check;regression")
        endif()

        add_executable(spin-faster-unknown-seed-sweep
            tests/full_rank_unknown_seed_sweep.cpp
        )
        target_compile_features(spin-faster-unknown-seed-sweep PRIVATE cxx_std_17)
        target_compile_definitions(spin-faster-unknown-seed-sweep PRIVATE
            integrateFullRankBoundaryMoments=integrateCanonicalFullRankBoundaryMoments
            evaluateFullRankMaterialMoments=canonicalFullRankMaterialImplementation
            makeFullRankMomentEvaluator=canonicalFullRankAdapterImplementation
            solveFullRankQuadricHollowing=canonicalFullRankSolverImplementation
        )
        target_link_libraries(spin-faster-unknown-seed-sweep
            PRIVATE make_it_stand::faster_fullrank
        )
        add_test(NAME spin-faster-unknown-seed-sweep
            COMMAND "$<TARGET_FILE:spin-faster-unknown-seed-sweep>"
        )
        set_tests_properties(spin-faster-unknown-seed-sweep PROPERTIES
            LABELS "spin-it-faster;canonical;unknown-answer-seeds;synthetic-regression"
            TIMEOUT 3600
        )

        add_executable(
            planar-fidget-dogleg-regression
            tests/planar_fidget_dogleg_regression.cpp
        )
        target_compile_features(planar-fidget-dogleg-regression PRIVATE cxx_std_17)
        set_target_properties(planar-fidget-dogleg-regression PROPERTIES CXX_EXTENSIONS OFF)
        target_link_libraries(planar-fidget-dogleg-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-planar-fidget-dogleg
                 COMMAND "$<TARGET_FILE:planar-fidget-dogleg-regression>")
        set_tests_properties(spin-faster-planar-fidget-dogleg
            PROPERTIES LABELS "spin-it;spin-it-faster;planar;regression")

        add_executable(
            planar-disk-intersection-regression
            tests/planar_disk_intersection_regression.cpp
        )
        target_compile_features(planar-disk-intersection-regression PRIVATE cxx_std_17)
        set_target_properties(planar-disk-intersection-regression PROPERTIES CXX_EXTENSIONS OFF)
        target_link_libraries(planar-disk-intersection-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-planar-disk-intersection
                 COMMAND "$<TARGET_FILE:planar-disk-intersection-regression>")
        set_tests_properties(spin-faster-planar-disk-intersection
            PROPERTIES LABELS "spin-it;spin-it-faster;planar;regression")

        add_executable(
            planar-spin-moment-problem-regression
            tests/planar_spin_moment_problem_regression.cpp
        )
        target_compile_features(planar-spin-moment-problem-regression PRIVATE cxx_std_17)
        set_target_properties(planar-spin-moment-problem-regression PROPERTIES CXX_EXTENSIONS OFF)
        target_link_libraries(planar-spin-moment-problem-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-planar-moment-problem
                 COMMAND "$<TARGET_FILE:planar-spin-moment-problem-regression>")
        set_tests_properties(spin-faster-planar-moment-problem
            PROPERTIES LABELS "spin-it;spin-it-faster;planar;regression")

        add_executable(
            planar-mass-moments-regression
            tests/planar_mass_moments_regression.cpp
        )
        target_compile_features(planar-mass-moments-regression PRIVATE cxx_std_17)
        set_target_properties(planar-mass-moments-regression PROPERTIES CXX_EXTENSIONS OFF)
        target_link_libraries(planar-mass-moments-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-planar-mass-moments
                 COMMAND "$<TARGET_FILE:planar-mass-moments-regression>")
        set_tests_properties(spin-faster-planar-mass-moments
            PROPERTIES LABELS "spin-it;spin-it-faster;planar;regression")

        add_executable(
            planar-extrusion-regression
            tests/planar_extrusion_regression.cpp
        )
        target_compile_features(planar-extrusion-regression PRIVATE cxx_std_17)
        set_target_properties(planar-extrusion-regression PROPERTIES CXX_EXTENSIONS OFF)
        target_link_libraries(planar-extrusion-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-planar-extrusion
                 COMMAND "$<TARGET_FILE:planar-extrusion-regression>")
        set_tests_properties(spin-faster-planar-extrusion
            PROPERTIES LABELS "spin-it;spin-it-faster;planar;regression")

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
            diagnose-boolean-arrangement-pwn
            experiments/spin_it_reproduction/strict/diagnose_boolean_arrangement_pwn.cpp
        )
        target_compile_features(diagnose-boolean-arrangement-pwn PRIVATE cxx_std_17)
        set_target_properties(
            diagnose-boolean-arrangement-pwn
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            diagnose-boolean-arrangement-pwn
            PRIVATE make_it_stand::core
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
        if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
            add_test(NAME spin-faster-rational-source-triangulator
                COMMAND "$<TARGET_FILE:quadric-triangulator-regression>" --rational-source)
            set_tests_properties(spin-faster-rational-source-triangulator PROPERTIES
                LABELS "spin-it-faster;geometry;exact;precision;regression" TIMEOUT 180)
            add_test(NAME spin-faster-rational-hyperboloid-source
                COMMAND "$<TARGET_FILE:quadric-triangulator-regression>" --rational-hyperboloid)
            set_tests_properties(spin-faster-rational-hyperboloid-source PROPERTIES
                LABELS "spin-it-faster;geometry;exact;precision;regression" TIMEOUT 180)
        endif()
        add_test(
            NAME spin-faster-adaptive-cap-refinement
            COMMAND "$<TARGET_FILE:quadric-triangulator-regression>"
                    --adaptive-cap-refinement-only
        )
        set_tests_properties(
            spin-faster-quadric-triangulator spin-faster-adaptive-cap-refinement
            PROPERTIES LABELS "spin-it;spin-it-faster;geometry;regression"
        )

        add_executable(
            quadric-inverse-gauss-map-regression
            tests/quadric_inverse_gauss_map_regression.cpp
        )
        target_compile_features(
            quadric-inverse-gauss-map-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            quadric-inverse-gauss-map-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            quadric-inverse-gauss-map-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-inverse-gauss-map
            COMMAND "$<TARGET_FILE:quadric-inverse-gauss-map-regression>"
        )
        set_tests_properties(
            spin-faster-inverse-gauss-map
            PROPERTIES LABELS "spin-it;spin-it-faster;geometry;regression"
        )

        foreach(authority_component IN ITEMS
                coefficient_authority authority_source_point authority_truncation
                authority_triangulation authority_material authority_finite_differential
                authority_root_differential authority_source_motion authority_transport
                authority_production_derivative authority_linear_solve)
            set(authority_target "quadric-${authority_component}-regression")
            add_executable(${authority_target}
                "tests/quadric_${authority_component}_regression.cpp")
            target_compile_features(${authority_target} PRIVATE cxx_std_17)
            set_target_properties(${authority_target} PROPERTIES CXX_EXTENSIONS OFF)
            target_link_libraries(${authority_target} PRIVATE make_it_stand::core)
            add_test(NAME "spin-faster-${authority_component}" COMMAND ${authority_target})
            set_tests_properties("spin-faster-${authority_component}" PROPERTIES
                LABELS "spin-it-faster;precision;local-kernel;regression" TIMEOUT 180)
        endforeach()

        add_executable(
            quadric-source-motion-regression
            tests/quadric_source_motion_regression.cpp
        )
        add_executable(quadric-projective-edge-regression tests/quadric_projective_edge_regression.cpp)
        target_compile_features(quadric-projective-edge-regression PRIVATE cxx_std_17)
        set_target_properties(quadric-projective-edge-regression PROPERTIES CXX_EXTENSIONS OFF)
        target_link_libraries(quadric-projective-edge-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-projective-edge COMMAND "$<TARGET_FILE:quadric-projective-edge-regression>")
        set_tests_properties(spin-faster-projective-edge PROPERTIES LABELS "spin-it-faster;geometry;derivative;regression")
        if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
            add_executable(rank-two-algebraic-backend-probe tests/rank_two_algebraic_backend_probe.cpp)
            target_compile_features(rank-two-algebraic-backend-probe PRIVATE cxx_std_17)
            target_link_libraries(rank-two-algebraic-backend-probe PRIVATE CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-backend COMMAND rank-two-algebraic-backend-probe)
            add_executable(rank-two-source-input-regression tests/rank_two_source_input_regression.cpp)
            target_sources(rank-two-source-input-regression PRIVATE src/spin/RankTwoSource.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
            target_include_directories(rank-two-source-input-regression PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-source-input-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-source-input-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-source-input COMMAND rank-two-source-input-regression)
            add_executable(rank-two-algebraic-regression tests/rank_two_algebraic_regression.cpp)
            target_sources(rank-two-algebraic-regression PRIVATE src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoSource.cpp src/spin/QuadricCoefficients.cpp
                src/manufacturing/SlicerValidation.cpp)
            target_include_directories(rank-two-algebraic-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-algebraic-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-algebraic-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            if(MAKE_IT_STAND_RANK_TWO_FLINT_ROOT)
                target_compile_definitions(rank-two-algebraic-regression PRIVATE
                    RANK_TWO_DEGREE144_FIXTURE="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/rank_two_degree144_root.txt")
            endif()
            option(MAKE_IT_STAND_BUILD_RANK_TWO_DIAGNOSTIC_PROFILE
                "Build the captured-root certification timing diagnostic" OFF)
            if(MAKE_IT_STAND_BUILD_RANK_TWO_DIAGNOSTIC_PROFILE)
                if(NOT MAKE_IT_STAND_RANK_TWO_FLINT_ROOT)
                    message(FATAL_ERROR "Certification timing requires the explicitly audited factor backend")
                endif()
                add_executable(rank-two-certification-profile tests/rank_two_certification_profile.cpp)
                target_sources(rank-two-certification-profile PRIVATE src/spin/RankTwoAlgebraic.cpp
                    src/spin/RankTwoSource.cpp src/spin/QuadricCoefficients.cpp
                    src/manufacturing/SlicerValidation.cpp)
                target_include_directories(rank-two-certification-profile PRIVATE
                    "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
                target_compile_features(rank-two-certification-profile PRIVATE cxx_std_17)
                target_link_libraries(rank-two-certification-profile PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
                add_test(NAME spin-faster-rank-two-certification-profile COMMAND rank-two-certification-profile
                    "${CMAKE_CURRENT_BINARY_DIR}/rank-two-certification-profile-input.txt")
                set_tests_properties(spin-faster-rank-two-certification-profile PROPERTIES
                    LABELS "spin-it-faster;diagnostic;exact;rank-two")
                foreach(fiber_mode IN ITEMS native decomposed optimized factor-resultant factor-resultant-reverse)
                    set(fiber_target "rank-two-fiber-${fiber_mode}-profile")
                    add_executable(${fiber_target} tests/rank_two_fiber_profile.cpp)
                    target_sources(${fiber_target} PRIVATE src/spin/RankTwoAlgebraic.cpp
                        src/spin/RankTwoSource.cpp src/spin/QuadricCoefficients.cpp
                        src/manufacturing/SlicerValidation.cpp)
                    target_include_directories(${fiber_target} PRIVATE
                        "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
                    target_compile_features(${fiber_target} PRIVATE cxx_std_17)
                    target_link_libraries(${fiber_target} PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
                    if(fiber_mode STREQUAL "decomposed")
                        target_compile_definitions(${fiber_target} PRIVATE RANK_TWO_DECOMPOSE_NATIVE_FIBER=1)
                    elseif(fiber_mode STREQUAL "optimized")
                        target_compile_definitions(${fiber_target} PRIVATE RANK_TWO_OPTIMIZED_NATIVE_FIBER=1)
                    elseif(fiber_mode MATCHES "^factor-resultant")
                        target_compile_definitions(${fiber_target} PRIVATE RANK_TWO_FACTOR_RESULTANT_PROBE=1)
                        if(fiber_mode STREQUAL "factor-resultant-reverse")
                            target_compile_definitions(${fiber_target} PRIVATE RANK_TWO_FACTOR_RESULTANT_REVERSE=1)
                        endif()
                    endif()
                    add_test(NAME "spin-faster-rank-two-fiber-${fiber_mode}-profile" COMMAND ${fiber_target}
                        "${CMAKE_CURRENT_BINARY_DIR}/rank-two-fiber-profile-input.txt")
                    set_tests_properties("spin-faster-rank-two-fiber-${fiber_mode}-profile" PROPERTIES
                        LABELS "spin-it-faster;diagnostic;exact;rank-two")
                endforeach()
                foreach(weighted_target IN ITEMS rank-two-weighted-profile rank-two-weighted-optimized-profile)
                add_executable(${weighted_target} tests/rank_two_weighted_profile.cpp)
                target_sources(${weighted_target} PRIVATE src/spin/RankTwoAlgebraic.cpp
                    src/spin/RankTwoCarrier.cpp src/spin/RankTwoSource.cpp
                    src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
                target_include_directories(${weighted_target} PRIVATE
                    "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
                target_compile_features(${weighted_target} PRIVATE cxx_std_17)
                target_link_libraries(${weighted_target} PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
                if(weighted_target STREQUAL "rank-two-weighted-optimized-profile")
                    target_compile_definitions(${weighted_target} PRIVATE RANK_TWO_EXPECT_AUTHENTICATED_WEIGHTED_ROOTS=1)
                endif()
                add_test(NAME spin-faster-${weighted_target} COMMAND ${weighted_target}
                    "${CMAKE_CURRENT_BINARY_DIR}/rank-two-weighted-profile-input.txt")
                set_tests_properties(spin-faster-${weighted_target} PROPERTIES
                    LABELS "spin-it-faster;diagnostic;exact;rank-two")
                endforeach()
            endif()
            add_test(NAME spin-faster-rank-two-algebraic COMMAND rank-two-algebraic-regression)
            add_test(NAME spin-faster-rank-two-algebraic-degree40
                COMMAND rank-two-algebraic-regression --degree40)
            set_tests_properties(spin-faster-rank-two-algebraic-degree40 PROPERTIES
                LABELS "spin-it-faster;geometry;regression")
            add_executable(rank-two-carrier-regression tests/rank_two_carrier_regression.cpp)
            target_sources(rank-two-carrier-regression PRIVATE src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoSource.cpp src/spin/QuadricCoefficients.cpp
                src/manufacturing/SlicerValidation.cpp)
            target_sources(rank-two-carrier-regression PRIVATE src/spin/RankTwoCarrier.cpp)
            target_include_directories(rank-two-carrier-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-carrier-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-carrier-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-carrier-kernels COMMAND rank-two-carrier-regression)
            if(TARGET rank_two_native_factor)
            get_target_property(_rank_two_carrier_sources rank-two-carrier-regression SOURCES)
            add_executable(rank-two-actual-parent-affine-regression ${_rank_two_carrier_sources})
            target_compile_definitions(rank-two-actual-parent-affine-regression PRIVATE RANK_TWO_ACTUAL_PARENT_AFFINE=1)
            target_include_directories(rank-two-actual-parent-affine-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-actual-parent-affine-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-actual-parent-affine-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-actual-parent-affine COMMAND rank-two-actual-parent-affine-regression
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/rank_two_r188_failed_affine.txt")
            set_tests_properties(spin-faster-rank-two-actual-parent-affine PROPERTIES
                TIMEOUT 500 LABELS "spin-it-faster;exact;rank-two;regression")
            endif()
            add_executable(rank-two-arc-strip-regression tests/rank_two_arc_strip_regression.cpp)
            add_executable(rank-two-strip-memo-regression tests/rank_two_strip_memo_regression.cpp)
            target_sources(rank-two-strip-memo-regression PRIVATE src/spin/RankTwoSource.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp src/spin/RankTwoArcStrip.cpp)
            target_include_directories(rank-two-strip-memo-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-strip-memo-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-strip-memo-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-strip-memo COMMAND rank-two-strip-memo-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-strip-memo-evidence.txt")
            set_tests_properties(spin-faster-rank-two-strip-memo PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            target_sources(rank-two-arc-strip-regression PRIVATE src/spin/RankTwoSource.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
            target_sources(rank-two-arc-strip-regression PRIVATE src/spin/RankTwoArcStrip.cpp)
            target_include_directories(rank-two-arc-strip-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-arc-strip-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-arc-strip-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-arc-strip COMMAND rank-two-arc-strip-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-arc-strip-evidence.txt")
            set_tests_properties(spin-faster-rank-two-arc-strip PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-endpoint-arc-regression tests/rank_two_endpoint_arc_regression.cpp)
            target_sources(rank-two-endpoint-arc-regression PRIVATE src/spin/RankTwoCarrier.cpp
                src/spin/RankTwoArrangement.cpp src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp
                src/spin/RankTwoArcStrip.cpp src/spin/QuadricCoefficients.cpp
                src/manufacturing/SlicerValidation.cpp)
            target_sources(rank-two-endpoint-arc-regression PRIVATE src/spin/RankTwoEndpointArc.cpp)
            target_include_directories(rank-two-endpoint-arc-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-endpoint-arc-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-endpoint-arc-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-endpoint-arc COMMAND rank-two-endpoint-arc-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-endpoint-arc-evidence.txt")
            set_tests_properties(spin-faster-rank-two-endpoint-arc PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-arc-sampling-regression tests/rank_two_arc_sampling_regression.cpp)
            add_executable(rank-two-arc-registry-regression tests/rank_two_arc_registry_regression.cpp)
            add_executable(rank-two-registry-edge-regression tests/rank_two_registry_edge_regression.cpp)
            add_executable(rank-two-chord-incidence-regression tests/rank_two_chord_incidence_regression.cpp)
            add_executable(rank-two-registry-rotation-regression tests/rank_two_registry_rotation_regression.cpp)
            add_executable(rank-two-regular-registry-regression tests/rank_two_regular_registry_regression.cpp)
            target_sources(rank-two-regular-registry-regression PRIVATE
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp
                src/spin/RankTwoCarrier.cpp src/spin/RankTwoArrangement.cpp src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp src/spin/RankTwoArcSampling.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp
                src/spin/RankTwoArcRegistry.cpp src/spin/RankTwoArcStrip.cpp
                src/spin/RankTwoEndpointArc.cpp src/spin/RankTwoEndpointFan.cpp src/spin/RankTwoEndpointNeighborhood.cpp)
            target_include_directories(rank-two-regular-registry-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-regular-registry-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-regular-registry-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-regular-registry COMMAND rank-two-regular-registry-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-regular-registry-evidence.txt")
            set_tests_properties(spin-faster-rank-two-regular-registry PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            target_sources(rank-two-registry-rotation-regression PRIVATE
                src/spin/RankTwoRegistryRotation.cpp
                src/spin/RankTwoChordIncidence.cpp
                src/spin/RankTwoCarrier.cpp src/spin/RankTwoArrangement.cpp src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp src/spin/RankTwoArcSampling.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp
                src/spin/RankTwoArcRegistry.cpp src/spin/RankTwoArcStrip.cpp
                src/spin/RankTwoEndpointArc.cpp src/spin/RankTwoEndpointFan.cpp src/spin/RankTwoEndpointNeighborhood.cpp)
            target_include_directories(rank-two-registry-rotation-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-registry-rotation-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-registry-rotation-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-registry-rotation COMMAND rank-two-registry-rotation-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-registry-rotation-evidence.txt")
            set_tests_properties(spin-faster-rank-two-registry-rotation PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            target_sources(rank-two-chord-incidence-regression PRIVATE
                src/spin/RankTwoChordIncidence.cpp
                src/spin/RankTwoCarrier.cpp src/spin/RankTwoArrangement.cpp src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp src/spin/RankTwoArcSampling.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp
                src/spin/RankTwoArcRegistry.cpp src/spin/RankTwoArcStrip.cpp
                src/spin/RankTwoEndpointArc.cpp src/spin/RankTwoEndpointFan.cpp src/spin/RankTwoEndpointNeighborhood.cpp)
            target_include_directories(rank-two-chord-incidence-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-chord-incidence-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-chord-incidence-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-chord-incidence COMMAND rank-two-chord-incidence-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-chord-incidence-evidence.txt")
            set_tests_properties(spin-faster-rank-two-chord-incidence PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-arc-refinement-regression tests/rank_two_arc_refinement_regression.cpp)
            target_sources(rank-two-arc-refinement-regression PRIVATE
                src/spin/RankTwoArcRefinement.cpp
                src/spin/RankTwoRegistryEdge.cpp
                src/spin/RankTwoCarrier.cpp src/spin/RankTwoArrangement.cpp src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp src/spin/RankTwoArcSampling.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp
                src/spin/RankTwoArcRegistry.cpp src/spin/RankTwoArcStrip.cpp
                src/spin/RankTwoEndpointArc.cpp src/spin/RankTwoEndpointFan.cpp src/spin/RankTwoEndpointNeighborhood.cpp)
            target_include_directories(rank-two-arc-refinement-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-arc-refinement-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-arc-refinement-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-arc-refinement COMMAND rank-two-arc-refinement-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-arc-refinement-evidence.txt")
            set_tests_properties(spin-faster-rank-two-arc-refinement PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            target_sources(rank-two-registry-edge-regression PRIVATE
                src/spin/RankTwoRegistryEdge.cpp
                src/spin/RankTwoCarrier.cpp src/spin/RankTwoArrangement.cpp src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp src/spin/RankTwoArcSampling.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp
                src/spin/RankTwoArcRegistry.cpp src/spin/RankTwoArcStrip.cpp
                src/spin/RankTwoEndpointArc.cpp src/spin/RankTwoEndpointFan.cpp src/spin/RankTwoEndpointNeighborhood.cpp)
            target_include_directories(rank-two-registry-edge-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-registry-edge-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-registry-edge-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-registry-edge COMMAND rank-two-registry-edge-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-registry-edge-evidence.txt")
            set_tests_properties(spin-faster-rank-two-registry-edge PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            # These standalone regressions used identical compilation flags for
            # these thirteen sources. Compile them once per FRESH build, not once
            # per test executable; each test still runs in its own process.
            set(rank_two_registry_shared_sources
                src/spin/RankTwoCarrier.cpp src/spin/RankTwoArrangement.cpp src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp src/spin/RankTwoArcSampling.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp
                src/spin/RankTwoArcRegistry.cpp src/spin/RankTwoArcStrip.cpp
                src/spin/RankTwoEndpointArc.cpp src/spin/RankTwoEndpointFan.cpp src/spin/RankTwoEndpointNeighborhood.cpp)
            add_library(rank-two-registry-test-kernels OBJECT ${rank_two_registry_shared_sources})
            target_include_directories(rank-two-registry-test-kernels PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-registry-test-kernels PRIVATE cxx_std_17)
            target_link_libraries(rank-two-registry-test-kernels PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_executable(rank-two-endpoint-neighborhood-regression tests/rank_two_endpoint_neighborhood_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-endpoint-neighborhood-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-endpoint-neighborhood-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-endpoint-neighborhood-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-endpoint-neighborhood COMMAND rank-two-endpoint-neighborhood-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-endpoint-neighborhood-evidence.txt")
            set_tests_properties(spin-faster-rank-two-endpoint-neighborhood PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-neighborhood-registry-regression tests/rank_two_neighborhood_registry_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-neighborhood-registry-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-neighborhood-registry-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-neighborhood-registry-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-neighborhood-registry COMMAND rank-two-neighborhood-registry-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-neighborhood-registry-evidence.txt")
            set_tests_properties(spin-faster-rank-two-neighborhood-registry PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-registry-subdivision-regression tests/rank_two_registry_subdivision_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-registry-subdivision-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-registry-subdivision-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-registry-subdivision-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-registry-subdivision COMMAND rank-two-registry-subdivision-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-registry-subdivision-evidence.txt")
            set_tests_properties(spin-faster-rank-two-registry-subdivision PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-chord-sweep-regression tests/rank_two_chord_sweep_regression.cpp
                src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-chord-sweep-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-chord-sweep-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-chord-sweep-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-chord-sweep COMMAND rank-two-chord-sweep-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-chord-sweep-evidence.txt")
            set_tests_properties(spin-faster-rank-two-chord-sweep PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-mixed-arc-refinement-regression tests/rank_two_mixed_arc_refinement_regression.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-mixed-arc-refinement-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-mixed-arc-refinement-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-mixed-arc-refinement-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-mixed-arc-refinement COMMAND rank-two-mixed-arc-refinement-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-mixed-arc-refinement-evidence.txt")
            set_tests_properties(spin-faster-rank-two-mixed-arc-refinement PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-mixed-intervals-regression tests/rank_two_mixed_intervals_regression.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-mixed-intervals-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-mixed-intervals-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-mixed-intervals-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-mixed-intervals COMMAND rank-two-mixed-intervals-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-mixed-intervals-evidence.txt")
            set_tests_properties(spin-faster-rank-two-mixed-intervals PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-affine-dependency-regression tests/rank_two_affine_dependency_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-affine-dependency-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-affine-dependency-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-affine-dependency-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-affine-dependency COMMAND rank-two-affine-dependency-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-affine-dependency-evidence.txt")
            set_tests_properties(spin-faster-rank-two-affine-dependency PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-finite-arc-coverage-regression tests/rank_two_finite_arc_coverage_regression.cpp
                src/spin/RankTwoFiniteArcCoverage.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-finite-arc-coverage-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-finite-arc-coverage-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-finite-arc-coverage-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-finite-arc-coverage COMMAND rank-two-finite-arc-coverage-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-finite-arc-coverage-evidence.txt")
            set_tests_properties(spin-faster-rank-two-finite-arc-coverage PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-finite-coverage-continuation-regression tests/rank_two_finite_coverage_continuation_regression.cpp
                src/spin/RankTwoFiniteArcCoverage.cpp src/spin/RankTwoArcRefinement.cpp
                src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-finite-coverage-continuation-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-finite-coverage-continuation-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-finite-coverage-continuation-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-finite-coverage-continuation COMMAND rank-two-finite-coverage-continuation-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-finite-coverage-continuation-evidence.txt")
            set_tests_properties(spin-faster-rank-two-finite-coverage-continuation PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-finite-partition-regression tests/rank_two_finite_partition_regression.cpp
                src/spin/RankTwoFinitePartition.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-finite-partition-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-finite-partition-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-finite-partition-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-finite-partition COMMAND rank-two-finite-partition-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-finite-partition-evidence.txt")
            set_tests_properties(spin-faster-rank-two-finite-partition PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-finite-mesh-scope-regression tests/rank_two_finite_mesh_scope_regression.cpp
                src/spin/RankTwoFinitePartition.cpp src/spin/RankTwoFiniteTriangulation.cpp
                src/spin/RankTwoFiniteMeshScope.cpp
                src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-finite-mesh-scope-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-finite-mesh-scope-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-finite-mesh-scope-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-finite-mesh-scope COMMAND rank-two-finite-mesh-scope-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-finite-mesh-scope-evidence.txt")
            set_tests_properties(spin-faster-rank-two-finite-mesh-scope PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-covered-finite-mesh-regression tests/rank_two_covered_finite_mesh_regression.cpp
                src/spin/RankTwoFiniteArcCoverage.cpp src/spin/RankTwoArcRefinement.cpp
                src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoFiniteTriangulation.cpp src/spin/RankTwoChordIncidence.cpp
                src/spin/RankTwoRegistryRotation.cpp $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoCoveredFiniteMesh.cpp")
                target_sources(rank-two-covered-finite-mesh-regression PRIVATE src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp)
            endif()
            target_include_directories(rank-two-covered-finite-mesh-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-covered-finite-mesh-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-covered-finite-mesh-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-covered-finite-mesh COMMAND rank-two-covered-finite-mesh-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-covered-finite-mesh-evidence.txt")
            set_tests_properties(spin-faster-rank-two-covered-finite-mesh PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-complete-endpoint-import-regression tests/rank_two_complete_endpoint_import_regression.cpp
                src/spin/RankTwoSharedFiniteSeams.cpp src/spin/RankTwoArcRegistrySeams.cpp
                src/spin/RankTwoSharedEndpointSeams.cpp src/spin/RankTwoFinitePartition.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-complete-endpoint-import-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-complete-endpoint-import-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-complete-endpoint-import-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-complete-endpoint-import COMMAND rank-two-complete-endpoint-import-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-complete-endpoint-import-evidence.txt")
            set_tests_properties(spin-faster-rank-two-complete-endpoint-import PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-complete-endpoint-binding-regression tests/rank_two_complete_endpoint_import_regression.cpp
                src/spin/RankTwoSharedFiniteSeams.cpp src/spin/RankTwoArcRegistrySeams.cpp
                src/spin/RankTwoSharedEndpointSeams.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_compile_definitions(rank-two-complete-endpoint-binding-regression PRIVATE RANK_TWO_VERIFY_ENDPOINT_CONSUMPTION=1)
            target_include_directories(rank-two-complete-endpoint-binding-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-complete-endpoint-binding-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-complete-endpoint-binding-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-complete-endpoint-binding COMMAND rank-two-complete-endpoint-binding-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-complete-endpoint-binding-evidence.txt")
            set_tests_properties(spin-faster-rank-two-complete-endpoint-binding PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-refinement-reuse-regression tests/rank_two_finite_coverage_continuation_regression.cpp
                src/spin/RankTwoFiniteArcCoverage.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_compile_definitions(rank-two-refinement-reuse-regression PRIVATE RANK_TWO_VERIFY_REFINEMENT_REUSE=1)
            target_include_directories(rank-two-refinement-reuse-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-refinement-reuse-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-refinement-reuse-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-refinement-reuse COMMAND rank-two-refinement-reuse-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-refinement-reuse-evidence.txt")
            set_tests_properties(spin-faster-rank-two-refinement-reuse PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-registry-reuse-regression tests/rank_two_registry_reuse_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-registry-reuse-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-registry-reuse-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-registry-reuse-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-registry-reuse COMMAND rank-two-registry-reuse-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-registry-reuse-evidence.txt")
            set_tests_properties(spin-faster-rank-two-registry-reuse PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-complete-finite-coverage-regression tests/rank_two_complete_endpoint_import_regression.cpp
                src/spin/RankTwoSharedFiniteSeams.cpp src/spin/RankTwoArcRegistrySeams.cpp
                src/spin/RankTwoSharedEndpointSeams.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoFiniteArcCoverage.cpp src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_compile_definitions(rank-two-complete-finite-coverage-regression PRIVATE RANK_TWO_VERIFY_COMPLETE_FINITE_COVERAGE=1)
            target_include_directories(rank-two-complete-finite-coverage-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-complete-finite-coverage-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-complete-finite-coverage-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-complete-finite-coverage COMMAND rank-two-complete-finite-coverage-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-complete-finite-coverage-evidence.txt")
            set_tests_properties(spin-faster-rank-two-complete-finite-coverage PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-complete-finite-mesh-regression tests/rank_two_complete_endpoint_import_regression.cpp
                src/spin/RankTwoSharedFiniteSeams.cpp src/spin/RankTwoArcRegistrySeams.cpp
                src/spin/RankTwoSharedEndpointSeams.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoFiniteArcCoverage.cpp src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp
                src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp src/spin/RankTwoFiniteTriangulation.cpp
                src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_compile_definitions(rank-two-complete-finite-mesh-regression PRIVATE RANK_TWO_VERIFY_COMPLETE_FINITE_MESH=1)
            target_include_directories(rank-two-complete-finite-mesh-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-complete-finite-mesh-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-complete-finite-mesh-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-complete-finite-mesh COMMAND rank-two-complete-finite-mesh-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-complete-finite-mesh-evidence.txt")
            set_tests_properties(spin-faster-rank-two-complete-finite-mesh PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;integration" TIMEOUT 3000 RUN_SERIAL TRUE)
            get_target_property(_rank_two_complete_mesh_sources rank-two-complete-finite-mesh-regression SOURCES)
            add_executable(rank-two-complete-owner-material-regression ${_rank_two_complete_mesh_sources}
                src/spin/RankTwoWholeInteriorInventory.cpp src/spin/RankTwoInteriorInventory.cpp
                src/spin/RankTwoAutomaticMaterialMesh.cpp src/spin/RankTwoStrictMaterialMesh.cpp
                src/spin/RankTwoMaterialLifts.cpp src/spin/RankTwoMovingMeshEdge.cpp
                src/spin/RankTwoInteriorMesh.cpp src/spin/RankTwoMovingInterior.cpp)
            target_compile_definitions(rank-two-complete-owner-material-regression PRIVATE RANK_TWO_VERIFY_COMPLETE_MATERIAL_MESH=1)
            target_include_directories(rank-two-complete-owner-material-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-complete-owner-material-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-complete-owner-material-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-complete-owner-material COMMAND rank-two-complete-owner-material-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-complete-owner-material-evidence.txt")
            set_tests_properties(spin-faster-rank-two-complete-owner-material PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;material-thickness;integration" TIMEOUT 3000 RUN_SERIAL TRUE)
            option(MAKE_IT_STAND_BUILD_RANK_TWO_COVERAGE_PROFILE
                "Build test-only wall/thread-CPU whole coverage observation" OFF)
            if(MAKE_IT_STAND_BUILD_RANK_TWO_COVERAGE_PROFILE)
                add_executable(rank-two-phase-profile-regression tests/rank_two_phase_profile_regression.cpp)
                target_compile_features(rank-two-phase-profile-regression PRIVATE cxx_std_17)
                add_test(NAME spin-faster-rank-two-phase-profile COMMAND rank-two-phase-profile-regression)
                set_tests_properties(spin-faster-rank-two-phase-profile PROPERTIES LABELS "spin-it-faster;diagnostic;rank-two")
                foreach(profile IN ITEMS endpoint-timing complete-finite-coverage-timing)
                    add_executable(rank-two-${profile}-profile tests/rank_two_complete_endpoint_import_regression.cpp
                        src/spin/RankTwoSharedFiniteSeams.cpp src/spin/RankTwoArcRegistrySeams.cpp
                        src/spin/RankTwoSharedEndpointSeams.cpp src/spin/RankTwoFinitePartition.cpp
                        src/spin/RankTwoFiniteArcCoverage.cpp src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp
                        $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
                    target_compile_definitions(rank-two-${profile}-profile PRIVATE RANK_TWO_VERIFY_COVERAGE_TIMING=1)
                    if(profile STREQUAL "endpoint-timing")
                        target_compile_definitions(rank-two-${profile}-profile PRIVATE RANK_TWO_VERIFY_ENDPOINT_CONSUMPTION=1)
                    else()
                        target_compile_definitions(rank-two-${profile}-profile PRIVATE RANK_TWO_VERIFY_COMPLETE_FINITE_COVERAGE=1)
                    endif()
                    target_include_directories(rank-two-${profile}-profile PRIVATE
                        "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
                    target_compile_features(rank-two-${profile}-profile PRIVATE cxx_std_17)
                    target_link_libraries(rank-two-${profile}-profile PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
                    add_test(NAME spin-faster-rank-two-${profile}-profile COMMAND rank-two-${profile}-profile
                        "${CMAKE_CURRENT_BINARY_DIR}/rank-two-${profile}-profile-evidence.txt")
                    set_tests_properties(spin-faster-rank-two-${profile}-profile PROPERTIES
                        LABELS "spin-it-faster;diagnostic;exact;rank-two")
                endforeach()
            endif()
            add_executable(rank-two-projection-endpoint-import-regression tests/rank_two_projection_endpoint_import_regression.cpp
                src/spin/RankTwoSharedFiniteSeams.cpp src/spin/RankTwoArcRegistrySeams.cpp
                src/spin/RankTwoSharedEndpointSeams.cpp src/spin/RankTwoFinitePartition.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-projection-endpoint-import-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-projection-endpoint-import-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-projection-endpoint-import-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-projection-endpoint-import COMMAND rank-two-projection-endpoint-import-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-projection-endpoint-import-evidence.txt")
            set_tests_properties(spin-faster-rank-two-projection-endpoint-import PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-finite-endpoint-inventory-regression tests/rank_two_finite_endpoint_inventory_regression.cpp
                src/spin/RankTwoSharedFiniteSeams.cpp src/spin/RankTwoArcRegistrySeams.cpp
                src/spin/RankTwoSharedEndpointSeams.cpp src/spin/RankTwoFinitePartition.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-finite-endpoint-inventory-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-finite-endpoint-inventory-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-finite-endpoint-inventory-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-finite-endpoint-inventory COMMAND rank-two-finite-endpoint-inventory-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-finite-endpoint-inventory-evidence.txt")
            set_tests_properties(spin-faster-rank-two-finite-endpoint-inventory PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-all-seam-endpoints-regression tests/rank_two_all_seam_endpoints_regression.cpp
                src/spin/RankTwoSharedFiniteSeams.cpp src/spin/RankTwoArcRegistrySeams.cpp
                src/spin/RankTwoSharedEndpointSeams.cpp src/spin/RankTwoFinitePartition.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-all-seam-endpoints-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-all-seam-endpoints-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-all-seam-endpoints-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-all-seam-endpoints COMMAND rank-two-all-seam-endpoints-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-all-seam-endpoints-evidence.txt")
            set_tests_properties(spin-faster-rank-two-all-seam-endpoints PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-seam-endpoint-regression tests/rank_two_seam_endpoint_regression.cpp
                src/spin/RankTwoSharedFiniteSeams.cpp src/spin/RankTwoArcRegistrySeams.cpp src/spin/RankTwoFinitePartition.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoSharedEndpointSeams.cpp")
                target_sources(rank-two-seam-endpoint-regression PRIVATE src/spin/RankTwoSharedEndpointSeams.cpp)
            endif()
            target_include_directories(rank-two-seam-endpoint-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-seam-endpoint-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-seam-endpoint-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-seam-endpoint COMMAND rank-two-seam-endpoint-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-seam-endpoint-evidence.txt")
            set_tests_properties(spin-faster-rank-two-seam-endpoint PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-seam-refinement-regression tests/rank_two_seam_refinement_regression.cpp
                src/spin/RankTwoArcRegistrySeams.cpp src/spin/RankTwoSharedFiniteSeams.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-seam-refinement-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-seam-refinement-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-seam-refinement-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-seam-refinement COMMAND rank-two-seam-refinement-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-seam-refinement-evidence.txt")
            set_tests_properties(spin-faster-rank-two-seam-refinement PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-seam-registry-regression tests/rank_two_seam_registry_regression.cpp
                src/spin/RankTwoArcRegistrySeams.cpp src/spin/RankTwoSharedFiniteSeams.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoFiniteTriangulation.cpp src/spin/RankTwoChordIncidence.cpp
                src/spin/RankTwoRegistryRotation.cpp $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            target_include_directories(rank-two-seam-registry-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-seam-registry-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-seam-registry-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-seam-registry COMMAND rank-two-seam-registry-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-seam-registry-evidence.txt")
            set_tests_properties(spin-faster-rank-two-seam-registry PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-shared-finite-seams-regression tests/rank_two_shared_finite_seams_regression.cpp
                src/spin/RankTwoCarrier.cpp src/spin/RankTwoArrangement.cpp
                src/spin/RankTwoAlgebraic.cpp src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoSharedFiniteSeams.cpp")
                target_sources(rank-two-shared-finite-seams-regression PRIVATE src/spin/RankTwoSharedFiniteSeams.cpp)
            endif()
            target_include_directories(rank-two-shared-finite-seams-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-shared-finite-seams-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-shared-finite-seams-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-shared-finite-seams COMMAND rank-two-shared-finite-seams-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-shared-finite-seams-evidence.txt")
            set_tests_properties(spin-faster-rank-two-shared-finite-seams PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;shared-registry;regression")
            add_executable(rank-two-finite-cell-topology-regression tests/rank_two_finite_cell_topology_regression.cpp
                src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp src/spin/RankTwoFiniteArcCoverage.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp
                src/spin/RankTwoFinitePartition.cpp src/spin/RankTwoFiniteTriangulation.cpp
                src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoFiniteCellTopology.cpp")
                target_sources(rank-two-finite-cell-topology-regression PRIVATE src/spin/RankTwoFiniteCellTopology.cpp)
            endif()
            target_include_directories(rank-two-finite-cell-topology-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-finite-cell-topology-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-finite-cell-topology-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-finite-cell-topology COMMAND rank-two-finite-cell-topology-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-finite-cell-topology-evidence.txt")
            set_tests_properties(spin-faster-rank-two-finite-cell-topology PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-moving-interior-regression tests/rank_two_moving_interior_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoMovingInterior.cpp")
                target_sources(rank-two-moving-interior-regression PRIVATE src/spin/RankTwoMovingInterior.cpp)
            endif()
            target_include_directories(rank-two-moving-interior-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-moving-interior-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-moving-interior-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-moving-interior COMMAND rank-two-moving-interior-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-moving-interior-evidence.txt")
            set_tests_properties(spin-faster-rank-two-moving-interior PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-interior-mesh-regression tests/rank_two_interior_mesh_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoInteriorMesh.cpp")
                target_sources(rank-two-interior-mesh-regression PRIVATE src/spin/RankTwoInteriorMesh.cpp
                    src/spin/RankTwoMovingInterior.cpp src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp
                    src/spin/RankTwoFiniteArcCoverage.cpp src/spin/RankTwoArcRefinement.cpp
                    src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                    src/spin/RankTwoFiniteTriangulation.cpp src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp)
                if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoMovingMeshEdge.cpp")
                    target_sources(rank-two-interior-mesh-regression PRIVATE src/spin/RankTwoMovingMeshEdge.cpp)
                endif()
            endif()
            target_include_directories(rank-two-interior-mesh-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-interior-mesh-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-interior-mesh-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-interior-mesh COMMAND rank-two-interior-mesh-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-interior-mesh-evidence.txt")
            set_tests_properties(spin-faster-rank-two-interior-mesh PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-finite-triangulation-regression tests/rank_two_finite_triangulation_regression.cpp
                src/spin/RankTwoFiniteTriangulation.cpp
                src/spin/RankTwoFinitePartition.cpp src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            add_executable(rank-two-moving-mesh-edge-regression tests/rank_two_moving_mesh_edge_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoMovingMeshEdge.cpp")
                target_sources(rank-two-moving-mesh-edge-regression PRIVATE src/spin/RankTwoMovingMeshEdge.cpp
                    src/spin/RankTwoInteriorMesh.cpp src/spin/RankTwoMovingInterior.cpp
                    src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp src/spin/RankTwoFiniteArcCoverage.cpp
                    src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                    src/spin/RankTwoFiniteTriangulation.cpp src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp)
            endif()
            target_include_directories(rank-two-moving-mesh-edge-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-moving-mesh-edge-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-moving-mesh-edge-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-moving-mesh-edge COMMAND rank-two-moving-mesh-edge-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-moving-mesh-edge-evidence.txt")
            set_tests_properties(spin-faster-rank-two-moving-mesh-edge PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-strict-material-mesh-regression tests/rank_two_strict_material_mesh_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>
                src/spin/RankTwoMaterialLifts.cpp src/spin/RankTwoMovingMeshEdge.cpp
                src/spin/RankTwoInteriorMesh.cpp src/spin/RankTwoMovingInterior.cpp
                src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp src/spin/RankTwoFiniteArcCoverage.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoFiniteTriangulation.cpp src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoStrictMaterialMesh.cpp")
                target_sources(rank-two-strict-material-mesh-regression PRIVATE src/spin/RankTwoStrictMaterialMesh.cpp)
            endif()
            target_include_directories(rank-two-strict-material-mesh-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-strict-material-mesh-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-strict-material-mesh-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-strict-material-mesh COMMAND rank-two-strict-material-mesh-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-strict-material-mesh-evidence.txt")
            add_test(NAME spin-faster-rank-two-strict-material-annulus COMMAND rank-two-strict-material-mesh-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-strict-material-annulus-evidence.txt" annulus)
            set_tests_properties(spin-faster-rank-two-strict-material-mesh spin-faster-rank-two-strict-material-annulus PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;material-thickness;regression" RUN_SERIAL TRUE TIMEOUT 500)
            add_executable(rank-two-interior-inventory-regression tests/rank_two_interior_inventory_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>
                src/spin/RankTwoMaterialLifts.cpp src/spin/RankTwoMovingMeshEdge.cpp
                src/spin/RankTwoStrictMaterialMesh.cpp src/spin/RankTwoInteriorMesh.cpp src/spin/RankTwoMovingInterior.cpp
                src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp src/spin/RankTwoFiniteArcCoverage.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoFiniteTriangulation.cpp src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoInteriorInventory.cpp")
                target_sources(rank-two-interior-inventory-regression PRIVATE src/spin/RankTwoInteriorInventory.cpp)
            endif()
            target_include_directories(rank-two-interior-inventory-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-interior-inventory-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-interior-inventory-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-interior-inventory COMMAND rank-two-interior-inventory-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-interior-inventory-evidence.txt")
            set_tests_properties(spin-faster-rank-two-interior-inventory PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;native-inventory;regression" RUN_SERIAL TRUE TIMEOUT 240)
            add_executable(rank-two-whole-interior-inventory-regression tests/rank_two_whole_interior_inventory_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>
                src/spin/RankTwoInteriorInventory.cpp src/spin/RankTwoMovingInterior.cpp)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoWholeInteriorInventory.cpp")
                target_sources(rank-two-whole-interior-inventory-regression PRIVATE src/spin/RankTwoWholeInteriorInventory.cpp)
            endif()
            target_include_directories(rank-two-whole-interior-inventory-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-whole-interior-inventory-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-whole-interior-inventory-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-whole-interior-inventory COMMAND rank-two-whole-interior-inventory-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-whole-interior-inventory-evidence.txt")
            set_tests_properties(spin-faster-rank-two-whole-interior-inventory PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;native-inventory;normal-owners;regression" RUN_SERIAL TRUE TIMEOUT 500)
            add_executable(rank-two-automatic-material-mesh-regression tests/rank_two_automatic_material_mesh_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>
                src/spin/RankTwoInteriorInventory.cpp src/spin/RankTwoStrictMaterialMesh.cpp
                src/spin/RankTwoMaterialLifts.cpp src/spin/RankTwoMovingMeshEdge.cpp
                src/spin/RankTwoInteriorMesh.cpp src/spin/RankTwoMovingInterior.cpp
                src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp src/spin/RankTwoFiniteArcCoverage.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoFiniteTriangulation.cpp src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoAutomaticMaterialMesh.cpp")
                target_sources(rank-two-automatic-material-mesh-regression PRIVATE src/spin/RankTwoAutomaticMaterialMesh.cpp)
            endif()
            target_include_directories(rank-two-automatic-material-mesh-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-automatic-material-mesh-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-automatic-material-mesh-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-automatic-material-mesh COMMAND rank-two-automatic-material-mesh-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-automatic-material-mesh-evidence.txt")
            set_tests_properties(spin-faster-rank-two-automatic-material-mesh PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;native-inventory;material-thickness;regression" RUN_SERIAL TRUE TIMEOUT 500)
            add_executable(rank-two-inventory-material-mesh-regression tests/rank_two_inventory_material_mesh_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>
                src/spin/RankTwoAutomaticMaterialMesh.cpp src/spin/RankTwoInteriorInventory.cpp src/spin/RankTwoStrictMaterialMesh.cpp
                src/spin/RankTwoMaterialLifts.cpp src/spin/RankTwoMovingMeshEdge.cpp
                src/spin/RankTwoInteriorMesh.cpp src/spin/RankTwoMovingInterior.cpp
                src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp src/spin/RankTwoFiniteArcCoverage.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoFiniteTriangulation.cpp src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp)
            target_include_directories(rank-two-inventory-material-mesh-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-inventory-material-mesh-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-inventory-material-mesh-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-inventory-material-mesh COMMAND rank-two-inventory-material-mesh-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-inventory-material-mesh-evidence.txt")
            set_tests_properties(spin-faster-rank-two-inventory-material-mesh PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;native-inventory;material-thickness;regression" RUN_SERIAL TRUE TIMEOUT 500)
            add_executable(rank-two-material-lifts-regression tests/rank_two_material_lifts_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/spin/RankTwoMaterialLifts.cpp")
                target_sources(rank-two-material-lifts-regression PRIVATE src/spin/RankTwoMaterialLifts.cpp
                    src/spin/RankTwoMovingMeshEdge.cpp src/spin/RankTwoInteriorMesh.cpp src/spin/RankTwoMovingInterior.cpp
                    src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp src/spin/RankTwoFiniteArcCoverage.cpp
                    src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                    src/spin/RankTwoFiniteTriangulation.cpp src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp)
            endif()
            target_include_directories(rank-two-material-lifts-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-material-lifts-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-material-lifts-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-material-lifts COMMAND rank-two-material-lifts-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-material-lifts-evidence.txt")
            set_tests_properties(spin-faster-rank-two-material-lifts PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-annulus-edge-repair-regression tests/rank_two_annulus_edge_repair_regression.cpp
                $<TARGET_OBJECTS:rank-two-registry-test-kernels>
                src/spin/RankTwoMaterialLifts.cpp src/spin/RankTwoMovingMeshEdge.cpp
                src/spin/RankTwoInteriorMesh.cpp src/spin/RankTwoMovingInterior.cpp
                src/spin/RankTwoCoveredFiniteMesh.cpp src/spin/RankTwoFiniteMeshScope.cpp src/spin/RankTwoFiniteArcCoverage.cpp
                src/spin/RankTwoArcRefinement.cpp src/spin/RankTwoRegistryEdge.cpp src/spin/RankTwoFinitePartition.cpp
                src/spin/RankTwoFiniteTriangulation.cpp src/spin/RankTwoChordIncidence.cpp src/spin/RankTwoRegistryRotation.cpp)
            target_include_directories(rank-two-annulus-edge-repair-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-annulus-edge-repair-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-annulus-edge-repair-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-annulus-edge-repair COMMAND rank-two-annulus-edge-repair-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-annulus-edge-repair-evidence.txt")
            set_tests_properties(spin-faster-rank-two-annulus-edge-repair PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_test(NAME spin-faster-rank-two-annulus-asymmetric-edge COMMAND rank-two-annulus-edge-repair-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-annulus-asymmetric-edge-evidence.txt" asymmetric)
            set_tests_properties(spin-faster-rank-two-annulus-asymmetric-edge PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_test(NAME spin-faster-rank-two-annulus-all-edges COMMAND rank-two-annulus-edge-repair-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-annulus-all-edges-evidence.txt" all)
            set_tests_properties(spin-faster-rank-two-annulus-all-edges PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            target_include_directories(rank-two-finite-triangulation-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-finite-triangulation-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-finite-triangulation-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-finite-triangulation COMMAND rank-two-finite-triangulation-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-finite-triangulation-evidence.txt")
            set_tests_properties(spin-faster-rank-two-finite-triangulation PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-finite-predicate-regression tests/rank_two_finite_predicate_regression.cpp
                src/spin/RankTwoAlgebraic.cpp src/spin/RankTwoSource.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
            target_include_directories(rank-two-finite-predicate-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-finite-predicate-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-finite-predicate-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-finite-predicate COMMAND rank-two-finite-predicate-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-finite-predicate-evidence.txt")
            set_tests_properties(spin-faster-rank-two-finite-predicate PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            foreach(registry_test IN ITEMS rank-two-registry-edge-regression rank-two-registry-rotation-regression
                    rank-two-arc-refinement-regression rank-two-regular-registry-regression)
                get_target_property(registry_test_sources ${registry_test} SOURCES)
                foreach(shared_source IN LISTS rank_two_registry_shared_sources)
                    if(NOT shared_source IN_LIST registry_test_sources)
                        message(FATAL_ERROR "${registry_test} lacks expected shared source ${shared_source}")
                    endif()
                    list(REMOVE_ITEM registry_test_sources "${shared_source}")
                endforeach()
                set_property(TARGET ${registry_test} PROPERTY SOURCES "${registry_test_sources}")
                target_sources(${registry_test} PRIVATE $<TARGET_OBJECTS:rank-two-registry-test-kernels>)
            endforeach()
            target_sources(rank-two-arc-registry-regression PRIVATE src/spin/RankTwoCarrier.cpp
                src/spin/RankTwoArrangement.cpp src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp
                src/spin/RankTwoArcSampling.cpp src/spin/QuadricCoefficients.cpp
                src/manufacturing/SlicerValidation.cpp)
            target_sources(rank-two-arc-registry-regression PRIVATE src/spin/RankTwoArcRegistry.cpp)
            target_sources(rank-two-arc-registry-regression PRIVATE src/spin/RankTwoArcStrip.cpp
                src/spin/RankTwoEndpointArc.cpp src/spin/RankTwoEndpointFan.cpp src/spin/RankTwoEndpointNeighborhood.cpp)
            target_include_directories(rank-two-arc-registry-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-arc-registry-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-arc-registry-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-arc-registry COMMAND rank-two-arc-registry-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-arc-registry-evidence.txt")
            add_test(NAME spin-faster-rank-two-arc-registry-profile COMMAND rank-two-arc-registry-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-arc-registry-profile-evidence.txt" --trace)
            set_tests_properties(spin-faster-rank-two-arc-registry-profile PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;diagnostic")
            set_tests_properties(spin-faster-rank-two-arc-registry PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-endpoint-fan-regression tests/rank_two_endpoint_fan_regression.cpp)
            target_sources(rank-two-endpoint-fan-regression PRIVATE src/spin/RankTwoCarrier.cpp
                src/spin/RankTwoArrangement.cpp src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp
                src/spin/RankTwoArcStrip.cpp src/spin/RankTwoEndpointArc.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
            target_sources(rank-two-endpoint-fan-regression PRIVATE src/spin/RankTwoEndpointFan.cpp)
            target_include_directories(rank-two-endpoint-fan-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-endpoint-fan-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-endpoint-fan-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-endpoint-fan COMMAND rank-two-endpoint-fan-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-endpoint-fan-evidence.txt")
            set_tests_properties(spin-faster-rank-two-endpoint-fan PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            target_sources(rank-two-arc-sampling-regression PRIVATE src/spin/RankTwoCarrier.cpp
                src/spin/RankTwoArrangement.cpp src/spin/RankTwoAlgebraic.cpp
                src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
            target_sources(rank-two-arc-sampling-regression PRIVATE src/spin/RankTwoArcSampling.cpp
                src/spin/RankTwoArcStrip.cpp)
            target_include_directories(rank-two-arc-sampling-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-arc-sampling-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-arc-sampling-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-arc-sampling COMMAND rank-two-arc-sampling-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-arc-sampling-evidence.txt")
            add_test(NAME spin-faster-rank-two-arc-sampling-profile COMMAND rank-two-arc-sampling-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-arc-sampling-profile-evidence.txt" --trace)
            add_test(NAME spin-faster-rank-two-arc-strip-mapping-profile COMMAND rank-two-arc-sampling-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-arc-strip-mapping-profile-evidence.txt" --strip-profile)
            set_tests_properties(spin-faster-rank-two-arc-strip-mapping-profile PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;diagnostic")
            set_tests_properties(spin-faster-rank-two-arc-sampling-profile PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;diagnostic")
            set_tests_properties(spin-faster-rank-two-arc-sampling PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            set_tests_properties(spin-faster-rank-two-carrier-kernels PROPERTIES
                LABELS "spin-it-faster;exact;rank-two;local-kernel;regression")
            add_executable(rank-two-arrangement-regression tests/rank_two_arrangement_regression.cpp)
            target_sources(rank-two-arrangement-regression PRIVATE src/spin/RankTwoArrangement.cpp
                src/spin/RankTwoAlgebraic.cpp src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
            target_include_directories(rank-two-arrangement-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-arrangement-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-arrangement-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-arrangement COMMAND rank-two-arrangement-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-global-arrangement-evidence.txt")
            add_executable(rank-two-normal-mesh-regression tests/rank_two_normal_mesh_regression.cpp)
            target_sources(rank-two-normal-mesh-regression PRIVATE src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
            target_include_directories(rank-two-normal-mesh-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-normal-mesh-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-normal-mesh-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-normal-mesh COMMAND rank-two-normal-mesh-regression
                "${CMAKE_CURRENT_SOURCE_DIR}/.superpowers/sdd/2026-09-29-seven-day-faster-sprint/rank-two-normal-boundary/diagnostic-01/bundle/normals.tsv"
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-normal-mesh-evidence.txt")
            add_executable(rank-two-whole-arrangement-regression tests/rank_two_whole_arrangement_regression.cpp)
            target_sources(rank-two-whole-arrangement-regression PRIVATE src/spin/RankTwoArrangement.cpp
                src/spin/RankTwoSettings.cpp
                src/spin/RankTwoAlgebraic.cpp src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
            target_include_directories(rank-two-whole-arrangement-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-whole-arrangement-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-whole-arrangement-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            add_test(NAME spin-faster-rank-two-whole-arrangement COMMAND rank-two-whole-arrangement-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-whole-arrangement-evidence.txt")
            add_test(NAME spin-faster-rank-two-whole-arrangement-profile COMMAND rank-two-whole-arrangement-regression
                "${CMAKE_CURRENT_BINARY_DIR}/rank-two-whole-arrangement-profile-evidence.txt" --trace)
            set_tests_properties(spin-faster-rank-two-whole-arrangement-profile PROPERTIES LABELS "spin-it-faster;exact;rank-two;diagnostic")
            set_tests_properties(spin-faster-rank-two-backend spin-faster-rank-two-source-input
                spin-faster-rank-two-algebraic spin-faster-rank-two-arrangement spin-faster-rank-two-normal-mesh
                spin-faster-rank-two-whole-arrangement
                PROPERTIES LABELS "spin-it-faster;exact;rank-two;regression")
            add_executable(rank-two-frozen-edges-regression tests/rank_two_frozen_edges_regression.cpp)
            target_sources(rank-two-frozen-edges-regression PRIVATE src/spin/RankTwoArrangement.cpp
                src/spin/RankTwoAlgebraic.cpp src/spin/RankTwoNormalMesh.cpp src/spin/RankTwoSource.cpp
                src/spin/QuadricCoefficients.cpp src/manufacturing/SlicerValidation.cpp)
            target_include_directories(rank-two-frozen-edges-regression PRIVATE
                "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include")
            target_compile_features(rank-two-frozen-edges-regression PRIVATE cxx_std_17)
            target_link_libraries(rank-two-frozen-edges-regression PRIVATE igl::core CGAL::CGAL CGAL::CGAL_Core)
            foreach(rank_two_scope first all stages scopes state-0 state-1 state-2 state-18 state-56 state-57 state-58 state-59)
                add_test(NAME spin-faster-rank-two-frozen-edges-${rank_two_scope}
                    COMMAND rank-two-frozen-edges-regression
                    "${CMAKE_CURRENT_SOURCE_DIR}/.superpowers/sdd/2026-09-29-seven-day-faster-sprint/rank-two-normal-boundary/diagnostic-01/bundle/normals.tsv"
                    "${CMAKE_CURRENT_SOURCE_DIR}/.superpowers/sdd/2026-09-29-seven-day-faster-sprint/rank-two-edge-certificate/green-03/bundle/native.tsv"
                    "${CMAKE_CURRENT_BINARY_DIR}/rank-two-frozen-edges-${rank_two_scope}.txt" ${rank_two_scope})
                set_tests_properties(spin-faster-rank-two-frozen-edges-${rank_two_scope}
                    PROPERTIES LABELS "spin-it-faster;exact;rank-two;frozen")
            endforeach()
            add_executable(quadric-projective-surface-regression tests/quadric_projective_surface_regression.cpp)
            target_link_libraries(quadric-projective-surface-regression PRIVATE make_it_stand::core)
            target_compile_features(quadric-projective-surface-regression PRIVATE cxx_std_17)
            add_test(NAME spin-faster-projective-surface COMMAND quadric-projective-surface-regression)
            add_executable(quadric-fixed-normal-paraboloid-regression tests/quadric_fixed_normal_paraboloid_regression.cpp)
            target_link_libraries(quadric-fixed-normal-paraboloid-regression PRIVATE make_it_stand::core)
            target_compile_features(quadric-fixed-normal-paraboloid-regression PRIVATE cxx_std_17)
            add_test(NAME spin-faster-fixed-normal-paraboloid COMMAND quadric-fixed-normal-paraboloid-regression)
            set_tests_properties(spin-faster-projective-surface PROPERTIES LABELS "spin-it-faster;geometry;derivative;regression")
        endif()
        target_compile_features(quadric-source-motion-regression PRIVATE cxx_std_17)
        set_target_properties(
            quadric-source-motion-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            quadric-source-motion-regression
            PRIVATE make_it_stand::faster_fullrank
        )
        add_test(
            NAME spin-faster-source-motion
            COMMAND "$<TARGET_FILE:quadric-source-motion-regression>"
        )
        set_tests_properties(
            spin-faster-source-motion
            PROPERTIES LABELS "spin-it;spin-it-faster;geometry;derivative;regression"
        )
        if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
            add_test(NAME spin-faster-rational-source-motion-reuse
                COMMAND "$<TARGET_FILE:quadric-source-motion-regression>" --reuse-rational-source)
            set_tests_properties(spin-faster-rational-source-motion-reuse PROPERTIES
                LABELS "spin-it-faster;geometry;derivative;regression" TIMEOUT 180)
        endif()

        add_executable(
            truncation-root-certificate-regression
            tests/truncation_root_certificate_regression.cpp
        )
        target_compile_features(
            truncation-root-certificate-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            truncation-root-certificate-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            truncation-root-certificate-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-truncation-root-certificate
            COMMAND "$<TARGET_FILE:truncation-root-certificate-regression>"
        )
        set_tests_properties(
            spin-faster-truncation-root-certificate
            PROPERTIES LABELS "spin-it;spin-it-faster;geometry;derivative;regression"
        )
        if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
            add_test(NAME spin-faster-rational-truncation-point
                COMMAND "$<TARGET_FILE:truncation-root-certificate-regression>" --rational-point)
            set_tests_properties(spin-faster-rational-truncation-point PROPERTIES
                LABELS "spin-it;spin-it-faster;geometry;derivative;regression" TIMEOUT 180)
            add_test(NAME spin-faster-rational-truncation-near-cone
                COMMAND "$<TARGET_FILE:truncation-root-certificate-regression>" --near-cone-rational-point)
            set_tests_properties(spin-faster-rational-truncation-near-cone PROPERTIES
                LABELS "spin-it;spin-it-faster;geometry;derivative;regression" TIMEOUT 180)
        endif()

        add_executable(
            quadric-moment-evaluator-regression
            tests/quadric_moment_evaluator_regression.cpp
        )
        add_executable(quadric-fixed-normal-request-hash-regression
            tests/quadric_fixed_normal_request_hash_regression.cpp)
        target_compile_features(quadric-fixed-normal-request-hash-regression
            PRIVATE cxx_std_17)
        set_target_properties(quadric-fixed-normal-request-hash-regression
            PROPERTIES CXX_EXTENSIONS OFF)
        target_link_libraries(quadric-fixed-normal-request-hash-regression
            PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-fixed-normal-request-hash
            COMMAND "$<TARGET_FILE:quadric-fixed-normal-request-hash-regression>")
        set_tests_properties(spin-faster-fixed-normal-request-hash PROPERTIES
            LABELS "spin-it;spin-it-faster;identity;regression")
        add_executable(quadric-discrete-moment-regression
            tests/quadric_discrete_moment_regression.cpp)
        target_compile_features(quadric-discrete-moment-regression PRIVATE cxx_std_17)
        target_link_libraries(quadric-discrete-moment-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-discrete-moments
            COMMAND "$<TARGET_FILE:quadric-discrete-moment-regression>")
        set_tests_properties(spin-faster-discrete-moments PROPERTIES
            LABELS "spin-it;spin-it-faster;geometry;regression")
        add_executable(quadric-stratum-tangent-regression
            tests/quadric_stratum_tangent_regression.cpp)
        target_compile_features(quadric-stratum-tangent-regression PRIVATE cxx_std_17)
        target_link_libraries(quadric-stratum-tangent-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-stratum-tangent
            COMMAND "$<TARGET_FILE:quadric-stratum-tangent-regression>")
        set_tests_properties(spin-faster-stratum-tangent PROPERTIES
            LABELS "spin-it;spin-it-faster;geometry;derivative;regression")
        add_executable(plane-family-regression
            tests/plane_family_regression.cpp)
        target_compile_features(plane-family-regression PRIVATE cxx_std_17)
        target_link_libraries(plane-family-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-plane-family
            COMMAND "$<TARGET_FILE:plane-family-regression>")
        set_tests_properties(spin-faster-plane-family PROPERTIES
            LABELS "spin-it;spin-it-faster;plane;geometry;derivative;regression")
        add_executable(plane-static-balance-regression
            tests/plane_static_balance_regression.cpp)
        target_compile_features(plane-static-balance-regression PRIVATE cxx_std_17)
        target_link_libraries(plane-static-balance-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-static-plane-balance
            COMMAND "$<TARGET_FILE:plane-static-balance-regression>")
        set_tests_properties(spin-faster-static-plane-balance PROPERTIES
            LABELS "spin-it;spin-it-faster;plane;solver;synthetic-regression" TIMEOUT 3600)
        add_executable(buoyant-plane-balance-regression
            tests/buoyant_plane_balance_regression.cpp)
        target_compile_features(buoyant-plane-balance-regression PRIVATE cxx_std_17)
        target_link_libraries(buoyant-plane-balance-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-buoyant-plane-balance
            COMMAND "$<TARGET_FILE:buoyant-plane-balance-regression>")
        set_tests_properties(spin-faster-buoyant-plane-balance PROPERTIES
            LABELS "spin-it;spin-it-faster;plane;buoyancy;solver;synthetic-regression"
            TIMEOUT 3600)
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
            PRIVATE make_it_stand::faster_fullrank
        )
        if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
            add_test(NAME spin-faster-rational-source-evaluator
                COMMAND "$<TARGET_FILE:quadric-moment-evaluator-regression>" --rational-source)
            set_tests_properties(spin-faster-rational-source-evaluator PROPERTIES
                LABELS "spin-it-faster;geometry;exact;precision;derivative;regression" TIMEOUT 180)
            add_test(NAME spin-faster-rational-hyperboloid-evaluator
                COMMAND "$<TARGET_FILE:quadric-moment-evaluator-regression>" --rational-hyperboloid)
            set_tests_properties(spin-faster-rational-hyperboloid-evaluator PROPERTIES
                LABELS "spin-it-faster;geometry;exact;precision;derivative;regression" TIMEOUT 180)
        endif()
        add_test(
            NAME spin-faster-quadric-moments
            COMMAND "$<TARGET_FILE:quadric-moment-evaluator-regression>"
        )
        set_tests_properties(
            spin-faster-quadric-moments
            PROPERTIES LABELS "spin-it;spin-it-faster;geometry;regression"
        )
        add_test(
            NAME spin-faster-quadric-degenerate-boolean-witness
            COMMAND "$<TARGET_FILE:quadric-moment-evaluator-regression>"
                --degenerate-boolean-witness
        )
        set_tests_properties(
            spin-faster-quadric-degenerate-boolean-witness
            PROPERTIES LABELS "spin-it-faster;geometry;diagnostic;regression"
        )
        add_test(
            NAME spin-faster-quadric-source-winding
            COMMAND "$<TARGET_FILE:quadric-moment-evaluator-regression>"
                --source-winding-regression
        )
        set_tests_properties(
            spin-faster-quadric-source-winding
            PROPERTIES LABELS "spin-it-faster;geometry;boolean;regression"
                TIMEOUT 600
        )
        add_test(
            NAME spin-faster-quadric-source-interior-side-witness
            COMMAND "$<TARGET_FILE:quadric-moment-evaluator-regression>"
                --source-interior-side-witness
        )
        set_tests_properties(
            spin-faster-quadric-source-interior-side-witness
            PROPERTIES LABELS "spin-it-faster;geometry;boolean;regression"
                TIMEOUT 600
        )

        add_executable(
            faster-theorem-conditions-regression
            tests/faster_theorem_conditions_regression.cpp
        )
        target_compile_features(
            faster-theorem-conditions-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            faster-theorem-conditions-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            faster-theorem-conditions-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-theorem-conditions
            COMMAND "$<TARGET_FILE:faster-theorem-conditions-regression>"
        )
        set_tests_properties(
            spin-faster-theorem-conditions
            PROPERTIES LABELS "spin-it;spin-it-faster;theorem;regression"
        )

        add_executable(
            faster-static-stability-contract-regression
            tests/faster_static_stability_contract_regression.cpp
        )
        target_compile_features(
            faster-static-stability-contract-regression
            PRIVATE cxx_std_17
        )
        set_target_properties(
            faster-static-stability-contract-regression
            PROPERTIES CXX_EXTENSIONS OFF
        )
        target_link_libraries(
            faster-static-stability-contract-regression
            PRIVATE make_it_stand::core
        )
        add_test(
            NAME spin-faster-static-stability-contract
            COMMAND "$<TARGET_FILE:faster-static-stability-contract-regression>"
        )
        set_tests_properties(
            spin-faster-static-stability-contract
            PROPERTIES LABELS "spin-it;spin-it-faster;static;application;regression"
        )

        if(MAKE_IT_STAND_BUILD_SPIN_IT_2014)
            add_executable(spin2014-authority-resource-admission-regression
                tests/replay_resource_admission_test.cpp)
            target_compile_features(spin2014-authority-resource-admission-regression PRIVATE cxx_std_17)
            add_test(NAME spin2014-authority-resource-admission
                COMMAND "$<TARGET_FILE:spin2014-authority-resource-admission-regression>")
            set_tests_properties(spin2014-authority-resource-admission PROPERTIES
                LABELS "spin2014;fixed-shape;precision;regression" TIMEOUT 30)
            add_executable(spin2014-authority-saved-step-replay
                tests/spin2014_authority_saved_step_replay.cpp)
            target_compile_features(spin2014-authority-saved-step-replay PRIVATE cxx_std_17)
            target_link_libraries(spin2014-authority-saved-step-replay PRIVATE make_it_stand::core)
            add_executable(spin2014-authority-objective-regression
                tests/spin2014_authority_objective_regression.cpp)
            target_compile_features(spin2014-authority-objective-regression PRIVATE cxx_std_17)
            target_link_libraries(spin2014-authority-objective-regression PRIVATE make_it_stand::core)
            add_test(NAME spin2014-authority-objective
                COMMAND "$<TARGET_FILE:spin2014-authority-objective-regression>")
            set_tests_properties(spin2014-authority-objective PROPERTIES
                LABELS "spin2014;spin-it-faster;fixed-shape;precision;local-kernel;regression" TIMEOUT 180)
            add_executable(spin2014-authority-newton-state-regression
                tests/spin2014_authority_newton_state_regression.cpp)
            target_compile_features(spin2014-authority-newton-state-regression PRIVATE cxx_std_17)
            target_link_libraries(spin2014-authority-newton-state-regression PRIVATE make_it_stand::core)
            add_test(NAME spin2014-authority-newton-state
                COMMAND "$<TARGET_FILE:spin2014-authority-newton-state-regression>")
            set_tests_properties(spin2014-authority-newton-state PROPERTIES
                LABELS "spin2014;spin-it-faster;fixed-shape;precision;local-kernel;regression" TIMEOUT 180)
            add_executable(spin2014-authority-newton-direction-regression
                tests/spin2014_authority_newton_direction_regression.cpp)
            target_compile_features(spin2014-authority-newton-direction-regression PRIVATE cxx_std_17)
            target_link_libraries(spin2014-authority-newton-direction-regression PRIVATE make_it_stand::core)
            add_test(NAME spin2014-authority-newton-direction
                COMMAND "$<TARGET_FILE:spin2014-authority-newton-direction-regression>")
            set_tests_properties(spin2014-authority-newton-direction PROPERTIES
                LABELS "spin2014;spin-it-faster;fixed-shape;precision;local-kernel;regression" TIMEOUT 180)
            add_executable(spin2014-authority-state-decision-regression
                tests/spin2014_authority_state_decision_regression.cpp)
            target_compile_features(spin2014-authority-state-decision-regression PRIVATE cxx_std_17)
            target_link_libraries(spin2014-authority-state-decision-regression PRIVATE make_it_stand::core)
            add_test(NAME spin2014-authority-state-decision
                COMMAND "$<TARGET_FILE:spin2014-authority-state-decision-regression>")
            set_tests_properties(spin2014-authority-state-decision PROPERTIES
                LABELS "spin2014;spin-it-faster;fixed-shape;precision;local-kernel;regression" TIMEOUT 180)
            if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS AND MAKE_IT_STAND_BUILD_PYTHON)
                add_test(NAME spin2014-authority-python-bindings
                    COMMAND "${Python3_EXECUTABLE}" -B -m unittest discover
                        -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p test_authority_python_bindings.py -v)
                set_tests_properties(spin2014-authority-python-bindings PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;precision;python;local-kernel"
                    TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;AUTHORITY_BINDING_TEST_MODULE=_make_it_stand;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
                add_test(NAME spin2014-authority-durable-events
                    COMMAND "${Python3_EXECUTABLE}" -B -m unittest discover
                        -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p test_authority_durable_events.py -v)
                set_tests_properties(spin2014-authority-durable-events PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;precision;python;durable-events"
                    TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
                add_test(NAME spin2014-authority-independent-equations
                    COMMAND "${Python3_EXECUTABLE}" -B -m unittest discover
                        -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p test_spin2014_authority_numerical_checker.py -v)
                set_tests_properties(spin2014-authority-independent-equations PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;precision;python;independent-checker"
                    TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
                add_test(NAME spin2014-authority-independent-hessians
                    COMMAND "${Python3_EXECUTABLE}" -B -m unittest discover
                        -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p test_authority_independent_hessians.py -v)
                set_tests_properties(spin2014-authority-independent-hessians PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;precision;python;independent-checker"
                    TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
                add_test(NAME spin2014-authority-exact-integration
                    COMMAND "${Python3_EXECUTABLE}" -B -m unittest discover
                        -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p test_independent_mass_moments.py -v)
                set_tests_properties(spin2014-authority-exact-integration PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;precision;python;independent-checker"
                    TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
                add_test(NAME spin2014-authority-event-journal-checker
                    COMMAND "${Python3_EXECUTABLE}" -B -m unittest discover
                        -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p test_authority_event_journal.py -v)
                set_tests_properties(spin2014-authority-event-journal-checker PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;precision;python;independent-checker"
                    TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
                add_test(NAME spin2014-authority-function-runner
                    COMMAND "${Python3_EXECUTABLE}" -B -m unittest discover
                        -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p test_authority_function_runner.py -v)
                set_tests_properties(spin2014-authority-function-runner PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;precision;python;functional-runner"
                    TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV};PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
            endif()
            add_executable(spin2014-quadric-moment-adapter-regression
                tests/spin2014_quadric_moment_adapter_regression.cpp)
            target_compile_features(spin2014-quadric-moment-adapter-regression PRIVATE cxx_std_17)
            target_link_libraries(spin2014-quadric-moment-adapter-regression PRIVATE make_it_stand::core)
            foreach(contract formulas kernel)
                add_test(NAME spin2014-quadric-moment-${contract}
                    COMMAND "$<TARGET_FILE:spin2014-quadric-moment-adapter-regression>" --${contract})
                set_tests_properties(spin2014-quadric-moment-${contract} PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;moment-adapter;regression")
            endforeach()
            if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
                add_executable(spin2014-authority-newton-trial-regression
                    tests/spin2014_authority_newton_trial_regression.cpp)
                target_compile_features(spin2014-authority-newton-trial-regression PRIVATE cxx_std_17)
                target_link_libraries(spin2014-authority-newton-trial-regression PRIVATE make_it_stand::core)
                add_test(NAME spin2014-authority-newton-trial
                    COMMAND "$<TARGET_FILE:spin2014-authority-newton-trial-regression>")
                set_tests_properties(spin2014-authority-newton-trial PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;precision;actual-material;regression" TIMEOUT 180)
                add_executable(quadric-authority-material-context-regression
                    tests/quadric_authority_material_context_regression.cpp)
                target_compile_features(quadric-authority-material-context-regression PRIVATE cxx_std_17)
                target_link_libraries(quadric-authority-material-context-regression PRIVATE make_it_stand::core)
                add_test(NAME quadric-authority-material-context
                    COMMAND "$<TARGET_FILE:quadric-authority-material-context-regression>")
                set_tests_properties(quadric-authority-material-context PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;precision;exact-material;regression" TIMEOUT 180)
                if(MAKE_IT_STAND_BUILD_PYTHON)
                    add_test(NAME spin2014-quadric-python-interface
                        COMMAND "${Python_EXECUTABLE}" -m unittest discover
                            -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                            -p test_spin2014_quadric_interface.py -v)
                    set_tests_properties(spin2014-quadric-python-interface PROPERTIES
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                        LABELS "spin2014;spin-it-faster;fixed-shape;geometry-adapter;python"
                        TIMEOUT 600)
                    add_test(NAME spin2014-quadric-python-rational-authority
                        COMMAND "${Python_EXECUTABLE}" -m unittest discover
                            -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                            -p test_rational_source_production.py -v)
                    set_tests_properties(spin2014-quadric-python-rational-authority PROPERTIES
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                        LABELS "spin2014;spin-it-faster;fixed-shape;exact;python;independent-check"
                        TIMEOUT 180)
                    add_test(NAME faster-exact-input-admission
                        COMMAND "${Python_EXECUTABLE}" -m unittest discover
                            -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                            -p test_faster_exact_input_admission.py -v)
                    set_tests_properties(faster-exact-input-admission PROPERTIES
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                        LABELS "spin-it-faster;exact;python;regression"
                        TIMEOUT 180)
                    add_test(NAME spin2014-quadric-native-stage-observer
                        COMMAND "${Python_EXECUTABLE}" -m unittest discover
                            -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                            -p test_quadric_stage_observer.py -v)
                    set_tests_properties(spin2014-quadric-native-stage-observer PROPERTIES
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                        LABELS "spin2014;spin-it-faster;fixed-shape;diagnostics;python"
                        TIMEOUT 180)
                    add_test(NAME spin2014-quadric-five-point-protocol-arithmetic
                        COMMAND "${Python_EXECUTABLE}" -m unittest discover
                            -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                            -p test_quadric_five_point_protocol.py -v)
                    add_test(NAME spin2014-quadric-point-job-manifest
                        COMMAND "${Python_EXECUTABLE}" -m unittest discover
                            -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                            -p test_quadric_point_job_manifest.py -v)
                    set_tests_properties(spin2014-quadric-point-job-manifest PROPERTIES
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                        LABELS "spin2014;spin-it-faster;fixed-shape;independent-check;accounting"
                        TIMEOUT 30)
                    add_test(NAME spin2014-quadric-independent-analytic-branch
                        COMMAND "${Python_EXECUTABLE}" -m unittest discover
                            -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                            -p test_exact_quadric_branch.py -v)
                    set_tests_properties(spin2014-quadric-independent-analytic-branch PROPERTIES
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                        LABELS "spin2014;spin-it-faster;fixed-shape;independent-check;analytic"
                        TIMEOUT 30)
                    add_test(NAME spin2014-quadric-exact-coplanar-patch
                        COMMAND "${Python_EXECUTABLE}" -m unittest discover
                            -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                            -p test_exact_coplanar_patch.py -v)
                    set_tests_properties(spin2014-quadric-exact-coplanar-patch PROPERTIES
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                        LABELS "spin2014;spin-it-faster;fixed-shape;exact;independent-check"
                        TIMEOUT 30)
                    add_test(NAME spin2014-quadric-compact-rational-source
                        COMMAND "${Python_EXECUTABLE}" -m unittest discover
                            -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                            -p test_compact_rational_source.py -v)
                    set_tests_properties(spin2014-quadric-compact-rational-source PROPERTIES
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                        LABELS "spin2014;spin-it-faster;fixed-shape;exact;precision;python"
                        TIMEOUT 180)
                    set_tests_properties(spin2014-quadric-five-point-protocol-arithmetic PROPERTIES
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                        LABELS "spin2014;spin-it-faster;fixed-shape;independent-check;arithmetic"
                        TIMEOUT 30)
                    add_test(NAME spin2014-quadric-python-rational-hyperboloid-persistence
                        COMMAND "${Python_EXECUTABLE}" -m unittest discover
                            -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
                            -p test_rational_hyperboloid_persistence.py -v)
                    set_tests_properties(spin2014-quadric-python-rational-hyperboloid-persistence PROPERTIES
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                        LABELS "spin2014;spin-it-faster;fixed-shape;exact;python;independent-check"
                        TIMEOUT 600)
                    foreach(module test_spin2014_quadric_frozen_batch test_spin2014_static_saved_material
                                   test_spin2014_static_material_output)
                        add_test(NAME spin2014-quadric-python-${module}
                            COMMAND "${Python_EXECUTABLE}" -m unittest discover
                                -s "${CMAKE_CURRENT_SOURCE_DIR}/tests/python" -p "${module}.py" -v)
                        set_tests_properties(spin2014-quadric-python-${module} PROPERTIES
                            ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                            LABELS "spin2014;spin-it-faster;fixed-shape;static;python;independent-check"
                            TIMEOUT 180)
                    endforeach()
                endif()
                add_executable(spin2014-quadric-request-regression tests/spin2014_quadric_request_regression.cpp)
                target_compile_features(spin2014-quadric-request-regression PRIVATE cxx_std_17)
                target_link_libraries(spin2014-quadric-request-regression PRIVATE make_it_stand::spin2014_quadric)
                foreach(contract admission real-call rational-complement)
                    add_test(NAME spin2014-quadric-request-${contract}
                        COMMAND "$<TARGET_FILE:spin2014-quadric-request-regression>" --${contract})
                    set_tests_properties(spin2014-quadric-request-${contract} PROPERTIES
                        LABELS "spin2014;spin-it-faster;fixed-shape;geometry-adapter;regression" TIMEOUT 120)
                endforeach()
                add_executable(spin2014-quadric-static-regression tests/spin2014_quadric_static_regression.cpp)
                target_compile_features(spin2014-quadric-static-regression PRIVATE cxx_std_17)
                target_link_libraries(spin2014-quadric-static-regression PRIVATE make_it_stand::spin2014_quadric)
                add_test(NAME spin2014-quadric-static-square-gate COMMAND "$<TARGET_FILE:spin2014-quadric-static-regression>")
                set_tests_properties(spin2014-quadric-static-square-gate PROPERTIES
                    LABELS "spin2014;spin-it-faster;fixed-shape;static;regression" TIMEOUT 120)
            endif()
        endif()

        add_executable(legacy-prediction-transaction-regression
            tests/legacy_prediction_transaction_regression.cpp)
        target_compile_features(legacy-prediction-transaction-regression PRIVATE cxx_std_17)
        target_link_libraries(legacy-prediction-transaction-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-legacy-transaction-contract
            COMMAND "$<TARGET_FILE:legacy-prediction-transaction-regression>" --contract)
    set_tests_properties(spin-faster-legacy-transaction-contract PROPERTIES
        LABELS "spin-it;spin-it-faster;legacy-transaction;regression")
    add_test(NAME spin-faster-legacy-transaction-disabled
        COMMAND "$<TARGET_FILE:legacy-prediction-transaction-regression>" --disabled)
    set_tests_properties(spin-faster-legacy-transaction-disabled PROPERTIES
        LABELS "spin-it;spin-it-faster;legacy-transaction;regression")

        add_executable(
            quadric-hollowing-solver-regression
            tests/quadric_hollowing_solver_regression.cpp
            tests/dual_monomial_newton_regression.cpp
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
        foreach(contract dual-pairing raw-stopping dual-coordinates
                         strict-positive-mu-restart
                         strict-nonpositive-stationary-restart)
            add_test(NAME spin-faster-newton-${contract}
                COMMAND "$<TARGET_FILE:quadric-hollowing-solver-regression>" --${contract})
            set_tests_properties(spin-faster-newton-${contract} PROPERTIES
                LABELS "spin-it;spin-it-faster;solver;regression")
        endforeach()
        set_tests_properties(
            spin-faster-newton
            PROPERTIES LABELS "spin-it;spin-it-faster;solver;regression"
        )

        add_executable(spin-moment-subspace-regression
            tests/spin_moment_subspace_regression.cpp)
        target_compile_features(spin-moment-subspace-regression PRIVATE cxx_std_17)
        target_link_libraries(spin-moment-subspace-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-subspace-linearization
            COMMAND "$<TARGET_FILE:spin-moment-subspace-regression>")
        set_tests_properties(spin-faster-subspace-linearization PROPERTIES
            LABELS "spin-it;spin-it-faster;solver;regression")

        add_executable(
            hydrostatic-geometry-regression
            tests/hydrostatic_geometry_regression.cpp
        )
        target_compile_features(hydrostatic-geometry-regression PRIVATE cxx_std_17)
        set_target_properties(hydrostatic-geometry-regression PROPERTIES CXX_EXTENSIONS OFF)
        target_link_libraries(hydrostatic-geometry-regression PRIVATE make_it_stand::core)
        add_test(NAME spin-faster-hydrostatic-geometry
                 COMMAND "$<TARGET_FILE:hydrostatic-geometry-regression>")
        set_tests_properties(spin-faster-hydrostatic-geometry
            PROPERTIES LABELS "spin-it;spin-it-faster;hydrostatics;regression")

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
            "-DMAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT=${CMAKE_CURRENT_BINARY_DIR}"
            "-DMAKE_IT_STAND_DEPENDENCY_TEST_HIGHS_SOURCE_DIR=${FETCHCONTENT_SOURCE_DIR_HIGHS}"
            "-DMAKE_IT_STAND_DEPENDENCY_TEST_HIGHS_REVISION=${MAKE_IT_STAND_HIGHS_GIT_TAG}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/assert_spin_it_2014_dependency.cmake"
    )
    add_test(
        NAME spin-2014-dependency-isolation
        COMMAND "${CMAKE_COMMAND}"
            "-DTEST_BINARY_ROOT=${CMAKE_CURRENT_BINARY_DIR}"
            "-DDEPENDENCY_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/assert_spin_it_2014_dependency.cmake"
            "-DALLOCATOR_HELPER=${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/allocate_retained_test_directory.cmake"
            "-DTEST_HIGHS_SOURCE_DIR=${FETCHCONTENT_SOURCE_DIR_HIGHS}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/assert_spin_it_2014_dependency_isolation.cmake"
    )
    set_tests_properties(
        spin-faster-dependency-contract
        PROPERTIES
            LABELS "spin-it;spin-it-faster;build-system"
            TIMEOUT 180
    )
    set_tests_properties(
        spin-2014-dependency-contract spin-2014-dependency-isolation
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
                TIMEOUT 300
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
        )

        add_test(
            NAME faster-section-6-6-exact-contract
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_faster_section_6_6_exact
                -v
        )
        set_tests_properties(
            faster-section-6-6-exact-contract
            PROPERTIES
                LABELS "paper;python;spin-it;spin-it-faster;exact-rational"
                TIMEOUT 30
                ENVIRONMENT
                    "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
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
            NAME independent-exact-feasibility-contract
            COMMAND "${Python3_EXECUTABLE}" -m unittest test_exact_feasibility -v
        )
        set_tests_properties(
            independent-exact-feasibility-contract
            PROPERTIES
                LABELS "paper;python;reachability;independent-verification"
                TIMEOUT 30
                ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
        )

        add_test(
            NAME feature-thickness-source-map-contract
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_feature_thickness_attachment -v
        )
        set_tests_properties(
            feature-thickness-source-map-contract
            PROPERTIES
                LABELS "paper;python;spin-it-2014;verification;feature-thickness;synthetic-regression"
                TIMEOUT 30
                ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/tests/python"
        )

        add_test(
            NAME run-evidence-v2-contract
            COMMAND "${Python3_EXECUTABLE}" -m unittest
                test_run_evidence_v2
                test_native_build_receipt
                test_native_build_receipt_linkage
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

        add_test(
            NAME source-archive-member-audit-contract
            COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_source_archive_member_audit.py"
        )
        set_tests_properties(
            source-archive-member-audit-contract
            PROPERTIES
                LABELS "paper;python;spin-it;verification;archive-member-provenance"
                TIMEOUT 60
                ENVIRONMENT
                    "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python:${CMAKE_CURRENT_SOURCE_DIR}/tests/python;PYTEST_DISABLE_PLUGIN_AUTOLOAD=1"
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
                add_test(
                    NAME spin-it-2014-supplement-algorithms-contract
                    COMMAND "${Python3_EXECUTABLE}" -m unittest
                        test_spin2014_supplement_algorithms -v
                )
                set_tests_properties(
                    spin-it-2014-supplement-algorithms-contract
                    PROPERTIES
                        LABELS "paper;python;native;spin-it-2014;independent-verification;supplement"
                        TIMEOUT 60
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                )
                add_test(
                    NAME spin-it-2014-givens-oracle-contract
                    COMMAND "${Python3_EXECUTABLE}" -m unittest
                        test_spin2014_givens_oracle -v
                )
                set_tests_properties(
                    spin-it-2014-givens-oracle-contract
                    PROPERTIES
                        LABELS "paper;python;native;spin-it-2014;independent-verification;givens"
                        TIMEOUT 60
                        ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
                )
                add_test(
                    NAME spin-it-2014-stability-threshold-contract
                    COMMAND "${Python3_EXECUTABLE}" -m unittest
                        test_spin2014_stability_threshold -v
                )
                set_tests_properties(
                    spin-it-2014-stability-threshold-contract
                    PROPERTIES
                        LABELS "paper;python;native;spin-it-2014;independent-verification;stability"
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
            COMMAND "${Python3_EXECUTABLE}" -c "import pybullet, trimesh, numpy, scipy"
            RESULT_VARIABLE MAKE_IT_STAND_PYBULLET_IMPORT_STATUS
            ERROR_VARIABLE MAKE_IT_STAND_PYBULLET_IMPORT_ERROR
        )
        if(NOT MAKE_IT_STAND_PYBULLET_IMPORT_STATUS EQUAL 0)
            message(FATAL_ERROR
                "PyBullet tests require pybullet, trimesh, numpy and scipy: "
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
        add_test(
            NAME pybullet-common-frame-contracts
            COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_bullet_common_frame.py" -v
        )
        set_tests_properties(
            pybullet-common-frame-contracts
            PROPERTIES
                LABELS "simulation;python;pybullet;material-contract;common-frame"
                TIMEOUT 60
                ENVIRONMENT "CONSTRUCTOR_MODULE=bullet_contact_common_frame;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1"
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
                NAME spin-faster-cherchi-face-provenance
                COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_cherchi_face_provenance.py"
            )
            set_tests_properties(
                spin-faster-cherchi-face-provenance
                PROPERTIES
                    LABELS "spin-it;spin-it-faster;boolean;provenance;python;regression"
                    TIMEOUT 60
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
            add_test(
                NAME spin-faster-cherchi-motion-differential
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                    test_cherchi_motion_differential -v
            )
            set_tests_properties(
                spin-faster-cherchi-motion-differential
                PROPERTIES
                    LABELS "spin-it-faster;boolean;differential;python;regression"
                    TIMEOUT 120
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
            add_test(
                NAME spin-faster-cherchi-newton-runner-config
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                    test_faster_cherchi_newton_integration_runner -v
            )
            set_tests_properties(
                spin-faster-cherchi-newton-runner-config
                PROPERTIES
                    LABELS "spin-it-faster;boolean;newton;python;regression"
                    TIMEOUT 60
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
            add_test(
                NAME spin-faster-newton-system-audit-unit
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                    test_faster_newton_system_audit -v
            )
            set_tests_properties(
                spin-faster-newton-system-audit-unit
                PROPERTIES
                    LABELS "spin-it-faster;newton;diagnostics;python;regression"
                    TIMEOUT 60
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
            add_test(
                NAME spin-faster-actual-boolean-jacobian-audit-unit
                COMMAND "${Python3_EXECUTABLE}" -m pytest -q
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/python/test_faster_actual_boolean_jacobian_audit.py"
            )
            set_tests_properties(
                spin-faster-actual-boolean-jacobian-audit-unit
                PROPERTIES
                    LABELS "spin-it-faster;boolean;differential;python;regression"
                    TIMEOUT 60
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
            if(MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS)
                add_test(NAME spin-faster-cherchi-cad-adapter
                    COMMAND "${Python3_EXECUTABLE}" -m unittest test_cherchi_cad -v)
                set_tests_properties(spin-faster-cherchi-cad-adapter PROPERTIES
                    LABELS "spin-it-faster;cherchi;cad;verification;python"
                    TIMEOUT 120 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
                add_test(NAME spin-faster-material-bridge-python
                    COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_material_bridge -v)
                set_tests_properties(spin-faster-material-bridge-python PROPERTIES
                    LABELS "spin-it;spin-it-faster;verification;exact-geometry;python;native"
                    TIMEOUT 300 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            endif()
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
                NAME spin-faster-planar-evaluation-ledger
                COMMAND "${Python3_EXECUTABLE}" -m unittest
                    test_planar_evaluation_ledger -v
            )
            set_tests_properties(spin-faster-planar-evaluation-ledger PROPERTIES
                LABELS "spin-it-faster;planar;python;bindings;regression"
                TIMEOUT 120
                ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
            add_test(NAME spin-faster-planar-saved-material
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_planar_saved_material -v)
            set_tests_properties(spin-faster-planar-saved-material PROPERTIES
                LABELS "spin-it-faster;planar;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-static-table-dispatch
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_static_table_dispatch test_static_full_domain_boundary test_static_volume_objective -v)
            set_tests_properties(spin-faster-static-table-dispatch PROPERTIES
                LABELS "spin-it-faster;static;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-cherchi-full-basis
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_cherchi_basis -v)
            set_tests_properties(spin-faster-cherchi-full-basis PROPERTIES
                LABELS "spin-it-faster;cherchi;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-cherchi-table-dispatch
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_cherchi_table_dispatch -v)
            set_tests_properties(spin-faster-cherchi-table-dispatch PROPERTIES
                LABELS "spin-it-faster;cherchi;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-algebraic-map-integration
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_algebraic_map_integration -v)
            set_tests_properties(spin-faster-algebraic-map-integration PROPERTIES
                LABELS "spin-it-faster;cherchi;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-cherchi-paraboloid-basis
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_cherchi_paraboloid_basis -v)
            add_test(NAME spin-faster-cherchi-fixed-normal-paraboloid
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_fixed_normal_paraboloid test_faster_fixed_normal_runner -v)
            set_tests_properties(spin-faster-cherchi-fixed-normal-paraboloid PROPERTIES
                LABELS "spin-it-faster;cherchi;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-cherchi-plane-basis
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_cherchi_plane_basis -v)
            set_tests_properties(spin-faster-cherchi-plane-basis PROPERTIES
                LABELS "spin-it-faster;cherchi;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-cherchi-pareto-dispatch
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_cherchi_pareto_dispatch -v)
            set_tests_properties(spin-faster-cherchi-pareto-dispatch PROPERTIES
                LABELS "spin-it-faster;cherchi;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-top-lowering-endpoint
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_top_lowering_endpoint -v)
            set_tests_properties(spin-faster-top-lowering-endpoint PROPERTIES
                LABELS "spin-it-faster;python;verification;regression"
                TIMEOUT 120 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-top-endpoint-transport
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_faster_top_endpoint_transport -v)
            set_tests_properties(spin-faster-top-endpoint-transport PROPERTIES
                LABELS "spin-it-faster;cherchi;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            set_tests_properties(spin-faster-cherchi-paraboloid-basis PROPERTIES
                LABELS "spin-it-faster;cherchi;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-buoyant-table-dispatch
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_buoyant_table_dispatch -v)
            set_tests_properties(spin-faster-buoyant-table-dispatch PROPERTIES
                LABELS "spin-it-faster;buoyant;python;verification;regression"
                TIMEOUT 180 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
            add_test(NAME spin-faster-hydrostatic-determinism
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_hydrostatic_determinism -v)
            set_tests_properties(spin-faster-hydrostatic-determinism PROPERTIES
                LABELS "spin-it-faster;buoyant;python;verification;regression"
                TIMEOUT 120 ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}")
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
                NAME independent-exact-feasibility-native
                COMMAND "${Python3_EXECUTABLE}" -m unittest test_exact_feasibility_native -v
            )
            set_tests_properties(
                independent-exact-feasibility-native
                PROPERTIES
                    LABELS "spin-it;reachability;python;native;independent-verification"
                    TIMEOUT 120
                    ENVIRONMENT "${MAKE_IT_STAND_PYTHON_TEST_ENV}"
            )
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
