#pragma once

#include <Eigen/Core>

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "deformation/DeformationSystem.h"
#include "physics/MassProperties.h"
#include "voxel/VoxelTopology.h"


class VoxelGrid;


namespace manufacturing
{


/** Required signed-volume orientation for every connected boundary shell. */
enum class BoundaryOrientation
{
    AnyNonzero,
    OutwardPositive,
    InwardNegative
};


/** Deterministic self-intersection audit statistics. */
struct IntersectionAudit
{
    bool completed = false;
    bool intersectionFree = false;
    std::size_t broadPhasePairCount = 0;
    std::size_t narrowPhasePairCount = 0;
    std::size_t intersectingPairCount = 0;
    int firstFace = -1;
    int secondFace = -1;
};


/** Topological, orientation, and geometric checks for one triangle boundary. */
struct BoundaryAudit
{
    bool empty = true;
    bool matrixShapeValid = false;
    bool finiteVertices = false;
    bool faceIndicesValid = false;
    bool nondegenerateFaces = false;
    bool edgeManifold = false;
    bool vertexManifold = false;
    bool closed = false;
    bool consistentlyOriented = false;
    bool componentOrientationsValid = false;

    int vertexCount = 0;
    int faceCount = 0;
    int referencedVertexCount = 0;
    int unreferencedVertexCount = 0;
    int degenerateFaceCount = 0;
    int boundaryEdgeCount = 0;
    int nonManifoldEdgeCount = 0;
    int inconsistentEdgeCount = 0;
    int nonManifoldVertexCount = 0;
    int faceComponentCount = 0;

    double signedVolume = 0.0;
    std::vector<double> componentSignedVolumes;
    IntersectionAudit selfIntersections;

    /** True for a nonempty strict printable triangle boundary. */
    bool strictBoundaryValid() const noexcept;
};


/** Exact closest-feature clearance between outer and inner triangle shells. */
struct WallClearanceAudit
{
    bool applicable = false;
    bool completed = false;
    bool intersectionFree = false;
    bool meetsRequestedThickness = false;
    double requestedThickness = 0.0;
    double measuredMinimumClearance =
        std::numeric_limits<double>::infinity();
    int outerFace = -1;
    int innerFace = -1;
    std::size_t broadPhasePairCount = 0;
    std::size_t narrowPhasePairCount = 0;
};


/** Per-component winding-number containment of every cavity shell. */
struct ContainmentAudit
{
    bool applicable = false;
    bool completed = false;
    bool allCavitiesInsideOuter = false;
    int testedComponentCount = 0;
    int failedComponent = -1;
    double failedWindingNumber = 0.0;
    int nestedCavityPairCount = 0;
    int firstContainingComponent = -1;
    int firstNestedComponent = -1;
};


/** Conservative sampled local-injectivity gate for deformed active cells. */
struct DeformedCellAudit
{
    bool completed = false;
    bool allActiveCellsValid = false;
    int activeCellCount = 0;
    int invalidCellCount = 0;
    int firstInvalidCell = -1;
    int degenerateBoundaryTriangleCount = 0;
    int nonpositiveJacobianSampleCount = 0;
    double minimumJacobianDeterminant =
        std::numeric_limits<double>::infinity();
};


/** Compact printable boundary; inner face indices follow the outer range. */
struct CoherentBoundaryMesh
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    int outerVertexCount = 0;
    int outerFaceCount = 0;
    int innerVertexCount = 0;
    int innerFaceCount = 0;
};


/** Strict pre-slicer gate and explicit reasons for every rejection. */
struct ManufacturingReadinessReport
{
    BoundaryAudit outerBoundary;
    BoundaryAudit innerBoundary;
    BoundaryAudit coherentBoundary;
    WallClearanceAudit wallClearance;
    ContainmentAudit containment;
    DeformedCellAudit deformedCells;
    VoxelTopologyReport voxelTopology;

    bool synchronizedInnerGeometry = false;
    bool deformedCellsValid = false;
    bool singleMaterialComponent = false;
    bool positiveMaterialVolume = false;
    bool massPropertiesConsistent = false;
    bool ready = false;

    double signedMaterialVolume = 0.0;
    MassProperties surfaceMassProperties;
    MassProperties cellOracleMassProperties;
    Eigen::Vector3d centerOfMass = Eigen::Vector3d::Zero();
    Eigen::Vector3d cellOracleCenterOfMass = Eigen::Vector3d::Zero();
    double massRelativeError =
        std::numeric_limits<double>::infinity();
    double firstMomentRelativeError =
        std::numeric_limits<double>::infinity();
    double centerOfMassRelativeError =
        std::numeric_limits<double>::infinity();
    std::string deformedCellFailure;
    std::vector<std::string> failures;

    /** Stable machine-readable report without host paths or timestamps. */
    std::string toJson() const;
};


struct ManufacturingGeometry
{
    CoherentBoundaryMesh boundary;
    ManufacturingReadinessReport report;
};


struct ManufacturingOptions
{
    /** Zero selects a scale-aware tolerance. */
    double geometricTolerance = 0.0;

    /** Require exactly one six-neighbor material component. */
    bool requireSingleMaterialComponent = true;
};


/**
 * Audit one triangle boundary independently of voxel/deformation state.
 * Empty input is accepted only when allowEmpty is true.
 */
BoundaryAudit auditTriangleBoundary(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    BoundaryOrientation requiredOrientation,
    bool allowEmpty = false,
    double geometricTolerance = 0.0
);


/**
 * Construct and audit the final outer-plus-cavity manufacturing boundary.
 *
 * This is a deterministic geometry-readiness gate. It does not constitute a
 * slicer, structural/FEM, printer, material, or physical-balance certificate.
 */
ManufacturingGeometry prepareManufacturingGeometry(
    const deformation::DeformedGeometry& geometry,
    const VoxelGrid& voxelGrid,
    const ManufacturingOptions& options = ManufacturingOptions()
);


enum class ManufacturingMeshFormat
{
    Obj,
    AsciiStl,
    BinaryStl
};


/** Export only a boundary whose readiness report passed every strict gate. */
void exportReadyManufacturingMesh(
    const ManufacturingGeometry& geometry,
    const std::string& filename,
    ManufacturingMeshFormat format
);


/** Write the stable JSON readiness report, whether it passed or failed. */
void writeManufacturingReadinessReport(
    const ManufacturingReadinessReport& report,
    const std::string& filename
);


} // namespace manufacturing
