#pragma once

#include <Eigen/Core>

#include <cstddef>
#include <string>
#include <vector>


namespace manufacturing
{


/** One externally authored print part; R8 never generates the partition. */
struct DeclaredManufacturingPart
{
    std::string id;
    std::string meshFilename;
    std::string expectedMeshSha256;
    std::string geometryReportFilename;
    std::string expectedGeometryReportSha256;
    std::string slicerReportFilename;
    std::string expectedSlicerReportSha256;
    Eigen::Matrix4d assemblyTransform = Eigen::Matrix4d::Identity();
};


/**
 * Optional interface for a user-declared process-specific split.
 *
 * This manifest records hashes and assembly transforms only.  It does not
 * design cuts, connectors, pins, drainage holes, supports, or assembly paths.
 */
struct ExternalPartitionManifest
{
    std::string sourceBoundaryFilename;
    std::string sourceBoundarySha256;
    double modelUnitsToMillimeters = 0.0;
    std::string externalGenerator;
    std::string externalGeneratorVersion;
    std::vector<DeclaredManufacturingPart> parts;
};


struct PartitionManifestAudit
{
    bool completed = false;
    bool schemaValid = false;
    bool uniquePartIds = false;
    bool allFilesAvailable = false;
    bool allHashesMatch = false;
    bool allTransformsFinite = false;
    bool allTransformsRigid = false;
    bool passed = false;

    std::size_t partCount = 0;
    int firstFailedPart = -1;
    std::vector<std::string> failures;

    std::string toJson() const;
};


/** Validate an explicit partition declaration without altering any mesh. */
PartitionManifestAudit auditExternalPartitionManifest(
    const ExternalPartitionManifest& manifest,
    const std::string& manifestDirectory
);


} // namespace manufacturing
