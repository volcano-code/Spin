#include "manufacturing/PartitionManifest.h"
#include "manufacturing/SlicerValidation.h"

#include <Eigen/Geometry>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <unistd.h>


namespace
{


int fail(const std::string& reason)
{
    std::cerr << "[FAIL] partition-manifest reason=" << reason << '\n';
    return EXIT_FAILURE;
}


void writeFile(const std::filesystem::path& path, const std::string& value)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output << value;
    if (!output) throw std::runtime_error("partition fixture write failed");
}


manufacturing::DeclaredManufacturingPart makePart(
    const std::filesystem::path& directory,
    const std::string& id,
    const std::string& filename
)
{
    const std::filesystem::path mesh = directory / filename;
    manufacturing::DeclaredManufacturingPart part;
    part.id = id;
    part.meshFilename = filename;
    part.expectedMeshSha256 = manufacturing::computeFileSha256(mesh.string());
    part.geometryReportFilename = id + "-geometry.json";
    part.expectedGeometryReportSha256 = manufacturing::computeFileSha256(
        (directory / part.geometryReportFilename).string()
    );
    part.slicerReportFilename = id + "-slicer.json";
    part.expectedSlicerReportSha256 = manufacturing::computeFileSha256(
        (directory / part.slicerReportFilename).string()
    );
    return part;
}


} // namespace


int main()
{
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() /
        ("make_it_stand_partition_" +
         std::to_string(static_cast<long long>(::getpid())));
    const std::filesystem::path outside =
        std::filesystem::temp_directory_path() /
        ("make_it_stand_partition_outside_" +
         std::to_string(static_cast<long long>(::getpid())) + ".stl");
    std::filesystem::create_directories(directory);
    try
    {
        writeFile(outside, "outside fixture\n");
        writeFile(directory / "left.stl", "left fixture\n");
        writeFile(directory / "right.stl", "right fixture\n");
        writeFile(directory / "left-geometry.json", "left geometry\n");
        writeFile(directory / "right-geometry.json", "right geometry\n");
        writeFile(directory / "left-slicer.json", "left slicer\n");
        writeFile(directory / "right-slicer.json", "right slicer\n");
        writeFile(directory / "source.obj", "source boundary\n");

        manufacturing::ExternalPartitionManifest manifest;
        manifest.sourceBoundaryFilename = "source.obj";
        manifest.sourceBoundarySha256 = manufacturing::computeFileSha256(
            (directory / manifest.sourceBoundaryFilename).string()
        );
        manifest.modelUnitsToMillimeters = 25.4;
        manifest.externalGenerator = "manual-test-partition";
        manifest.externalGeneratorVersion = "1";
        manifest.parts.push_back(makePart(directory, "left", "left.stl"));
        manifest.parts.push_back(makePart(directory, "right", "right.stl"));
        manifest.parts[1].assemblyTransform.topRightCorner<3, 1>() =
            Eigen::Vector3d(1.0, 2.0, 3.0);

        const auto passed = manufacturing::auditExternalPartitionManifest(
            manifest,
            directory.string()
        );
        if (!passed.passed || passed.partCount != 2U ||
            !passed.manifestFingerprintAvailable ||
            passed.manifestFingerprint.size() != 64U ||
            passed.toJson().find("\"passed\": true") == std::string::npos)
        {
            throw std::runtime_error("valid manifest was rejected");
        }

        manufacturing::ExternalPartitionManifest changed = manifest;
        changed.externalGeneratorVersion = "2";
        const auto changedAudit =
            manufacturing::auditExternalPartitionManifest(
                changed, directory.string()
            );
        if (!changedAudit.passed ||
            changedAudit.manifestFingerprint == passed.manifestFingerprint)
        {
            throw std::runtime_error("manifest fingerprint was not binding");
        }

        manufacturing::ExternalPartitionManifest duplicate = manifest;
        duplicate.parts[1].id = duplicate.parts[0].id;
        if (manufacturing::auditExternalPartitionManifest(
                duplicate, directory.string()).passed)
        {
            throw std::runtime_error("duplicate part id was accepted");
        }

        manufacturing::ExternalPartitionManifest hash = manifest;
        hash.parts[0].expectedMeshSha256 = std::string(64, '0');
        if (manufacturing::auditExternalPartitionManifest(
                hash, directory.string()).passed)
        {
            throw std::runtime_error("wrong part hash was accepted");
        }

        manufacturing::ExternalPartitionManifest sourceHash = manifest;
        sourceHash.sourceBoundarySha256 = std::string(64, '0');
        if (manufacturing::auditExternalPartitionManifest(
                sourceHash, directory.string()).passed)
        {
            throw std::runtime_error("wrong source hash was accepted");
        }

        manufacturing::ExternalPartitionManifest reportHash = manifest;
        reportHash.parts[0].expectedSlicerReportSha256 = std::string(64, '0');
        if (manufacturing::auditExternalPartitionManifest(
                reportHash, directory.string()).passed)
        {
            throw std::runtime_error("wrong report hash was accepted");
        }

        manufacturing::ExternalPartitionManifest traversal = manifest;
        traversal.parts[0].meshFilename = "../outside.stl";
        if (manufacturing::auditExternalPartitionManifest(
                traversal, directory.string()).passed)
        {
            throw std::runtime_error("path traversal was accepted");
        }

        manufacturing::ExternalPartitionManifest missing = manifest;
        missing.parts[0].meshFilename = "missing.stl";
        const auto missingAudit =
            manufacturing::auditExternalPartitionManifest(
                missing, directory.string()
            );
        if (missingAudit.passed || missingAudit.allFilesAvailable ||
            missingAudit.allHashesMatch)
        {
            throw std::runtime_error(
                "missing artifact retained a positive hash aggregate"
            );
        }

        std::filesystem::create_symlink(outside, directory / "escaped.stl");
        manufacturing::ExternalPartitionManifest escaped = manifest;
        escaped.parts[0].meshFilename = "escaped.stl";
        escaped.parts[0].expectedMeshSha256 =
            manufacturing::computeFileSha256(outside.string());
        const auto escapedAudit =
            manufacturing::auditExternalPartitionManifest(
                escaped, directory.string()
            );
        if (escapedAudit.passed || escapedAudit.allFilesAvailable ||
            escapedAudit.allHashesMatch)
        {
            throw std::runtime_error(
                "symlink escape outside manifest root was accepted"
            );
        }

        manufacturing::ExternalPartitionManifest scaled = manifest;
        scaled.parts[0].assemblyTransform(0, 0) = 2.0;
        if (manufacturing::auditExternalPartitionManifest(
                scaled, directory.string()).passed)
        {
            throw std::runtime_error("non-rigid assembly transform was accepted");
        }
    }
    catch (const std::exception& error)
    {
        std::error_code ignored;
        std::filesystem::remove_all(directory, ignored);
        std::filesystem::remove(outside, ignored);
        return fail(error.what());
    }

    std::error_code ignored;
    std::filesystem::remove_all(directory, ignored);
    std::filesystem::remove(outside, ignored);
    std::cout << "[PASS] partition-manifest\n";
    return EXIT_SUCCESS;
}
