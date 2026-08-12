#include "manufacturing/PartitionManifest.h"

#include "manufacturing/SlicerValidation.h"

#include <Eigen/LU>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <locale>
#include <set>
#include <sstream>
#include <string>


namespace manufacturing
{
namespace
{


bool hashTokenValid(const std::string& value)
{
    return value.size() == 64U && std::all_of(
        value.begin(),
        value.end(),
        [](unsigned char character)
        {
            return std::isxdigit(character) != 0;
        }
    );
}


std::string lowercase(std::string value)
{
    std::transform(
        value.begin(), value.end(), value.begin(),
        [](unsigned char character)
        {
            return static_cast<char>(std::tolower(character));
        }
    );
    return value;
}


bool regularFile(const std::filesystem::path& path)
{
    std::error_code error;
    return std::filesystem::is_regular_file(path, error);
}


bool safeRelativePath(const std::filesystem::path& path)
{
    return !path.empty() && !path.is_absolute() && std::none_of(
        path.begin(), path.end(),
        [](const std::filesystem::path& component)
        {
            return component == "..";
        }
    );
}


bool rigidTransform(const Eigen::Matrix4d& transform)
{
    if (!transform.allFinite() ||
        !transform.row(3).isApprox(Eigen::RowVector4d(0.0, 0.0, 0.0, 1.0), 1.0e-12))
    {
        return false;
    }
    const Eigen::Matrix3d rotation = transform.topLeftCorner<3, 3>();
    return
        rotation.determinant() > 0.0 &&
        std::abs(rotation.determinant() - 1.0) <= 1.0e-10 &&
        (rotation.transpose() * rotation - Eigen::Matrix3d::Identity()).norm()
            <= 1.0e-10;
}


const char* jsonBool(bool value)
{
    return value ? "true" : "false";
}


} // namespace


PartitionManifestAudit auditExternalPartitionManifest(
    const ExternalPartitionManifest& manifest,
    const std::string& manifestDirectory
)
{
    PartitionManifestAudit audit;
    audit.completed = true;
    audit.partCount = manifest.parts.size();
    audit.schemaValid =
        hashTokenValid(manifest.sourceBoundarySha256) &&
        std::isfinite(manifest.modelUnitsToMillimeters) &&
        manifest.modelUnitsToMillimeters > 0.0 &&
        !manifest.externalGenerator.empty() &&
        !manifest.externalGeneratorVersion.empty() &&
        !manifest.parts.empty() &&
        std::filesystem::path(manifestDirectory).is_absolute();

    std::set<std::string> identifiers;
    audit.uniquePartIds = true;
    audit.allFilesAvailable = true;
    audit.allHashesMatch = true;
    audit.allTransformsFinite = true;
    audit.allTransformsRigid = true;

    const std::filesystem::path base(manifestDirectory);
    for (std::size_t index = 0; index < manifest.parts.size(); ++index)
    {
        const DeclaredManufacturingPart& part = manifest.parts[index];
        const auto failPart = [&audit, index]()
        {
            if (audit.firstFailedPart < 0)
            {
                audit.firstFailedPart = static_cast<int>(index);
            }
        };

        if (part.id.empty())
        {
            audit.schemaValid = false;
            audit.uniquePartIds = false;
            failPart();
        }
        else if (!identifiers.insert(part.id).second)
        {
            audit.uniquePartIds = false;
            failPart();
        }
        if (!part.assemblyTransform.allFinite())
        {
            audit.allTransformsFinite = false;
            failPart();
        }
        if (!rigidTransform(part.assemblyTransform))
        {
            audit.allTransformsRigid = false;
            failPart();
        }

        const std::array<std::pair<std::string, std::string>, 3> artifacts{{
            {part.meshFilename, part.expectedMeshSha256},
            {part.geometryReportFilename,
             part.expectedGeometryReportSha256},
            {part.slicerReportFilename, part.expectedSlicerReportSha256}
        }};
        for (const auto& artifact : artifacts)
        {
            const std::filesystem::path relative(artifact.first);
            const bool safeRelative = safeRelativePath(relative);
            const std::filesystem::path file = base / relative;
            const bool available = safeRelative && regularFile(file);
            const bool hashValid = hashTokenValid(artifact.second);
            if (!safeRelative || !hashValid)
            {
                audit.schemaValid = false;
                failPart();
            }
            if (!available)
            {
                audit.allFilesAvailable = false;
                failPart();
            }
            if (!hashValid ||
                (available && computeFileSha256(file.string()) !=
                    lowercase(artifact.second)))
            {
                audit.allHashesMatch = false;
                failPart();
            }
        }
    }

    if (!audit.schemaValid) audit.failures.push_back("partition_schema_invalid");
    if (!audit.uniquePartIds) audit.failures.push_back("partition_part_ids_not_unique");
    if (!audit.allFilesAvailable) audit.failures.push_back("partition_files_unavailable");
    if (!audit.allHashesMatch) audit.failures.push_back("partition_hash_mismatch");
    if (!audit.allTransformsFinite) audit.failures.push_back("partition_transform_nonfinite");
    if (!audit.allTransformsRigid) audit.failures.push_back("partition_transform_not_rigid");
    audit.passed = audit.failures.empty();
    return audit;
}


std::string PartitionManifestAudit::toJson() const
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << "{\n"
        << "  \"schema\": \"make-it-stand-r8-partition-audit-v1\",\n"
        << "  \"claim_ceiling\": \"external partition declaration only\",\n"
        << "  \"passed\": " << jsonBool(passed) << ",\n"
        << "  \"part_count\": " << partCount << ",\n"
        << "  \"schema_valid\": " << jsonBool(schemaValid) << ",\n"
        << "  \"unique_part_ids\": " << jsonBool(uniquePartIds) << ",\n"
        << "  \"all_files_available\": " << jsonBool(allFilesAvailable) << ",\n"
        << "  \"all_hashes_match\": " << jsonBool(allHashesMatch) << ",\n"
        << "  \"all_transforms_finite\": " << jsonBool(allTransformsFinite) << ",\n"
        << "  \"all_transforms_rigid\": " << jsonBool(allTransformsRigid) << ",\n"
        << "  \"first_failed_part\": " << firstFailedPart << ",\n"
        << "  \"failures\": [";
    for (std::size_t index = 0; index < failures.size(); ++index)
    {
        if (index > 0) output << ", ";
        output << '"' << failures[index] << '"';
    }
    output << "]\n}\n";
    return output.str();
}


} // namespace manufacturing
