#include "viewer/ResultView.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace result_view
{
namespace
{

void validateMesh(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    const char* name
)
{
    if (vertices.rows() == 0 || vertices.cols() != 3 || !vertices.allFinite())
    {
        throw std::invalid_argument(std::string(name) +
            " vertices must be a nonempty finite n-by-3 matrix");
    }
    if (faces.rows() == 0 || faces.cols() != 3)
    {
        throw std::invalid_argument(std::string(name) +
            " faces must be a nonempty m-by-3 matrix");
    }
    if (faces.minCoeff() < 0 || faces.maxCoeff() >= vertices.rows())
    {
        throw std::invalid_argument(std::string(name) +
            " faces reference a vertex outside the mesh");
    }
}


std::string readText(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input)
    {
        throw std::runtime_error("failed to open Spin-It diagnostics: " + path);
    }
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    );
}


std::string stringField(const std::string& text, const std::string& key)
{
    const std::string marker = "\"" + key + "\"";
    std::size_t cursor = text.find(marker);
    if (cursor == std::string::npos) throw std::invalid_argument("missing " + key);
    cursor = text.find(':', cursor + marker.size());
    cursor = text.find('"', cursor + 1);
    const std::size_t end = text.find('"', cursor + 1);
    if (cursor == std::string::npos || end == std::string::npos)
    {
        throw std::invalid_argument("invalid string field " + key);
    }
    return text.substr(cursor + 1, end - cursor - 1);
}


std::vector<double> numericArrayField(
    const std::string& text,
    const std::string& key
)
{
    const std::string marker = "\"" + key + "\"";
    std::size_t begin = text.find(marker);
    if (begin == std::string::npos) throw std::invalid_argument("missing " + key);
    begin = text.find('[', begin + marker.size());
    if (begin == std::string::npos) throw std::invalid_argument("invalid " + key);
    int depth = 0;
    std::size_t end = begin;
    for (; end < text.size(); ++end)
    {
        if (text[end] == '[') ++depth;
        if (text[end] == ']' && --depth == 0) break;
    }
    if (end == text.size()) throw std::invalid_argument("unterminated " + key);
    std::vector<double> values;
    const char* cursor = text.c_str() + begin + 1;
    const char* limit = text.c_str() + end;
    while (cursor < limit)
    {
        char* parsedEnd = nullptr;
        errno = 0;
        const double value = std::strtod(cursor, &parsedEnd);
        if (parsedEnd != cursor)
        {
            if (errno == ERANGE || !std::isfinite(value))
            {
                throw std::invalid_argument("non-finite numeric field " + key);
            }
            values.push_back(value);
            cursor = parsedEnd;
        }
        else
        {
            ++cursor;
        }
    }
    return values;
}


double numericField(const std::string& text, const std::string& key)
{
    const std::string marker = "\"" + key + "\"";
    std::size_t cursor = text.find(marker);
    if (cursor == std::string::npos) throw std::invalid_argument("missing " + key);
    cursor = text.find(':', cursor + marker.size());
    if (cursor == std::string::npos) throw std::invalid_argument("invalid " + key);
    char* parsedEnd = nullptr;
    errno = 0;
    const double value = std::strtod(text.c_str() + cursor + 1, &parsedEnd);
    if (parsedEnd == text.c_str() + cursor + 1 || errno == ERANGE ||
        !std::isfinite(value))
    {
        throw std::invalid_argument("non-finite numeric field " + key);
    }
    return value;
}

} // namespace

ResultViewArguments parseResultViewArguments(
    int argc,
    const char* const argv[]
)
{
    ResultViewArguments result;
    if (argc < 2 || std::string(argv[1]) != "--result")
    {
        if (argc < 2 || std::string(argv[1]) != "--spin-result")
        {
            return result;
        }
        if (argc < 7 || argc > 8)
        {
            throw std::invalid_argument(
                "usage: make-it-stand --spin-result OUTER DESIGN MATERIAL "
                "CAVITY DIAGNOSTICS [WINDOW_TITLE]"
            );
        }
        result.enabled = true;
        result.spinEnabled = true;
        result.outerPath = argv[2];
        result.designPath = argv[3];
        result.materialPath = argv[4];
        result.innerPath = argv[5];
        result.diagnosticsPath = argv[6];
        if (argc == 8 && argv[7] != nullptr && !std::string(argv[7]).empty())
        {
            result.title = argv[7];
        }
        return result;
    }
    if (argc < 4 || argc > 5 || argv[2] == nullptr || argv[3] == nullptr ||
        std::string(argv[2]).empty() || std::string(argv[3]).empty())
    {
        throw std::invalid_argument(
            "usage: make-it-stand --result OUTER INNER [WINDOW_TITLE]"
        );
    }
    result.enabled = true;
    result.outerPath = argv[2];
    result.innerPath = argv[3];
    if (argc == 5 && argv[4] != nullptr && !std::string(argv[4]).empty())
    {
        result.title = argv[4];
    }
    return result;
}


SpinResultViewMetadata loadSpinResultViewMetadata(
    const std::string& diagnosticsPath
)
{
    const std::string text = readText(diagnosticsPath);
    const std::string schema = stringField(text, "schema");
    if ((schema != "spin-it-faster-artifact-v1" &&
         schema != "spin-it-2014-artifact-v1") ||
        text.find("\"geometry_sealed\": true") == std::string::npos)
    {
        throw std::invalid_argument(
            "Spin-It viewer requires a sealed artifact-v1 diagnostics file"
        );
    }
    SpinResultViewMetadata result;
    if (schema == "spin-it-2014-artifact-v1")
    {
        result.artifactFamily = "Spin-It 2014";
        result.cavityLayerName = "Cavity / adaptive voxel field";
    }
    result.mode = stringField(text, "mode");
    if (result.mode == "BUOYANT_STABILITY")
    {
        result.buoyancyAvailable = true;
        result.planeLayerName = "Waterline / buoyancy plane";
        result.displacedLiquidVolume =
            numericField(text, "displaced_liquid_volume");
        result.buoyantMetacenterHeight =
            numericField(text, "buoyant_metacenter_height");
        result.buoyantStabilityMargin =
            numericField(text, "buoyant_stability_margin");
    }
    const std::vector<double> center =
        numericArrayField(text, "center_of_mass_model");
    const std::vector<double> origin =
        numericArrayField(text, "target_origin_model");
    const std::vector<double> target =
        numericArrayField(text, "target_axes_model");
    const std::vector<double> principal =
        numericArrayField(text, "principal_axes_model");
    if (center.size() != 3U || origin.size() != 3U ||
        target.size() != 9U || principal.size() != 9U)
    {
        throw std::invalid_argument("Spin-It diagnostics vector shape is invalid");
    }
    for (int index = 0; index < 3; ++index)
    {
        result.centerOfMass(index) = center[static_cast<std::size_t>(index)];
        result.targetOrigin(index) = origin[static_cast<std::size_t>(index)];
        for (int column = 0; column < 3; ++column)
        {
            result.targetAxes(index, column) =
                target[static_cast<std::size_t>(3 * index + column)];
            result.principalAxes(index, column) =
                principal[static_cast<std::size_t>(3 * index + column)];
        }
    }
    return result;
}


SpinResultViewGeometry prepareSpinResultViewGeometry(
    const Eigen::MatrixXd& outerVertices,
    const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& designVertices,
    const Eigen::MatrixXi& designFaces,
    const Eigen::MatrixXd& materialVertices,
    const Eigen::MatrixXi& materialFaces,
    const Eigen::MatrixXd& cavityVertices,
    const Eigen::MatrixXi& cavityFaces,
    const SpinResultViewMetadata& metadata,
    double displayScaleFraction
)
{
    validateMesh(outerVertices, outerFaces, "outer");
    validateMesh(designVertices, designFaces, "design");
    validateMesh(materialVertices, materialFaces, "material");
    validateMesh(cavityVertices, cavityFaces, "cavity");
    if (!metadata.centerOfMass.allFinite() ||
        !metadata.targetOrigin.allFinite() ||
        !metadata.targetAxes.allFinite() ||
        !metadata.principalAxes.allFinite() ||
        (metadata.buoyancyAvailable &&
         (!std::isfinite(metadata.displacedLiquidVolume) ||
          metadata.displacedLiquidVolume <= 0.0 ||
          !std::isfinite(metadata.buoyantMetacenterHeight) ||
          !std::isfinite(metadata.buoyantStabilityMargin))) ||
        !std::isfinite(displayScaleFraction) || displayScaleFraction <= 0.0)
    {
        throw std::invalid_argument("Spin-It view metadata must be finite");
    }
    const Eigen::Matrix3d targetGram =
        metadata.targetAxes.transpose() * metadata.targetAxes;
    const Eigen::Matrix3d principalGram =
        metadata.principalAxes.transpose() * metadata.principalAxes;
    if (!targetGram.isApprox(Eigen::Matrix3d::Identity(), 1e-10) ||
        !principalGram.isApprox(Eigen::Matrix3d::Identity(), 1e-10))
    {
        throw std::invalid_argument("Spin-It view axes must be orthonormal");
    }

    SpinResultViewGeometry result;
    result.outerVertices = outerVertices;
    result.outerFaces = outerFaces;
    result.designVertices = designVertices;
    result.designFaces = designFaces;
    result.materialVertices = materialVertices;
    result.materialFaces = materialFaces;
    result.cavityVertices = cavityVertices;
    result.cavityFaces = cavityFaces;
    result.metadata = metadata;
    const Eigen::Vector3d extent =
        outerVertices.colwise().maxCoeff().transpose() -
        outerVertices.colwise().minCoeff().transpose();
    result.axisDisplayLength = displayScaleFraction * extent.maxCoeff();
    if (!(result.axisDisplayLength > 0.0)) result.axisDisplayLength = 1.0;

    const Eigen::Vector3d e1 = metadata.targetAxes.col(0);
    const Eigen::Vector3d e2 = metadata.targetAxes.col(1);
    const double radius = result.axisDisplayLength;
    result.targetPlaneVertices.resize(4, 3);
    result.targetPlaneVertices.row(0) =
        (metadata.targetOrigin - radius * e1 - radius * e2).transpose();
    result.targetPlaneVertices.row(1) =
        (metadata.targetOrigin + radius * e1 - radius * e2).transpose();
    result.targetPlaneVertices.row(2) =
        (metadata.targetOrigin - radius * e1 + radius * e2).transpose();
    result.targetPlaneVertices.row(3) =
        (metadata.targetOrigin + radius * e1 + radius * e2).transpose();
    result.targetPlaneFaces.resize(2, 3);
    result.targetPlaneFaces << 0, 1, 2, 2, 1, 3;
    return result;
}

ResultViewGeometry prepareResultViewGeometry(
    const Eigen::MatrixXd& outerVertices,
    const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& innerVertices,
    const Eigen::MatrixXi& innerFaces,
    double groundMarginFraction
)
{
    validateMesh(outerVertices, outerFaces, "outer");
    validateMesh(innerVertices, innerFaces, "inner");
    if (!std::isfinite(groundMarginFraction) || groundMarginFraction < 0.0)
    {
        throw std::invalid_argument("ground margin fraction must be finite and nonnegative");
    }

    ResultViewGeometry result;
    // The sealed artifact has already been transformed into its experiment
    // frame. In particular, applying alignToAxis here destroys the requested
    // standing posture and de-registers the support evidence.
    result.outerVertices = outerVertices;
    result.outerFaces = outerFaces;
    result.innerVertices = innerVertices;
    result.innerFaces = innerFaces;
    result.groundHeight = outerVertices.col(1).minCoeff();

    const double minimumX = outerVertices.col(0).minCoeff();
    const double maximumX = outerVertices.col(0).maxCoeff();
    const double minimumZ = outerVertices.col(2).minCoeff();
    const double maximumZ = outerVertices.col(2).maxCoeff();
    const double largestExtent = std::max(maximumX - minimumX, maximumZ - minimumZ);
    const double margin = groundMarginFraction *
        (largestExtent > 0.0 ? largestExtent : 1.0);

    result.groundVertices.resize(4, 3);
    result.groundVertices <<
        minimumX - margin, result.groundHeight, minimumZ - margin,
        maximumX + margin, result.groundHeight, minimumZ - margin,
        minimumX - margin, result.groundHeight, maximumZ + margin,
        maximumX + margin, result.groundHeight, maximumZ + margin;
    result.groundFaces.resize(2, 3);
    result.groundFaces << 0, 2, 1, 2, 3, 1;
    return result;
}

} // namespace result_view
