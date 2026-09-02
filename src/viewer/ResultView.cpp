#include "viewer/ResultView.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

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

} // namespace

ResultViewArguments parseResultViewArguments(
    int argc,
    const char* const argv[]
)
{
    ResultViewArguments result;
    if (argc < 2 || std::string(argv[1]) != "--result")
    {
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
