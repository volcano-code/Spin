#include "verification/ExactLibiglBoolean.h"
#include "verification/ExactRationalCoordinates.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gmpq.h>
#include <CGAL/number_utils.h>
#include <igl/copyleft/cgal/mesh_boolean.h>

#include <locale>
#include <sstream>
#include <stdexcept>

namespace verification {
namespace {
template<class Matrix>
ExactLibiglBooleanConstruction construct(
    const Matrix& designVertices, const Eigen::MatrixXi& designFaces,
    const Matrix& cutterVertices, const Eigen::MatrixXi& cutterFaces,
    ExactMaterialBooleanOperation operation) {
    igl::MeshBooleanType type;
    if (operation == ExactMaterialBooleanOperation::Intersection)
        type = igl::MESH_BOOLEAN_TYPE_INTERSECT;
    else if (operation == ExactMaterialBooleanOperation::Difference)
        type = igl::MESH_BOOLEAN_TYPE_MINUS;
    else
        throw std::invalid_argument("unsupported exact libigl Boolean operation");

    using Exact = CGAL::Epeck::FT;
    Eigen::Matrix<Exact, Eigen::Dynamic, Eigen::Dynamic> exactVertices;
    ExactLibiglBooleanConstruction result;
    result.succeeded = igl::copyleft::cgal::mesh_boolean(
        designVertices, designFaces, cutterVertices, cutterFaces, type,
        exactVertices, result.faces, result.birthFaces);
    if (!result.succeeded)
        return result;

    result.vertices.reserve(static_cast<std::size_t>(exactVertices.rows()));
    for (Eigen::Index row = 0; row < exactVertices.rows(); ++row) {
        std::array<std::string, 3> coordinates;
        for (Eigen::Index axis = 0; axis < 3; ++axis) {
            const CGAL::Gmpq rational = CGAL::exact(exactVertices(row, axis));
            std::ostringstream text;
            text.imbue(std::locale::classic());
            text << rational;
            coordinates[static_cast<std::size_t>(axis)] = text.str();
        }
        result.vertices.push_back(std::move(coordinates));
    }
    return result;
}

using RationalMatrix = Eigen::Matrix<CGAL::Epeck::FT, Eigen::Dynamic, Eigen::Dynamic>;

RationalMatrix rationalOperand(
    const std::vector<std::array<std::string, 3>>& coordinates,
    const Eigen::MatrixXi& faces) {
    if (coordinates.empty() || faces.rows() == 0 || faces.cols() != 3)
        throw std::invalid_argument("exact Boolean requires nonempty triangulated operands");
    for (Eigen::Index row = 0; row < faces.rows(); ++row)
        for (Eigen::Index axis = 0; axis < 3; ++axis)
            if (faces(row, axis) < 0 || static_cast<std::size_t>(faces(row, axis)) >= coordinates.size())
                throw std::invalid_argument("exact Boolean triangle index outside vertex range");
    RationalMatrix result(coordinates.size(), 3);
    for (std::size_t row = 0; row < coordinates.size(); ++row)
        for (std::size_t axis = 0; axis < 3; ++axis)
            result(row, axis) = CGAL::Epeck::FT(parseExactRationalCoordinate(coordinates[row][axis]));
    return result;
}
} // namespace

ExactLibiglBooleanConstruction computeExactLibiglBooleanConstruction(
    const Eigen::MatrixXd& designVertices, const Eigen::MatrixXi& designFaces,
    const Eigen::MatrixXd& cutterVertices, const Eigen::MatrixXi& cutterFaces,
    ExactMaterialBooleanOperation operation) {
    return construct(designVertices, designFaces, cutterVertices, cutterFaces, operation);
}

ExactLibiglBooleanConstruction computeExactLibiglBooleanConstructionFromRationals(
    const std::vector<std::array<std::string, 3>>& designVertices,
    const Eigen::MatrixXi& designFaces,
    const std::vector<std::array<std::string, 3>>& cutterVertices,
    const Eigen::MatrixXi& cutterFaces, ExactMaterialBooleanOperation operation) {
    const auto design = rationalOperand(designVertices, designFaces);
    const auto cutter = rationalOperand(cutterVertices, cutterFaces);
    return construct(design, designFaces, cutter, cutterFaces, operation);
}

} // namespace verification
