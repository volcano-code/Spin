#pragma once

#include "verification/ExactMaterialRegion.h"

#include <Eigen/Core>

#include <array>
#include <string>
#include <vector>

namespace verification {

// Diagnostic only: run the pinned libigl Boolean with exact output scalars,
// before its usual assignment to a binary64 material mesh. Inputs remain
// binary64 polyhedra; this is not a continuous-quadric or solver output.
struct ExactLibiglBooleanConstruction {
    bool succeeded = false;
    std::vector<std::array<std::string, 3>> vertices;
    Eigen::MatrixXi faces;
    Eigen::VectorXi birthFaces;
};

ExactLibiglBooleanConstruction computeExactLibiglBooleanConstruction(
    const Eigen::MatrixXd& designVertices, const Eigen::MatrixXi& designFaces,
    const Eigen::MatrixXd& cutterVertices, const Eigen::MatrixXi& cutterFaces,
    ExactMaterialBooleanOperation operation);

// Exact polyhedral operands stay rational throughout construction. This is
// not a continuous-quadric claim or an automatic material acceptance gate.
ExactLibiglBooleanConstruction computeExactLibiglBooleanConstructionFromRationals(
    const std::vector<std::array<std::string, 3>>& designVertices,
    const Eigen::MatrixXi& designFaces,
    const std::vector<std::array<std::string, 3>>& cutterVertices,
    const Eigen::MatrixXi& cutterFaces, ExactMaterialBooleanOperation operation);

} // namespace verification
