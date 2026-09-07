#pragma once
#include "spin/QuadricTriangulator.h"
#include <string>
#include <vector>

namespace spin {
/** Fixed input-feature identities, not coordinate-sorted Boolean vertex IDs. */
struct QuadricDiscreteMomentDifferential {
    bool defined = false;
    std::string diagnostic;
    Eigen::Matrix<double,10,10> jacobian = Eigen::Matrix<double,10,10>::Zero();
    std::vector<std::string> vertexProvenance;
    int constructedIntersectionVertexCount = 0;
    double maximumConstraintResidual = 0;
    double minimumConstraintSingularRatio = 1;
};

/** Differentiate a supported, fixed-provenance Boolean boundary in target coordinates.
 * No finite differences; the inverse-Gauss map already includes coefficient
 * factor-two and normalization, so no second tangent projection is applied. */
QuadricDiscreteMomentDifferential differentiateDiscreteEllipsoidMoments(
    const TriangleMesh& design, const QuadricTriangulation& triangulation,
    const QuadricCoefficients10& quadric, const TriangleMesh& material,
    const Eigen::VectorXi& birthFaces, double density);
}
