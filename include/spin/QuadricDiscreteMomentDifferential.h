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

/** Discrete Boolean-moment differential restricted to one rank-deficient
 * quadric stratum. Currently supports planes, simple-spectrum elliptic,
 * hyperbolic, and parabolic cylinders, and simple-spectrum paraboloids. The Jacobian has ten moment rows
 * and d coefficient columns; `coefficientBasis` maps those d coordinates into
 * the normalized ten-vector.
 */
struct QuadricStratumDiscreteMomentDifferential {
    bool defined=false;
    std::string diagnostic;
    Eigen::MatrixXd coefficientBasis;
    Eigen::MatrixXd jacobian;
    std::vector<std::string> vertexProvenance;
    int constructedIntersectionVertexCount=0;
    double maximumConstraintResidual=0;
    double minimumConstraintSingularRatio=1;
};

/** Differentiate a supported, fixed-provenance Boolean boundary in target coordinates.
 * No finite differences; the inverse-Gauss map already includes coefficient
 * factor-two and normalization, so no second tangent projection is applied. */
QuadricDiscreteMomentDifferential differentiateDiscreteEllipsoidMoments(
    const TriangleMesh& design, const QuadricTriangulation& triangulation,
    const QuadricCoefficients10& quadric, const TriangleMesh& material,
    const Eigen::VectorXi& birthFaces, double density);

/** Differentiate the actual fixed-topology Boolean moment map for supported
 * full-rank quadrics. Hyperboloids include the motion of truncation roots and
 * cap centroids. The requested triangulation must reproduce the exact source
 * mesh used to construct the Boolean result; topology/radius-branch changes
 * and rank-deficient quadric strata fail closed. */
QuadricDiscreteMomentDifferential differentiateDiscreteQuadricMoments(
    const TriangleMesh& design, const QuadricTriangulation& triangulation,
    const QuadricCoefficients10& quadric, const TriangleMesh& material,
    const Eigen::VectorXi& birthFaces, double density,
    const QuadricTriangulationSettings& triangulationSettings);

/** Differentiate the same fixed-topology Boolean-to-moment map on a supported
 * rank-deficient stratum (plane, elliptic/hyperbolic/parabolic cylinder, or paraboloid).
 * This returns a reduced Jacobian and is not an Equation-14 full-rank solver
 * interface.
 */
QuadricStratumDiscreteMomentDifferential differentiateDiscreteQuadricStratumMoments(
    const TriangleMesh& design, const QuadricTriangulation& triangulation,
    const QuadricCoefficients10& quadric, const TriangleMesh& material,
    const Eigen::VectorXi& birthFaces, double density,
    const QuadricTriangulationSettings& triangulationSettings);
}
