#pragma once
#include "spin/QuadricMomentEvaluator.h"
#include "spin/PolyhedralMomentDifferential.h"
#include <optional>

namespace spin {
struct FullRankMomentEvaluation {
    physics::RawMassMoments10 variableMomentsTarget;
    physics::RawMassMoments10 totalMomentsTarget;
    TriangleMesh materialBoundaryTarget;
    ExactBooleanVertexCoordinates materialExactVerticesTarget;
    ExactBooleanVertexCoordinates designExactVerticesTarget;
    ExactBooleanVertexCoordinates sourceExactVerticesTarget;
    ExactBooleanVertexCoordinates sourceCoordinateLowerTarget;
    ExactBooleanVertexCoordinates sourceCoordinateUpperTarget;
    std::string sourceConstructionIdentity;
    /** Same-operation branch metadata, excluding moving coordinates. Present
     * only for explicit rational production. A change is NOT proof that a
     * transition is differentiable; consumers must reject/re-audit it. */
    std::string rationalBranchIdentity;
    Eigen::MatrixXi sourceFaces;
    Eigen::VectorXi materialBirthFaces;
    std::vector<std::string> variableMomentsExactTarget;
    std::optional<Eigen::Matrix<double,10,10>> coefficientJacobian;
    std::string geometryHash;
    /**
     * Provenance-aware Boolean topology identity for the material boundary.
     * This is deliberately distinct from geometryHash: a coordinate-sorted
     * mesh can have identical connectivity while its birth features change.
     */
    std::string topologyHash;
    /** Full-rank Newton always consumes this same discrete Boolean derivative,
     * even when a legacy value-side setting is supplied for compatibility. */
    std::string effectiveDerivativeMethod;
    std::string booleanBackendIdentity;
    std::string sourceScientificHash;
    std::string scientificHash;
    std::string diagnostic;
};
/** Exact binary64 comparison of oriented triangle sets, independent of mesh
 * vertex/face storage order. This is triangulation-sensitive by design: two
 * boundaries with equal mass moments but different faces are different maps. */
bool sameOrientedTriangleBoundary(const TriangleMesh&, const TriangleMesh&);
/** Opt-in returned-native-moment differential. Old evaluator/solver unchanged.
 * Derivative columns perturb unit q.values() and renormalize once.
 * Unsupported source strata and ambiguous Boolean motion return absent J,
 * including the native scale-aware 1e-10 selector exclusion unless an exact
 * strict INACTIVE-centroid/uniform-interface margin is established, followed
 * by the separately required matching regular source and Boolean motion.
 * never a zero placeholder. Geometry failures may throw. Not a solver/seal.
 * Values, returned boundary, derivative and provenance are copied from one
 * native evaluator call; the adapter must not replay the Boolean.
 */
FullRankMomentEvaluation evaluateFullRankMaterialMoments(
    const SpinDesignRequest&, const QuadricCoefficients10&,
    const QuadricMomentEvaluatorSettings& settings = {});
}
