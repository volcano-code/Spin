#pragma once
#include "spin/PolyhedralMomentDifferential.h"
#include "spin/QuadricInverseGaussMap.h"

namespace spin {
struct QuadricSourceMotion {
    bool defined=false;
    QuadricTriangulation triangulation;
    std::vector<SourceVertexJacobian> velocity;
    std::vector<std::string> provenance;
    std::string diagnostic;
    std::string scientificHash;
};

/** Versioned rational derivative attachment for an actual authoritative source.
 * Keeps rational velocities and their proof bounds; never recovers authority
 * from binary64 views. Local analytic construction derivatives are not the
 * derivative of coordinate quantization or a Boolean/moment certificate. */
struct QuadricRationalSourceMotion {
    using Matrix=QuadricRationalFiniteDifferentialResult::Matrix;
    bool defined=false;
    QuadricTriangulation triangulation;
    std::vector<Matrix> coefficientJacobian,lower,upper;
    std::vector<std::array<std::string,3>> pointProofLower,pointProofUpper;
    std::vector<std::string> provenance;
    std::string diagnostic,scientificHash;
};

/** Restricted local source motion on a rank-deficient quadric stratum.
 * `coefficientBasis` is 10 x d and each vertex velocity is 3 x d. It is
 * intentionally not represented as SourceVertexJacobian (3 x 10): these
 * derivatives are only defined for tangent directions preserving the fixed
 * surface stratum. Supported strata are planes, simple-spectrum elliptic,
 * hyperbolic, and parabolic cylinders, and simple-spectrum paraboloids; unsupported or
 * branch-unstable cases fail closed.
 */
struct QuadricStratumSourceMotion {
    bool defined=false;
    QuadricTriangulation triangulation;
    Eigen::MatrixXd coefficientBasis;
    std::vector<Eigen::MatrixXd> vertexVelocity;
    std::vector<std::string> provenance;
    std::string diagnostic;
    std::string scientificHash;
};
/** Full-rank native finite construction, with unit coefficient tangent chain.
 * Local real-arithmetic derivative, not the derivative of 80 floating-point
 * bisection decisions. Requires a unique regular truncation root agreeing with
 * the returned point, certified by exact-rational root isolation. Hyperboloid
 * adaptive radius doubling is supported by binding the source differential to
 * the returned effective radius. Derivatives assume unchanged finite topology
 * and radius-selection branch; this function does not certify a full parameter
 * neighborhood. Invalid/native geometry construction may throw.
 * Rank-deficient paraboloid, elliptic/hyperbolic/parabolic-cylinder, and plane strata are rejected
 * with an explicit `RANK_DEFICIENT_STRATUM_DERIVATIVE_UNDEFINED` diagnostic:
 * their tangent spaces are lower-dimensional than the solver's ten coefficient
 * coordinates and must not be represented by a fabricated full Jacobian.
 * No manufacture, geometric-equivalence or solver-convergence certification.
 */
QuadricSourceMotion evaluateQuadricSourceMotion(
    const QuadricCoefficients10&, const Eigen::MatrixXd& designVertices,
    const QuadricTriangulationSettings& settings = {});

/** Same-operation counterpart for the explicit rational source profile.
 * Consumes the ACTUAL native Boolean operand, verifies its coefficient and
 * complete content/provenance binding, then retains all existing source/root
 * differential checks. No new source geometry is constructed. Invalid or
 * stale input returns undefined; the historical entry remains compatible. */
QuadricSourceMotion evaluateQuadricSourceMotion(
    const QuadricCoefficients10&, const Eigen::MatrixXd& designVertices,
    const QuadricTriangulationSettings&, const QuadricTriangulation& actualSource);

/** Explicit authority path for same-operation full-rank rational sources.
 * Verifies actual finite/root coordinates and cap loops and preserves all
 * rational derivatives. Stale identity, branch or uncertified motion fails
 * closed; legacy represented-coefficient entries remain unchanged. */
QuadricRationalSourceMotion evaluateQuadricSourceMotion(
    const QuadricCoefficientAuthority10&,const Eigen::MatrixXd& designVertices,
    const QuadricTriangulationSettings&,const QuadricTriangulation& actualSource);

/** Exact local derivative of supported rank-deficient native finite
 * constructions in their reduced coefficient tangent coordinates. This does
 * not provide a full-rank Eq. (14) Jacobian or a solver result.
 */
QuadricStratumSourceMotion evaluateQuadricStratumSourceMotion(
    const QuadricCoefficients10&, const Eigen::MatrixXd& designVertices,
    const QuadricTriangulationSettings& settings = {});
}
