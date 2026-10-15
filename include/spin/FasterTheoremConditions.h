#pragma once

#include "spin/QuadricCoefficients.h"
#include "spin/SpinMomentProblem.h"

#include <Eigen/Core>

#include <string>


namespace spin
{


/** Outcome of a finite-precision audit of the numerical parts of Eq. (7).
 *
 * This is intentionally not a theorem prover.  In particular, it cannot
 * establish the infinite-dimensional local-optimality premise of Theorem 2,
 * nor can it replace geometric proof that a material characteristic function
 * agrees with beta > 0 except on a measure-zero set.
 */
enum class FasterTheoremConditionStatus
{
    SatisfiedNumerically,
    Violated,
    Unknown
};


const char* fasterTheoremConditionStatusToken(
    FasterTheoremConditionStatus status
) noexcept;


struct FasterTheoremOneNumericalAudit
{
    FasterTheoremConditionStatus status = FasterTheoremConditionStatus::Unknown;
    double unitNormError = 0.0;
    double constraintResidual = 0.0;
    double stationarityResidual = 0.0;
    double multiplierMu = 0.0;
    std::string diagnostic;
    /** This literal scope is evidence-visible and must not be widened. */
    std::string proofScope =
        "NUMERICAL_EQ7_AND_EQUALITY_CONDITIONS_ONLY_NOT_A_THEOREM_PROOF";
};


/** A finite-partition certificate for the sign argument in Theorem 1.
 *
 * For a partition with strictly positive cell measures this records the
 * discrete analogue of Eqs. (6), (8), and the final sign argument in the
 * proof of Theorem 1.  It is intentionally unable to represent a
 * measure-zero beta=0 set: a zero-beta cell has positive partition measure,
 * so the result is Unknown rather than silently being treated as harmless.
 *
 * This is a checkable numerical corollary for one supplied variation, not a
 * proof of the infinite-dimensional theorem or a replacement for the
 * theorem's KKT/regularity premises.
 */
struct FasterTheoremOneVariationAudit
{
    FasterTheoremConditionStatus status = FasterTheoremConditionStatus::Unknown;
    double minimumAbsoluteBeta = 0.0;
    double weightedVariationL1 = 0.0;
    double integralVariationTimesBeta = 0.0;
    double objectiveIncrease = 0.0;
    double strictIncreaseLowerBound = 0.0;
    std::string diagnostic;
    /** This literal scope is evidence-visible and must not be widened. */
    std::string proofScope =
        "FINITE_PARTITION_SIGN_CERTIFICATE_ONLY_NOT_A_THEOREM_PROOF";
};


/**
 * Audit the finite-dimensional numerical hypotheses occurring in Faster
 * Theorem 1, Eq. (7): p has unit norm, equalities vanish, mu is positive,
 * and grad_r f = -mu p + Jac_r(g)^T lambda.
 *
 * `geometryQuadric` is converted to the literal monomial coordinate p used
 * by the paper: cross-term coefficients are doubled before normalizing.  No
 * caller-supplied boolean can declare a theorem passed.  Malformed numerical
 * data returns Unknown; finite data that violates a stated numeric condition
 * returns Violated.
 */
FasterTheoremOneNumericalAudit auditFasterTheoremOneNumericalConditions(
    const QuadricCoefficients10& geometryQuadric,
    const SpinMomentEvaluation& momentEvaluation,
    double multiplierMu,
    const Eigen::VectorXd& equalityMultipliers,
    double tolerance = 1e-12
);


/**
 * Check the discrete strict-increase argument from Theorem 1's proof.
 *
 * ``beta`` stores the literal quadratic function on a finite partition,
 * ``variation`` stores an admissible first variation, and ``cellMeasures``
 * gives the positive measure of every cell.  For beta > 0 (where Eq. (6)
 * selects material), the variation must be non-positive; for beta < 0 it
 * must be non-negative.  A strictly positive result certifies
 *
 *   -mu * sum_i(cellMeasures_i * variation_i * beta_i) > 0.
 *
 * A beta value indistinguishable from zero on a finite-measure cell is
 * Unknown rather than being mistaken for the paper's measure-zero zero set.
 */
FasterTheoremOneVariationAudit auditFasterTheoremOneDiscreteVariation(
    const Eigen::VectorXd& beta,
    const Eigen::VectorXd& variation,
    const Eigen::VectorXd& cellMeasures,
    double multiplierMu,
    double tolerance = 1e-12
);


} // namespace spin
