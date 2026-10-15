#include "spin/FasterTheoremConditions.h"

#include "spin/DualMonomialCoordinates.h"

#include <algorithm>
#include <cmath>
#include <limits>


namespace
{


bool finiteNonnegative(const double value)
{
    return std::isfinite(value) && value >= 0.0;
}


} // namespace


namespace spin
{


const char* fasterTheoremConditionStatusToken(
    const FasterTheoremConditionStatus status
) noexcept
{
    switch (status)
    {
    case FasterTheoremConditionStatus::SatisfiedNumerically:
        return "SATISFIED_NUMERICALLY";
    case FasterTheoremConditionStatus::Violated:
        return "VIOLATED";
    case FasterTheoremConditionStatus::Unknown:
        return "UNKNOWN";
    }
    return "UNKNOWN";
}


FasterTheoremOneNumericalAudit auditFasterTheoremOneNumericalConditions(
    const QuadricCoefficients10& geometryQuadric,
    const SpinMomentEvaluation& momentEvaluation,
    const double multiplierMu,
    const Eigen::VectorXd& equalityMultipliers,
    const double tolerance
)
{
    FasterTheoremOneNumericalAudit result;
    result.multiplierMu = multiplierMu;
    if (!finiteNonnegative(tolerance) || tolerance == 0.0)
    {
        result.diagnostic = "numerical theorem-audit tolerance must be finite and positive";
        return result;
    }
    if (
        !std::isfinite(multiplierMu) ||
        !equalityMultipliers.allFinite() ||
        !momentEvaluation.objectiveGradient.allFinite() ||
        !momentEvaluation.constraints.allFinite() ||
        !momentEvaluation.constraintJacobian.allFinite() ||
        momentEvaluation.constraintJacobian.rows() != momentEvaluation.constraints.size() ||
        momentEvaluation.constraintJacobian.cols() !=
            momentEvaluation.objectiveGradient.size() ||
        equalityMultipliers.size() != momentEvaluation.constraints.size()
    )
    {
        result.diagnostic = "Eq. 7 numerical inputs are non-finite or dimensionally inconsistent";
        return result;
    }

    const auto dual = DualMonomialCoordinates10::fromGeometry(geometryQuadric);
    const auto& p = dual.values();
    result.unitNormError = std::abs(p.norm() - 1.0);
    result.constraintResidual = momentEvaluation.constraints.norm();
    const auto stationarity = momentEvaluation.objectiveGradient +
        multiplierMu * p -
        momentEvaluation.constraintJacobian.transpose() * equalityMultipliers;
    result.stationarityResidual = stationarity.norm();
    if (
        !finiteNonnegative(result.unitNormError) ||
        !finiteNonnegative(result.constraintResidual) ||
        !finiteNonnegative(result.stationarityResidual)
    )
    {
        result.diagnostic = "Eq. 7 numerical residual is non-finite";
        return result;
    }

    const double stationarityScale = std::max({
        1.0,
        momentEvaluation.objectiveGradient.norm(),
        std::abs(multiplierMu),
        momentEvaluation.constraintJacobian.norm() * equalityMultipliers.norm()
    });
    if (
        result.unitNormError > tolerance ||
        result.constraintResidual > tolerance ||
        !(multiplierMu > 0.0) ||
        result.stationarityResidual > tolerance * stationarityScale
    )
    {
        result.status = FasterTheoremConditionStatus::Violated;
        result.diagnostic =
            "finite numerical data violates a required Eq. 7 or equality condition";
        return result;
    }

    result.status = FasterTheoremConditionStatus::SatisfiedNumerically;
    result.diagnostic =
        "Eq. 7, equality residuals, unit p, and positive mu pass finite-precision audit; "
        "this is not a proof of Theorem 1 or Theorem 2";
    return result;
}


FasterTheoremOneVariationAudit auditFasterTheoremOneDiscreteVariation(
    const Eigen::VectorXd& beta,
    const Eigen::VectorXd& variation,
    const Eigen::VectorXd& cellMeasures,
    const double multiplierMu,
    const double tolerance
)
{
    FasterTheoremOneVariationAudit result;
    if (!finiteNonnegative(tolerance) || tolerance == 0.0)
    {
        result.diagnostic = "discrete variation tolerance must be finite and positive";
        return result;
    }
    if (
        beta.size() == 0 || beta.size() != variation.size() ||
        beta.size() != cellMeasures.size() || !beta.allFinite() ||
        !variation.allFinite() || !cellMeasures.allFinite() ||
        !std::isfinite(multiplierMu)
    )
    {
        result.diagnostic = "discrete variation inputs are non-finite or dimensionally inconsistent";
        return result;
    }
    if (!(multiplierMu > 0.0) || (cellMeasures.array() <= 0.0).any())
    {
        result.status = FasterTheoremConditionStatus::Violated;
        result.diagnostic = "Theorem 1 finite-partition certificate requires positive mu and cell measures";
        return result;
    }

    result.minimumAbsoluteBeta = beta.cwiseAbs().minCoeff();
    if (!std::isfinite(result.minimumAbsoluteBeta) || result.minimumAbsoluteBeta <= tolerance)
    {
        result.diagnostic =
            "a finite-measure partition cell intersects beta=0 or is numerically unresolved; "
            "the measure-zero premise is not certified";
        return result;
    }

    for (Eigen::Index index = 0; index < beta.size(); ++index)
    {
        const double signedVariation = variation(index);
        if ((beta(index) > 0.0 && signedVariation > tolerance) ||
            (beta(index) < 0.0 && signedVariation < -tolerance))
        {
            result.status = FasterTheoremConditionStatus::Violated;
            result.diagnostic = "variation violates Theorem 1 occupied/empty-region sign condition";
            return result;
        }
    }

    result.weightedVariationL1 = (cellMeasures.array() * variation.array().abs()).sum();
    if (!std::isfinite(result.weightedVariationL1) || result.weightedVariationL1 <= tolerance)
    {
        result.status = FasterTheoremConditionStatus::Violated;
        result.diagnostic = "Theorem 1 strict-increase conclusion requires a non-vanishing variation";
        return result;
    }

    result.integralVariationTimesBeta =
        (cellMeasures.array() * variation.array() * beta.array()).sum();
    result.objectiveIncrease = -multiplierMu * result.integralVariationTimesBeta;
    result.strictIncreaseLowerBound =
        multiplierMu * result.minimumAbsoluteBeta * result.weightedVariationL1;
    if (
        !std::isfinite(result.integralVariationTimesBeta) ||
        !std::isfinite(result.objectiveIncrease) ||
        !std::isfinite(result.strictIncreaseLowerBound) ||
        !(result.integralVariationTimesBeta < -tolerance) ||
        !(result.objectiveIncrease + tolerance >= result.strictIncreaseLowerBound)
    )
    {
        result.status = FasterTheoremConditionStatus::Violated;
        result.diagnostic = "finite sign certificate does not establish strict objective increase";
        return result;
    }

    result.status = FasterTheoremConditionStatus::SatisfiedNumerically;
    result.diagnostic =
        "finite partition satisfies Theorem 1 sign argument for this non-vanishing variation; "
        "this is not a proof of Theorem 1 or Theorem 2";
    return result;
}


} // namespace spin
