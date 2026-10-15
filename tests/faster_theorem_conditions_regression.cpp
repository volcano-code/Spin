#include "spin/FasterTheoremConditions.h"
#include "spin/SpinMomentProblem.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace
{

int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=faster-theorem-conditions " << message << '\n';
    return EXIT_FAILURE;
}


spin::QuadricCoefficients10 constantPositiveQuadric()
{
    spin::QuadricCoefficients10::Vector coefficients =
        spin::QuadricCoefficients10::Vector::Zero();
    coefficients(0) = 1.0;
    return spin::QuadricCoefficients10::fromUnnormalized(coefficients);
}


} // namespace


int main()
{
    const spin::SpinMomentProblem problem = spin::SpinMomentProblem::staticBalance(
        spin::StaticObjective::MaximumMass
    );
    spin::SpinMomentProblem::Vector10 raw = spin::SpinMomentProblem::Vector10::Zero();
    raw(0) = 1.0;
    const spin::SpinMomentEvaluation evaluated = problem.evaluate(raw);

    const spin::FasterTheoremOneNumericalAudit accepted =
        spin::auditFasterTheoremOneNumericalConditions(
            constantPositiveQuadric(), evaluated, 1.0,
            Eigen::VectorXd::Zero(evaluated.constraints.size())
        );
    if (
        accepted.status != spin::FasterTheoremConditionStatus::SatisfiedNumerically ||
        accepted.unitNormError != 0.0 ||
        accepted.constraintResidual != 0.0 ||
        accepted.stationarityResidual != 0.0 ||
        accepted.proofScope !=
            "NUMERICAL_EQ7_AND_EQUALITY_CONDITIONS_ONLY_NOT_A_THEOREM_PROOF"
    )
    {
        return fail("exact Eq. 7 and equality fixture was not accepted in its narrow scope");
    }

    spin::SpinMomentProblem::Vector10 infeasibleRaw = raw;
    infeasibleRaw(1) = 1e-4;
    const spin::FasterTheoremOneNumericalAudit infeasible =
        spin::auditFasterTheoremOneNumericalConditions(
            constantPositiveQuadric(), problem.evaluate(infeasibleRaw), 1.0,
            Eigen::VectorXd::Zero(2)
        );
    if (
        infeasible.status != spin::FasterTheoremConditionStatus::Violated ||
        infeasible.constraintResidual <= 1e-12
    )
    {
        return fail("nonzero equality constraint was accepted");
    }

    const spin::FasterTheoremOneNumericalAudit nonpositiveMu =
        spin::auditFasterTheoremOneNumericalConditions(
            constantPositiveQuadric(), evaluated, 0.0,
            Eigen::VectorXd::Zero(evaluated.constraints.size())
        );
    if (nonpositiveMu.status != spin::FasterTheoremConditionStatus::Violated)
    {
        return fail("mu <= 0 was accepted");
    }

    const spin::FasterTheoremOneNumericalAudit malformed =
        spin::auditFasterTheoremOneNumericalConditions(
            constantPositiveQuadric(), evaluated, 1.0, Eigen::VectorXd::Zero(1)
        );
    if (malformed.status != spin::FasterTheoremConditionStatus::Unknown)
    {
        return fail("wrong multiplier dimension was not rejected as unknown");
    }

    Eigen::VectorXd beta(3);
    beta << 2.0, -3.0, 1.0;
    Eigen::VectorXd admissibleVariation(3);
    admissibleVariation << -0.5, 0.2, -0.25;
    Eigen::VectorXd cellMeasures(3);
    cellMeasures << 0.5, 2.0, 1.5;
    const spin::FasterTheoremOneVariationAudit strictVariation =
        spin::auditFasterTheoremOneDiscreteVariation(
            beta, admissibleVariation, cellMeasures, 4.0
        );
    if (
        strictVariation.status != spin::FasterTheoremConditionStatus::SatisfiedNumerically ||
        std::abs(strictVariation.weightedVariationL1 - 1.025) > 1e-14 ||
        std::abs(strictVariation.integralVariationTimesBeta + 2.075) > 1e-14 ||
        std::abs(strictVariation.objectiveIncrease - 8.3) > 1e-14 ||
        std::abs(strictVariation.strictIncreaseLowerBound - 4.1) > 1e-14 ||
        strictVariation.proofScope !=
            "FINITE_PARTITION_SIGN_CERTIFICATE_ONLY_NOT_A_THEOREM_PROOF"
    )
    {
        return fail("admissible finite partition variation did not certify strict increase");
    }

    Eigen::VectorXd forbiddenVariation(3);
    forbiddenVariation << 0.1, 0.0, 0.0;
    const spin::FasterTheoremOneVariationAudit forbidden =
        spin::auditFasterTheoremOneDiscreteVariation(
            beta, forbiddenVariation, cellMeasures, 4.0
        );
    if (forbidden.status != spin::FasterTheoremConditionStatus::Violated)
    {
        return fail("variation against occupied-region sign condition was accepted");
    }

    beta(1) = 0.0;
    const spin::FasterTheoremOneVariationAudit zeroSetCell =
        spin::auditFasterTheoremOneDiscreteVariation(
            beta, admissibleVariation, cellMeasures, 4.0
        );
    if (zeroSetCell.status != spin::FasterTheoremConditionStatus::Unknown)
    {
        return fail("positive-measure discrete beta zero cell was not reported as unknown");
    }

    std::cout << "[PASS] Eq. 7 numerical hypotheses accept exact fixture\n";
    std::cout << "[PASS] Eq. 7 numerical hypotheses reject infeasible and nonpositive-mu fixtures\n";
    std::cout << "[PASS] finite signed variation certificate proves a strict discrete increase\n";
    std::cout << "[PASS] finite beta-zero cells never stand in for a measure-zero proof\n";
    std::cout << "[PASS] theorem audit scope remains numerical only\n";
    return EXIT_SUCCESS;
}
