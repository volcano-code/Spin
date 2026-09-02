#include "spin/SpinMomentProblem.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


using Vector10 = Eigen::Matrix<double, 10, 1>;


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-moment-problem " << message << '\n';
    return EXIT_FAILURE;
}


bool throwsDomain(const std::function<void()>& operation)
{
    try
    {
        operation();
    }
    catch (const std::domain_error&)
    {
        return true;
    }
    return false;
}


Vector10 fixture()
{
    Vector10 r;
    r << 8.0, 0.4, -0.3, 1.0, 3.0, 3.5, 1.0, 0.07, -0.04, 0.1;
    return r;
}


bool derivativesAgree(
    const spin::SpinMomentProblem& problem,
    const Vector10& r,
    std::string* diagnostic
)
{
    const spin::SpinMomentEvaluation analytic = problem.evaluate(r);
    constexpr double gradientStep = 1e-6;
    constexpr double hessianStep = 2e-5;
    Eigen::Matrix<double, 10, 1> numericGradient;
    Eigen::Matrix<double, 10, 10> numericHessian;
    Eigen::MatrixXd numericJacobian(
        analytic.constraints.size(),
        10
    );
    std::vector<Eigen::Matrix<double, 10, 10>> numericConstraintHessians(
        static_cast<std::size_t>(analytic.constraints.size()),
        Eigen::Matrix<double, 10, 10>::Zero()
    );
    for (int column = 0; column < 10; ++column)
    {
        Vector10 plus = r;
        Vector10 minus = r;
        plus(column) += gradientStep;
        minus(column) -= gradientStep;
        const spin::SpinMomentEvaluation plusEvaluation =
            problem.evaluate(plus);
        const spin::SpinMomentEvaluation minusEvaluation =
            problem.evaluate(minus);
        numericGradient(column) =
            (plusEvaluation.objective - minusEvaluation.objective) /
            (2.0 * gradientStep);
        if (analytic.constraints.size() > 0)
        {
            numericJacobian.col(column) =
                (plusEvaluation.constraints - minusEvaluation.constraints) /
                (2.0 * gradientStep);
        }

        plus = r;
        minus = r;
        plus(column) += hessianStep;
        minus(column) -= hessianStep;
        numericHessian.col(column) =
            (problem.evaluate(plus).objectiveGradient -
             problem.evaluate(minus).objectiveGradient) /
            (2.0 * hessianStep);
        const spin::SpinMomentEvaluation plusConstraintEvaluation =
            problem.evaluate(plus);
        const spin::SpinMomentEvaluation minusConstraintEvaluation =
            problem.evaluate(minus);
        for (Eigen::Index row = 0; row < analytic.constraints.size(); ++row)
        {
            numericConstraintHessians[static_cast<std::size_t>(row)].col(
                column
            ) = (
                plusConstraintEvaluation.constraintJacobian.row(row) -
                minusConstraintEvaluation.constraintJacobian.row(row)
            ).transpose() / (2.0 * hessianStep);
        }
    }
    if (!analytic.objectiveGradient.isApprox(numericGradient, 2e-7))
    {
        *diagnostic = "objective gradient finite-difference mismatch";
        return false;
    }
    if (!analytic.objectiveHessian.isApprox(numericHessian, 2e-6))
    {
        *diagnostic = "objective Hessian finite-difference mismatch";
        return false;
    }
    if (!analytic.constraintJacobian.isApprox(numericJacobian, 2e-7))
    {
        *diagnostic = "constraint Jacobian finite-difference mismatch";
        return false;
    }
    if (analytic.constraintHessians.size() !=
        numericConstraintHessians.size())
    {
        *diagnostic = "constraint Hessian count mismatch";
        return false;
    }
    for (std::size_t row = 0; row < analytic.constraintHessians.size(); ++row)
    {
        if (!analytic.constraintHessians[row].isApprox(
                numericConstraintHessians[row],
                2e-6
            ))
        {
            *diagnostic = "constraint Hessian finite-difference mismatch";
            return false;
        }
    }
    return true;
}


} // namespace


int main()
{
    Vector10 r = fixture();

    const spin::SpinMomentProblem staticMass =
        spin::SpinMomentProblem::staticBalance(
            spin::StaticObjective::MaximumMass
        );
    spin::SpinMomentEvaluation evaluation = staticMass.evaluate(r);
    if (
        std::abs(evaluation.objective + 8.0) > 1e-14 ||
        evaluation.constraints.size() != 2 ||
        std::abs(evaluation.constraints(0) - r(1)) > 1e-14 ||
        std::abs(evaluation.constraints(1) - r(2)) > 1e-14
    )
    {
        return fail("static maximum-mass formula is incorrect");
    }

    const spin::SpinMomentProblem staticLow =
        spin::SpinMomentProblem::staticBalance(
            spin::StaticObjective::LowerCenterOfMass
        );
    evaluation = staticLow.evaluate(r);
    if (std::abs(evaluation.objective - 0.125) > 1e-14)
    {
        return fail("static lower-CoM formula is incorrect");
    }

    const spin::SpinMomentProblem fidget =
        spin::SpinMomentProblem::fidgetSpinner();
    evaluation = fidget.evaluate(r);
    if (std::abs(evaluation.objective + 6.5) > 1e-14)
    {
        return fail("fidget-spinner formula is incorrect");
    }

    Vector10 symmetric;
    symmetric << 8.0, 0.0, 0.0, 0.0, 3.0, 3.0, 1.0,
        0.0, 0.0, 0.0;
    const spin::SpinMomentProblem yoyo = spin::SpinMomentProblem::yoYo();
    evaluation = yoyo.evaluate(symmetric);
    if (
        std::abs(evaluation.objective - 8.0 / 9.0) > 1e-14 ||
        evaluation.constraints.size() != 5 ||
        evaluation.constraints.norm() > 1e-14
    )
    {
        return fail("yo-yo objective or constraints are incorrect");
    }

    Vector10 topMoments = symmetric;
    topMoments(3) = 1.0;
    const spin::SpinMomentProblem top =
        spin::SpinMomentProblem::spinningTop(2.0);
    evaluation = top.evaluate(topMoments);
    const double topExpected =
        2.0 * 3.875 * 3.875 / 36.0 + 2.0 * 0.125;
    if (
        std::abs(evaluation.objective - topExpected) > 1e-14 ||
        evaluation.constraints.size() != 4 ||
        evaluation.constraints.norm() > 1e-14
    )
    {
        return fail("spinning-top objective or constraints are incorrect");
    }

    std::string diagnostic;
    const std::vector<spin::SpinMomentProblem> problems = {
        staticMass,
        staticLow,
        fidget,
        yoyo,
        top
    };
    for (const spin::SpinMomentProblem& problem : problems)
    {
        if (!derivativesAgree(problem, r, &diagnostic))
        {
            return fail(
                std::string(spin::spinModeToken(problem.mode())) + ": " +
                diagnostic
            );
        }
    }

    Vector10 invalid = r;
    invalid(0) = 0.0;
    if (!throwsDomain([&] { yoyo.evaluate(invalid); }))
    {
        return fail("zero mass was accepted for a central-moment problem");
    }
    invalid = r;
    invalid(4) = 0.0;
    invalid(5) = 0.0;
    invalid(1) = 0.0;
    invalid(2) = 0.0;
    if (!throwsDomain([&] { yoyo.evaluate(invalid); }))
    {
        return fail("zero axial inertia denominator was accepted");
    }

    std::cout << "[PASS] scenario=spin-moment-problem\n";
    return EXIT_SUCCESS;
}
