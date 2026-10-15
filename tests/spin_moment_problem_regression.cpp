#include "spin/FasterTheoremConditions.h"
#include "spin/QuadricCoefficients.h"
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


spin::SpinMomentEvaluation independentBuoyantEvaluation(
    const Vector10& raw,
    const double targetVolume,
    const Eigen::Vector2d& horizontalCenterOfBuoyancy
)
{
    spin::SpinMomentEvaluation result;
    const double volume = raw(0);
    result.objective = raw(3) / volume;
    result.objectiveGradient(0) = -raw(3) / (volume * volume);
    result.objectiveGradient(3) = 1.0 / volume;
    result.objectiveHessian(0, 0) = 2.0 * raw(3) /
        (volume * volume * volume);
    result.objectiveHessian(0, 3) = -1.0 / (volume * volume);
    result.objectiveHessian(3, 0) = result.objectiveHessian(0, 3);
    result.constraints.resize(3);
    result.constraints << raw(1) - horizontalCenterOfBuoyancy.x() * volume,
        raw(2) - horizontalCenterOfBuoyancy.y() * volume,
        volume - targetVolume;
    result.constraintJacobian = Eigen::MatrixXd::Zero(3, 10);
    result.constraintJacobian(0, 0) = -horizontalCenterOfBuoyancy.x();
    result.constraintJacobian(0, 1) = 1.0;
    result.constraintJacobian(1, 0) = -horizontalCenterOfBuoyancy.y();
    result.constraintJacobian(1, 2) = 1.0;
    result.constraintJacobian(2, 0) = 1.0;
    result.constraintHessians = {
        Eigen::Matrix<double, 10, 10>::Zero(),
        Eigen::Matrix<double, 10, 10>::Zero(),
        Eigen::Matrix<double, 10, 10>::Zero(),
    };
    return result;
}


spin::SpinMomentEvaluation independentYoYoEvaluation(const Vector10& raw)
{
    const double i11 = raw(5) + raw(6);
    const double i22 = raw(4) + raw(6);
    const double i33 = raw(4) + raw(5);
    const double i12 = -raw(9);
    const double numerator = i11 * i11 + i22 * i22 + 2.0 * i12 * i12;
    spin::SpinMomentEvaluation result;
    result.objective = numerator / (i33 * i33);
    result.objectiveGradient(4) =
        2.0 * i22 / (i33 * i33) - 2.0 * numerator / (i33 * i33 * i33);
    result.objectiveGradient(5) =
        2.0 * i11 / (i33 * i33) - 2.0 * numerator / (i33 * i33 * i33);
    result.objectiveGradient(6) = 2.0 * (i11 + i22) / (i33 * i33);
    result.objectiveGradient(9) = -4.0 * i12 / (i33 * i33);
    result.constraints.resize(5);
    result.constraints << raw(1), raw(2), raw(3), -raw(8), -raw(7);
    result.constraintJacobian = Eigen::MatrixXd::Zero(5, 10);
    result.constraintJacobian(0, 1) = 1.0;
    result.constraintJacobian(1, 2) = 1.0;
    result.constraintJacobian(2, 3) = 1.0;
    result.constraintJacobian(3, 8) = -1.0;
    result.constraintJacobian(4, 7) = -1.0;
    result.constraintHessians = {
        Eigen::Matrix<double, 10, 10>::Zero(),
        Eigen::Matrix<double, 10, 10>::Zero(),
        Eigen::Matrix<double, 10, 10>::Zero(),
        Eigen::Matrix<double, 10, 10>::Zero(),
        Eigen::Matrix<double, 10, 10>::Zero(),
    };
    return result;
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

    const spin::SpinMomentProblem staticPareto =
        spin::SpinMomentProblem::staticBalancePareto(0.25);
    evaluation = staticPareto.evaluate(r);
    if (std::abs(evaluation.objective - (-0.25 * 8.0 + 0.75 * 0.125)) > 1e-14)
    {
        return fail("static Pareto objective is incorrect");
    }

    const spin::SpinMomentProblem fidget =
        spin::SpinMomentProblem::fidgetSpinner();
    // rho=2, M=8, H=1, kappa=1/4: volume=4, center=1/8.
    // Expectations are hand-derived, independent of the AD implementation.
    const auto volumePareto=spin::SpinMomentProblem::staticBalancePareto(.25,.5);
    const auto volumeEvaluation=volumePareto.evaluate(r);
    if(std::abs(volumeEvaluation.objective+.90625)>1e-14 ||
       std::abs(volumeEvaluation.objectiveGradient[0]+.13671875)>1e-14 ||
       std::abs(volumeEvaluation.objectiveGradient[3]-.09375)>1e-14 ||
       std::abs(volumeEvaluation.objectiveHessian(0,0)-.0029296875)>1e-14 ||
       std::abs(volumeEvaluation.objectiveHessian(0,3)+.01171875)>1e-14)
        return fail("volume Pareto objective/gradient/Hessian lost density conversion");
    for(double step : {1e-3,1e-4,1e-5}) {
        for(int coordinate : {0,3}) {
            Vector10 plus=r,minus=r;plus[coordinate]+=step;minus[coordinate]-=step;
            const auto a=volumePareto.evaluate(plus),b=volumePareto.evaluate(minus);
            const double derivative=(a.objective-b.objective)/(2*step);
            if(std::abs(derivative-volumeEvaluation.objectiveGradient[coordinate])>
               1e-6*std::abs(volumeEvaluation.objectiveGradient[coordinate]))
                return fail("volume Pareto three-step objective derivative mismatch");
            const Vector10 hessianColumn=(a.objectiveGradient-b.objectiveGradient)/(2*step);
            if((hessianColumn-volumeEvaluation.objectiveHessian.col(coordinate)).norm()>
               1e-6*volumeEvaluation.objectiveHessian.col(coordinate).norm())
                return fail("volume Pareto three-step gradient derivative mismatch");
        }
    }
    if(std::abs(spin::SpinMomentProblem::staticBalance(spin::StaticObjective::MaximumMass,.5).evaluate(r).objective+4.)>1e-14)
        return fail("volume pure endpoint lost inverse density");
    evaluation = fidget.evaluate(r);
    if (std::abs(evaluation.objective + 6.5) > 1e-14)
    {
        return fail("fidget-spinner formula is incorrect");
    }

    const spin::SpinMomentProblem buoyant =
        spin::SpinMomentProblem::buoyantStability(
            8.0,
            Eigen::Vector2d(0.05, -0.0375)
        );
    evaluation = buoyant.evaluate(r);
    if (
        std::abs(evaluation.objective - 0.125) > 1e-14 ||
        evaluation.constraints.size() != 3 ||
        evaluation.constraints.norm() > 1e-14
    )
    {
        return fail("buoyant-stability objective or constraints are incorrect");
    }

    spin::BuoyantStabilityData waterline;
    waterline.displacedLiquidVolume = 2.0;
    waterline.centerOfBuoyancyTarget = Eigen::Vector3d(0.0, 0.0, -0.1);
    waterline.waterlineCentralSecondMoment =
        (Eigen::Matrix2d() << 0.5, 0.0, 0.0, 0.25).finished();
    const double metacenter = spin::buoyantMetacenterHeight(waterline);
    if (std::abs(metacenter - 0.025) > 1e-14)
    {
        return fail("buoyant metacenter formula is incorrect");
    }
    Vector10 floating = r;
    floating(3) = 0.16;
    if (
        std::abs(spin::buoyantStabilityMargin(floating, waterline) - 0.005) >
        1e-14
    )
    {
        return fail("buoyant stability margin is incorrect");
    }

    // Faster Section 4.3: fixed shell and liquid densities determine the
    // variable-domain volume target, while the waterline geometry determines
    // the metacenter.  This deliberately recomputes that contract without
    // using the production objective/constraint implementation.
    constexpr double materialDensity = 2.0;
    constexpr double liquidDensity = 3.0;
    constexpr double displacedLiquidVolume = 10.0;
    constexpr double fixedShellVolume = 4.0;
    const double prescribedVariableVolume =
        liquidDensity / materialDensity * displacedLiquidVolume - fixedShellVolume;
    Vector10 buoyantRaw = Vector10::Zero();
    buoyantRaw(0) = prescribedVariableVolume;
    buoyantRaw(3) = -4.4;
    const Eigen::Vector2d buoyancyHorizontal(0.0, 0.0);
    const spin::SpinMomentEvaluation buoyantReference =
        independentBuoyantEvaluation(
            buoyantRaw, prescribedVariableVolume, buoyancyHorizontal
        );
    const spin::SpinMomentEvaluation buoyantProduction =
        spin::SpinMomentProblem::buoyantStability(
            prescribedVariableVolume, buoyancyHorizontal
        ).evaluate(buoyantRaw);
    if (
        std::abs(prescribedVariableVolume - 11.0) > 1e-14 ||
        std::abs(buoyantReference.objective - buoyantProduction.objective) > 1e-14 ||
        !buoyantReference.objectiveGradient.isApprox(
            buoyantProduction.objectiveGradient, 1e-14
        ) ||
        !buoyantReference.objectiveHessian.isApprox(
            buoyantProduction.objectiveHessian, 1e-14
        ) ||
        !buoyantReference.constraints.isApprox(buoyantProduction.constraints, 1e-14) ||
        !buoyantReference.constraintJacobian.isApprox(
            buoyantProduction.constraintJacobian, 1e-14
        )
    )
    {
        return fail("buoyant density balance or independently derived objective/constraints are incorrect");
    }
    spin::BuoyantStabilityData buoyantWaterline;
    buoyantWaterline.displacedLiquidVolume = displacedLiquidVolume;
    buoyantWaterline.centerOfBuoyancyTarget = Eigen::Vector3d(0.0, 0.0, -0.2);
    buoyantWaterline.waterlineCentralSecondMoment =
        (Eigen::Matrix2d() << 4.0, 0.0, 0.0, 6.0).finished();
    const double independentlyDerivedMetacenter = -0.2 + 4.0 / displacedLiquidVolume;
    const double independentlyDerivedCenterHeight = buoyantRaw(3) / buoyantRaw(0);
    const double independentlyDerivedMargin =
        independentlyDerivedMetacenter - independentlyDerivedCenterHeight;
    if (
        std::abs(spin::buoyantMetacenterHeight(buoyantWaterline) -
                 independentlyDerivedMetacenter) > 1e-14 ||
        std::abs(spin::buoyantStabilityMargin(buoyantRaw, buoyantWaterline) -
                 independentlyDerivedMargin) > 1e-14 ||
        !(independentlyDerivedMargin > 0.0)
    )
    {
        return fail("buoyant metacenter construction or stable-feasibility conclusion is incorrect");
    }
    Vector10 buoyantLiteral = Vector10::Zero();
    constexpr double lambdaOne = 0.07;
    constexpr double lambdaTwo = -0.04;
    constexpr double lambdaThree = 0.2;
    buoyantLiteral << independentlyDerivedCenterHeight / prescribedVariableVolume + lambdaThree,
        lambdaOne, lambdaTwo, -1.0 / prescribedVariableVolume,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    const spin::QuadricCoefficients10 buoyantQuadric =
        spin::QuadricCoefficients10::fromUnnormalized(buoyantLiteral);
    Eigen::VectorXd buoyantMultipliers(3);
    buoyantMultipliers << lambdaOne, lambdaTwo, lambdaThree;
    const spin::FasterTheoremOneNumericalAudit buoyantAudit =
        spin::auditFasterTheoremOneNumericalConditions(
            buoyantQuadric, buoyantReference, buoyantLiteral.norm(), buoyantMultipliers
        );
    // Unlike the unconstrained static lowering problem, the prescribed volume
    // introduces lambda_3, so the buoyant interface need not pass through c.
    const double buoyantPlaneHeight = -buoyantQuadric.constant() /
        buoyantQuadric.linearVector().z();
    const Eigen::Vector3d buoyantPlanePoint(0.0, 0.0, buoyantPlaneHeight);
    const Eigen::Vector3d buoyantBelow(0.0, 0.0, buoyantPlaneHeight - 0.1);
    const Eigen::Vector3d buoyantAbove(0.0, 0.0, buoyantPlaneHeight + 0.1);
    if (
        buoyantAudit.status != spin::FasterTheoremConditionStatus::SatisfiedNumerically ||
        !(buoyantQuadric.values()(3) < 0.0) ||
        buoyantQuadric.values().segment<3>(4).norm() != 0.0 ||
        std::abs(buoyantQuadric.evaluate(buoyantPlanePoint)) > 1e-14 ||
        !(buoyantQuadric.evaluate(buoyantBelow) > 0.0) ||
        !(buoyantQuadric.evaluate(buoyantAbove) < 0.0)
    )
    {
        return fail("buoyant Eq. 13 variant lacks the published nonvertical solid-below-plane condition");
    }
    std::cout << "[PASS] buoyant density/shell balance, target-frame equilibrium, and independently derived objective/constraints agree\n";
    std::cout << "[PASS] buoyant metacenter stability, Eq. 13 variant, p!=0, and nonvertical solid-below-plane coefficients agree\n";
    std::cout << "[PASS] buoyant application scope is synthetic method validation only\n";

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

    // Faster Section 4.5: e3 is the requested spin axis (not gravity).  At
    // l=0 the central and origin inertia agree, making this a direct,
    // independent raw-moment derivation of the published yo-yo contract.
    const spin::SpinMomentEvaluation yoyoReference =
        independentYoYoEvaluation(symmetric);
    const spin::SpinMomentEvaluation yoyoProduction = yoyo.evaluate(symmetric);
    if (
        std::abs(yoyoReference.objective - yoyoProduction.objective) > 1e-14 ||
        !yoyoReference.objectiveGradient.isApprox(
            yoyoProduction.objectiveGradient, 1e-14
        ) ||
        !yoyoReference.constraints.isApprox(yoyoProduction.constraints, 1e-14) ||
        !yoyoReference.constraintJacobian.isApprox(
            yoyoProduction.constraintJacobian, 1e-14
        )
    )
    {
        return fail("yo-yo target-frame objective or inertia-block constraints differ from the direct formula");
    }
    Vector10 yoyoLiteral = Vector10::Zero();
    // With I11=I22=4 and I33=6, -grad f has the paper's a=0 and
    // A33<0 nondegenerate signature.  No multiplier is needed for this
    // symmetric, exactly feasible fixture.
    yoyoLiteral(4) = -yoyoReference.objectiveGradient(4);
    yoyoLiteral(5) = -yoyoReference.objectiveGradient(5);
    yoyoLiteral(6) = -yoyoReference.objectiveGradient(6);
    const spin::QuadricCoefficients10 yoyoQuadric =
        spin::QuadricCoefficients10::fromUnnormalized(yoyoLiteral);
    const spin::FasterTheoremOneNumericalAudit yoyoAudit =
        spin::auditFasterTheoremOneNumericalConditions(
            yoyoQuadric, yoyoReference, yoyoLiteral.norm(), Eigen::VectorXd::Zero(5)
        );
    if (
        yoyoAudit.status != spin::FasterTheoremConditionStatus::SatisfiedNumerically ||
        std::abs(yoyoQuadric.constant()) > 1e-14 ||
        !(yoyoQuadric.values()(6) < 0.0) ||
        !(yoyoQuadric.values().norm() > 0.0)
    )
    {
        return fail("yo-yo Eq. 7 witness lacks the published a=0, A33<0, p!=0 condition");
    }
    std::cout << "[PASS] yo-yo spin-frame objective and center/inertia-block constraints agree independently\n";
    std::cout << "[PASS] yo-yo Eq. 7 coefficient condition a=0, A33<0, and p!=0 agree independently\n";
    std::cout << "[PASS] yo-yo application scope is synthetic method validation only\n";

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
    // Faster Section 4.6 differs from the yo-yo contract: the spin axis is
    // gravity-aligned and c3 is deliberately free.  For this fixture c3=1/8,
    // I^c_11=I^c_22=31/8, and I^c_33=6, which gives the direct expression
    // below without invoking SpinMomentProblem's objective implementation.
    const double topDirectInertia = 31.0 / 8.0;
    const double topDirectObjective =
        2.0 * topDirectInertia * topDirectInertia / 36.0 + 2.0 / 8.0;
    if (
        std::abs(evaluation.objective - topDirectObjective) > 1e-14 ||
        evaluation.constraints.size() != 4 ||
        std::abs(evaluation.constraints(0)) > 1e-14 ||
        std::abs(evaluation.constraints(1)) > 1e-14 ||
        std::abs(evaluation.constraints(2)) > 1e-14 ||
        std::abs(evaluation.constraints(3)) > 1e-14
    )
    {
        return fail("spinning-top gravity-frame objective or free-vertical-center constraints differ from the direct formula");
    }
    Vector10 topLiteral = -evaluation.objectiveGradient;
    const spin::QuadricCoefficients10 topQuadric =
        spin::QuadricCoefficients10::fromUnnormalized(topLiteral);
    const spin::FasterTheoremOneNumericalAudit topAudit =
        spin::auditFasterTheoremOneNumericalConditions(
            topQuadric, evaluation, topLiteral.norm(), Eigen::VectorXd::Zero(4)
        );
    if (
        topAudit.status != spin::FasterTheoremConditionStatus::SatisfiedNumerically ||
        !(topQuadric.values()(6) < 0.0) ||
        !(topQuadric.values().norm() > 0.0)
    )
    {
        return fail("spinning-top Eq. 7 witness lacks the published A33<0 and p!=0 condition");
    }
    std::cout << "[PASS] spinning-top gravity-axis objective, free vertical center, and inertia-block constraints agree independently\n";
    std::cout << "[PASS] spinning-top Eq. 7 coefficient condition A33<0 and p!=0 agrees independently\n";
    std::cout << "[PASS] spinning-top application scope is synthetic method validation only\n";


    const auto topLower = spin::SpinMomentProblem::spinningTop(
        0.0, spin::SpinningTopObjective::LowerCenterOfMass);
    const auto lower = topLower.evaluate(r);
    Vector10 lowerGradient=Vector10::Zero();
    lowerGradient(0)=-1.0/64.0; lowerGradient(3)=1.0/8.0;
    Eigen::Matrix<double,10,10> lowerHessian=Eigen::Matrix<double,10,10>::Zero();
    lowerHessian(0,0)=1.0/256.0;
    lowerHessian(0,3)=lowerHessian(3,0)=-1.0/64.0;
    if(std::abs(lower.objective-0.125)>1e-14 ||
       (lower.objectiveGradient-lowerGradient).norm()>1e-14 ||
       (lower.objectiveHessian-lowerHessian).norm()>1e-14 ||
       lower.constraints.size()!=4 || std::abs(lower.constraints(2)-0.09)>1e-14 ||
       std::abs(lower.constraints(3)+0.1075)>1e-14)
        return fail("pure-height top endpoint must retain both cross-inertia constraints and exact derivatives");
    double previousError=1.0;
    for(double h : {1e-2,3e-3,1e-3}) {
        double error=0.0;
        for(int j=0;j<10;++j) {
            Vector10 plus=r,minus=r; plus(j)+=h; minus(j)-=h;
            const auto a=topLower.evaluate(plus),b=topLower.evaluate(minus);
            error=std::max(error,std::abs((a.objective-b.objective)/(2*h)-lowerGradient(j)));
            error=std::max(error,((a.objectiveGradient-b.objectiveGradient)/(2*h)-lowerHessian.col(j)).norm());
        }
        if(error>1e-6 || !(error<previousError))return fail("pure-height three-step differential convergence");
        previousError=error;
    }
    Vector10 zeroAxial=r;zeroAxial(1)=zeroAxial(2)=zeroAxial(4)=zeroAxial(5)=0;
    if(!std::isfinite(topLower.evaluate(zeroAxial).objective))
        return fail("removed inertia-ratio objective still divides by axial inertia");
    bool rejectedExtraWeight=false;
    try { (void)spin::SpinMomentProblem::spinningTop(1.0,spin::SpinningTopObjective::LowerCenterOfMass); }
    catch(const std::invalid_argument&) { rejectedExtraWeight=true; }
    if(!rejectedExtraWeight)return fail("pure endpoint silently accepted an extra weight");
    std::cout << "[PASS] pure-height top endpoint, four constraints and three-step derivatives\n";

    std::string derivativeDiagnostic;
    if (!derivativesAgree(staticPareto, r, &derivativeDiagnostic))
    {
        return fail("static Pareto " + derivativeDiagnostic);
    }
    const spin::SpinMomentEvaluation paretoLower =
        spin::SpinMomentProblem::staticBalancePareto(0.0).evaluate(r);
    const spin::SpinMomentEvaluation paretoMass =
        spin::SpinMomentProblem::staticBalancePareto(1.0).evaluate(r);
    if (
        std::abs(paretoLower.objective - staticLow.evaluate(r).objective) > 1e-14 ||
        std::abs(paretoMass.objective - staticMass.evaluate(r).objective) > 1e-14
    )
    {
        return fail("static Pareto endpoints do not match published objectives");
    }

    std::string diagnostic;
    const std::vector<spin::SpinMomentProblem> problems = {
        staticMass,
        staticLow,
        buoyant,
        fidget,
        yoyo,
        top,
        topLower
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
