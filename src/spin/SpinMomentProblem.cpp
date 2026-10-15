#include "spin/SpinMomentProblem.h"
#include "spin/MomentInertiaAlgebra.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>


namespace
{


constexpr int dimension = 10;
using Gradient = Eigen::Matrix<double, dimension, 1>;
using Hessian = Eigen::Matrix<double, dimension, dimension>;


struct SecondOrderScalar
{
    double value = 0.0;
    Gradient gradient = Gradient::Zero();
    Hessian hessian = Hessian::Zero();

    static SecondOrderScalar variable(double value, int coordinate)
    {
        SecondOrderScalar result;
        result.value = value;
        result.gradient(coordinate) = 1.0;
        return result;
    }
};


SecondOrderScalar operator+(
    const SecondOrderScalar& left,
    const SecondOrderScalar& right
)
{
    SecondOrderScalar result;
    result.value = left.value + right.value;
    result.gradient = left.gradient + right.gradient;
    result.hessian = left.hessian + right.hessian;
    return result;
}


SecondOrderScalar operator-(
    const SecondOrderScalar& left,
    const SecondOrderScalar& right
)
{
    SecondOrderScalar result;
    result.value = left.value - right.value;
    result.gradient = left.gradient - right.gradient;
    result.hessian = left.hessian - right.hessian;
    return result;
}


SecondOrderScalar operator-(const SecondOrderScalar& value)
{
    SecondOrderScalar result;
    result.value = -value.value;
    result.gradient = -value.gradient;
    result.hessian = -value.hessian;
    return result;
}


SecondOrderScalar operator*(
    const SecondOrderScalar& left,
    const SecondOrderScalar& right
)
{
    SecondOrderScalar result;
    result.value = left.value * right.value;
    result.gradient =
        left.gradient * right.value + right.gradient * left.value;
    result.hessian =
        left.hessian * right.value + right.hessian * left.value +
        left.gradient * right.gradient.transpose() +
        right.gradient * left.gradient.transpose();
    return result;
}


SecondOrderScalar operator*(double scale, const SecondOrderScalar& value)
{
    SecondOrderScalar result;
    result.value = scale * value.value;
    result.gradient = scale * value.gradient;
    result.hessian = scale * value.hessian;
    return result;
}


SecondOrderScalar reciprocal(const SecondOrderScalar& value)
{
    const double inverse = 1.0 / value.value;
    SecondOrderScalar result;
    result.value = inverse;
    result.gradient = -inverse * inverse * value.gradient;
    result.hessian =
        2.0 * inverse * inverse * inverse *
            (value.gradient * value.gradient.transpose()) -
        inverse * inverse * value.hessian;
    return result;
}


SecondOrderScalar operator/(
    const SecondOrderScalar& numerator,
    const SecondOrderScalar& denominator
)
{
    return numerator * reciprocal(denominator);
}


SecondOrderScalar squared(const SecondOrderScalar& value)
{
    return value * value;
}


using CentralInertia = spin::detail::CentralInertia10<SecondOrderScalar>;


CentralInertia centralInertia(
    const std::vector<SecondOrderScalar>& r
)
{
    return spin::detail::centralInertia10(r);
}


SecondOrderScalar yoyoObjective(const CentralInertia& inertia)
{
    return spin::detail::transverseInertiaRatio(inertia);
}


spin::SpinMomentEvaluation package(
    const SecondOrderScalar& objective,
    const std::vector<SecondOrderScalar>& constraints
)
{
    spin::SpinMomentEvaluation result;
    result.objective = objective.value;
    result.objectiveGradient = objective.gradient;
    result.objectiveHessian = objective.hessian;
    result.constraints.resize(
        static_cast<Eigen::Index>(constraints.size())
    );
    result.constraintJacobian.resize(
        static_cast<Eigen::Index>(constraints.size()),
        dimension
    );
    result.constraintHessians.reserve(constraints.size());
    for (std::size_t row = 0; row < constraints.size(); ++row)
    {
        result.constraints(static_cast<Eigen::Index>(row)) =
            constraints[row].value;
        result.constraintJacobian.row(static_cast<Eigen::Index>(row)) =
            constraints[row].gradient.transpose();
        result.constraintHessians.push_back(constraints[row].hessian);
    }
    return result;
}


} // namespace


namespace spin
{


SpinMomentProblem::SpinMomentProblem(
    SpinMode mode,
    StaticObjective staticObjective,
    double centerOfMassWeight,
    double staticParetoWeight,
    double targetMass,
    const Eigen::Vector2d& horizontalCenterOfBuoyancy
)
    : mode_(mode),
      staticObjective_(staticObjective),
      centerOfMassWeight_(centerOfMassWeight),
      staticParetoWeight_(staticParetoWeight),
      targetMass_(targetMass),
      horizontalCenterOfBuoyancy_(horizontalCenterOfBuoyancy)
{
}


SpinMomentProblem SpinMomentProblem::buoyantStability(
    double targetMass,
    const Eigen::Vector2d& horizontalCenterOfBuoyancy
)
{
    if (
        !std::isfinite(targetMass) || targetMass <= 0.0 ||
        !horizontalCenterOfBuoyancy.allFinite()
    )
    {
        throw std::invalid_argument(
            "buoyant target mass must be positive and center must be finite"
        );
    }
    return SpinMomentProblem(
        SpinMode::BuoyantStability,
        StaticObjective::MaximumMass,
        0.0,
        0.0,
        targetMass,
        horizontalCenterOfBuoyancy
    );
}


SpinMomentProblem SpinMomentProblem::staticBalance(StaticObjective objective, double massScale)
{
    if (objective == StaticObjective::ParetoTradeoff)
    {
        throw std::invalid_argument(
            "Pareto static balance requires an explicit mass weight"
        );
    }
    if(!std::isfinite(massScale)||massScale<=0.0)
        throw std::invalid_argument("static objective mass scale must be finite and positive");
    auto result=SpinMomentProblem(SpinMode::StaticBalance, objective, 0.0, 0.0);
    result.staticMassScale_=massScale;
    return result;
}


SpinMomentProblem SpinMomentProblem::staticBalancePareto(double massWeight, double massScale)
{
    if (!std::isfinite(massWeight) || massWeight < 0.0 || massWeight > 1.0)
    {
        throw std::invalid_argument(
            "static Pareto mass weight must lie in [0,1]"
        );
    }
    if(!std::isfinite(massScale)||massScale<=0.0)
        throw std::invalid_argument("static objective mass scale must be finite and positive");
    auto result=SpinMomentProblem(
        SpinMode::StaticBalance,
        StaticObjective::ParetoTradeoff,
        0.0,
        massWeight
    );
    result.staticMassScale_=massScale;
    return result;
}


SpinMomentProblem SpinMomentProblem::fidgetSpinner()
{
    return SpinMomentProblem(
        SpinMode::FidgetSpinner,
        StaticObjective::MaximumMass,
        0.0,
        0.0
    );
}


SpinMomentProblem SpinMomentProblem::yoYo()
{
    return SpinMomentProblem(
        SpinMode::YoYo,
        StaticObjective::MaximumMass,
        0.0,
        0.0
    );
}


SpinMomentProblem SpinMomentProblem::spinningTop(
    double centerOfMassWeight
)
{
    return spinningTop(centerOfMassWeight, SpinningTopObjective::InertiaAndLowering);
}


SpinMomentProblem SpinMomentProblem::spinningTop(
    double centerOfMassWeight, SpinningTopObjective objective
)
{
    if (!std::isfinite(centerOfMassWeight) || centerOfMassWeight < 0.0)
    {
        throw std::invalid_argument(
            "spinning-top center-of-mass weight must be finite and nonnegative"
        );
    }
    if (objective != SpinningTopObjective::InertiaAndLowering &&
        objective != SpinningTopObjective::LowerCenterOfMass)
        throw std::invalid_argument("unsupported spinning-top objective");
    if (objective == SpinningTopObjective::LowerCenterOfMass && centerOfMassWeight != 0.0)
        throw std::invalid_argument("pure-height endpoint requires zero extra center-of-mass weight");
    auto result = SpinMomentProblem(
        SpinMode::SpinningTop,
        StaticObjective::MaximumMass,
        centerOfMassWeight,
        0.0
    );
    result.spinningTopObjective_ = objective;
    return result;
}


SpinMode SpinMomentProblem::mode() const noexcept
{
    return mode_;
}


StaticObjective SpinMomentProblem::staticObjective() const noexcept
{
    return staticObjective_;
}


double SpinMomentProblem::centerOfMassWeight() const noexcept
{
    return centerOfMassWeight_;
}


double SpinMomentProblem::staticParetoWeight() const noexcept
{
    return staticParetoWeight_;
}


SpinMomentEvaluation SpinMomentProblem::evaluate(
    const Vector10& rawMoments
) const
{
    if (!rawMoments.allFinite())
    {
        throw std::domain_error("raw moments must be finite");
    }
    if (rawMoments(0) <= 0.0)
    {
        throw std::domain_error("raw moments require positive mass/volume");
    }

    std::vector<SecondOrderScalar> r;
    r.reserve(dimension);
    for (int coordinate = 0; coordinate < dimension; ++coordinate)
    {
        r.push_back(
            SecondOrderScalar::variable(rawMoments(coordinate), coordinate)
        );
    }

    std::vector<SecondOrderScalar> constraints;
    SecondOrderScalar objective;
    switch (mode_)
    {
    case SpinMode::StaticBalance:
        constraints = {r[1], r[2]};
        if (staticObjective_ == StaticObjective::MaximumMass)
        {
            objective = -staticMassScale_ * r[0];
        }
        else
        if (staticObjective_ == StaticObjective::LowerCenterOfMass)
        {
            objective = r[3] / r[0];
        }
        else
        {
            objective = -(staticParetoWeight_ * staticMassScale_) * r[0] +
                (1.0 - staticParetoWeight_) * (r[3] / r[0]);
        }
        break;

    case SpinMode::BuoyantStability:
    {
        SecondOrderScalar prescribedMass;
        prescribedMass.value = targetMass_;
        objective = r[3] / r[0];
        constraints = {
            r[1] - horizontalCenterOfBuoyancy_.x() * r[0],
            r[2] - horizontalCenterOfBuoyancy_.y() * r[0],
            r[0] - prescribedMass
        };
        break;
    }

    case SpinMode::FidgetSpinner:
        constraints = {r[1], r[2]};
        objective = -(r[4] + r[5]);
        break;

    case SpinMode::YoYo:
    {
        const CentralInertia inertia = centralInertia(r);
        if (std::abs(inertia.i33.value) <= 1e-14)
        {
            throw std::domain_error(
                "yo-yo axial central inertia denominator is zero"
            );
        }
        objective = yoyoObjective(inertia);
        constraints = {r[1], r[2], r[3], inertia.i13, inertia.i23};
        break;
    }

    case SpinMode::SpinningTop:
    {
        const CentralInertia inertia = centralInertia(r);
        constraints = {r[1], r[2], inertia.i13, inertia.i23};
        if (spinningTopObjective_ == SpinningTopObjective::LowerCenterOfMass)
        {
            // Figure 12 endpoint: do not evaluate the removed inertia ratio,
            // or replace this four-constraint problem with static balance.
            objective = r[3] / r[0];
            break;
        }
        if (std::abs(inertia.i33.value) <= 1e-14)
        {
            throw std::domain_error(
                "spinning-top axial central inertia denominator is zero"
            );
        }
        objective = yoyoObjective(inertia) +
            centerOfMassWeight_ * (r[3] / r[0]);
        break;
    }
    }
    return package(objective, constraints);
}


double buoyantMetacenterHeight(const BuoyantStabilityData& buoyancy)
{
    if (
        !std::isfinite(buoyancy.displacedLiquidVolume) ||
        buoyancy.displacedLiquidVolume <= 0.0 ||
        !buoyancy.centerOfBuoyancyTarget.allFinite() ||
        !buoyancy.waterlineCentralSecondMoment.allFinite() ||
        !buoyancy.waterlineCentralSecondMoment.isApprox(
            buoyancy.waterlineCentralSecondMoment.transpose(),
            1e-12)
    )
    {
        throw std::invalid_argument("invalid buoyant metacenter inputs");
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigen(
        buoyancy.waterlineCentralSecondMoment
    );
    if (
        eigen.info() != Eigen::Success ||
        eigen.eigenvalues().minCoeff() <= 0.0
    )
    {
        throw std::invalid_argument(
            "waterline central second moment must be positive definite"
        );
    }
    return buoyancy.centerOfBuoyancyTarget.z() +
        eigen.eigenvalues().minCoeff() /
        buoyancy.displacedLiquidVolume;
}


double buoyantStabilityMargin(
    const Eigen::Matrix<double, 10, 1>& rawMoments,
    const BuoyantStabilityData& buoyancy
)
{
    if (!rawMoments.allFinite() || rawMoments(0) <= 0.0)
    {
        throw std::invalid_argument(
            "buoyant stability requires positive finite mass moments"
        );
    }
    return buoyantMetacenterHeight(buoyancy) - rawMoments(3) / rawMoments(0);
}


} // namespace spin
