#include "spin2014/QuadricMomentProblem.h"
#include "manufacturing/SlicerValidation.h"

#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace spin2014
{
namespace
{
SpinIt2014MomentProblem declaredProblem(const QuadricObjectiveSettings& settings)
{
    return SpinIt2014MomentProblem::fromObjective(
        settings.objective, settings.inertiaWeight, settings.centerWeight);
}

Eigen::Matrix<double,10,1> momentScales(const QuadricObjectiveSettings& settings)
{
    Eigen::Matrix<double,10,1> d;
    const double m=settings.solverMassUnitInDeclaredUnits;
    const double l=settings.solverLengthUnitInDeclaredUnits;
    d(0)=m; d.segment<3>(1).setConstant(m*l);
    d.tail<6>().setConstant(m*l*l);
    return d;
}
}

const char* spin2014QuadricObjectiveProfile() noexcept
{
    return "SPIN2014_OBJECTIVES_QUADRIC_V1";
}

void validateSpin2014QuadricObjectiveSettings(const QuadricObjectiveSettings& settings)
{
    (void)declaredProblem(settings); // Reuse the original selector/weight validation.
    if (settings.declaredLengthUnit.empty() || settings.declaredMassUnit.empty() ||
        !std::isfinite(settings.solverMassUnitInDeclaredUnits) ||
        !std::isfinite(settings.solverLengthUnitInDeclaredUnits) ||
        settings.solverMassUnitInDeclaredUnits<=0 || settings.solverLengthUnitInDeclaredUnits<=0)
        throw std::invalid_argument("declared moment units and positive finite unit factors are required");
    const auto d=momentScales(settings);
    if (!d.allFinite() || (d.array()<=0).any())
        throw std::invalid_argument("moment unit factors overflow or underflow");
}

std::string spin2014QuadricObjectiveHash(const QuadricObjectiveSettings& settings)
{
    validateSpin2014QuadricObjectiveSettings(settings);
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << spin2014QuadricObjectiveProfile() << ';'
        << spin2014ObjectiveToken(settings.objective) << ';'
        << spin::spinModeToken(spin2014ObjectiveMode(settings.objective)) << ';'
        << std::hexfloat << settings.inertiaWeight << ';' << settings.centerWeight << ';'
        << settings.solverMassUnitInDeclaredUnits << ';' << settings.solverLengthUnitInDeclaredUnits << ';'
        << settings.declaredMassUnit.size() << ':' << settings.declaredMassUnit << ';'
        << settings.declaredLengthUnit.size() << ':' << settings.declaredLengthUnit;
    return manufacturing::computeContentSha256(out.str());
}

spin::MassMomentProblem makeSpin2014QuadricMomentProblem(const QuadricObjectiveSettings& settings)
{
    validateSpin2014QuadricObjectiveSettings(settings);
    const auto problem=declaredProblem(settings);
    const Eigen::Matrix<double,10,10> d=momentScales(settings).asDiagonal();
    return spin::MassMomentProblem([problem,d](const spin::MassMomentProblem::Vector10& moments) {
        auto result=problem.evaluate(d*moments);
        result.objectiveGradient=(d*result.objectiveGradient).eval();
        result.objectiveHessian=(d*result.objectiveHessian*d).eval();
        result.constraintJacobian=(result.constraintJacobian*d).eval();
        for(auto& h:result.constraintHessians) h=(d*h*d).eval();
        return result;
    });
}

} // namespace spin2014
