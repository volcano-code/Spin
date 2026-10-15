#pragma once

#include "spin/MassMomentProblem.h"
#include "spin2014/SpinIt2014MomentProblem.h"

#include <string>

namespace spin2014
{

/** Fixed-shape reimplementation, not the 2014 voxel/Givens trajectory.
 * Unit factors convert ALREADY density-weighted solver moments to declared
 * units; no density is applied by this objective adapter. The ordinary mesh
 * evaluator currently emits target-frame moments in its input units, hence
 * both factors default to one rather than an implicit geometry normalization.
 */
struct QuadricObjectiveSettings
{
    Spin2014Objective objective = Spin2014Objective::TopMassAndLowering;
    double inertiaWeight = 1.0;
    double centerWeight = 1.0;
    std::string declaredLengthUnit = "model_length";
    std::string declaredMassUnit = "model_mass";
    double solverMassUnitInDeclaredUnits = 1.0;
    double solverLengthUnitInDeclaredUnits = 1.0;
};

const char* spin2014QuadricObjectiveProfile() noexcept;
void validateSpin2014QuadricObjectiveSettings(const QuadricObjectiveSettings& settings);
std::string spin2014QuadricObjectiveHash(const QuadricObjectiveSettings& settings);
spin::MassMomentProblem makeSpin2014QuadricMomentProblem(const QuadricObjectiveSettings& settings);

} // namespace spin2014
