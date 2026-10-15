#pragma once
#include "spin/QuadricHollowingSolver.h"
#include <string>
#include <vector>

namespace spin {
/** Attachment for a protected double candidate checked against its owned,
 * canonical rational equation. Not a nonlinear convergence certificate.
 * The absolute limit is explicit: callers must bind it to their profile.
 */
struct SpinNewtonAuthorityLinearSolveResult {
    bool solved = false;
    unsigned precisionBits = 1024;
    std::vector<std::vector<std::string>> matrixAuthority;
    std::vector<std::string> rightHandSideAuthority;
    std::string absoluteResidualLimit;
    SpinNewtonLinearSolveResult candidate;
    std::string residualInfinityLower;
    std::string residualInfinityUpper;
    std::string backwardErrorUpper;
    std::string diagnostic;
};

SpinNewtonAuthorityLinearSolveResult solveSpinNewtonAuthorityLinearSystem(
    const std::vector<std::vector<std::string>>& matrix,
    const std::vector<std::string>& rightHandSide,
    const std::string& absoluteResidualLimit);
} // namespace spin
