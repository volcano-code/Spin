#pragma once
#include "spin2014/QuadricMomentProblem.h"
#include <array>
#include <string>
#include <vector>

namespace spin2014 {
using AuthorityMomentVector10 = std::array<std::string,10>;
using AuthorityMomentMatrix10 = std::array<AuthorityMomentVector10,10>;
/** MPFR1024 represented results, not exact analytic objective values. This
 * attachment never derives authoritative moments from compatibility doubles.
 * Inputs are already density-weighted, in the declared solver coordinates.
 */
struct QuadricAuthorityMomentEvaluation {
    unsigned precisionBits=1024;
    std::string objectiveProfile;
    std::string objectiveIdentity;
    AuthorityMomentVector10 rawMomentsAuthority;
    std::string objective;
    AuthorityMomentVector10 objectiveGradient;
    AuthorityMomentMatrix10 objectiveHessian;
    std::vector<std::string> constraints;
    std::vector<AuthorityMomentVector10> constraintJacobian;
    std::vector<AuthorityMomentMatrix10> constraintHessians;
};
QuadricAuthorityMomentEvaluation evaluateSpin2014AuthorityMomentProblem(
    const QuadricObjectiveSettings& settings,
    const AuthorityMomentVector10& rawMoments);
} // namespace spin2014
