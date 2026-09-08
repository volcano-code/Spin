#pragma once
#include "spin2014/SpinIt2014MomentProblem.h"
#include <vector>

namespace spin2014 {
enum class Spin2014Formulation { LegacyAngleEliminated, PaperGivens };
const char* spin2014FormulationToken(Spin2014Formulation formulation) noexcept;

struct Spin2014GivensEvaluation {
    double objective = 0.;
    Eigen::Matrix<double,11,1> objectiveGradient;
    Eigen::Matrix<double,11,11> objectiveHessian;
    Eigen::VectorXd constraints;
    Eigen::MatrixXd constraintJacobian;
    std::vector<Eigen::Matrix<double,11,11>> constraintHessians;
};

/** Paper Eqs. 5--8 in raw moment order, with finite unbounded phi at index 10.
 * Static uses no angle: the unified evaluator has zero phi derivatives.
 */
class Spin2014GivensProblem {
public:
    using Vector11 = Eigen::Matrix<double,11,1>;
    Spin2014GivensProblem(Spin2014Objective objective, double inertiaWeight=1., double centerWeight=1.);
    Spin2014GivensEvaluation evaluate(const Vector11& state) const;
    bool hasAngle() const noexcept { return objective_ != Spin2014Objective::StaticLowering; }
private:
    Spin2014Objective objective_;
    double inertiaWeight_, centerWeight_;
};
}
