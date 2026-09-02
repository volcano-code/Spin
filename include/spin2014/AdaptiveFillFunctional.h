#pragma once

#include "physics/RawMassMoments10.h"
#include "spin2014/AdaptiveVoxelDomain.h"
#include "spin2014/SpinIt2014MomentProblem.h"

#include <Eigen/Core>


namespace spin2014
{


struct AdaptiveFillEvaluation
{
    double objective = 0.0;
    double paperObjective = 0.0;
    double regularizationEnergy = 0.0;
    Eigen::VectorXd gradient;
    Eigen::VectorXd constraints;
    Eigen::MatrixXd constraintJacobian;
    physics::RawMassMoments10 finalMoments;
};


/** Section 5 relaxed fill functional and its analytic first derivatives. */
class AdaptiveFillFunctional
{
public:
    AdaptiveFillFunctional(
        const AdaptiveVoxelDomain& domain,
        SpinIt2014MomentProblem problem,
        double laplacianWeight
    );

    AdaptiveFillEvaluation evaluate(const Eigen::VectorXd& beta) const;
    double laplacianEnergy(const Eigen::VectorXd& beta) const;
    const Eigen::MatrixXd& momentColumns() const noexcept;
    const Eigen::MatrixXd& laplacian() const noexcept;

private:
    const AdaptiveVoxelDomain* domain_ = nullptr;
    SpinIt2014MomentProblem problem_ =
        SpinIt2014MomentProblem::staticBalance();
    double laplacianWeight_ = 0.0;
    Eigen::MatrixXd momentColumns_;
    Eigen::MatrixXd laplacian_;
};


} // namespace spin2014
