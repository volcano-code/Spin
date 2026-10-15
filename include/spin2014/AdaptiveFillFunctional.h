#pragma once

#include "physics/RawMassMoments10.h"
#include "spin/SpinDesignRequest.h"
#include "spin2014/AdaptiveVoxelDomain.h"
#include "spin2014/SpinIt2014MomentProblem.h"
#include "spin2014/Spin2014GivensProblem.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>


namespace spin2014
{
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
class SourceAdaptiveDomain;
#endif


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
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    AdaptiveFillFunctional(const SourceAdaptiveDomain& domain,
        SpinIt2014MomentProblem problem,double laplacianWeight,double density=1.0,
        const spin::SpinTargetFrame& targetFrame = {});
    AdaptiveFillFunctional(const SourceAdaptiveDomain& domain,
        SpinIt2014MomentProblem problem,double laplacianWeight,
        double boundaryDensity,double interiorDensity,
        const spin::SpinTargetFrame& targetFrame = {});
#endif
    AdaptiveFillFunctional(
        const AdaptiveVoxelDomain& domain,
        SpinIt2014MomentProblem problem,
        double laplacianWeight
    );
    AdaptiveFillFunctional(
        const AdaptiveVoxelDomain& domain,
        SpinIt2014MomentProblem problem,
        double laplacianWeight,
        double boundaryDensity,
        double interiorDensity,
        const spin::SpinTargetFrame& targetFrame = {}
    );

    AdaptiveFillEvaluation evaluate(const Eigen::VectorXd& beta) const;
    AdaptiveFillEvaluation evaluatePaper(const Eigen::VectorXd& beta, double phi) const;
    Eigen::Index betaCount() const noexcept { return momentColumns_.cols(); }
    bool paperHasAngle() const noexcept { return problem_.objective() != Spin2014Objective::StaticLowering; }
    double laplacianEnergy(const Eigen::VectorXd& beta) const;
    const Eigen::MatrixXd& momentColumns() const noexcept;
    const Eigen::SparseMatrix<double>& laplacian() const noexcept;

private:
    physics::RawMassMoments10 assembleMoments(const Eigen::VectorXd& beta) const;
    physics::RawMassMoments10 transformToTarget(
        const physics::RawMassMoments10& moments
    ) const;
    void setTargetFrame(const spin::SpinTargetFrame& targetFrame);
    const AdaptiveVoxelDomain* domain_ = nullptr;
    SpinIt2014MomentProblem problem_ =
        SpinIt2014MomentProblem::staticBalance();
    double laplacianWeight_ = 0.0;
    physics::RawMassMoments10 baseMoments_;
    Eigen::MatrixXd momentColumns_;
    Eigen::MatrixXd targetMomentColumns_;
    Eigen::SparseMatrix<double> laplacian_;
    Eigen::Matrix3d targetLinear_ = Eigen::Matrix3d::Identity();
    Eigen::Vector3d targetTranslation_ = Eigen::Vector3d::Zero();
};


} // namespace spin2014
