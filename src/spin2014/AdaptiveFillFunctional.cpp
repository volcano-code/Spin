#include "spin2014/AdaptiveFillFunctional.h"

#include <cmath>
#include <stdexcept>
#include <utility>


namespace spin2014
{


AdaptiveFillFunctional::AdaptiveFillFunctional(
    const AdaptiveVoxelDomain& domain,
    SpinIt2014MomentProblem problem,
    double laplacianWeight
)
    : domain_(&domain),
      problem_(std::move(problem)),
      laplacianWeight_(laplacianWeight),
      momentColumns_(
          static_cast<Eigen::Index>(physics::RawMassMoments10::Count),
          static_cast<Eigen::Index>(domain.leaves().size())
      ),
      laplacian_(domain.uniformLaplacian())
{
    if (!std::isfinite(laplacianWeight) || laplacianWeight < 0.0)
    {
        throw std::invalid_argument(
            "Spin-It 2014 Laplacian weight must be finite and nonnegative"
        );
    }
    for (std::size_t index = 0; index < domain.leaves().size(); ++index)
    {
        momentColumns_.col(static_cast<Eigen::Index>(index)) =
            domain.leaves()[index].moments.values();
    }
}


AdaptiveFillEvaluation AdaptiveFillFunctional::evaluate(
    const Eigen::VectorXd& beta
) const
{
    if (
        beta.size() != momentColumns_.cols() || !beta.allFinite() ||
        (beta.array() < 0.0).any() || (beta.array() > 1.0).any()
    )
    {
        throw std::invalid_argument(
            "Spin-It 2014 fill values must be finite and lie in [0,1]"
        );
    }
    const Eigen::Matrix<double, 10, 1> raw =
        domain_->solidMoments().values() - momentColumns_ * beta;
    const physics::RawMassMoments10 finalMoments =
        physics::RawMassMoments10::fromMassMoments(
            [&raw]()
            {
                MassMoments result;
                result.mass = raw(0);
                result.firstMoment = raw.segment<3>(1);
                result.secondMoment <<
                    raw(4), raw(9), raw(8),
                    raw(9), raw(5), raw(7),
                    raw(8), raw(7), raw(6);
                return result;
            }()
        );
    const spin::SpinMomentEvaluation paper = problem_.evaluate(raw);

    AdaptiveFillEvaluation result;
    result.paperObjective = paper.objective;
    result.regularizationEnergy = laplacianEnergy(beta);
    result.objective = result.paperObjective + result.regularizationEnergy;
    result.gradient =
        -momentColumns_.transpose() * paper.objectiveGradient +
        laplacianWeight_ * laplacian_ * beta;
    result.constraints = paper.constraints;
    result.constraintJacobian =
        -paper.constraintJacobian * momentColumns_;
    result.finalMoments = finalMoments;
    return result;
}


double AdaptiveFillFunctional::laplacianEnergy(
    const Eigen::VectorXd& beta
) const
{
    if (beta.size() != laplacian_.cols() || !beta.allFinite())
    {
        throw std::invalid_argument("invalid Laplacian fill field");
    }
    return 0.5 * laplacianWeight_ * beta.dot(laplacian_ * beta);
}


const Eigen::MatrixXd& AdaptiveFillFunctional::momentColumns() const noexcept
{
    return momentColumns_;
}


const Eigen::MatrixXd& AdaptiveFillFunctional::laplacian() const noexcept
{
    return laplacian_;
}


} // namespace spin2014
