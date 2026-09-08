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
    : AdaptiveFillFunctional(
          domain,
          std::move(problem),
          laplacianWeight,
          1.0,
          1.0
      )
{
}


AdaptiveFillFunctional::AdaptiveFillFunctional(
    const AdaptiveVoxelDomain& domain,
    SpinIt2014MomentProblem problem,
    double laplacianWeight,
    double boundaryDensity,
    double interiorDensity
)
    : domain_(&domain),
      problem_(std::move(problem)),
      laplacianWeight_(laplacianWeight),
      baseMoments_(
          (domain.solidMoments() - domain.variableMoments()).scaledDensity(
              boundaryDensity
          ) + domain.variableMoments().scaledDensity(interiorDensity)
      ),
      momentColumns_(
          static_cast<Eigen::Index>(physics::RawMassMoments10::Count),
          static_cast<Eigen::Index>(domain.leaves().size())
      ),
      laplacian_(domain.uniformLaplacian())
{
    if (
        !std::isfinite(laplacianWeight) || laplacianWeight < 0.0 ||
        !std::isfinite(boundaryDensity) || boundaryDensity <= 0.0 ||
        !std::isfinite(interiorDensity) || interiorDensity <= 0.0
    )
    {
        throw std::invalid_argument(
            "Spin-It 2014 Laplacian weight must be finite and nonnegative"
        );
    }
    for (std::size_t index = 0; index < domain.leaves().size(); ++index)
    {
        momentColumns_.col(static_cast<Eigen::Index>(index)) =
            interiorDensity * domain.leaves()[index].moments.values();
    }
}


physics::RawMassMoments10 AdaptiveFillFunctional::assembleMoments(
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
        baseMoments_.values() - momentColumns_ * beta;
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
    return finalMoments;
}

AdaptiveFillEvaluation AdaptiveFillFunctional::evaluate(const Eigen::VectorXd& beta) const
{
    const auto finalMoments = assembleMoments(beta);
    const spin::SpinMomentEvaluation paper = problem_.evaluate(finalMoments.values());

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

AdaptiveFillEvaluation AdaptiveFillFunctional::evaluatePaper(const Eigen::VectorXd& beta, double phi) const
{
    if (!std::isfinite(phi) || (!paperHasAngle() && phi != 0.))
        throw std::invalid_argument("static paper functional requires phi=0; angle must be finite");
    const auto finalMoments = assembleMoments(beta);
    Spin2014GivensProblem::Vector11 state;
    state.head<10>() = finalMoments.values(); state(10) = phi;
    const auto paper = Spin2014GivensProblem(problem_.objective(),problem_.inertiaWeight(),problem_.centerWeight()).evaluate(state);
    AdaptiveFillEvaluation result;
    result.finalMoments=finalMoments;result.paperObjective=paper.objective;
    result.regularizationEnergy=laplacianEnergy(beta);
    result.objective=paper.objective+result.regularizationEnergy;
    const Eigen::Index n=betaCount();
    result.gradient.resize(n+(paperHasAngle()?1:0));
    result.gradient.head(n)=-momentColumns_.transpose()*paper.objectiveGradient.head<10>()+laplacianWeight_*laplacian_*beta;
    result.constraints=paper.constraints;
    result.constraintJacobian.resize(paper.constraints.size(),result.gradient.size());
    result.constraintJacobian.leftCols(n)=-paper.constraintJacobian.leftCols(10)*momentColumns_;
    if(paperHasAngle()) {result.gradient(n)=paper.objectiveGradient(10);result.constraintJacobian.col(n)=paper.constraintJacobian.col(10);}
    if(!std::isfinite(result.objective)||!result.gradient.allFinite()||!result.constraintJacobian.allFinite())
        throw std::invalid_argument("nonfinite paper fill evaluation");
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


const Eigen::SparseMatrix<double>& AdaptiveFillFunctional::laplacian() const noexcept
{
    return laplacian_;
}


} // namespace spin2014
