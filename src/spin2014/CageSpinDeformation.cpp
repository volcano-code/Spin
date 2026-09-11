#include "spin2014/CageSpinDeformation.h"

#include "spin2014/SurfaceRawMomentDifferential.h"
#include "spin2014/Spin2014GivensProblem.h"

#include <Eigen/Core>

#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>


namespace
{


Eigen::MatrixXd uniformMeshLaplacian(
    Eigen::Index vertexCount,
    const Eigen::MatrixXi& faces
)
{
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(vertexCount, vertexCount);
    std::set<std::pair<int, int>> edges;
    for (Eigen::Index row = 0; row < faces.rows(); ++row)
    {
        for (int local = 0; local < 3; ++local)
        {
            int first = faces(row, local);
            int second = faces(row, (local + 1) % 3);
            if (first > second) std::swap(first, second);
            edges.emplace(first, second);
        }
    }
    for (const auto& edge : edges)
    {
        result(edge.first, edge.second) = -1.0;
        result(edge.second, edge.first) = -1.0;
        result(edge.first, edge.first) += 1.0;
        result(edge.second, edge.second) += 1.0;
    }
    return result;
}


physics::RawMassMoments10 fromVector(
    const Eigen::Matrix<double, 10, 1>& raw
)
{
    MassMoments result;
    result.mass = raw(0);
    result.firstMoment = raw.segment<3>(1);
    result.secondMoment <<
        raw(4), raw(9), raw(8),
        raw(9), raw(5), raw(7),
        raw(8), raw(7), raw(6);
    return physics::RawMassMoments10::fromMassMoments(result);
}


void validateWeights(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& weights,
    const deformation::HandleSet& handles
)
{
    deformation::validateSkinningWeights(
        vertices,
        weights,
        handles.size()
    );
}

void chainSpin(
    spin2014::CageSpinDeformationEvaluation& result,
    double objective,
    const Eigen::Matrix<double, 10, 1>& gradient,
    const Eigen::VectorXd& constraints,
    const Eigen::MatrixXd& jacobian
)
{
    result.spinObjective = objective;
    result.objective += objective;
    result.fillGradient += result.momentFillJacobian.transpose() * gradient;
    result.handleGradient += result.momentHandleJacobian.transpose() * gradient;
    result.constraints = constraints;
    result.fillConstraintJacobian = jacobian * result.momentFillJacobian;
    result.handleConstraintJacobian = jacobian * result.momentHandleJacobian;
    if (!std::isfinite(result.objective) || !std::isfinite(result.spinObjective) ||
        !std::isfinite(result.fillRegularizationEnergy) ||
        !std::isfinite(result.shapeEnergy) || !std::isfinite(result.displacementEnergy) ||
        !result.fillGradient.allFinite() || !result.handleGradient.allFinite() ||
        !result.constraints.allFinite() || !result.fillConstraintJacobian.allFinite() ||
        !result.handleConstraintJacobian.allFinite() ||
        !result.momentFillJacobian.allFinite() || !result.momentHandleJacobian.allFinite() ||
        !result.finalMoments.values().allFinite() || !result.deformedOuterVertices.allFinite())
        throw std::invalid_argument("nonfinite cage functional output");
    for (const auto& vertices : result.deformedVoidCellVertices)
        if (!vertices.allFinite()) throw std::invalid_argument("nonfinite cage cell geometry");
}


} // namespace


namespace spin2014
{


struct CageSpinDeformationFunctional::CellModel
{
    deformation::FreeParameterLbsModel lbs;

    CellModel(
        const spin::TriangleMesh& mesh,
        const Eigen::MatrixXd& weights,
        const deformation::HandleSet& handles
    ) : lbs(mesh.vertices, weights, handles)
    {
    }
};


CageSpinDeformationFunctional::CageSpinDeformationFunctional(
    const CageSpinDeformationRequest& request
)
    : outerSurface_(request.outerSurface),
      voidCells_(request.voidCells),
      handles_(request.handles),
      outerModel_(
          request.outerSurface.vertices,
          request.outerWeights,
          handles_
      ),
      problem_(request.problem),
      fillLaplacian_(request.fillLaplacian),
      outerLaplacian_(uniformMeshLaplacian(
          request.outerSurface.vertices.rows(), request.outerSurface.faces
      )),
      fillLaplacianWeight_(request.fillLaplacianWeight),
      biharmonicWeight_(request.biharmonicWeight),
      displacementWeight_(request.displacementWeight),
      minimumShellThickness_(request.minimumShellThickness)
{
    if (
        outerSurface_.vertices.rows() == 0 ||
        outerSurface_.faces.rows() == 0 ||
        voidCells_.empty() ||
        voidCells_.size() != request.voidCellWeights.size() ||
        fillLaplacian_.rows() != static_cast<Eigen::Index>(voidCells_.size()) ||
        fillLaplacian_.cols() != static_cast<Eigen::Index>(voidCells_.size()) ||
        !fillLaplacian_.allFinite() ||
        !fillLaplacian_.isApprox(fillLaplacian_.transpose(), 1e-12) ||
        !std::isfinite(fillLaplacianWeight_) || fillLaplacianWeight_ < 0.0 ||
        !std::isfinite(biharmonicWeight_) || biharmonicWeight_ < 0.0 ||
        !std::isfinite(displacementWeight_) || displacementWeight_ < 0.0 ||
        !std::isfinite(minimumShellThickness_) ||
        minimumShellThickness_ < 0.0
    )
    {
        throw std::invalid_argument("invalid Spin-It 2014 cage request");
    }
    validateWeights(outerSurface_.vertices, request.outerWeights, handles_);
    cellModels_.reserve(voidCells_.size());
    for (std::size_t index = 0; index < voidCells_.size(); ++index)
    {
        validateWeights(
            voidCells_[index].vertices,
            request.voidCellWeights[index],
            handles_
        );
        cellModels_.push_back(std::make_shared<const CellModel>(
            voidCells_[index], request.voidCellWeights[index], handles_
        ));
    }
}


Eigen::VectorXd CageSpinDeformationFunctional::initialHandleParameters() const
{
    return outerModel_.initialParameters();
}

Eigen::Index CageSpinDeformationFunctional::betaCount() const noexcept
{ return static_cast<Eigen::Index>(voidCells_.size()); }
bool CageSpinDeformationFunctional::paperHasAngle() const noexcept
{ return problem_.objective() != Spin2014Objective::StaticLowering; }
Eigen::Index CageSpinDeformationFunctional::paperHandleOffset() const noexcept
{ return betaCount() + (paperHasAngle() ? 1 : 0); }
Eigen::Index CageSpinDeformationFunctional::jointParameterCount() const noexcept
{ return paperHandleOffset() + handles_.optimizableParameterCount(); }
CageSpinPaperEvaluation CageSpinDeformationFunctional::evaluatePaper(
    const Eigen::VectorXd& beta, double phi,
    const Eigen::VectorXd& handleParameters) const
{
    if (!std::isfinite(phi) || (!paperHasAngle() && phi != 0.0))
        throw std::invalid_argument("invalid explicit cage angle");
    CageSpinPaperEvaluation result;
    result.deformation = evaluateGeometryRegularization(beta, handleParameters);
    auto& d = result.deformation;
    Spin2014GivensProblem::Vector11 state;
    state.head<10>() = d.finalMoments.values();
    state(10) = phi;
    const auto spinValue = Spin2014GivensProblem(
        problem_.objective(), problem_.inertiaWeight(), problem_.centerWeight()
    ).evaluate(state);
    chainSpin(d, spinValue.objective, spinValue.objectiveGradient.head<10>(),
              spinValue.constraints, spinValue.constraintJacobian.leftCols(10));
    result.angleGradient = spinValue.objectiveGradient(10);
    result.jointGradient = Eigen::VectorXd::Zero(jointParameterCount());
    result.jointGradient.head(betaCount()) = d.fillGradient;
    result.jointGradient.tail(handleParameters.size()) = d.handleGradient;
    result.jointConstraintJacobian = Eigen::MatrixXd::Zero(d.constraints.size(), jointParameterCount());
    result.jointConstraintJacobian.leftCols(betaCount()) = d.fillConstraintJacobian;
    result.jointConstraintJacobian.rightCols(handleParameters.size()) = d.handleConstraintJacobian;
    result.rawMomentJacobian = Eigen::MatrixXd::Zero(10, jointParameterCount());
    result.rawMomentJacobian.leftCols(betaCount()) = d.momentFillJacobian;
    result.rawMomentJacobian.rightCols(handleParameters.size()) = d.momentHandleJacobian;
    if (paperHasAngle())
    {
        result.jointGradient(betaCount()) = result.angleGradient;
        result.jointConstraintJacobian.col(betaCount()) = spinValue.constraintJacobian.col(10);
    }
    return result;
}
CageSpinPaperEvaluation CageSpinDeformationFunctional::evaluatePaperJoint(
    const Eigen::VectorXd& state) const
{
    if (state.size() != jointParameterCount() || !state.allFinite())
        throw std::invalid_argument("invalid assembled cage state");
    return evaluatePaper(state.head(betaCount()), paperHasAngle() ? state(betaCount()) : 0.,
                         state.tail(handles_.optimizableParameterCount()));
}


const Eigen::MatrixXi& CageSpinDeformationFunctional::outerFaces() const noexcept
{
    return outerSurface_.faces;
}


const std::vector<spin::TriangleMesh>&
CageSpinDeformationFunctional::voidCells() const noexcept
{
    return voidCells_;
}


double CageSpinDeformationFunctional::minimumShellThickness() const noexcept
{
    return minimumShellThickness_;
}


CageSpinDeformationEvaluation CageSpinDeformationFunctional::evaluate(
    const Eigen::VectorXd& beta,
    const Eigen::VectorXd& handleParameters
) const
{
    auto result = evaluateGeometryRegularization(beta, handleParameters);
    const auto spinValue = problem_.evaluate(result.finalMoments.values());
    chainSpin(result, spinValue.objective, spinValue.objectiveGradient,
              spinValue.constraints, spinValue.constraintJacobian);
    return result;
}

CageSpinDeformationEvaluation CageSpinDeformationFunctional::evaluateGeometryRegularization(
    const Eigen::VectorXd& beta,
    const Eigen::VectorXd& handleParameters
) const
{
    if (
        beta.size() != static_cast<Eigen::Index>(voidCells_.size()) ||
        !beta.allFinite() || (beta.array() < 0.0).any() ||
        (beta.array() > 1.0).any() ||
        handleParameters.size() != outerModel_.initialParameters().size() ||
        !handleParameters.allFinite()
    )
    {
        throw std::invalid_argument("invalid joint cage optimization state");
    }

    CageSpinDeformationEvaluation result;
    result.deformedOuterVertices = outerModel_.evaluate(handleParameters);
    const SurfaceRawMomentDifferential outer = differentiateSurfaceRawMoments(
        result.deformedOuterVertices, outerSurface_.faces
    );
    Eigen::Matrix<double, 10, 1> raw = outer.moments.values();
    result.momentHandleJacobian =
        outer.jacobian * outerModel_.jacobian();
    result.momentFillJacobian.resize(10, beta.size());
    result.deformedVoidCellVertices.reserve(voidCells_.size());

    for (std::size_t index = 0; index < voidCells_.size(); ++index)
    {
        const Eigen::MatrixXd deformed =
            cellModels_[index]->lbs.evaluate(handleParameters);
        result.deformedVoidCellVertices.push_back(deformed);
        const SurfaceRawMomentDifferential cell =
            differentiateSurfaceRawMoments(deformed, voidCells_[index].faces);
        const Eigen::Index column = static_cast<Eigen::Index>(index);
        raw -= beta(column) * cell.moments.values();
        result.momentFillJacobian.col(column) = -cell.moments.values();
        result.momentHandleJacobian -= beta(column) *
            cell.jacobian * cellModels_[index]->lbs.jacobian();
    }
    result.finalMoments = fromVector(raw);
    result.fillRegularizationEnergy = 0.5 * fillLaplacianWeight_ *
        beta.dot(fillLaplacian_ * beta);

    const Eigen::MatrixXd laplacianPositions =
        outerLaplacian_ * result.deformedOuterVertices;
    const Eigen::MatrixXd displacement =
        result.deformedOuterVertices - outerSurface_.vertices;
    result.shapeEnergy = biharmonicWeight_ * laplacianPositions.squaredNorm();
    result.displacementEnergy =
        displacementWeight_ * displacement.squaredNorm();
    result.objective = result.fillRegularizationEnergy +
        result.shapeEnergy + result.displacementEnergy;

    result.fillGradient =
        fillLaplacianWeight_ * fillLaplacian_ * beta;
    const Eigen::MatrixXd vertexShapeGradient =
        2.0 * biharmonicWeight_ *
            outerLaplacian_.transpose() * laplacianPositions +
        2.0 * displacementWeight_ * displacement;
    Eigen::VectorXd stackedShapeGradient(3 * vertexShapeGradient.rows());
    for (Eigen::Index row = 0; row < vertexShapeGradient.rows(); ++row)
    {
        stackedShapeGradient.segment<3>(3 * row) =
            vertexShapeGradient.row(row).transpose();
    }
    result.handleGradient =
        outerModel_.jacobian().transpose() * stackedShapeGradient;
    return result;
}


} // namespace spin2014
