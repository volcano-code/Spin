#pragma once

#include "deformation/HandleSet.h"
#include "deformation/LinearBlendSkinning.h"
#include "physics/RawMassMoments10.h"
#include "spin/SpinDesignRequest.h"
#include "spin2014/SpinIt2014MomentProblem.h"

#include <Eigen/Core>

#include <memory>
#include <vector>


namespace spin2014
{


struct CageSpinDeformationRequest
{
    // Legacy dense mesh adapter: both outer and void material use unit density.
    spin::TriangleMesh outerSurface;
    std::vector<spin::TriangleMesh> voidCells;
    Eigen::MatrixXd outerWeights;
    std::vector<Eigen::MatrixXd> voidCellWeights;
    deformation::HandleVector handles;
    SpinIt2014MomentProblem problem =
        SpinIt2014MomentProblem::spinningTop(1.0);
    Eigen::MatrixXd fillLaplacian;
    double fillLaplacianWeight = 1e-3;
    double biharmonicWeight = 1e-3;
    double displacementWeight = 1e-3;
    double minimumShellThickness = 0.0;
};


struct CageSpinDeformationEvaluation
{
    double objective = 0.0;
    double spinObjective = 0.0;
    double fillRegularizationEnergy = 0.0;
    double shapeEnergy = 0.0;
    double displacementEnergy = 0.0;
    Eigen::VectorXd fillGradient;
    Eigen::VectorXd handleGradient;
    Eigen::VectorXd constraints;
    Eigen::MatrixXd fillConstraintJacobian;
    Eigen::MatrixXd handleConstraintJacobian;
    // Raw order: mass, x, y, z, xx, yy, zz, yz, xz, xy.
    // Signed derivatives of finalMoments, including the negative void term.
    Eigen::MatrixXd momentFillJacobian;
    Eigen::MatrixXd momentHandleJacobian;
    physics::RawMassMoments10 finalMoments;
    Eigen::MatrixXd deformedOuterVertices;
    std::vector<Eigen::MatrixXd> deformedVoidCellVertices;
};

// Total gradients include all regularizers. Raw moment phi column is zero.
struct CageSpinPaperEvaluation
{
    CageSpinDeformationEvaluation deformation;
    double angleGradient = 0.0;
    Eigen::VectorXd jointGradient;
    Eigen::MatrixXd jointConstraintJacobian;
    Eigen::MatrixXd rawMomentJacobian;
};

/**
 * Section 6 joint fill/deformation functional with analytic derivatives.
 * Handle rotations remain fixed; optimized parameters are translations and
 * uniform scales through the existing affine LBS map.
 * The uniform outer graph Laplacian is an IMPLEMENTATION_DEFAULT, not the
 * paper's rotation-aware Sorkine operator or recovered author discretization.
 */
class CageSpinDeformationFunctional
{
public:
    explicit CageSpinDeformationFunctional(
        const CageSpinDeformationRequest& request
    );

    Eigen::VectorXd initialHandleParameters() const;
    Eigen::Index betaCount() const noexcept;
    bool paperHasAngle() const noexcept;
    Eigen::Index paperHandleOffset() const noexcept;
    Eigen::Index jointParameterCount() const noexcept;
    // Layout: [beta, phi if paperHasAngle(), existing free handle parameters].
    // StaticLowering has no angle column and evaluatePaper requires phi == 0.
    CageSpinPaperEvaluation evaluatePaper(
        const Eigen::VectorXd& beta, double phi,
        const Eigen::VectorXd& handleParameters
    ) const;
    CageSpinPaperEvaluation evaluatePaperJoint(const Eigen::VectorXd& state) const;
    const Eigen::MatrixXi& outerFaces() const noexcept;
    const std::vector<spin::TriangleMesh>& voidCells() const noexcept;
    double minimumShellThickness() const noexcept;
    CageSpinDeformationEvaluation evaluate(
        const Eigen::VectorXd& beta,
        const Eigen::VectorXd& handleParameters
    ) const;

private:
    struct CellModel;
    CageSpinDeformationEvaluation evaluateGeometryRegularization(
        const Eigen::VectorXd& beta, const Eigen::VectorXd& handleParameters
    ) const;

    spin::TriangleMesh outerSurface_;
    std::vector<spin::TriangleMesh> voidCells_;
    deformation::HandleSet handles_;
    deformation::FreeParameterLbsModel outerModel_;
    std::vector<std::shared_ptr<const CellModel>> cellModels_;
    SpinIt2014MomentProblem problem_;
    Eigen::MatrixXd fillLaplacian_;
    Eigen::MatrixXd outerLaplacian_;
    double fillLaplacianWeight_ = 0.0;
    double biharmonicWeight_ = 0.0;
    double displacementWeight_ = 0.0;
    double minimumShellThickness_ = 0.0;
};


} // namespace spin2014
