#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>


namespace deformation
{


/**
 * Paper Section 5.1 implicit-similarity Laplacian editing energy.
 *
 * The fixed residual operator M is built from the rest mesh using the local
 * seven-parameter linear similarity model of Sorkine et al. The unscaled
 * shape term is E_M(x)=0.5*||M*x||^2, where coordinates are stacked per
 * vertex as [x0,y0,z0,x1,...]. The paper's mesh-specific lambda is applied
 * by DeformationObjective rather than stored here.
 */
class LaplacianEditingEnergy
{
public:
    LaplacianEditingEnergy(
        const Eigen::MatrixXd& restVertices,
        const Eigen::MatrixXi& faces
    );

    const Eigen::SparseMatrix<double>&
    residualOperator() const noexcept;

    /** Explicit Q=M^T M from paper Eq. (12), primarily for diagnostics. */
    Eigen::SparseMatrix<double> quadraticMatrix() const;

    double energy(
        const Eigen::MatrixXd& vertices
    ) const;

    /** Gradient with respect to interleaved vertex coordinates. */
    Eigen::VectorXd vertexGradient(
        const Eigen::MatrixXd& vertices
    ) const;

private:
    Eigen::Index vertexCount_;
    Eigen::SparseMatrix<double> residualOperator_;
};


} // namespace deformation
