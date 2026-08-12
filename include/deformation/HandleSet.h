#pragma once

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <cstddef>
#include <vector>


namespace deformation
{

class GridDeformationDomain;


enum class HandleRole
{
    Support,
    User
};


/**
 * One volumetric deformation handle in its rest binding and current pose.
 *
 * fixedRotation is deliberately not part of the optimizer parameter vector:
 * the paper lets a user choose rotations, then optimizes only translation and
 * uniform scale. constrainedNodes contains indices into the weight domain.
 */
struct Handle
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Eigen::Vector3d restPivot = Eigen::Vector3d::Zero();
    Eigen::Vector3d translation = Eigen::Vector3d::Zero();
    Eigen::Matrix3d fixedRotation = Eigen::Matrix3d::Identity();
    double uniformScale = 1.0;
    HandleRole role = HandleRole::User;
    bool locked = false;
    std::vector<int> constrainedNodes;
};


using HandleVector =
    std::vector<Handle, Eigen::aligned_allocator<Handle>>;


/**
 * Bind support contacts to complete active grid cells.
 *
 * Every corner of every deterministically embedded contact cell is assigned
 * to one locked identity handle.  Consequently the contact cells, rather
 * than only their sampled points, remain fixed under interpolation.
 */
Handle makeFixedSupportHandle(
    const GridDeformationDomain& domain,
    const Eigen::MatrixXd& contactPoints
);


/**
 * Validated deformation handles with the R5 optimizer parameter contract.
 *
 * A support handle is immutable and must have the identity pose. Every
 * unlocked user handle contributes [tx, ty, tz, uniformScale] to the
 * optimization vector. setOptimizationParameters validates the complete
 * input before changing any handle.
 */
class HandleSet
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    explicit HandleSet(HandleVector handles);

    std::size_t size() const noexcept;
    Eigen::Index fullParameterCount() const noexcept;
    Eigen::Index optimizableParameterCount() const noexcept;

    const Handle& at(std::size_t index) const;
    const Handle& operator[](std::size_t index) const;
    const HandleVector& handles() const noexcept;

    std::vector<std::vector<int>> constraintNodeSets() const;

    Eigen::VectorXd optimizationParameters() const;

    void setOptimizationParameters(
        const Eigen::VectorXd& parameters
    );

    /** Compare only the rest binding, not the current handle pose. */
    bool sameBindingAs(
        const HandleSet& other,
        double tolerance = 1e-12
    ) const;

    /** Binding equality plus the fixed rotations baked into Equation (10). */
    bool sameLinearBasisAs(
        const HandleSet& other,
        double tolerance = 1e-12
    ) const;

private:
    HandleVector handles_;
};


} // namespace deformation
