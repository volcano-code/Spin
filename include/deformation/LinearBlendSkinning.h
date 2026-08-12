#pragma once

#include "deformation/HandleSet.h"

#include <Eigen/Core>


namespace deformation
{


/** Affine map from free [translation, scale] values to full Eq. (10) H. */
struct StackedParameterReduction
{
    Eigen::MatrixXd linear;
    Eigen::VectorXd offset;
};


/** Validate the common finite/bounded/partition-of-unity LBS contract. */
void validateSkinningWeights(
    const Eigen::MatrixXd& restPositions,
    const Eigen::MatrixXd& weights,
    std::size_t handleCount,
    double weightTolerance = 1e-10
);


/**
 * Evaluate Equation (9) from rest positions:
 *
 *   x_i' = sum_j w_ij [p_j + t_j + s_j R_j (x_i - p_j)].
 *
 * Positions are rows of restPositions and weights(i,j) binds vertex i to
 * handle j. Inputs must be finite; every weight must be bounded in [0,1],
 * and each row must sum to one within weightTolerance.
 */
Eigen::MatrixXd deformLinearBlendSkinning(
    const Eigen::MatrixXd& restPositions,
    const Eigen::MatrixXd& weights,
    const HandleSet& handles,
    double weightTolerance = 1e-10
);


/**
 * Construct the fixed Equation (10) matrix M, with shape 3n x 4N.
 *
 * Handle j occupies columns [4j,4j+3]. Its first three parameters are
 * the transformed pivot p_j+t_j and the fourth is uniform scale s_j.
 * Rotations, rest pivots, rest positions, and weights are fixed in M.
 */
Eigen::MatrixXd buildReducedLbsMatrix(
    const Eigen::MatrixXd& restPositions,
    const Eigen::MatrixXd& weights,
    const HandleSet& handles,
    double weightTolerance = 1e-10
);


/** Stack every handle as [p_j+t_j.x/y/z, s_j], including locked handles. */
Eigen::VectorXd stackedHandleParameters(
    const HandleSet& handles
);


/** Build H = offset + linear * q for HandleSet::optimizationParameters(). */
StackedParameterReduction buildStackedParameterReduction(
    const HandleSet& handles
);


/** Multiply a validated 3n x 4N Equation (10) matrix by 4N parameters. */
Eigen::MatrixXd deformWithReducedLbsMatrix(
    const Eigen::MatrixXd& reducedMatrix,
    const Eigen::VectorXd& stackedParameters
);


/**
 * Bound Equation (10) operator that rejects a stale pivot/rotation basis.
 */
class ReducedLbsModel
{
public:
    ReducedLbsModel(
        const Eigen::MatrixXd& restPositions,
        const Eigen::MatrixXd& weights,
        const HandleSet& basisHandles,
        double weightTolerance = 1e-10
    );

    const Eigen::MatrixXd& matrix() const noexcept;
    const HandleSet& basisHandles() const noexcept;

    Eigen::MatrixXd evaluate(
        const HandleSet& currentHandles
    ) const;

private:
    HandleSet basisHandles_;
    Eigen::MatrixXd matrix_;
};


/**
 * Affine LBS map from only the unlocked-user parameters q to positions.
 *
 * The fixed rotations and all locked-handle poses are frozen when this model
 * is built. Positions therefore satisfy x(q) = offset + jacobian * q, with
 * q ordered exactly as HandleSet::optimizationParameters().
 */
class FreeParameterLbsModel
{
public:
    FreeParameterLbsModel(
        const Eigen::MatrixXd& restPositions,
        const Eigen::MatrixXd& weights,
        const HandleSet& basisHandles,
        double weightTolerance = 1e-10
    );

    const Eigen::MatrixXd& jacobian() const noexcept;
    const Eigen::VectorXd& offset() const noexcept;
    const HandleSet& basisHandles() const noexcept;
    Eigen::VectorXd initialParameters() const;

    Eigen::MatrixXd evaluate(
        const Eigen::VectorXd& freeParameters
    ) const;

    Eigen::MatrixXd evaluate(
        const HandleSet& currentHandles
    ) const;

private:
    HandleSet basisHandles_;
    StackedParameterReduction reduction_;
    Eigen::MatrixXd jacobian_;
    Eigen::VectorXd offset_;
};


} // namespace deformation
