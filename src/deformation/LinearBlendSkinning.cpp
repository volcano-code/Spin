#include "deformation/LinearBlendSkinning.h"

#include <cmath>
#include <stdexcept>


namespace deformation
{
namespace
{


void validateTolerance(double tolerance)
{
    if (!std::isfinite(tolerance) || tolerance < 0.0)
    {
        throw std::invalid_argument(
            "weight tolerance must be finite and non-negative"
        );
    }
}


void validateWeightContract(
    const Eigen::MatrixXd& restPositions,
    const Eigen::MatrixXd& weights,
    std::size_t handleCount,
    double weightTolerance
)
{
    validateTolerance(weightTolerance);

    if (restPositions.cols() != 3)
    {
        throw std::invalid_argument(
            "LBS rest positions must have exactly three columns"
        );
    }

    if (
        weights.rows() != restPositions.rows() ||
        weights.cols() != static_cast<Eigen::Index>(handleCount)
    )
    {
        throw std::invalid_argument(
            "LBS weight dimensions do not match positions and handles"
        );
    }

    if (!restPositions.allFinite() || !weights.allFinite())
    {
        throw std::invalid_argument(
            "LBS rest positions and weights must be finite"
        );
    }

    for (Eigen::Index row = 0; row < weights.rows(); ++row)
    {
        double sum = 0.0;

        for (Eigen::Index column = 0;
             column < weights.cols();
             ++column)
        {
            const double weight = weights(row, column);

            if (weight < 0.0 || weight > 1.0)
            {
                throw std::invalid_argument(
                    "LBS weights must lie in [0,1]"
                );
            }

            sum += weight;
        }

        if (std::abs(sum - 1.0) > weightTolerance)
        {
            throw std::invalid_argument(
                "every LBS weight row must sum to one"
            );
        }
    }
}


void validateSkinningInputs(
    const Eigen::MatrixXd& restPositions,
    const Eigen::MatrixXd& weights,
    const HandleSet& handles,
    double weightTolerance
)
{
    validateWeightContract(
        restPositions,
        weights,
        handles.size(),
        weightTolerance
    );
}


void validateReducedInputs(
    const Eigen::MatrixXd& reducedMatrix,
    const Eigen::VectorXd& parameters
)
{
    if (
        reducedMatrix.rows() % 3 != 0 ||
        reducedMatrix.cols() == 0 ||
        reducedMatrix.cols() % 4 != 0 ||
        reducedMatrix.cols() != parameters.size()
    )
    {
        throw std::invalid_argument(
            "reduced LBS dimensions must be 3n by 4N"
        );
    }

    if (!reducedMatrix.allFinite() || !parameters.allFinite())
    {
        throw std::invalid_argument(
            "reduced LBS matrix and parameters must be finite"
        );
    }

    for (Eigen::Index scale = 3;
         scale < parameters.size();
         scale += 4)
    {
        if (parameters(scale) <= 0.0)
        {
            throw std::invalid_argument(
                "reduced LBS scales must be positive"
            );
        }
    }
}


} // namespace


void validateSkinningWeights(
    const Eigen::MatrixXd& restPositions,
    const Eigen::MatrixXd& weights,
    std::size_t handleCount,
    double weightTolerance
)
{
    if (handleCount == 0U)
    {
        throw std::invalid_argument(
            "LBS validation requires at least one handle"
        );
    }

    validateWeightContract(
        restPositions,
        weights,
        handleCount,
        weightTolerance
    );
}


Eigen::MatrixXd deformLinearBlendSkinning(
    const Eigen::MatrixXd& restPositions,
    const Eigen::MatrixXd& weights,
    const HandleSet& handles,
    double weightTolerance
)
{
    validateSkinningInputs(
        restPositions,
        weights,
        handles,
        weightTolerance
    );

    Eigen::MatrixXd deformed =
        Eigen::MatrixXd::Zero(restPositions.rows(), 3);

    for (Eigen::Index vertex = 0;
         vertex < restPositions.rows();
         ++vertex)
    {
        const Eigen::Vector3d rest =
            restPositions.row(vertex).transpose();

        Eigen::Vector3d position = Eigen::Vector3d::Zero();

        for (std::size_t handleIndex = 0;
             handleIndex < handles.size();
             ++handleIndex)
        {
            const Handle& handle = handles[handleIndex];
            const double weight = weights(
                vertex,
                static_cast<Eigen::Index>(handleIndex)
            );

            const Eigen::Vector3d transformed =
                handle.restPivot +
                handle.translation +
                handle.uniformScale *
                    handle.fixedRotation *
                    (rest - handle.restPivot);

            position += weight * transformed;
        }

        deformed.row(vertex) = position.transpose();
    }

    if (!deformed.allFinite())
    {
        throw std::domain_error(
            "LBS evaluation produced non-finite positions"
        );
    }

    return deformed;
}


Eigen::MatrixXd buildReducedLbsMatrix(
    const Eigen::MatrixXd& restPositions,
    const Eigen::MatrixXd& weights,
    const HandleSet& handles,
    double weightTolerance
)
{
    validateSkinningInputs(
        restPositions,
        weights,
        handles,
        weightTolerance
    );

    Eigen::MatrixXd reduced = Eigen::MatrixXd::Zero(
        3 * restPositions.rows(),
        handles.fullParameterCount()
    );

    for (Eigen::Index vertex = 0;
         vertex < restPositions.rows();
         ++vertex)
    {
        const Eigen::Vector3d rest =
            restPositions.row(vertex).transpose();
        const Eigen::Index row = 3 * vertex;

        for (std::size_t handleIndex = 0;
             handleIndex < handles.size();
             ++handleIndex)
        {
            const Handle& handle = handles[handleIndex];
            const double weight = weights(
                vertex,
                static_cast<Eigen::Index>(handleIndex)
            );
            const Eigen::Index column =
                4 * static_cast<Eigen::Index>(handleIndex);

            reduced.block<3, 3>(row, column) =
                weight * Eigen::Matrix3d::Identity();

            reduced.block<3, 1>(row, column + 3) =
                weight * handle.fixedRotation *
                    (rest - handle.restPivot);
        }
    }

    return reduced;
}


Eigen::VectorXd stackedHandleParameters(
    const HandleSet& handles
)
{
    Eigen::VectorXd parameters(handles.fullParameterCount());

    for (std::size_t index = 0;
         index < handles.size();
         ++index)
    {
        const Handle& handle = handles[index];
        const Eigen::Index offset =
            4 * static_cast<Eigen::Index>(index);

        parameters.segment<3>(offset) =
            handle.restPivot + handle.translation;
        parameters(offset + 3) = handle.uniformScale;
    }

    return parameters;
}


StackedParameterReduction buildStackedParameterReduction(
    const HandleSet& handles
)
{
    StackedParameterReduction reduction;
    reduction.linear = Eigen::MatrixXd::Zero(
        handles.fullParameterCount(),
        handles.optimizableParameterCount()
    );
    reduction.offset = Eigen::VectorXd::Zero(
        handles.fullParameterCount()
    );

    Eigen::Index freeOffset = 0;

    for (std::size_t index = 0;
         index < handles.size();
         ++index)
    {
        const Handle& handle = handles[index];
        const Eigen::Index fullOffset =
            4 * static_cast<Eigen::Index>(index);

        if (handle.role == HandleRole::User && !handle.locked)
        {
            reduction.offset.segment<3>(fullOffset) = handle.restPivot;
            reduction.linear.block<4, 4>(fullOffset, freeOffset).setIdentity();
            freeOffset += 4;
        }
        else
        {
            reduction.offset.segment<3>(fullOffset) =
                handle.restPivot + handle.translation;
            reduction.offset(fullOffset + 3) = handle.uniformScale;
        }
    }

    return reduction;
}


Eigen::MatrixXd deformWithReducedLbsMatrix(
    const Eigen::MatrixXd& reducedMatrix,
    const Eigen::VectorXd& stackedParametersValue
)
{
    validateReducedInputs(
        reducedMatrix,
        stackedParametersValue
    );

    const Eigen::VectorXd stackedPositions =
        reducedMatrix * stackedParametersValue;

    Eigen::MatrixXd positions(
        reducedMatrix.rows() / 3,
        3
    );

    for (Eigen::Index vertex = 0;
         vertex < positions.rows();
         ++vertex)
    {
        positions.row(vertex) =
            stackedPositions.segment<3>(3 * vertex).transpose();
    }

    if (!positions.allFinite())
    {
        throw std::domain_error(
            "reduced LBS evaluation produced non-finite positions"
        );
    }

    return positions;
}


ReducedLbsModel::ReducedLbsModel(
    const Eigen::MatrixXd& restPositions,
    const Eigen::MatrixXd& weights,
    const HandleSet& basisHandles,
    double weightTolerance
)
    :
    basisHandles_(basisHandles),
    matrix_(buildReducedLbsMatrix(
        restPositions,
        weights,
        basisHandles,
        weightTolerance
    ))
{
}


const Eigen::MatrixXd& ReducedLbsModel::matrix() const noexcept
{
    return matrix_;
}


const HandleSet& ReducedLbsModel::basisHandles() const noexcept
{
    return basisHandles_;
}


Eigen::MatrixXd ReducedLbsModel::evaluate(
    const HandleSet& currentHandles
) const
{
    if (!basisHandles_.sameLinearBasisAs(currentHandles))
    {
        throw std::invalid_argument(
            "current handles do not match the reduced LBS basis"
        );
    }

    return deformWithReducedLbsMatrix(
        matrix_,
        stackedHandleParameters(currentHandles)
    );
}


} // namespace deformation
