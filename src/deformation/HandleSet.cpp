#include "deformation/HandleSet.h"

#include "deformation/GridDeformationDomain.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <set>
#include <utility>


namespace deformation
{
namespace
{


constexpr double ROTATION_TOLERANCE = 1e-10;


bool isProperRotation(const Eigen::Matrix3d& rotation)
{
    if (!rotation.allFinite())
    {
        return false;
    }

    const Eigen::Matrix3d gram =
        rotation.transpose() * rotation;

    return
        gram.isApprox(
            Eigen::Matrix3d::Identity(),
            ROTATION_TOLERANCE
        ) &&
        std::abs(rotation.determinant() - 1.0) <=
            ROTATION_TOLERANCE;
}


void canonicalizeConstraints(Handle& handle)
{
    std::sort(
        handle.constrainedNodes.begin(),
        handle.constrainedNodes.end()
    );

    handle.constrainedNodes.erase(
        std::unique(
            handle.constrainedNodes.begin(),
            handle.constrainedNodes.end()
        ),
        handle.constrainedNodes.end()
    );
}


void validateHandle(const Handle& handle)
{
    if (!handle.restPivot.allFinite())
    {
        throw std::invalid_argument(
            "handle rest pivot must be finite"
        );
    }

    if (!handle.translation.allFinite())
    {
        throw std::invalid_argument(
            "handle translation must be finite"
        );
    }

    if (!isProperRotation(handle.fixedRotation))
    {
        throw std::invalid_argument(
            "handle rotation must be finite, orthonormal, and proper"
        );
    }

    if (
        !std::isfinite(handle.uniformScale) ||
        handle.uniformScale <= 0.0
    )
    {
        throw std::invalid_argument(
            "handle uniform scale must be finite and positive"
        );
    }

    if (handle.constrainedNodes.empty())
    {
        throw std::invalid_argument(
            "each handle must constrain at least one domain node"
        );
    }

    if (
        std::any_of(
            handle.constrainedNodes.begin(),
            handle.constrainedNodes.end(),
            [](int node) { return node < 0; }
        )
    )
    {
        throw std::invalid_argument(
            "handle constraint indices must be non-negative"
        );
    }

    if (handle.role == HandleRole::Support)
    {
        if (!handle.locked)
        {
            throw std::invalid_argument(
                "the support handle must be locked"
            );
        }

        if (
            !handle.translation.isZero(0.0) ||
            !handle.fixedRotation.isIdentity(0.0) ||
            handle.uniformScale != 1.0
        )
        {
            throw std::invalid_argument(
                "the support handle must retain the identity pose"
            );
        }
    }
}


bool isOptimizable(const Handle& handle)
{
    return
        handle.role == HandleRole::User &&
        !handle.locked;
}


} // namespace


Handle makeFixedSupportHandle(
    const GridDeformationDomain& domain,
    const Eigen::MatrixXd& contactPoints
)
{
    if (
        contactPoints.rows() == 0 ||
        contactPoints.cols() != 3 ||
        !contactPoints.allFinite()
    )
    {
        throw std::invalid_argument(
            "support contacts must be a nonempty finite N-by-3 matrix"
        );
    }

    std::set<int> constrainedNodes;

    for (Eigen::Index row = 0; row < contactPoints.rows(); ++row)
    {
        const TrilinearEmbedding embedding =
            domain.embed(contactPoints.row(row).transpose());

        constrainedNodes.insert(
            embedding.nodes.begin(),
            embedding.nodes.end()
        );
    }

    Handle support;
    support.restPivot = contactPoints.colwise().mean().transpose();
    support.role = HandleRole::Support;
    support.locked = true;
    support.constrainedNodes.assign(
        constrainedNodes.begin(),
        constrainedNodes.end()
    );
    return support;
}


HandleSet::HandleSet(HandleVector handles)
    :
    handles_(std::move(handles))
{
    if (handles_.empty())
    {
        throw std::invalid_argument(
            "a deformation binding requires at least one handle"
        );
    }

    for (Handle& handle : handles_)
    {
        canonicalizeConstraints(handle);
        validateHandle(handle);
    }
}


std::size_t HandleSet::size() const noexcept
{
    return handles_.size();
}


Eigen::Index HandleSet::fullParameterCount() const noexcept
{
    return 4 * static_cast<Eigen::Index>(handles_.size());
}


Eigen::Index HandleSet::optimizableParameterCount() const noexcept
{
    Eigen::Index count = 0;

    for (const Handle& handle : handles_)
    {
        if (isOptimizable(handle))
        {
            count += 4;
        }
    }

    return count;
}


const Handle& HandleSet::at(std::size_t index) const
{
    return handles_.at(index);
}


const Handle& HandleSet::operator[](std::size_t index) const
{
    return handles_.at(index);
}


const HandleVector& HandleSet::handles() const noexcept
{
    return handles_;
}


std::vector<std::vector<int>> HandleSet::constraintNodeSets() const
{
    std::vector<std::vector<int>> constraints;
    constraints.reserve(handles_.size());

    for (const Handle& handle : handles_)
    {
        constraints.push_back(handle.constrainedNodes);
    }

    return constraints;
}


Eigen::VectorXd HandleSet::optimizationParameters() const
{
    Eigen::VectorXd parameters(
        optimizableParameterCount()
    );

    Eigen::Index offset = 0;

    for (const Handle& handle : handles_)
    {
        if (!isOptimizable(handle))
        {
            continue;
        }

        parameters.segment<3>(offset) = handle.translation;
        parameters(offset + 3) = handle.uniformScale;
        offset += 4;
    }

    return parameters;
}


void HandleSet::setOptimizationParameters(
    const Eigen::VectorXd& parameters
)
{
    if (parameters.size() != optimizableParameterCount())
    {
        throw std::invalid_argument(
            "optimizer parameter vector has the wrong size"
        );
    }

    if (!parameters.allFinite())
    {
        throw std::invalid_argument(
            "optimizer parameters must be finite"
        );
    }

    for (Eigen::Index offset = 0;
         offset < parameters.size();
         offset += 4)
    {
        if (parameters(offset + 3) <= 0.0)
        {
            throw std::invalid_argument(
                "optimized handle scales must be positive"
            );
        }
    }

    // Commit only after every value has passed validation.
    Eigen::Index offset = 0;

    for (Handle& handle : handles_)
    {
        if (!isOptimizable(handle))
        {
            continue;
        }

        handle.translation = parameters.segment<3>(offset);
        handle.uniformScale = parameters(offset + 3);
        offset += 4;
    }
}


bool HandleSet::sameBindingAs(
    const HandleSet& other,
    double tolerance
) const
{
    if (!std::isfinite(tolerance) || tolerance < 0.0)
    {
        throw std::invalid_argument(
            "binding comparison tolerance must be finite and non-negative"
        );
    }

    if (handles_.size() != other.handles_.size())
    {
        return false;
    }

    for (std::size_t index = 0;
         index < handles_.size();
         ++index)
    {
        const Handle& first = handles_[index];
        const Handle& second = other.handles_[index];

        if (
            first.role != second.role ||
            first.locked != second.locked ||
            first.constrainedNodes != second.constrainedNodes ||
            (first.restPivot - second.restPivot).norm() >
                tolerance * std::max(
                    1.0,
                    std::max(
                        first.restPivot.norm(),
                        second.restPivot.norm()
                    )
                )
        )
        {
            return false;
        }
    }

    return true;
}


bool HandleSet::sameLinearBasisAs(
    const HandleSet& other,
    double tolerance
) const
{
    if (!sameBindingAs(other, tolerance))
    {
        return false;
    }

    for (std::size_t index = 0;
         index < handles_.size();
         ++index)
    {
        if (
            !handles_[index].fixedRotation.isApprox(
                other.handles_[index].fixedRotation,
                tolerance
            )
        )
        {
            return false;
        }
    }

    return true;
}


} // namespace deformation
