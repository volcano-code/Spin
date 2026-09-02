#include "spin2014/AdaptiveVoxelDomain.h"

#include "manufacturing/SlicerValidation.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>


namespace
{


physics::RawMassMoments10 cuboidMoments(
    const Eigen::Vector3d& minimum,
    const Eigen::Vector3d& maximum
)
{
    const Eigen::Vector3d length = maximum - minimum;
    const double mass = length.prod();
    MassMoments result;
    result.mass = mass;
    result.firstMoment = mass * 0.5 * (minimum + maximum);

    Eigen::Vector3d firstIntegral;
    Eigen::Vector3d secondIntegral;
    for (int axis = 0; axis < 3; ++axis)
    {
        firstIntegral(axis) = 0.5 * (
            maximum(axis) * maximum(axis) -
            minimum(axis) * minimum(axis)
        );
        secondIntegral(axis) = (
            maximum(axis) * maximum(axis) * maximum(axis) -
            minimum(axis) * minimum(axis) * minimum(axis)
        ) / 3.0;
    }
    result.secondMoment.setZero();
    for (int axis = 0; axis < 3; ++axis)
    {
        const int otherA = (axis + 1) % 3;
        const int otherB = (axis + 2) % 3;
        result.secondMoment(axis, axis) =
            secondIntegral(axis) * length(otherA) * length(otherB);
    }
    result.secondMoment(0, 1) =
        firstIntegral(0) * firstIntegral(1) * length(2);
    result.secondMoment(1, 0) = result.secondMoment(0, 1);
    result.secondMoment(0, 2) =
        firstIntegral(0) * firstIntegral(2) * length(1);
    result.secondMoment(2, 0) = result.secondMoment(0, 2);
    result.secondMoment(1, 2) =
        firstIntegral(1) * firstIntegral(2) * length(0);
    result.secondMoment(2, 1) = result.secondMoment(1, 2);
    return physics::RawMassMoments10::fromMassMoments(result);
}


int powerOfTwo(int exponent)
{
    if (exponent < 0 || exponent >= 30)
    {
        throw std::invalid_argument("octree level must lie in [0,29]");
    }
    return 1 << exponent;
}


bool keyLess(
    const spin2014::AdaptiveVoxelCell& left,
    const spin2014::AdaptiveVoxelCell& right
)
{
    return std::tie(
        left.level,
        left.coordinate(2),
        left.coordinate(1),
        left.coordinate(0)
    ) < std::tie(
        right.level,
        right.coordinate(2),
        right.coordinate(1),
        right.coordinate(0)
    );
}


struct Bounds
{
    Eigen::Vector3i minimum;
    Eigen::Vector3i maximum;
};


Bounds bounds(const spin2014::AdaptiveVoxelCell& cell, int maximumLevel)
{
    const int span = powerOfTwo(maximumLevel - cell.level);
    return {span * cell.coordinate, span * (cell.coordinate.array() + 1).matrix()};
}


bool faceNeighbors(
    const spin2014::AdaptiveVoxelCell& left,
    const spin2014::AdaptiveVoxelCell& right,
    int maximumLevel
)
{
    const Bounds a = bounds(left, maximumLevel);
    const Bounds b = bounds(right, maximumLevel);
    int touchingAxis = -1;
    for (int axis = 0; axis < 3; ++axis)
    {
        const bool touches =
            a.maximum(axis) == b.minimum(axis) ||
            b.maximum(axis) == a.minimum(axis);
        if (touches)
        {
            if (touchingAxis >= 0)
            {
                return false;
            }
            touchingAxis = axis;
        }
        else if (
            std::max(a.minimum(axis), b.minimum(axis)) >=
            std::min(a.maximum(axis), b.maximum(axis))
        )
        {
            return false;
        }
    }
    return touchingAxis >= 0;
}


void appendInt(std::string& bytes, int value)
{
    const std::uint32_t bits = static_cast<std::uint32_t>(value);
    for (int shift = 0; shift < 32; shift += 8)
    {
        bytes.push_back(static_cast<char>((bits >> shift) & 0xffU));
    }
}


void appendDouble(std::string& bytes, double value)
{
    std::uint64_t bits = 0U;
    static_assert(sizeof(bits) == sizeof(value));
    std::memcpy(&bits, &value, sizeof(bits));
    for (int shift = 0; shift < 64; shift += 8)
    {
        bytes.push_back(static_cast<char>((bits >> shift) & 0xffU));
    }
}


} // namespace


namespace spin2014
{


AdaptiveVoxelDomain::AdaptiveVoxelDomain(
    std::shared_ptr<const voxel::VoxelDomain> source,
    int initialLevel,
    int maximumLevel
)
    : source_(std::move(source)),
      initialLevel_(initialLevel),
      maximumLevel_(maximumLevel),
      rootSide_(powerOfTwo(maximumLevel))
{
    if (!source_)
    {
        throw std::invalid_argument("adaptive voxel source must not be null");
    }
    if (initialLevel < 0 || initialLevel > maximumLevel)
    {
        throw std::invalid_argument(
            "initial octree level must not exceed maximum level"
        );
    }
    if ((source_->grid().dimensions.array() > rootSide_).any())
    {
        throw std::invalid_argument(
            "maximum octree level cannot contain the source lattice"
        );
    }

    solidMoments_ = momentsForCell(0, Eigen::Vector3i::Zero(), false);
    variableMoments_ = momentsForCell(0, Eigen::Vector3i::Zero(), true);
    const int cellsPerAxis = powerOfTwo(initialLevel_);
    for (int z = 0; z < cellsPerAxis; ++z)
    {
        for (int y = 0; y < cellsPerAxis; ++y)
        {
            for (int x = 0; x < cellsPerAxis; ++x)
            {
                AdaptiveVoxelCell cell;
                cell.level = initialLevel_;
                cell.coordinate = Eigen::Vector3i(x, y, z);
                cell.moments = momentsForCell(
                    cell.level,
                    cell.coordinate,
                    true
                );
                if (cell.moments.mass() > 0.0)
                {
                    leaves_.push_back(cell);
                }
            }
        }
    }
    if (leaves_.empty())
    {
        throw std::invalid_argument("adaptive domain contains no carveable mass");
    }
    sortLeaves();
}


const std::shared_ptr<const voxel::VoxelDomain>&
AdaptiveVoxelDomain::source() const noexcept
{
    return source_;
}


int AdaptiveVoxelDomain::initialLevel() const noexcept { return initialLevel_; }
int AdaptiveVoxelDomain::maximumLevel() const noexcept { return maximumLevel_; }


const std::vector<AdaptiveVoxelCell>&
AdaptiveVoxelDomain::leaves() const noexcept
{
    return leaves_;
}


const physics::RawMassMoments10&
AdaptiveVoxelDomain::solidMoments() const noexcept
{
    return solidMoments_;
}


const physics::RawMassMoments10&
AdaptiveVoxelDomain::variableMoments() const noexcept
{
    return variableMoments_;
}


physics::RawMassMoments10 AdaptiveVoxelDomain::momentsForCell(
    int level,
    const Eigen::Vector3i& coordinate,
    bool carveableOnly
) const
{
    const int span = powerOfTwo(maximumLevel_ - level);
    const Eigen::Vector3i lower = span * coordinate;
    const Eigen::Vector3i upper = span * (coordinate.array() + 1).matrix();
    const Eigen::Vector3i dimensions = source_->grid().dimensions;
    physics::RawMassMoments10 result;
    for (int z = lower(2); z < std::min(upper(2), dimensions(2)); ++z)
    {
        for (int y = lower(1); y < std::min(upper(1), dimensions(1)); ++y)
        {
            for (int x = lower(0); x < std::min(upper(0), dimensions(0)); ++x)
            {
                const std::size_t index = source_->cellIndex(x, y, z);
                if (
                    !source_->material(index) ||
                    (carveableOnly && !source_->carveable(index))
                )
                {
                    continue;
                }
                const Eigen::Vector3d minimum =
                    source_->grid().minimum +
                    source_->spacing().cwiseProduct(
                        Eigen::Vector3d(x, y, z)
                    );
                result += cuboidMoments(
                    minimum,
                    minimum + source_->spacing()
                );
            }
        }
    }
    return result;
}


void AdaptiveVoxelDomain::sortLeaves()
{
    std::sort(leaves_.begin(), leaves_.end(), keyLess);
}


void AdaptiveVoxelDomain::splitLeaf(std::size_t leafIndex)
{
    if (leafIndex >= leaves_.size())
    {
        throw std::out_of_range("adaptive leaf index is outside the partition");
    }
    const AdaptiveVoxelCell parent = leaves_[leafIndex];
    if (parent.level >= maximumLevel_)
    {
        throw std::invalid_argument("maximum-resolution leaf cannot be split");
    }
    leaves_.erase(leaves_.begin() + static_cast<std::ptrdiff_t>(leafIndex));
    for (int z = 0; z < 2; ++z)
    {
        for (int y = 0; y < 2; ++y)
        {
            for (int x = 0; x < 2; ++x)
            {
                AdaptiveVoxelCell child;
                child.level = parent.level + 1;
                child.coordinate = 2 * parent.coordinate +
                    Eigen::Vector3i(x, y, z);
                child.moments = momentsForCell(
                    child.level,
                    child.coordinate,
                    true
                );
                if (child.moments.mass() > 0.0)
                {
                    leaves_.push_back(child);
                }
            }
        }
    }
    sortLeaves();
}


std::size_t AdaptiveVoxelDomain::splitFractionalLeaves(
    Eigen::VectorXd& fill,
    double epsilon
)
{
    if (
        fill.size() != static_cast<Eigen::Index>(leaves_.size()) ||
        !fill.allFinite() || !std::isfinite(epsilon) || epsilon < 0.0 ||
        epsilon >= 0.5
    )
    {
        throw std::invalid_argument("invalid adaptive split state");
    }
    std::vector<std::pair<AdaptiveVoxelCell, double>> refined;
    std::size_t splitCount = 0U;
    for (std::size_t index = 0; index < leaves_.size(); ++index)
    {
        const AdaptiveVoxelCell& parent = leaves_[index];
        const double parentFill = fill(static_cast<Eigen::Index>(index));
        const bool fractional =
            parentFill >= epsilon && parentFill <= 1.0 - epsilon;
        if (!fractional || parent.level >= maximumLevel_)
        {
            refined.emplace_back(parent, parentFill);
            continue;
        }
        ++splitCount;
        for (int z = 0; z < 2; ++z)
        {
            for (int y = 0; y < 2; ++y)
            {
                for (int x = 0; x < 2; ++x)
                {
                    AdaptiveVoxelCell child;
                    child.level = parent.level + 1;
                    child.coordinate = 2 * parent.coordinate +
                        Eigen::Vector3i(x, y, z);
                    child.moments = momentsForCell(
                        child.level,
                        child.coordinate,
                        true
                    );
                    if (child.moments.mass() > 0.0)
                    {
                        refined.emplace_back(child, parentFill);
                    }
                }
            }
        }
    }
    if (splitCount == 0U)
    {
        return 0U;
    }
    std::sort(refined.begin(), refined.end(), [](const auto& left, const auto& right)
        {
            return keyLess(left.first, right.first);
        });
    leaves_.clear();
    fill.resize(static_cast<Eigen::Index>(refined.size()));
    for (std::size_t index = 0; index < refined.size(); ++index)
    {
        leaves_.push_back(refined[index].first);
        fill(static_cast<Eigen::Index>(index)) = refined[index].second;
    }
    return splitCount;
}


std::size_t AdaptiveVoxelDomain::mergeBinarySiblings(
    Eigen::VectorXd& fill,
    double epsilon
)
{
    if (
        fill.size() != static_cast<Eigen::Index>(leaves_.size()) ||
        !fill.allFinite() || !std::isfinite(epsilon) || epsilon < 0.0 ||
        epsilon >= 0.5
    )
    {
        throw std::invalid_argument("invalid adaptive merge state");
    }
    using ParentKey = std::tuple<int, int, int, int>;
    std::map<ParentKey, std::vector<std::size_t>> groups;
    for (std::size_t index = 0; index < leaves_.size(); ++index)
    {
        const AdaptiveVoxelCell& child = leaves_[index];
        if (child.level <= initialLevel_)
        {
            continue;
        }
        const Eigen::Vector3i parent = child.coordinate / 2;
        groups[{child.level - 1, parent(0), parent(1), parent(2)}].push_back(index);
    }

    std::vector<unsigned char> remove(leaves_.size(), 0U);
    std::vector<std::pair<AdaptiveVoxelCell, double>> additions;
    for (const auto& entry : groups)
    {
        const std::vector<std::size_t>& indices = entry.second;
        bool allZero = true;
        bool allOne = true;
        for (std::size_t index : indices)
        {
            allZero = allZero && fill(static_cast<Eigen::Index>(index)) <= epsilon;
            allOne = allOne && fill(static_cast<Eigen::Index>(index)) >= 1.0 - epsilon;
        }
        if (!allZero && !allOne)
        {
            continue;
        }
        AdaptiveVoxelCell parent;
        parent.level = std::get<0>(entry.first);
        parent.coordinate = Eigen::Vector3i(
            std::get<1>(entry.first),
            std::get<2>(entry.first),
            std::get<3>(entry.first)
        );
        parent.moments = momentsForCell(parent.level, parent.coordinate, true);
        for (std::size_t index : indices)
        {
            remove[index] = 1U;
        }
        additions.emplace_back(parent, allOne ? 1.0 : 0.0);
    }
    if (additions.empty())
    {
        return 0U;
    }

    std::vector<std::pair<AdaptiveVoxelCell, double>> combined;
    combined.reserve(leaves_.size() + additions.size());
    for (std::size_t index = 0; index < leaves_.size(); ++index)
    {
        if (remove[index] == 0U)
        {
            combined.emplace_back(
                leaves_[index],
                fill(static_cast<Eigen::Index>(index))
            );
        }
    }
    combined.insert(combined.end(), additions.begin(), additions.end());
    std::sort(combined.begin(), combined.end(), [](const auto& left, const auto& right)
        {
            return keyLess(left.first, right.first);
        });
    leaves_.clear();
    fill.resize(static_cast<Eigen::Index>(combined.size()));
    for (std::size_t index = 0; index < combined.size(); ++index)
    {
        leaves_.push_back(combined[index].first);
        fill(static_cast<Eigen::Index>(index)) = combined[index].second;
    }
    return additions.size();
}


std::vector<unsigned char> AdaptiveVoxelDomain::binaryVoidMask(
    const Eigen::VectorXd& fill,
    double threshold
) const
{
    if (
        fill.size() != static_cast<Eigen::Index>(leaves_.size()) ||
        !fill.allFinite() || !std::isfinite(threshold) ||
        threshold < 0.0 || threshold > 1.0
    )
    {
        throw std::invalid_argument("invalid binary void-field request");
    }
    std::vector<unsigned char> result(source_->cellCount(), 0U);
    const Eigen::Vector3i dimensions = source_->grid().dimensions;
    for (std::size_t leafIndex = 0; leafIndex < leaves_.size(); ++leafIndex)
    {
        if (fill(static_cast<Eigen::Index>(leafIndex)) < threshold)
        {
            continue;
        }
        const Bounds cellBounds = bounds(leaves_[leafIndex], maximumLevel_);
        for (int z = cellBounds.minimum(2);
             z < std::min(cellBounds.maximum(2), dimensions(2)); ++z)
        {
            for (int y = cellBounds.minimum(1);
                 y < std::min(cellBounds.maximum(1), dimensions(1)); ++y)
            {
                for (int x = cellBounds.minimum(0);
                     x < std::min(cellBounds.maximum(0), dimensions(0)); ++x)
                {
                    const std::size_t index = source_->cellIndex(x, y, z);
                    if (source_->carveable(index))
                    {
                        result[index] = 1U;
                    }
                }
            }
        }
    }
    return result;
}


Eigen::MatrixXd AdaptiveVoxelDomain::uniformLaplacian() const
{
    const Eigen::Index count = static_cast<Eigen::Index>(leaves_.size());
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(count, count);
    for (Eigen::Index left = 0; left < count; ++left)
    {
        for (Eigen::Index right = left + 1; right < count; ++right)
        {
            if (!faceNeighbors(
                    leaves_[static_cast<std::size_t>(left)],
                    leaves_[static_cast<std::size_t>(right)],
                    maximumLevel_
                ))
            {
                continue;
            }
            result(left, right) = -1.0;
            result(right, left) = -1.0;
            result(left, left) += 1.0;
            result(right, right) += 1.0;
        }
    }
    return result;
}


bool AdaptiveVoxelDomain::validatePartition() const
{
    const Eigen::Vector3i dimensions = source_->grid().dimensions;
    std::vector<unsigned char> coverage(source_->cellCount(), 0U);
    for (const AdaptiveVoxelCell& cell : leaves_)
    {
        const Bounds cellBounds = bounds(cell, maximumLevel_);
        for (int z = cellBounds.minimum(2);
             z < std::min(cellBounds.maximum(2), dimensions(2)); ++z)
        {
            for (int y = cellBounds.minimum(1);
                 y < std::min(cellBounds.maximum(1), dimensions(1)); ++y)
            {
                for (int x = cellBounds.minimum(0);
                     x < std::min(cellBounds.maximum(0), dimensions(0)); ++x)
                {
                    const std::size_t index = source_->cellIndex(x, y, z);
                    if (source_->carveable(index))
                    {
                        if (++coverage[index] != 1U)
                        {
                            return false;
                        }
                    }
                }
            }
        }
    }
    for (std::size_t index = 0; index < coverage.size(); ++index)
    {
        if ((coverage[index] == 1U) != source_->carveable(index))
        {
            return false;
        }
    }
    physics::RawMassMoments10 sum;
    for (const AdaptiveVoxelCell& cell : leaves_)
    {
        sum += cell.moments;
    }
    return sum.values().isApprox(variableMoments_.values(), 1e-10);
}


std::string AdaptiveVoxelDomain::partitionHash() const
{
    std::string bytes("spin-it-2014-adaptive-partition-v1");
    bytes += source_->artifactId();
    appendInt(bytes, initialLevel_);
    appendInt(bytes, maximumLevel_);
    for (const AdaptiveVoxelCell& cell : leaves_)
    {
        appendInt(bytes, cell.level);
        for (int axis = 0; axis < 3; ++axis)
        {
            appendInt(bytes, cell.coordinate(axis));
        }
        for (int coordinate = 0;
             coordinate < physics::RawMassMoments10::Count;
             ++coordinate)
        {
            appendDouble(bytes, cell.moments.values()(coordinate));
        }
    }
    return manufacturing::computeContentSha256(bytes);
}


} // namespace spin2014
