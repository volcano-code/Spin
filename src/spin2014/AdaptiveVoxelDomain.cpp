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
#include <unordered_map>
#include <unordered_set>
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


struct FaceTile
{
    int axis = 0;
    int plane = 0;
    int first = 0;
    int second = 0;

    bool operator==(const FaceTile& other) const noexcept
    {
        return axis == other.axis && plane == other.plane &&
            first == other.first && second == other.second;
    }
};


struct FaceTileHash
{
    std::size_t operator()(const FaceTile& value) const noexcept
    {
        std::size_t seed = 0xcbf29ce484222325ULL;
        for (const int coordinate : {
                 value.axis, value.plane, value.first, value.second
             })
        {
            seed ^= std::hash<int>{}(coordinate) + 0x9e3779b9U +
                (seed << 6U) + (seed >> 2U);
        }
        return seed;
    }
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
    const int cellsPerAxis = powerOfTwo(level);
    const Eigen::Vector3d extent =
        source_->grid().maximum - source_->grid().minimum;
    const Eigen::Vector3d cellMinimum = source_->grid().minimum +
        extent.cwiseProduct(coordinate.cast<double>()) /
            static_cast<double>(cellsPerAxis);
    const Eigen::Vector3d cellMaximum = source_->grid().minimum +
        extent.cwiseProduct((coordinate.array() + 1).matrix().cast<double>()) /
            static_cast<double>(cellsPerAxis);
    const Eigen::Vector3i dimensions = source_->grid().dimensions;
    Eigen::Vector3i lower;
    Eigen::Vector3i upper;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double lo = (cellMinimum(axis) - source_->grid().minimum(axis)) /
            source_->spacing()(axis);
        const double hi = (cellMaximum(axis) - source_->grid().minimum(axis)) /
            source_->spacing()(axis);
        lower(axis) = std::clamp(
            static_cast<int>(std::floor(lo + 1e-12)), 0, dimensions(axis)
        );
        upper(axis) = std::clamp(
            static_cast<int>(std::ceil(hi - 1e-12)), 0, dimensions(axis)
        );
    }
    physics::RawMassMoments10 result;
    for (int z = lower(2); z < upper(2); ++z)
    {
        for (int y = lower(1); y < upper(1); ++y)
        {
            for (int x = lower(0); x < upper(0); ++x)
            {
                const std::size_t index = source_->cellIndex(x, y, z);
                if (
                    !source_->material(index) ||
                    (carveableOnly && !source_->carveable(index))
                )
                {
                    continue;
                }
                const Eigen::Vector3d sourceMinimum =
                    source_->grid().minimum +
                    source_->spacing().cwiseProduct(
                        Eigen::Vector3d(x, y, z)
                    );
                const Eigen::Vector3d overlapMinimum =
                    sourceMinimum.cwiseMax(cellMinimum);
                const Eigen::Vector3d overlapMaximum =
                    (sourceMinimum + source_->spacing()).cwiseMin(cellMaximum);
                if ((overlapMaximum.array() > overlapMinimum.array()).all())
                {
                    result += cuboidMoments(overlapMinimum, overlapMaximum);
                }
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
    const double exactBoundTolerance =
        64.0 * std::numeric_limits<double>::epsilon();
    for (const auto& entry : groups)
    {
        const std::vector<std::size_t>& indices = entry.second;
        bool allZero = true;
        bool allOne = true;
        for (std::size_t index : indices)
        {
            allZero = allZero &&
                fill(static_cast<Eigen::Index>(index)) <= exactBoundTolerance;
            allOne = allOne &&
                fill(static_cast<Eigen::Index>(index)) >=
                    1.0 - exactBoundTolerance;
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
        physics::RawMassMoments10 childMoments;
        for (std::size_t index : indices)
        {
            childMoments += leaves_[index].moments;
        }
        if (!childMoments.values().isApprox(parent.moments.values(), 1e-12))
        {
            // Some of the parent's region is represented by a coarser leaf,
            // so replacing this partial sibling set would overlap material.
            continue;
        }
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
    std::vector<double> voidVolume(source_->cellCount(), 0.0);
    const Eigen::Vector3i dimensions = source_->grid().dimensions;
    for (std::size_t leafIndex = 0; leafIndex < leaves_.size(); ++leafIndex)
    {
        if (fill(static_cast<Eigen::Index>(leafIndex)) < threshold)
        {
            continue;
        }
        const AdaptiveVoxelCell& leaf = leaves_[leafIndex];
        const int cellsPerAxis = powerOfTwo(leaf.level);
        const Eigen::Vector3d extent =
            source_->grid().maximum - source_->grid().minimum;
        const Eigen::Vector3d leafMinimum = source_->grid().minimum +
            extent.cwiseProduct(leaf.coordinate.cast<double>()) /
                static_cast<double>(cellsPerAxis);
        const Eigen::Vector3d leafMaximum = source_->grid().minimum +
            extent.cwiseProduct(
                (leaf.coordinate.array() + 1).matrix().cast<double>()
            ) / static_cast<double>(cellsPerAxis);
        Eigen::Vector3i lower;
        Eigen::Vector3i upper;
        for (int axis = 0; axis < 3; ++axis)
        {
            lower(axis) = std::clamp(static_cast<int>(std::floor(
                (leafMinimum(axis) - source_->grid().minimum(axis)) /
                source_->spacing()(axis) + 1e-12)), 0, dimensions(axis));
            upper(axis) = std::clamp(static_cast<int>(std::ceil(
                (leafMaximum(axis) - source_->grid().minimum(axis)) /
                source_->spacing()(axis) - 1e-12)), 0, dimensions(axis));
        }
        for (int z = lower(2); z < upper(2); ++z)
        {
            for (int y = lower(1); y < upper(1); ++y)
            {
                for (int x = lower(0); x < upper(0); ++x)
                {
                    const std::size_t index = source_->cellIndex(x, y, z);
                    if (source_->carveable(index))
                    {
                        const Eigen::Vector3d sourceMinimum =
                            source_->grid().minimum +
                            source_->spacing().cwiseProduct(Eigen::Vector3d(x,y,z));
                        const Eigen::Vector3d overlapMinimum =
                            sourceMinimum.cwiseMax(leafMinimum);
                        const Eigen::Vector3d overlapMaximum =
                            (sourceMinimum + source_->spacing()).cwiseMin(leafMaximum);
                        if ((overlapMaximum.array() > overlapMinimum.array()).all())
                        {
                            voidVolume[index] +=
                                (overlapMaximum - overlapMinimum).prod();
                        }
                    }
                }
            }
        }
    }
    const double sourceCellVolume = source_->spacing().prod();
    for (std::size_t index = 0; index < result.size(); ++index)
    {
        result[index] = source_->carveable(index) &&
            voidVolume[index] >= sourceCellVolume * (1.0 - 1e-10) ? 1U : 0U;
    }
    return result;
}


std::shared_ptr<const voxel::VoxelDomain>
AdaptiveVoxelDomain::maximumResolutionMaterialDomain(
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
        throw std::invalid_argument("invalid maximum-resolution material request");
    }
    const int side = powerOfTwo(maximumLevel_);
    const Eigen::Vector3i sourceDimensions = source_->grid().dimensions;
    bool divisible = true;
    for (int axis = 0; axis < 3; ++axis)
    {
        divisible = divisible && side % sourceDimensions(axis) == 0;
    }
    if (!divisible)
    {
        throw std::invalid_argument(
            "maximum-resolution grid must be divisible by source dimensions"
        );
    }
    const std::size_t cellCount = static_cast<std::size_t>(side) *
        static_cast<std::size_t>(side) * static_cast<std::size_t>(side);
    std::vector<unsigned char> material(cellCount, 0U);
    for (int z = 0; z < side; ++z)
    {
        const int sourceZ = z * sourceDimensions(2) / side;
        for (int y = 0; y < side; ++y)
        {
            const int sourceY = y * sourceDimensions(1) / side;
            for (int x = 0; x < side; ++x)
            {
                const int sourceX = x * sourceDimensions(0) / side;
                material[static_cast<std::size_t>(x + side * (y + side * z))] =
                    source_->material(source_->cellIndex(sourceX, sourceY, sourceZ))
                    ? 1U : 0U;
            }
        }
    }
    for (std::size_t leafIndex = 0; leafIndex < leaves_.size(); ++leafIndex)
    {
        if (fill(static_cast<Eigen::Index>(leafIndex)) < threshold) continue;
        const Bounds leafBounds = bounds(leaves_[leafIndex], maximumLevel_);
        for (int z = leafBounds.minimum(2); z < leafBounds.maximum(2); ++z)
            for (int y = leafBounds.minimum(1); y < leafBounds.maximum(1); ++y)
                for (int x = leafBounds.minimum(0); x < leafBounds.maximum(0); ++x)
                {
                    const int sourceX = x * sourceDimensions(0) / side;
                    const int sourceY = y * sourceDimensions(1) / side;
                    const int sourceZ = z * sourceDimensions(2) / side;
                    if (source_->carveable(source_->cellIndex(
                            sourceX, sourceY, sourceZ
                        )))
                    {
                        material[static_cast<std::size_t>(
                            x + side * (y + side * z)
                        )] = 0U;
                    }
                }
    }
    voxel::GridSpec grid = source_->grid();
    grid.dimensions = Eigen::Vector3i::Constant(side);
    return voxel::VoxelDomain::create(grid, material, material, material);
}


Eigen::SparseMatrix<double> AdaptiveVoxelDomain::uniformLaplacian() const
{
    const Eigen::Index count = static_cast<Eigen::Index>(leaves_.size());
    std::unordered_map<FaceTile, int, FaceTileHash> faceOwner;
    std::unordered_set<std::uint64_t> neighborPairs;
    const auto pairKey = [](int first, int second)
    {
        if (first > second) std::swap(first, second);
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(first)) << 32U) |
            static_cast<std::uint32_t>(second);
    };
    for (int index = 0; index < static_cast<int>(leaves_.size()); ++index)
    {
        const Bounds cell = bounds(leaves_[static_cast<std::size_t>(index)], maximumLevel_);
        for (int axis = 0; axis < 3; ++axis)
        {
            const int firstAxis = (axis + 1) % 3;
            const int secondAxis = (axis + 2) % 3;
            for (const int plane : {cell.minimum(axis), cell.maximum(axis)})
            {
                for (int first = cell.minimum(firstAxis);
                     first < cell.maximum(firstAxis); ++first)
                {
                    for (int second = cell.minimum(secondAxis);
                         second < cell.maximum(secondAxis); ++second)
                    {
                        const FaceTile key{axis, plane, first, second};
                        const auto inserted = faceOwner.emplace(key, index);
                        if (!inserted.second && inserted.first->second != index)
                        {
                            neighborPairs.insert(pairKey(inserted.first->second, index));
                        }
                    }
                }
            }
        }
    }
    std::vector<double> degree(static_cast<std::size_t>(count), 0.0);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(2U * neighborPairs.size() + static_cast<std::size_t>(count));
    for (const std::uint64_t encoded : neighborPairs)
    {
        const int first = static_cast<int>(encoded >> 32U);
        const int second = static_cast<int>(encoded & 0xffffffffU);
        triplets.emplace_back(first, second, -1.0);
        triplets.emplace_back(second, first, -1.0);
        degree[static_cast<std::size_t>(first)] += 1.0;
        degree[static_cast<std::size_t>(second)] += 1.0;
    }
    for (Eigen::Index index = 0; index < count; ++index)
    {
        if (degree[static_cast<std::size_t>(index)] != 0.0)
            triplets.emplace_back(index, index, degree[static_cast<std::size_t>(index)]);
    }
    Eigen::SparseMatrix<double> result(count, count);
    result.setFromTriplets(triplets.begin(), triplets.end());
    result.makeCompressed();
    return result;
}


bool AdaptiveVoxelDomain::validatePartition() const
{
    for (const AdaptiveVoxelCell& cell : leaves_)
    {
        const int side = powerOfTwo(cell.level);
        if ((cell.coordinate.array() < 0).any() ||
            (cell.coordinate.array() >= side).any())
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
