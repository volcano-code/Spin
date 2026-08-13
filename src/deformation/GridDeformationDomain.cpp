#include "deformation/GridDeformationDomain.h"

#include "voxel/Voxel.h"
#include "voxel/VoxelGrid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace deformation
{
namespace
{


int checkedProduct(const Eigen::Vector3i& dimensions)
{
    if ((dimensions.array() <= 0).any())
    {
        throw std::invalid_argument(
            "deformation domain requires positive grid dimensions"
        );
    }

    const long long product =
        static_cast<long long>(dimensions(0)) *
        static_cast<long long>(dimensions(1)) *
        static_cast<long long>(dimensions(2));

    if (product > std::numeric_limits<int>::max())
    {
        throw std::overflow_error(
            "deformation grid exceeds supported integer indexing"
        );
    }

    return static_cast<int>(product);
}


int cellIndex(
    int x,
    int y,
    int z,
    const Eigen::Vector3i& dimensions
)
{
    return
        (x * dimensions(1) + y) * dimensions(2) + z;
}


int cornerIndex(
    int x,
    int y,
    int z,
    const Eigen::Vector3i& cornerDimensions
)
{
    // igl::grid(), used by VoxelGrid, increments its first axis first.
    return
        x + cornerDimensions(0) *
        (y + cornerDimensions(1) * z);
}


std::array<int, 8> globalCellCorners(
    int x,
    int y,
    int z,
    const Eigen::Vector3i& cornerDimensions
)
{
    std::array<int, 8> result{};

    for (int bit = 0; bit < 8; ++bit)
    {
        result[static_cast<std::size_t>(bit)] = cornerIndex(
            x + ((bit & 1) != 0 ? 1 : 0),
            y + ((bit & 2) != 0 ? 1 : 0),
            z + ((bit & 4) != 0 ? 1 : 0),
            cornerDimensions
        );
    }

    return result;
}


bool finiteMatrix(const Eigen::MatrixXd& matrix)
{
    for (Eigen::Index row = 0; row < matrix.rows(); ++row)
    {
        for (Eigen::Index column = 0; column < matrix.cols(); ++column)
        {
            if (!std::isfinite(matrix(row, column)))
            {
                return false;
            }
        }
    }

    return true;
}


void validatePointMatrix(
    const Eigen::MatrixXd& points,
    const char* name
)
{
    if (points.cols() != 3 || !finiteMatrix(points))
    {
        throw std::invalid_argument(
            std::string(name) + " must be a finite N-by-3 matrix"
        );
    }
}


void validatePartitionWeights(
    const Eigen::MatrixXd& weights,
    Eigen::Index expectedRows
)
{
    constexpr double tolerance = 1.0e-8;

    if (
        weights.rows() != expectedRows ||
        weights.cols() <= 0 ||
        !finiteMatrix(weights)
    )
    {
        throw std::invalid_argument(
            "node weights must be a finite node-by-handle matrix"
        );
    }

    for (Eigen::Index row = 0; row < weights.rows(); ++row)
    {
        double sum = 0.0;

        for (Eigen::Index column = 0; column < weights.cols(); ++column)
        {
            const double weight = weights(row, column);
            if (weight < -tolerance || weight > 1.0 + tolerance)
            {
                throw std::invalid_argument(
                    "node weights must lie in the closed unit interval"
                );
            }
            sum += weight;
        }

        if (std::abs(sum - 1.0) > tolerance)
        {
            throw std::invalid_argument(
                "every node-weight row must sum to one"
            );
        }
    }
}


std::vector<int> axisCellCandidates(
    double coordinate,
    double minimum,
    double spacing,
    int cellCount,
    double tolerance
)
{
    const double normalized = (coordinate - minimum) / spacing;
    const double normalizedTolerance = tolerance / spacing;

    if (normalized <= normalizedTolerance)
    {
        return {0};
    }

    if (
        normalized >=
        static_cast<double>(cellCount) - normalizedTolerance
    )
    {
        return {cellCount - 1};
    }

    const double nearestPlane = std::round(normalized);
    if (
        nearestPlane > 0.0 &&
        nearestPlane < static_cast<double>(cellCount) &&
        std::abs(normalized - nearestPlane) <= normalizedTolerance
    )
    {
        const int upper = static_cast<int>(nearestPlane);
        return {upper - 1, upper};
    }

    const int index = std::clamp(
        static_cast<int>(std::floor(normalized)),
        0,
        cellCount - 1
    );
    return {index};
}


bool sameCornerSet(
    const std::set<int>& actual,
    const std::array<int, 8>& expected
)
{
    if (actual.size() != expected.size())
    {
        return false;
    }

    for (const int corner : expected)
    {
        if (actual.count(corner) == 0U)
        {
            return false;
        }
    }

    return true;
}


bool isExpectedCorner(
    const Eigen::Vector3d& actual,
    const Eigen::Vector3d& expected,
    const Eigen::Vector3d& voxelSize,
    double scaleTolerance
)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const double axisTolerance = std::max(
            scaleTolerance,
            1024.0 * std::numeric_limits<double>::epsilon() *
                std::max(voxelSize(axis), 1.0)
        );

        if (std::abs(actual(axis) - expected(axis)) > axisTolerance)
        {
            return false;
        }
    }

    return true;
}


} // namespace


GridDeformationDomain GridDeformationDomain::fromVoxelGrid(
    const VoxelGrid& voxelGrid
)
{
    GridDeformationDomain domain;
    domain.dimensions_ = voxelGrid.dimensions();
    domain.voxelSize_ = voxelGrid.voxelSize();
    domain.fullRestCorners_ = voxelGrid.corners();

    const int fullCellCount = checkedProduct(domain.dimensions_);
    const Eigen::Vector3i cornerDimensions =
        domain.dimensions_ + Eigen::Vector3i::Ones();
    const int fullCornerCount = checkedProduct(cornerDimensions);

    if (
        !domain.voxelSize_.allFinite() ||
        (domain.voxelSize_.array() <= 0.0).any()
    )
    {
        throw std::invalid_argument(
            "deformation domain requires finite positive voxel sizes"
        );
    }

    if (
        domain.fullRestCorners_.rows() != fullCornerCount ||
        domain.fullRestCorners_.cols() != 3 ||
        !finiteMatrix(domain.fullRestCorners_)
    )
    {
        throw std::invalid_argument(
            "voxel-grid corners do not match the declared dimensions"
        );
    }

    domain.gridMinimum_ =
        domain.fullRestCorners_.colwise().minCoeff().transpose();
    domain.gridMaximum_ =
        domain.fullRestCorners_.colwise().maxCoeff().transpose();

    const double coordinateScale = std::max(
        {
            1.0,
            domain.fullRestCorners_.cwiseAbs().maxCoeff(),
            domain.voxelSize_.norm(),
            (domain.gridMaximum_ - domain.gridMinimum_).norm()
        }
    );
    domain.geometricTolerance_ =
        512.0 * std::numeric_limits<double>::epsilon() * coordinateScale;

    domain.cellNodes_ = Eigen::MatrixXi::Constant(fullCellCount, 8, -1);
    domain.cellGlobalCorners_.resize(fullCellCount, 8);
    domain.activeCellMask_ = Eigen::VectorXi::Zero(fullCellCount);

    std::vector<unsigned char> activeGlobalCorners(
        static_cast<std::size_t>(fullCornerCount),
        0U
    );

    const double halfDiagonal = 0.5 * domain.voxelSize_.norm();

    for (int x = 0; x < domain.dimensions_(0); ++x)
    {
        for (int y = 0; y < domain.dimensions_(1); ++y)
        {
            for (int z = 0; z < domain.dimensions_(2); ++z)
            {
                const int linear = cellIndex(
                    x,
                    y,
                    z,
                    domain.dimensions_
                );
                const Voxel& voxel = voxelGrid(x, y, z);

                if (
                    !voxel.center.allFinite() ||
                    !std::isfinite(voxel.distanceToMesh) ||
                    voxel.xIdx != x ||
                    voxel.yIdx != y ||
                    voxel.zIdx != z
                )
                {
                    throw std::invalid_argument(
                        "voxel metadata is inconsistent with its grid index"
                    );
                }

                const Eigen::Vector3d expectedCenter =
                    domain.gridMinimum_ +
                    domain.voxelSize_.cwiseProduct(
                        Eigen::Vector3d(
                            static_cast<double>(x) + 0.5,
                            static_cast<double>(y) + 0.5,
                            static_cast<double>(z) + 0.5
                        )
                    );
                if (!isExpectedCorner(
                        voxel.center,
                        expectedCenter,
                        domain.voxelSize_,
                        8.0 * domain.geometricTolerance_
                    ))
                {
                    throw std::invalid_argument(
                        "voxel centers are inconsistent with the corner lattice"
                    );
                }

                const std::array<int, 8> globals = globalCellCorners(
                    x,
                    y,
                    z,
                    cornerDimensions
                );

                if (
                    voxel.cornerIndices.size() != 8U ||
                    std::any_of(
                        voxel.cornerIndices.begin(),
                        voxel.cornerIndices.end(),
                        [fullCornerCount](int corner)
                        {
                            return
                                corner < 0 ||
                                corner >= fullCornerCount;
                        }
                    )
                )
                {
                    throw std::invalid_argument(
                        "voxel corner metadata contains invalid indices"
                    );
                }

                if (!sameCornerSet(voxel.cornerIndices, globals))
                {
                    /*
                     * A non-cubic lattice can expose an indexing convention
                     * mismatch between a coordinate formula and igl::grid.
                     * Resolve it from the stored coordinates, then still
                     * require exact agreement with Voxel::cornerIndices.
                     */
                    std::array<int, 8> resolved{};
                    std::set<int> used;

                    for (int bit = 0; bit < 8; ++bit)
                    {
                        const Eigen::Vector3d expectedCorner =
                            domain.gridMinimum_ +
                            domain.voxelSize_.cwiseProduct(
                                Eigen::Vector3d(
                                    static_cast<double>(
                                        x + ((bit & 1) != 0)
                                    ),
                                    static_cast<double>(
                                        y + ((bit & 2) != 0)
                                    ),
                                    static_cast<double>(
                                        z + ((bit & 4) != 0)
                                    )
                                )
                            );

                        int match = -1;
                        double best =
                            std::numeric_limits<double>::infinity();

                        for (const int candidate : voxel.cornerIndices)
                        {
                            const Eigen::Vector3d actualCorner =
                                domain.fullRestCorners_.row(candidate).
                                    transpose();
                            const double error =
                                (actualCorner - expectedCorner).
                                    cwiseQuotient(domain.voxelSize_).
                                    cwiseAbs().maxCoeff();

                            if (
                                used.count(candidate) == 0U &&
                                error < best
                            )
                            {
                                best = error;
                                match = candidate;
                            }
                        }

                        if (
                            match < 0 ||
                            !isExpectedCorner(
                                domain.fullRestCorners_.row(match).
                                    transpose(),
                                expectedCorner,
                                domain.voxelSize_,
                                8.0 * domain.geometricTolerance_
                            )
                        )
                        {
                            throw std::invalid_argument(
                                "voxel corner metadata is inconsistent with the lattice"
                            );
                        }

                        resolved[static_cast<std::size_t>(bit)] = match;
                        used.insert(match);
                    }

                    if (!sameCornerSet(voxel.cornerIndices, resolved))
                    {
                        throw std::invalid_argument(
                            "voxel corner metadata does not form one regular cell"
                        );
                    }

                    for (int bit = 0; bit < 8; ++bit)
                    {
                        domain.cellGlobalCorners_(linear, bit) =
                            resolved[static_cast<std::size_t>(bit)];
                    }
                }
                else
                {
                    for (int bit = 0; bit < 8; ++bit)
                    {
                        domain.cellGlobalCorners_(linear, bit) =
                            globals[static_cast<std::size_t>(bit)];
                    }
                }

                for (int bit = 0; bit < 8; ++bit)
                {
                    const int global =
                        domain.cellGlobalCorners_(linear, bit);
                    const Eigen::Vector3d expectedCorner =
                        domain.gridMinimum_ +
                        domain.voxelSize_.cwiseProduct(
                            Eigen::Vector3d(
                                static_cast<double>(x + ((bit & 1) != 0)),
                                static_cast<double>(y + ((bit & 2) != 0)),
                                static_cast<double>(z + ((bit & 4) != 0))
                            )
                        );
                    if (!isExpectedCorner(
                            domain.fullRestCorners_.row(global).transpose(),
                            expectedCorner,
                            domain.voxelSize_,
                            8.0 * domain.geometricTolerance_
                        ))
                    {
                        throw std::invalid_argument(
                            "corner coordinates are inconsistent with grid indexing"
                        );
                    }
                }

                if (
                    voxel.distanceToMesh <=
                    halfDiagonal + domain.geometricTolerance_
                )
                {
                    domain.activeCellMask_(linear) = 1;
                    ++domain.activeCellCount_;

                    for (int bit = 0; bit < 8; ++bit)
                    {
                        const int global =
                            domain.cellGlobalCorners_(linear, bit);
                        activeGlobalCorners[
                            static_cast<std::size_t>(global)
                        ] = 1U;
                    }
                }
            }
        }
    }

    if (domain.activeCellCount_ == 0)
    {
        throw std::invalid_argument(
            "deformation domain contains no object-intersecting cells"
        );
    }

    int activeNodeCount = 0;
    for (const unsigned char active : activeGlobalCorners)
    {
        activeNodeCount += active != 0U ? 1 : 0;
    }

    domain.restNodes_.resize(activeNodeCount, 3);
    domain.localToGlobalCorner_.resize(activeNodeCount);
    domain.globalToLocalCorner_ =
        Eigen::VectorXi::Constant(fullCornerCount, -1);

    int nextLocal = 0;
    for (int global = 0; global < fullCornerCount; ++global)
    {
        if (activeGlobalCorners[static_cast<std::size_t>(global)] == 0U)
        {
            continue;
        }

        domain.globalToLocalCorner_(global) = nextLocal;
        domain.localToGlobalCorner_(nextLocal) = global;
        domain.restNodes_.row(nextLocal) =
            domain.fullRestCorners_.row(global);
        ++nextLocal;
    }

    std::set<std::pair<int, int>> uniqueEdges;
    constexpr std::array<std::pair<int, int>, 12> cubeEdges{{
        {0, 1}, {2, 3}, {4, 5}, {6, 7},
        {0, 2}, {1, 3}, {4, 6}, {5, 7},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    }};

    for (int cell = 0; cell < fullCellCount; ++cell)
    {
        if (domain.activeCellMask_(cell) == 0)
        {
            continue;
        }

        for (int bit = 0; bit < 8; ++bit)
        {
            const int global = domain.cellGlobalCorners_(cell, bit);
            const int local = domain.globalToLocalCorner_(global);
            if (local < 0)
            {
                throw std::logic_error(
                    "active cell references a missing deformation node"
                );
            }
            domain.cellNodes_(cell, bit) = local;
        }

        for (const auto& edge : cubeEdges)
        {
            int first = domain.cellNodes_(cell, edge.first);
            int second = domain.cellNodes_(cell, edge.second);
            if (first > second)
            {
                std::swap(first, second);
            }
            uniqueEdges.emplace(first, second);
        }
    }

    domain.edges_.resize(
        static_cast<Eigen::Index>(uniqueEdges.size()),
        2
    );
    Eigen::Index edgeRow = 0;
    for (const auto& edge : uniqueEdges)
    {
        domain.edges_(edgeRow, 0) = edge.first;
        domain.edges_(edgeRow, 1) = edge.second;
        ++edgeRow;
    }

    std::vector<std::vector<int>> adjacency(
        static_cast<std::size_t>(activeNodeCount)
    );
    for (Eigen::Index row = 0; row < domain.edges_.rows(); ++row)
    {
        const int first = domain.edges_(row, 0);
        const int second = domain.edges_(row, 1);
        adjacency[static_cast<std::size_t>(first)].push_back(second);
        adjacency[static_cast<std::size_t>(second)].push_back(first);
    }

    domain.nodeComponents_ = Eigen::VectorXi::Constant(activeNodeCount, -1);
    for (int start = 0; start < activeNodeCount; ++start)
    {
        if (domain.nodeComponents_(start) >= 0)
        {
            continue;
        }

        const int component = domain.componentCount_++;
        std::queue<int> pending;
        pending.push(start);
        domain.nodeComponents_(start) = component;

        while (!pending.empty())
        {
            const int node = pending.front();
            pending.pop();

            for (const int neighbor : adjacency[static_cast<std::size_t>(node)])
            {
                if (domain.nodeComponents_(neighbor) < 0)
                {
                    domain.nodeComponents_(neighbor) = component;
                    pending.push(neighbor);
                }
            }
        }
    }

    if (
        domain.componentCount_ <= 0 ||
        (domain.nodeComponents_.array() < 0).any()
    )
    {
        throw std::logic_error(
            "deformation-domain component metadata is incomplete"
        );
    }

    return domain;
}


Eigen::Index GridDeformationDomain::nodeCount() const noexcept
{
    return restNodes_.rows();
}


Eigen::Index GridDeformationDomain::edgeCount() const noexcept
{
    return edges_.rows();
}


Eigen::Index GridDeformationDomain::cellCount() const noexcept
{
    return cellNodes_.rows();
}


Eigen::Index GridDeformationDomain::activeCellCount() const noexcept
{
    return activeCellCount_;
}


const Eigen::MatrixXd& GridDeformationDomain::restNodes() const noexcept
{
    return restNodes_;
}


const Eigen::MatrixXi& GridDeformationDomain::edges() const noexcept
{
    return edges_;
}


const Eigen::VectorXi&
GridDeformationDomain::localToGlobalCorner() const noexcept
{
    return localToGlobalCorner_;
}


const Eigen::VectorXi& GridDeformationDomain::localToGlobal() const noexcept
{
    return localToGlobalCorner_;
}


const Eigen::VectorXi&
GridDeformationDomain::globalToLocalCorners() const noexcept
{
    return globalToLocalCorner_;
}


const Eigen::VectorXi& GridDeformationDomain::globalToLocal() const noexcept
{
    return globalToLocalCorner_;
}


const Eigen::VectorXi&
GridDeformationDomain::globalToLocalCorner() const noexcept
{
    return globalToLocalCorner_;
}


const Eigen::MatrixXi& GridDeformationDomain::cellNodes() const noexcept
{
    return cellNodes_;
}


const Eigen::VectorXi&
GridDeformationDomain::activeCellMask() const noexcept
{
    return activeCellMask_;
}


const Eigen::Vector3i& GridDeformationDomain::dimensions() const noexcept
{
    return dimensions_;
}


const Eigen::Vector3d& GridDeformationDomain::voxelSize() const noexcept
{
    return voxelSize_;
}


const Eigen::MatrixXd&
GridDeformationDomain::fullRestCorners() const noexcept
{
    return fullRestCorners_;
}


int GridDeformationDomain::componentCount() const noexcept
{
    return componentCount_;
}


const Eigen::VectorXi&
GridDeformationDomain::nodeComponents() const noexcept
{
    return nodeComponents_;
}


bool GridDeformationDomain::connected() const noexcept
{
    return componentCount_ == 1;
}


int GridDeformationDomain::globalToLocalCorner(int globalCorner) const
{
    if (
        globalCorner < 0 ||
        globalCorner >= globalToLocalCorner_.size()
    )
    {
        throw std::out_of_range("global corner index is outside the grid");
    }

    return globalToLocalCorner_(globalCorner);
}


int GridDeformationDomain::localToGlobalCorner(int localNode) const
{
    if (localNode < 0 || localNode >= localToGlobalCorner_.size())
    {
        throw std::out_of_range("local node index is outside the domain");
    }

    return localToGlobalCorner_(localNode);
}


int GridDeformationDomain::localNodeForGlobalCorner(
    int globalCorner
) const
{
    return globalToLocalCorner(globalCorner);
}


int GridDeformationDomain::globalCornerForLocalNode(int localNode) const
{
    return localToGlobalCorner(localNode);
}


int GridDeformationDomain::cellLinearIndex(int x, int y, int z) const
{
    if (
        x < 0 || x >= dimensions_(0) ||
        y < 0 || y >= dimensions_(1) ||
        z < 0 || z >= dimensions_(2)
    )
    {
        throw std::out_of_range("cell coordinate is outside the grid");
    }

    return cellIndex(x, y, z, dimensions_);
}


Eigen::Vector3i GridDeformationDomain::cellCoordinates(
    int linearIndex
) const
{
    if (linearIndex < 0 || linearIndex >= cellNodes_.rows())
    {
        throw std::out_of_range("linear cell index is outside the grid");
    }

    const int yz = dimensions_(1) * dimensions_(2);
    const int x = linearIndex / yz;
    const int remainder = linearIndex % yz;
    const int y = remainder / dimensions_(2);
    const int z = remainder % dimensions_(2);
    return Eigen::Vector3i(x, y, z);
}


bool GridDeformationDomain::isActiveCell(int linearIndex) const
{
    if (linearIndex < 0 || linearIndex >= activeCellMask_.size())
    {
        throw std::out_of_range("linear cell index is outside the grid");
    }

    return activeCellMask_(linearIndex) != 0;
}


std::array<int, 8> GridDeformationDomain::cellNodes(
    int linearIndex
) const
{
    if (linearIndex < 0 || linearIndex >= cellNodes_.rows())
    {
        throw std::out_of_range("linear cell index is outside the grid");
    }

    std::array<int, 8> result{};
    for (int bit = 0; bit < 8; ++bit)
    {
        result[static_cast<std::size_t>(bit)] =
            cellNodes_(linearIndex, bit);
    }
    return result;
}


std::array<int, 8> GridDeformationDomain::cellGlobalCorners(
    int linearIndex
) const
{
    if (linearIndex < 0 || linearIndex >= cellGlobalCorners_.rows())
    {
        throw std::out_of_range("linear cell index is outside the grid");
    }

    std::array<int, 8> result{};
    for (int bit = 0; bit < 8; ++bit)
    {
        result[static_cast<std::size_t>(bit)] =
            cellGlobalCorners_(linearIndex, bit);
    }
    return result;
}


TrilinearEmbedding GridDeformationDomain::embed(
    const Eigen::Vector3d& point
) const
{
    if (!point.allFinite())
    {
        throw std::invalid_argument("embedding point must be finite");
    }

    Eigen::Vector3d clamped = point;
    for (int axis = 0; axis < 3; ++axis)
    {
        if (
            point(axis) < gridMinimum_(axis) - geometricTolerance_ ||
            point(axis) > gridMaximum_(axis) + geometricTolerance_
        )
        {
            throw std::out_of_range(
                "embedding point lies outside the voxel lattice"
            );
        }

        clamped(axis) = std::clamp(
            point(axis),
            gridMinimum_(axis),
            gridMaximum_(axis)
        );
    }

    std::array<std::vector<int>, 3> candidates;
    for (int axis = 0; axis < 3; ++axis)
    {
        candidates[static_cast<std::size_t>(axis)] = axisCellCandidates(
            clamped(axis),
            gridMinimum_(axis),
            voxelSize_(axis),
            dimensions_(axis),
            geometricTolerance_
        );
    }

    int selectedCell = -1;
    Eigen::Vector3i selectedCoordinates = Eigen::Vector3i::Zero();

    for (const int x : candidates[0])
    {
        for (const int y : candidates[1])
        {
            for (const int z : candidates[2])
            {
                const int linear = cellIndex(x, y, z, dimensions_);
                if (
                    activeCellMask_(linear) != 0 &&
                    (selectedCell < 0 || linear < selectedCell)
                )
                {
                    selectedCell = linear;
                    selectedCoordinates = Eigen::Vector3i(x, y, z);
                }
            }
        }
    }

    if (selectedCell < 0)
    {
        throw std::out_of_range(
            "embedding point is not contained in an active deformation cell"
        );
    }

    Eigen::Vector3d local =
        (clamped - gridMinimum_).cwiseQuotient(voxelSize_) -
        selectedCoordinates.cast<double>();
    local = local.cwiseMax(0.0).cwiseMin(1.0);

    TrilinearEmbedding result;
    result.cellLinearIndex = selectedCell;

    for (int bit = 0; bit < 8; ++bit)
    {
        result.nodes[static_cast<std::size_t>(bit)] =
            cellNodes_(selectedCell, bit);

        double coefficient = 1.0;
        for (int axis = 0; axis < 3; ++axis)
        {
            const bool upper = (bit & (1 << axis)) != 0;
            coefficient *= upper
                ? local(axis)
                : 1.0 - local(axis);
        }
        result.coeffs[static_cast<std::size_t>(bit)] = coefficient;
    }

    return result;
}


int GridDeformationDomain::nearestActiveNode(
    const Eigen::Vector3d& point
) const
{
    if (!point.allFinite())
    {
        throw std::invalid_argument("nearest-node query must be finite");
    }

    if (restNodes_.rows() == 0)
    {
        throw std::logic_error("deformation domain has no active nodes");
    }

    int bestNode = 0;
    double bestSquaredDistance =
        (restNodes_.row(0).transpose() - point).squaredNorm();

    for (Eigen::Index node = 1; node < restNodes_.rows(); ++node)
    {
        const double squaredDistance =
            (restNodes_.row(node).transpose() - point).squaredNorm();
        if (squaredDistance < bestSquaredDistance)
        {
            bestNode = static_cast<int>(node);
            bestSquaredDistance = squaredDistance;
        }
    }

    return bestNode;
}


Eigen::MatrixXd GridDeformationDomain::interpolateWeights(
    const Eigen::MatrixXd& points,
    const Eigen::MatrixXd& nodeWeights
) const
{
    validatePointMatrix(points, "interpolation points");
    validatePartitionWeights(nodeWeights, restNodes_.rows());

    Eigen::MatrixXd result(points.rows(), nodeWeights.cols());
    result.setZero();

    for (Eigen::Index row = 0; row < points.rows(); ++row)
    {
        const TrilinearEmbedding embedding =
            embed(points.row(row).transpose());

        for (int bit = 0; bit < 8; ++bit)
        {
            result.row(row) +=
                embedding.coeffs[static_cast<std::size_t>(bit)] *
                nodeWeights.row(
                    embedding.nodes[static_cast<std::size_t>(bit)]
                );
        }

        // Each interpolated row is a convex combination mathematically, but
        // BBW values and trilinear coefficients can leave roundoff-sized
        // excursions outside [0,1].  Project only values already accepted by
        // validatePartitionWeights() and renormalize before the strict LBS
        // contract sees them.  Material excursions still fail closed above.
        for (Eigen::Index column = 0; column < result.cols(); ++column)
        {
            result(row, column) = std::clamp(
                result(row, column),
                0.0,
                1.0
            );
        }
        const double rowSum = result.row(row).sum();
        if (!(rowSum > 0.0) || !std::isfinite(rowSum))
        {
            throw std::logic_error(
                "interpolated weight row has no positive finite mass"
            );
        }
        result.row(row) /= rowSum;
    }

    validatePartitionWeights(result, points.rows());
    return result;
}


Eigen::MatrixXd GridDeformationDomain::deformedCellCenters(
    const Eigen::MatrixXd& fullDeformedCorners
) const
{
    validatePointMatrix(fullDeformedCorners, "deformed grid corners");

    if (fullDeformedCorners.rows() != fullRestCorners_.rows())
    {
        throw std::invalid_argument(
            "deformed corner array does not match the full grid"
        );
    }

    Eigen::MatrixXd centers(cellGlobalCorners_.rows(), 3);
    centers.setZero();

    for (Eigen::Index cell = 0; cell < cellGlobalCorners_.rows(); ++cell)
    {
        for (int bit = 0; bit < 8; ++bit)
        {
            centers.row(cell) += fullDeformedCorners.row(
                cellGlobalCorners_(cell, bit)
            );
        }
        centers.row(cell) /= 8.0;
    }

    return centers;
}


} // namespace deformation
