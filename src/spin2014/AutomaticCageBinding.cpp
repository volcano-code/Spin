#include "spin2014/AutomaticCageBinding.h"

#include "deformation/GridDeformationDomain.h"
#include "manufacturing/SlicerValidation.h"
#include "voxel/VoxelGrid.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


void appendInt(std::string& bytes, int value)
{
    bytes.append(reinterpret_cast<const char*>(&value), sizeof(value));
}


void appendDouble(std::string& bytes, double value)
{
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    for (int shift = 0; shift < 64; shift += 8)
    {
        bytes.push_back(static_cast<char>((bits >> shift) & 0xffU));
    }
}


void validatePoints(const Eigen::MatrixXd& points, const char* name)
{
    if (points.rows() == 0 || points.cols() != 3 || !points.allFinite())
    {
        throw std::invalid_argument(
            std::string(name) + " must be a nonempty finite N-by-3 matrix"
        );
    }
}


std::vector<int> boundaryNodes(
    const deformation::GridDeformationDomain& grid
)
{
    std::vector<int> incidence(
        static_cast<std::size_t>(grid.nodeCount()), 0
    );
    for (Eigen::Index cell = 0; cell < grid.cellCount(); ++cell)
    {
        if (!grid.isActiveCell(static_cast<int>(cell))) continue;
        for (const int node : grid.cellNodes(static_cast<int>(cell)))
        {
            ++incidence[static_cast<std::size_t>(node)];
        }
    }
    std::vector<int> result;
    for (int node = 0; node < static_cast<int>(incidence.size()); ++node)
    {
        if (incidence[static_cast<std::size_t>(node)] > 0 &&
            incidence[static_cast<std::size_t>(node)] < 8)
        {
            result.push_back(node);
        }
    }
    return result;
}


std::vector<int> farthestPointSample(
    const Eigen::MatrixXd& positions,
    const std::vector<int>& candidates,
    int count
)
{
    std::vector<int> selected;
    selected.reserve(static_cast<std::size_t>(count));
    selected.push_back(candidates.front());
    std::vector<double> minimumSquaredDistance(
        candidates.size(), std::numeric_limits<double>::infinity()
    );
    while (static_cast<int>(selected.size()) < count)
    {
        const int latest = selected.back();
        int best = -1;
        double bestDistance = -1.0;
        for (std::size_t index = 0; index < candidates.size(); ++index)
        {
            const int candidate = candidates[index];
            minimumSquaredDistance[index] = std::min(
                minimumSquaredDistance[index],
                (positions.row(candidate) - positions.row(latest)).squaredNorm()
            );
            const double distance = minimumSquaredDistance[index];
            if (distance > bestDistance + 1e-15 ||
                (std::abs(distance - bestDistance) <= 1e-15 &&
                 (best < 0 || candidate < best)))
            {
                bestDistance = distance;
                best = candidate;
            }
        }
        if (best < 0 || std::find(selected.begin(), selected.end(), best) !=
            selected.end())
        {
            break;
        }
        selected.push_back(best);
    }
    return selected;
}


} // namespace


namespace spin2014
{


AutomaticCageBinding buildAutomaticCageBinding(
    const AutomaticCageBindingRequest& request
)
{
    if (!request.domain)
    {
        throw std::invalid_argument("automatic cage requires a voxel domain");
    }
    validatePoints(request.outerVertices, "outer vertices");
    for (const Eigen::MatrixXd& vertices : request.voidCellVertices)
    {
        validatePoints(vertices, "void-cell vertices");
    }
    if (request.maximumHandleCount < 2)
    {
        throw std::invalid_argument(
            "automatic cage requires at least two requested handles"
        );
    }

    const VoxelGrid source = VoxelGrid::fromMaterialOccupancy(
        request.domain, 0.0, 0
    );
    const deformation::GridDeformationDomain grid =
        deformation::GridDeformationDomain::fromVoxelGrid(source);
    const std::vector<int> candidates = boundaryNodes(grid);
    if (candidates.size() < 2U)
    {
        throw std::invalid_argument(
            "automatic cage found fewer than two boundary nodes"
        );
    }
    const int handleCount = std::min(
        request.maximumHandleCount,
        static_cast<int>(candidates.size())
    );
    const std::vector<int> selected = farthestPointSample(
        grid.restNodes(), candidates, handleCount
    );

    AutomaticCageBinding result;
    result.boundaryCandidateCount = static_cast<int>(candidates.size());
    result.handles.reserve(selected.size());
    for (const int node : selected)
    {
        deformation::Handle handle;
        handle.restPivot = grid.restNodes().row(node).transpose();
        handle.constrainedNodes = {node};
        result.handles.push_back(handle);
    }
    const deformation::HandleSet handleSet(result.handles);
    const deformation::WeightSolveResult weights =
        deformation::solveBoundedBiharmonicWeights(
            grid, handleSet, request.weightOptions
        );
    if (!weights.diagnostics.converged)
    {
        throw std::runtime_error(
            "automatic cage bounded-biharmonic solve did not converge"
        );
    }
    result.nodeWeights = weights.weights;
    result.weightDiagnostics = weights.diagnostics;
    result.outerWeights = grid.interpolateWeights(
        request.outerVertices, result.nodeWeights
    );
    result.voidCellWeights.reserve(request.voidCellVertices.size());
    for (const Eigen::MatrixXd& vertices : request.voidCellVertices)
    {
        result.voidCellWeights.push_back(
            grid.interpolateWeights(vertices, result.nodeWeights)
        );
    }

    std::string bytes("spin-it-2014-automatic-cage-v1");
    bytes += request.domain->artifactId();
    appendInt(bytes, result.boundaryCandidateCount);
    for (const int node : selected) appendInt(bytes, node);
    for (Eigen::Index row = 0; row < result.nodeWeights.rows(); ++row)
    {
        for (Eigen::Index column = 0;
             column < result.nodeWeights.cols(); ++column)
        {
            appendDouble(bytes, result.nodeWeights(row, column));
        }
    }
    result.bindingHash = manufacturing::computeContentSha256(bytes);
    return result;
}


} // namespace spin2014
