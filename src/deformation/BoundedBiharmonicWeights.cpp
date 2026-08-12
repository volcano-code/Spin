#include "deformation/BoundedBiharmonicWeights.h"

#include "deformation/GridDeformationDomain.h"
#include "deformation/HandleSet.h"

#include <igl/SolverStatus.h>
#include <igl/active_set.h>
#include <igl/min_quad_with_fixed.h>

#include <Eigen/Sparse>

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>


namespace deformation
{
namespace
{

struct ValidatedGraph
{
    Eigen::SparseMatrix<double> quadraticEnergy;
    Eigen::VectorXi constrainedNodes;
    Eigen::MatrixXd boundaryValues;
};


void validateOptions(const WeightSolveOptions& options)
{
    if (options.maxIterations <= 0)
    {
        throw std::invalid_argument(
            "weight solve requires a positive iteration limit"
        );
    }

    if (
        !std::isfinite(options.tolerance) ||
        options.tolerance <= 0.0
    )
    {
        throw std::invalid_argument(
            "weight solve tolerance must be finite and positive"
        );
    }
}


void validateRestNodes(const Eigen::MatrixXd& restNodes)
{
    if (restNodes.rows() == 0 || restNodes.cols() != 3)
    {
        throw std::invalid_argument(
            "weight domain must contain three-dimensional rest nodes"
        );
    }

    for (Eigen::Index row = 0; row < restNodes.rows(); ++row)
    {
        for (Eigen::Index column = 0; column < restNodes.cols(); ++column)
        {
            if (!std::isfinite(restNodes(row, column)))
            {
                throw std::invalid_argument(
                    "weight-domain rest nodes must be finite"
                );
            }
        }
    }
}


ValidatedGraph validateAndBuildGraph(
    int nodeCount,
    const Eigen::MatrixXi& edges,
    const std::vector<std::vector<int>>& constraintSets
)
{
    if (edges.cols() != 2)
    {
        throw std::invalid_argument(
            "weight-domain edges must have exactly two columns"
        );
    }

    if (constraintSets.empty())
    {
        throw std::invalid_argument(
            "bounded biharmonic weights require at least one handle"
        );
    }

    std::set<std::pair<int, int>> uniqueEdges;
    std::vector<std::vector<int>> adjacency(
        static_cast<std::size_t>(nodeCount)
    );

    for (Eigen::Index row = 0; row < edges.rows(); ++row)
    {
        int first = edges(row, 0);
        int second = edges(row, 1);

        if (
            first < 0 || first >= nodeCount ||
            second < 0 || second >= nodeCount
        )
        {
            throw std::invalid_argument(
                "weight-domain edge endpoint is outside the node domain"
            );
        }

        if (first == second)
        {
            throw std::invalid_argument(
                "weight-domain graph must not contain self edges"
            );
        }

        if (first > second)
        {
            std::swap(first, second);
        }

        if (uniqueEdges.emplace(first, second).second)
        {
            adjacency[static_cast<std::size_t>(first)].push_back(second);
            adjacency[static_cast<std::size_t>(second)].push_back(first);
        }
    }

    std::vector<int> owner(static_cast<std::size_t>(nodeCount), -1);
    int constrainedNodeCount = 0;

    for (
        std::size_t handleIndex = 0;
        handleIndex < constraintSets.size();
        ++handleIndex
    )
    {
        if (constraintSets[handleIndex].empty())
        {
            throw std::invalid_argument(
                "every deformation handle must constrain at least one node"
            );
        }

        std::set<int> uniqueConstraints;
        for (const int node : constraintSets[handleIndex])
        {
            if (node < 0 || node >= nodeCount)
            {
                throw std::invalid_argument(
                    "handle constraint is outside the weight domain"
                );
            }

            if (!uniqueConstraints.emplace(node).second)
            {
                throw std::invalid_argument(
                    "a handle must not repeat a constrained node"
                );
            }

            int& nodeOwner = owner[static_cast<std::size_t>(node)];
            if (nodeOwner >= 0)
            {
                throw std::invalid_argument(
                    "a weight-domain node cannot belong to multiple handles"
                );
            }

            nodeOwner = static_cast<int>(handleIndex);
            ++constrainedNodeCount;
        }
    }

    std::vector<unsigned char> visited(
        static_cast<std::size_t>(nodeCount),
        0U
    );

    for (int start = 0; start < nodeCount; ++start)
    {
        if (visited[static_cast<std::size_t>(start)] != 0U)
        {
            continue;
        }

        bool componentHasConstraint = false;
        std::queue<int> pending;
        pending.push(start);
        visited[static_cast<std::size_t>(start)] = 1U;

        while (!pending.empty())
        {
            const int node = pending.front();
            pending.pop();

            componentHasConstraint =
                componentHasConstraint ||
                owner[static_cast<std::size_t>(node)] >= 0;

            for (const int neighbor : adjacency[static_cast<std::size_t>(node)])
            {
                if (visited[static_cast<std::size_t>(neighbor)] == 0U)
                {
                    visited[static_cast<std::size_t>(neighbor)] = 1U;
                    pending.push(neighbor);
                }
            }
        }

        if (!componentHasConstraint)
        {
            throw std::invalid_argument(
                "every graph component must contain a handle constraint"
            );
        }
    }

    std::vector<Eigen::Triplet<double>> laplacianEntries;
    laplacianEntries.reserve(
        static_cast<std::size_t>(nodeCount) +
        2U * uniqueEdges.size()
    );

    for (int node = 0; node < nodeCount; ++node)
    {
        laplacianEntries.emplace_back(
            node,
            node,
            static_cast<double>(
                adjacency[static_cast<std::size_t>(node)].size()
            )
        );
    }

    for (const auto& edge : uniqueEdges)
    {
        laplacianEntries.emplace_back(edge.first, edge.second, -1.0);
        laplacianEntries.emplace_back(edge.second, edge.first, -1.0);
    }

    Eigen::SparseMatrix<double> laplacian(nodeCount, nodeCount);
    laplacian.setFromTriplets(
        laplacianEntries.begin(),
        laplacianEntries.end()
    );

    ValidatedGraph result;
    result.quadraticEnergy = laplacian.transpose() * laplacian;
    result.quadraticEnergy.makeCompressed();
    result.constrainedNodes.resize(constrainedNodeCount);
    result.boundaryValues = Eigen::MatrixXd::Zero(
        constrainedNodeCount,
        static_cast<Eigen::Index>(constraintSets.size())
    );

    int boundaryRow = 0;
    for (
        std::size_t handleIndex = 0;
        handleIndex < constraintSets.size();
        ++handleIndex
    )
    {
        for (const int node : constraintSets[handleIndex])
        {
            result.constrainedNodes(boundaryRow) = node;
            result.boundaryValues(
                boundaryRow,
                static_cast<Eigen::Index>(handleIndex)
            ) = 1.0;
            ++boundaryRow;
        }
    }

    return result;
}


WeightSolverStatus convertStatus(igl::SolverStatus status)
{
    switch (status)
    {
        case igl::SOLVER_STATUS_CONVERGED:
            return WeightSolverStatus::Converged;
        case igl::SOLVER_STATUS_MAX_ITER:
            return WeightSolverStatus::MaximumIterationsReached;
        case igl::SOLVER_STATUS_ERROR:
        default:
            return WeightSolverStatus::Error;
    }
}


void clampToleranceViolations(
    Eigen::MatrixXd& weights,
    double tolerance
)
{
    for (Eigen::Index row = 0; row < weights.rows(); ++row)
    {
        for (Eigen::Index column = 0; column < weights.cols(); ++column)
        {
            double& value = weights(row, column);

            if (value < 0.0 && value >= -tolerance)
            {
                value = 0.0;
            }
            else if (value > 1.0 && value <= 1.0 + tolerance)
            {
                value = 1.0;
            }

            if (std::abs(value) <= tolerance)
            {
                value = 0.0;
            }
            else if (std::abs(value - 1.0) <= tolerance)
            {
                value = 1.0;
            }
        }
    }
}


void normalizeRows(Eigen::MatrixXd& weights, double tolerance)
{
    for (Eigen::Index row = 0; row < weights.rows(); ++row)
    {
        const double rowSum = weights.row(row).sum();

        if (!std::isfinite(rowSum) || rowSum <= tolerance)
        {
            throw std::runtime_error(
                "bounded biharmonic weights have a non-positive row sum"
            );
        }

        weights.row(row) /= rowSum;
    }
}


void rewriteConstraints(
    Eigen::MatrixXd& weights,
    const Eigen::VectorXi& constrainedNodes,
    const Eigen::MatrixXd& boundaryValues
)
{
    for (
        Eigen::Index boundaryRow = 0;
        boundaryRow < constrainedNodes.rows();
        ++boundaryRow
    )
    {
        weights.row(constrainedNodes(boundaryRow)) =
            boundaryValues.row(boundaryRow);
    }
}


void auditWeights(
    const Eigen::MatrixXd& weights,
    const Eigen::VectorXi& constrainedNodes,
    const Eigen::MatrixXd& boundaryValues,
    double tolerance,
    WeightSolveDiagnostics& diagnostics
)
{
    diagnostics.maximumBoundViolation = 0.0;
    diagnostics.maximumRowSumError = 0.0;
    diagnostics.maximumConstraintError = 0.0;

    for (Eigen::Index row = 0; row < weights.rows(); ++row)
    {
        double rowSum = 0.0;
        for (Eigen::Index column = 0; column < weights.cols(); ++column)
        {
            const double value = weights(row, column);
            if (!std::isfinite(value))
            {
                throw std::runtime_error(
                    "bounded biharmonic solve produced a non-finite weight"
                );
            }

            diagnostics.maximumBoundViolation = std::max(
                diagnostics.maximumBoundViolation,
                std::max(0.0, std::max(-value, value - 1.0))
            );
            rowSum += value;
        }

        diagnostics.maximumRowSumError = std::max(
            diagnostics.maximumRowSumError,
            std::abs(rowSum - 1.0)
        );
    }

    for (
        Eigen::Index boundaryRow = 0;
        boundaryRow < constrainedNodes.rows();
        ++boundaryRow
    )
    {
        const int node = constrainedNodes(boundaryRow);
        for (Eigen::Index column = 0; column < weights.cols(); ++column)
        {
            diagnostics.maximumConstraintError = std::max(
                diagnostics.maximumConstraintError,
                std::abs(
                    weights(node, column) -
                    boundaryValues(boundaryRow, column)
                )
            );
        }
    }

    if (
        diagnostics.maximumBoundViolation > tolerance ||
        diagnostics.maximumRowSumError > tolerance ||
        diagnostics.maximumConstraintError > tolerance
    )
    {
        throw std::runtime_error(
            "bounded biharmonic weights failed feasibility validation"
        );
    }
}

} // namespace


const char* weightSolverStatusName(WeightSolverStatus status)
{
    switch (status)
    {
        case WeightSolverStatus::Converged:
            return "converged";
        case WeightSolverStatus::MaximumIterationsReached:
            return "maximum-iterations-reached";
        case WeightSolverStatus::Error:
            return "error";
    }

    return "unknown";
}


WeightSolveResult solveBoundedBiharmonicWeights(
    const Eigen::MatrixXd& restNodes,
    const Eigen::MatrixXi& edges,
    const std::vector<std::vector<int>>& constraintSets,
    const WeightSolveOptions& options
)
{
    validateOptions(options);
    validateRestNodes(restNodes);

    const int nodeCount = static_cast<int>(restNodes.rows());
    const int handleCount = static_cast<int>(constraintSets.size());
    const ValidatedGraph graph = validateAndBuildGraph(
        nodeCount,
        edges,
        constraintSets
    );

    if (graph.constrainedNodes.size() == nodeCount)
    {
        WeightSolveResult result;
        result.weights = Eigen::MatrixXd::Zero(nodeCount, handleCount);
        rewriteConstraints(
            result.weights,
            graph.constrainedNodes,
            graph.boundaryValues
        );
        result.diagnostics.statusByHandle.assign(
            static_cast<std::size_t>(handleCount),
            WeightSolverStatus::Converged
        );
        result.diagnostics.converged = true;
        auditWeights(
            result.weights,
            graph.constrainedNodes,
            graph.boundaryValues,
            options.tolerance,
            result.diagnostics
        );
        return result;
    }

    WeightSolveResult result;
    result.diagnostics.statusByHandle.assign(
        static_cast<std::size_t>(handleCount),
        WeightSolverStatus::Converged
    );

    if (handleCount == 1)
    {
        result.weights = Eigen::MatrixXd::Ones(nodeCount, 1);
        result.diagnostics.converged = true;
        auditWeights(
            result.weights,
            graph.constrainedNodes,
            graph.boundaryValues,
            options.tolerance,
            result.diagnostics
        );
        return result;
    }

    const Eigen::MatrixXd linearTerm = Eigen::MatrixXd::Zero(
        nodeCount,
        handleCount
    );
    const Eigen::SparseMatrix<double> noEqualities(0, nodeCount);
    const Eigen::MatrixXd noEqualityValues(0, handleCount);

    if (!igl::min_quad_with_fixed(
        graph.quadraticEnergy,
        linearTerm,
        graph.constrainedNodes,
        graph.boundaryValues,
        noEqualities,
        noEqualityValues,
        true,
        result.weights
    ))
    {
        throw std::runtime_error(
            "unconstrained biharmonic initialization failed"
        );
    }

    const Eigen::VectorXd zeroLinearTerm = Eigen::VectorXd::Zero(nodeCount);
    const Eigen::SparseMatrix<double> noInequalities(0, nodeCount);
    const Eigen::VectorXd noEqualityVector(0);
    const Eigen::VectorXd noInequalityVector(0);
    const Eigen::VectorXd lowerBounds = Eigen::VectorXd::Zero(nodeCount);
    const Eigen::VectorXd upperBounds = Eigen::VectorXd::Ones(nodeCount);

    igl::active_set_params solverOptions;
    solverOptions.Auu_pd = true;
    solverOptions.max_iter = options.maxIterations;
    solverOptions.inactive_threshold = options.tolerance;
    solverOptions.constraint_threshold = options.tolerance;
    solverOptions.solution_diff_threshold =
        options.tolerance * options.tolerance;

    bool allConverged = true;

    for (int handle = 0; handle < handleCount; ++handle)
    {
        Eigen::VectorXd column = result.weights.col(handle);
        const Eigen::VectorXd boundaryColumn =
            graph.boundaryValues.col(handle);

        const igl::SolverStatus rawStatus = igl::active_set(
            graph.quadraticEnergy,
            zeroLinearTerm,
            graph.constrainedNodes,
            boundaryColumn,
            noEqualities,
            noEqualityVector,
            noInequalities,
            noInequalityVector,
            lowerBounds,
            upperBounds,
            solverOptions,
            column
        );

        const WeightSolverStatus status = convertStatus(rawStatus);
        result.diagnostics.statusByHandle[static_cast<std::size_t>(handle)] =
            status;

        if (status == WeightSolverStatus::Error)
        {
            throw std::runtime_error(
                "active-set bounded weight solve failed for handle " +
                std::to_string(handle)
            );
        }

        allConverged =
            allConverged && status == WeightSolverStatus::Converged;
        result.weights.col(handle) = column;
    }

    result.diagnostics.converged = allConverged;

    clampToleranceViolations(result.weights, options.tolerance);
    normalizeRows(result.weights, options.tolerance);
    rewriteConstraints(
        result.weights,
        graph.constrainedNodes,
        graph.boundaryValues
    );
    clampToleranceViolations(result.weights, options.tolerance);
    auditWeights(
        result.weights,
        graph.constrainedNodes,
        graph.boundaryValues,
        options.tolerance,
        result.diagnostics
    );

    return result;
}


WeightSolveResult solveBoundedBiharmonicWeights(
    const GridDeformationDomain& domain,
    const HandleSet& handles,
    const WeightSolveOptions& options
)
{
    return solveBoundedBiharmonicWeights(
        domain.restNodes(),
        domain.edges(),
        handles.constraintNodeSets(),
        options
    );
}

} // namespace deformation
