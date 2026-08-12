#pragma once

#include <Eigen/Core>

#include <vector>


namespace deformation
{

class GridDeformationDomain;
class HandleSet;


/** Options for the open, graph-based bounded-biharmonic weight solve. */
struct WeightSolveOptions
{
    int maxIterations = 200;
    double tolerance = 1.0e-8;
};


/** Terminal state reported by libigl's open active-set solver. */
enum class WeightSolverStatus
{
    Converged,
    MaximumIterationsReached,
    Error
};


const char* weightSolverStatusName(WeightSolverStatus status);


/** Auditable feasibility and convergence checks for one weight solve. */
struct WeightSolveDiagnostics
{
    std::vector<WeightSolverStatus> statusByHandle;
    bool converged = false;
    /** True when feasibility was audited without running the local solver. */
    bool verifiedImport = false;
    double maximumBoundViolation = 0.0;
    double maximumRowSumError = 0.0;
    double maximumConstraintError = 0.0;
};


struct WeightSolveResult
{
    /** One row per graph node and one column per handle. */
    Eigen::MatrixXd weights;
    WeightSolveDiagnostics diagnostics;
};


/** Canonical descriptive result name used by the synchronized R5 facade. */
using BoundedBiharmonicWeightResult = WeightSolveResult;


/**
 * Solve bounded biharmonic weights on an undirected, unweighted graph.
 *
 * The graph Laplacian has vertex valence on its diagonal and -1 for every
 * edge.  The quadratic objective is L^T L.  constraintSets[h] contains the
 * one-hot Dirichlet nodes for handle h.  Every connected component must be
 * anchored by at least one handle.
 */
WeightSolveResult solveBoundedBiharmonicWeights(
    const Eigen::MatrixXd& restNodes,
    const Eigen::MatrixXi& edges,
    const std::vector<std::vector<int>>& constraintSets,
    const WeightSolveOptions& options = WeightSolveOptions()
);


/** Solve weights using the canonical R5 grid-domain and handle containers. */
WeightSolveResult solveBoundedBiharmonicWeights(
    const GridDeformationDomain& domain,
    const HandleSet& handles,
    const WeightSolveOptions& options = WeightSolveOptions()
);

} // namespace deformation
