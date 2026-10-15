#pragma once

#include "spin2014/SlqpHistory.h"
#include <Eigen/Core>
#include <string>
#include <vector>

namespace spin2014 {
struct TangentQpSettings {
    // Primal V2 only: -1 means one possible working face per bounded variable.
    // Numeric >=0 budgets retain their legacy meaning; iteration/time/memory
    // limits remain separate. Other QP policies reject the -1 sentinel.
    int maximumFaces = 256;
    int maximumIterations = 1028;
    double activeTolerance = 1e-10;
    double residualTolerance = 1e-9;
};
struct TangentQpStep {
    std::vector<int> workingSet;
    std::string action;
    int changedIndex = -1;
    double equalityResidual = 0.;
    double stationarity = 0.;
};
struct TangentQpResult {
    bool attempted = false;
    bool checked = false;
    std::string status = "NOT_ATTEMPTED";
    Eigen::VectorXd direction;
    std::vector<int> workingSet;
    Eigen::VectorXd multipliers;
    double equalityResidual = 0.;
    double stationarity = 0.;
    std::vector<TangentQpStep> trace;
};

// Minimize g'd + .5 d'Bd subject to a*d+c=0 and the tangent cone
// of the [0,1] box on x.head(bounded). B is the inverse of history's
// inverse L-BFGS action. Working-set rows are +e_i, so their multiplier
// signs are nonpositive at lower faces and nonnegative at upper faces.
// Interior coordinates still need a separate finite-step box limit.
// Unknown results are neither directions certified for use nor proofs
// of infeasibility. Trace actions describe the failure or set transition.
TangentQpResult solveTangentQp(const Eigen::VectorXd& x,
                             const Eigen::VectorXd& g,
                             const Eigen::MatrixXd& a,
                             const Eigen::VectorXd& c,
                             Eigen::Index bounded,
                             const LbfgsHistory& history,
                             const TangentQpSettings& settings = {});

// Full displacement box: -x_i <= d_i <= 1-x_i. Same trace carrier,
// distinguished by BOX_QP_* status; face IDs are 2*i (lower), 2*i+1 (upper).
// Uses bounded Schur working sets and matrix-free BFGS, not a dense n*n B.
// Rank deficiency, cycling and budgets return UNKNOWN, never infeasibility.
TangentQpResult solveFiniteBoxQp(const Eigen::VectorXd& x,
                               const Eigen::VectorXd& g,
                               const Eigen::MatrixXd& a,
                               const Eigen::VectorXd& c,
                               Eigen::Index bounded,
                               const LbfgsHistory& history,
                               const TangentQpSettings& settings = {});

// V2: keep a feasible primal iterate, move only to the first blocking face,
// and release a wrong-sign multiplier only at a working-set minimizer.
// Uses the same bounds, metric and raw KKT tolerances as V1; UNKNOWN is not
// an infeasibility certificate. Zero-objective V1 supplies Phase I if needed.
TangentQpResult solvePrimalBoxQp(const Eigen::VectorXd& x,
                               const Eigen::VectorXd& g,
                               const Eigen::MatrixXd& a,
                               const Eigen::VectorXd& c,
                               Eigen::Index bounded,
                               const LbfgsHistory& history,
                               const TangentQpSettings& settings = {});
}
