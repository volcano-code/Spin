#pragma once
#include "spin2014/AdaptiveFillFunctional.h"
#include <limits>
#include <string>
#include <vector>

namespace spin2014 {
// All controls below are IMPLEMENTATION_DEFAULT, not recovered author values.
struct SlqpAuditSettings {
    double restorationTrustBound = 1.;
    int lpIterationLimit = 10000;
    double lpTimeLimit = 10.;
    double lpFeasibilityTolerance = 1e-9;
};
struct KktResidual {
    bool valid = false;
    std::string status = "INVALID_DATA";
    Eigen::VectorXd lambda, alpha, eta;
    double equalityPrimal = std::numeric_limits<double>::infinity();
    double boundPrimal = std::numeric_limits<double>::infinity();
    double stationarity = std::numeric_limits<double>::infinity();
    double dualSign = std::numeric_limits<double>::infinity();
    double complementarity = std::numeric_limits<double>::infinity();
    double full = std::numeric_limits<double>::infinity();
};
struct RestorationTrial {
    double step = 0.;
    double constraintL1 = std::numeric_limits<double>::infinity();
    bool valid = false;
    bool accepted = false;
    std::string reason;
};
struct RestorationEvidence {
    bool attempted = false;
    bool accepted = false;
    std::string status;
    std::string reason;
    Eigen::VectorXd direction;
    double predictedResidual = std::numeric_limits<double>::infinity();
    std::vector<RestorationTrial> trials;
};
bool validSlqpEvaluation(const AdaptiveFillEvaluation&, Eigen::Index variables);
void validateSlqpAuditSettings(const SlqpAuditSettings&);
KktResidual evaluateKkt(const AdaptiveFillEvaluation&, const Eigen::VectorXd& x,
    Eigen::Index betaCount, const Eigen::VectorXd& lambda,
    const Eigen::VectorXd& alpha, const Eigen::VectorXd& eta);
KktResidual computeKktWitness(const AdaptiveFillEvaluation&, const Eigen::VectorXd& x,
    Eigen::Index betaCount, double activeTolerance, const SlqpAuditSettings&);
RestorationEvidence solveRestoration(const AdaptiveFillEvaluation&, const Eigen::VectorXd& x,
    Eigen::Index betaCount, const SlqpAuditSettings&);
} // namespace spin2014
