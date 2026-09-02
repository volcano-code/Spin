#pragma once

#include <Eigen/Core>

#include <string>
#include <vector>

namespace reachability
{

struct Halfspaces2D
{
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
};

struct ReachabilityCell
{
    int cellIndex = -1;
    double mass = 0.0;
    Eigen::Vector2d firstMoment2D = Eigen::Vector2d::Zero();
};

struct CarvingReachabilityInput
{
    double fixedMass = 0.0;
    Eigen::Vector2d fixedFirstMoment2D = Eigen::Vector2d::Zero();
    std::vector<ReachabilityCell> variableCells;
    Eigen::MatrixXd safePolygon2D;
    double feasibilityTolerance = 1e-9;
};

enum class ReachabilityStatus
{
    AlreadySafe,
    ProvenUnreachable,
    RelaxationReachable,
    Indeterminate
};

struct ReachabilityCertificate
{
    Eigen::VectorXd weights;
    double minimumWeight = 0.0;
    double normalizationResidual = 0.0;
    double stationarityInfinityNorm = 0.0;
    double contradiction = 0.0;
    bool valid = false;
};

struct CarvingReachabilityResult
{
    ReachabilityStatus status = ReachabilityStatus::Indeterminate;
    int inputCellCount = 0;
    Eigen::VectorXd retainedFractions;
    double finalMass = 0.0;
    Eigen::Vector2d finalFirstMoment2D = Eigen::Vector2d::Zero();
    Eigen::Vector2d finalCenter2D = Eigen::Vector2d::Zero();
    double maximumPrimalViolation = 0.0;
    ReachabilityCertificate certificate;
    std::string primalSolverStatus;
    std::string certificateSolverStatus;
};

Halfspaces2D safePolygonHalfspaces(const Eigen::MatrixXd& polygon);

bool contains(
    const Halfspaces2D& halfspaces,
    const Eigen::Vector2d& point,
    double tolerance = 0.0
);

void validateInput(const CarvingReachabilityInput& input);

ReachabilityCertificate validateFarkasCertificate(
    const Eigen::MatrixXd& C,
    const Eigen::VectorXd& d,
    const Eigen::VectorXd& weights,
    double tolerance
);

const char* toString(ReachabilityStatus status) noexcept;

class CarvingReachabilityAnalyzer
{
public:
    CarvingReachabilityResult analyze(
        const CarvingReachabilityInput& input
    ) const;
};

} // namespace reachability
