#pragma once

#include "spin/SpinDesignRequest.h"

#include <Eigen/Core>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace spin {

using Vector10 = Eigen::Matrix<double, 10, 1>;
using Matrix10 = Eigen::Matrix<double, 10, 10>;

struct LegacyPredictionResponse;

// A predictor has no authority to certify geometry or a Newton derivative.
// contextIdentity binds the owned request, geometries, frame, units, densities,
// evaluator settings and profile. stateIdentity additionally binds this proposal.
struct LegacyPredictionRequest {
    int schemaVersion = 1;
    SpinEquationProfile equationProfile = SpinEquationProfile::LegacyNativeNormalizedV1;
    std::string contextIdentity, stateIdentity, checkpointIdentity;
    std::string predictorIdentity, previousTentativeIdentity;
    std::uint64_t stateVersion = 0, transactionId = 0;
    int tentativeStep = 0;
    Vector10 q = Vector10::Zero();
    // Optional engineering passthrough query, owned and billed by the host.
    // Valid only during this synchronous callback; never a persistent checkpoint.
    std::function<LegacyPredictionResponse()> queryTrueResponse;
};

struct LegacyPredictionResponse {
    int schemaVersion = 1;
    SpinEquationProfile equationProfile = SpinEquationProfile::LegacyNativeNormalizedV1;
    // Target-frame variable mass moments, density applied, fixed material excluded.
    Vector10 variableMoments = Vector10::Zero();
    // Operational direction only; NEVER a certified mathematical Jacobian.
    Matrix10 operationalDirectionMatrix = Matrix10::Zero();
    bool valid = false;
    std::string contextIdentity, stateIdentity, predictorIdentity;
    std::string errorType, diagnostic;
    double reportedSeconds = 0;
};

struct LegacyTransactionSettings {
    int schemaVersion = 1;
    int maximumTentativeSteps = 1;
    std::string predictorIdentity = "legacy-native-synchronous-predictor-v1";
};

struct LegacyTransactionEvent {
    std::uint64_t sequence = 0, callId = 0, parentCallId = 0, transactionId = 0;
    std::uint64_t trueCalls = 0, logicalTrials = 0, stateVersion = 0;
    std::string role, outcome, reason, errorType;
    std::string checkpointIdentity, stateIdentity;
    Vector10 coefficients = Vector10::Zero();
    std::optional<Vector10> tentativeMoments;
    bool geometryCertified = false;
    double wallSeconds = 0, exclusiveSeconds = 0;
    std::optional<double> reportedSeconds;
};

struct LegacyTransactionAttachment {
    int schemaVersion = 1;
    std::string contextIdentity, predictorIdentity, checkpointIdentity, scientificHash;
    std::uint64_t predictorCalls = 0, trueCalls = 0, logicalTrials = 0;
    std::uint64_t committedBlocks = 0, rolledBackBlocks = 0;
    int consecutiveFailedBlocks = 0;
    bool predictorDisabled = false;
    double rootWallSeconds = 0;
    std::optional<Vector10> lastCertifiedCoefficients;
    std::vector<LegacyTransactionEvent> events;
};

using LegacyPredictor = std::function<LegacyPredictionResponse(
    const LegacyPredictionRequest&)>;

void validateLegacyTransactionSettings(const LegacyTransactionSettings& settings);
std::string legacyPredictionStateIdentity(const LegacyPredictionRequest& request);
// Empty means the payload contract is valid, NOT that its geometry is certified.
std::string validateLegacyPrediction(const LegacyPredictionRequest& request,
                                     const LegacyPredictionResponse& response);

} // namespace spin
