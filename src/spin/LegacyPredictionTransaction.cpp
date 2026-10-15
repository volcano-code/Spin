#include "spin/LegacyPredictionTransaction.h"

#include "manufacturing/SlicerValidation.h"

#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace spin {
namespace {
void appendIdentity(std::ostream& out, const std::string& identity) {
    // Length prefixes prevent ambiguous concatenations, including embedded NULs.
    out << identity.size() << ':';
    out.write(identity.data(), static_cast<std::streamsize>(identity.size()));
    out << ';';
}
std::string requestError(const LegacyPredictionRequest& request) {
    if (request.schemaVersion != 1) return "REQUEST_SCHEMA_UNSUPPORTED";
    if (request.equationProfile != SpinEquationProfile::LegacyNativeNormalizedV1)
        return "REQUEST_PROFILE_NOT_LEGACY";
    if (request.contextIdentity.empty() || request.checkpointIdentity.empty() ||
        request.predictorIdentity.empty()) return "REQUEST_IDENTITY_MISSING";
    if (request.stateVersion == 0 || request.transactionId == 0 ||
        request.tentativeStep < 1 || request.tentativeStep > 3)
        return "REQUEST_STATE_INVALID";
    if ((request.tentativeStep == 1 && !request.previousTentativeIdentity.empty()) ||
        (request.tentativeStep > 1 && request.previousTentativeIdentity.empty()))
        return "REQUEST_CHAIN_INVALID";
    if (!request.q.allFinite() || std::abs(request.q.norm() - 1.0) > 1e-12)
        return "REQUEST_Q_INVALID";
    return {};
}
} // namespace

void validateLegacyTransactionSettings(const LegacyTransactionSettings& settings) {
    if (settings.schemaVersion != 1 || settings.maximumTentativeSteps < 1 ||
        settings.maximumTentativeSteps > 3 || settings.predictorIdentity.empty())
        throw std::invalid_argument("LEGACY transactions require schema 1 and one to three steps");
}

std::string legacyPredictionStateIdentity(const LegacyPredictionRequest& request) {
    const auto error = requestError(request);
    if (!error.empty()) throw std::invalid_argument(error);
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << "LegacyPredictionState-v1;" << request.schemaVersion << ';';
    appendIdentity(out, spinEquationProfileToken(request.equationProfile));
    appendIdentity(out, request.contextIdentity);
    appendIdentity(out, request.checkpointIdentity);
    appendIdentity(out, request.predictorIdentity);
    appendIdentity(out, request.previousTentativeIdentity);
    out << request.stateVersion << ';' << request.transactionId << ';'
        << request.tentativeStep << ';' << std::hexfloat;
    for (int i = 0; i < 10; ++i) out << request.q(i) << ';';
    return manufacturing::computeContentSha256(out.str());
}

std::string validateLegacyPrediction(const LegacyPredictionRequest& request,
                                     const LegacyPredictionResponse& response) {
    const auto error = requestError(request);
    if (!error.empty()) return error;
    if (request.stateIdentity.empty() ||
        request.stateIdentity != legacyPredictionStateIdentity(request))
        return "REQUEST_STATE_IDENTITY_MISMATCH";
    if (response.schemaVersion != 1) return "RESPONSE_SCHEMA_UNSUPPORTED";
    if (response.equationProfile != SpinEquationProfile::LegacyNativeNormalizedV1)
        return "RESPONSE_PROFILE_NOT_LEGACY";
    if (response.contextIdentity != request.contextIdentity)
        return "RESPONSE_CONTEXT_MISMATCH";
    if (response.stateIdentity != request.stateIdentity)
        return "RESPONSE_STATE_MISMATCH";
    if (response.predictorIdentity != request.predictorIdentity)
        return "RESPONSE_PREDICTOR_MISMATCH";
    if (!response.valid) return "RESPONSE_DECLARED_INVALID";
    if (!response.variableMoments.allFinite()) return "RESPONSE_MOMENTS_NONFINITE";
    if (!response.operationalDirectionMatrix.allFinite())
        return "RESPONSE_OPERATIONAL_MATRIX_NONFINITE";
    if (!std::isfinite(response.reportedSeconds) || response.reportedSeconds < 0)
        return "RESPONSE_TIMING_INVALID";
    return {};
}

} // namespace spin
