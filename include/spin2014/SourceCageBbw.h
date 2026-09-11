#pragma once
#include "spin2014/SourceCageFem.h"

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
enum class SourceCageBbwProtocol { COUPLED_PARTITION, INDEPENDENT_NORMALIZED };
enum class SourceCageBbwTerminal {
    VALIDATED, INVALID_CONFIG, RESOURCE_LIMIT, NUMERICAL_ERROR,
    BACKEND_ERROR, ITERATION_LIMIT, TIME_LIMIT, KKT_REJECTED, NORMALIZATION_REJECTED
};
struct SourceCageBbwConfig {
    SourceCageBbwProtocol protocol=SourceCageBbwProtocol::COUPLED_PARTITION;
    double primalTolerance=1e-9;
    double optimalityTolerance=1e-8;
    double backendTolerance=1e-10;
    double timeLimitSeconds=300; // Per backend call; no implicit retry.
    int iterationLimit=100000;
    std::size_t maximumQPairContributions=50000000;
    std::size_t maximumQNonzeros=10000000;
    std::size_t maximumVariables=1000000;
    // Counts the conservative simultaneously-live dense workspace bound,
    // including retained raw/final fields and FEM evaluation temporaries.
    std::size_t maximumDenseValues=20000000;
    std::size_t maximumBytes=1073741824;
};
struct SourceCageBbwAssessment {
    bool available=false, primalPassed=false, optimalityPassed=false;
    std::string reason;
    double energy=0, energyScale=0;
    double boundaryResidual=0, boundsResidual=0, partitionResidual=0;
    double stationarityResidual=0, dualResidual=0, complementarityResidual=0;
    double projectedGradientResidual=0;
    double unscaledStationarityResidual=0, unscaledDualResidual=0;
    // The unscaled projected residual is scale * ||w-P(w-(g+lambda)/scale)||,
    // i.e. the gradient mapping with step 1/scale, not a changed candidate.
    double unscaledComplementarityResidual=0, unscaledProjectedGradientResidual=0;
};
struct SourceCageBbwAttempt {
    int handle=-1; // -1 is the complete coupled problem.
    int passModelStatus=-1, returnStatus=-1, modelStatus=-1;
    long long qpIterations=0;
    bool runInvoked=false, candidateAvailable=false, diagnosticsTruncated=false;
    double startUnixSeconds=0, endUnixSeconds=0, elapsedSeconds=0;
    std::string terminalReason, diagnostics;
    SourceCageBbwAssessment assessment;
};
/** Evaluates the supplied field using the original FEM operator. No caller
 * status or backend dual is trusted. Independent mode assesses RAW columns,
 * without imposing partition; coupled mode includes all free-node row sums.
 * Invalid shape/nonfinite data yields available=false, never a valid report.
 */
SourceCageBbwAssessment assessSourceCageBbw(
    const SourceCageFem&, const Eigen::MatrixXd&, const SourceCageBbwConfig&);

/** Actual immutable run of a reconstructed BBW QP using pinned HiGHS qpasm.
 * No constructor accepts weights or a success flag. Boundary elimination is
 * exact. Independent normalization is explicit and only follows all-column
 * validation; final coupled optimality is assessed separately. Traces contain
 * actual backend calls/statuses/iteration counts, not internal iterates.
 * The per-call diagnostic text is capped at 65536 bytes with an explicit
 * truncation flag. Timing and diagnostic text are excluded from the scientific
 * digest; settings, statuses, iteration counts, assessments and fields bind it.
 * Independent QPs always include a redundant 0=0 row to avoid the pinned
 * backend's defective zero-row bounded-start path. This preserves the exact
 * objective/feasible set and is never a fallback or a post-solve repair.
 * Failure feasibility is UNKNOWN. No geometry/material/physical certificate.
 */
class SourceCageBbw {
public:
    SourceCageBbw(const SourceCageFem&,const SourceCageBbwConfig& = {});
    bool success() const;
    SourceCageBbwTerminal terminal() const;
    const std::string& terminalReason() const;
    const SourceCageBbwConfig& config() const;
    const std::string& backendSettings() const;
    const SourceCageFem& fem() const;
    // On incomplete solves raw may contain partial candidate columns; only
    // attempts with candidateAvailable identify observed backend columns.
    // Unsolved free entries are workspace zeros, never validated weights.
    const Eigen::MatrixXd& rawWeights() const;
    const Eigen::MatrixXd& finalWeights() const;
    const SourceCageBbwAssessment& rawAssessment() const;
    const SourceCageBbwAssessment& finalAssessment() const;
    const std::vector<SourceCageBbwAttempt>& attempts() const;
    double maximumNormalizationChange() const;
    double energyScale() const;
    std::size_t qPairContributions() const;
    std::size_t qNonzeros() const;
    std::size_t variablesPerCoupledProblem() const;
    std::size_t checkedStorageBytes() const;
    const std::string& configSha256() const;
    const std::string& outputSha256() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
};
}
#endif
