#pragma once
#include "spin2014/SourceAdaptiveCheckpoint.h"
#include <map>
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
// The persisted stage directory is zero based: stage-00000019 has completed
// twenty passes.  Keep this binding constant next to the read-only API so a
// caller cannot silently turn the G7 stage index into a pass count.
inline constexpr int spin2014G7PostSolvePrePassCount=20;
// Versioned metadata (<=64 MiB), streamed sealed histories only. This codec
// does no publication, authentication or wall-clock accounting. A caller must
// bind code/binary/input identities and durably publish these bytes together
// with immutable trace references (e.g. StageCheckpointStore).
std::string encodeSourceAdaptiveCheckpoint(const SourceAdaptiveCheckpoint&);
// Reconstruct source-derived moments and split/merge history, check evaluations
// and final KKT witnesses; do not solve earlier stages. No output file is made.
// Same-environment bitwise restart only; rejects changed numerical results.
std::shared_ptr<const SourceAdaptiveCheckpoint> decodeSourceAdaptiveCheckpoint(
    const SourceAdaptiveSpinRequest&,const std::string&);

// Read-only reconstruction of one already-completed SLQP pass. The checkpoint
// and trace are authenticated before any transition is replayed. This API
// cannot resume, solve, round, export geometry, or grant scientific status.
struct SourceAdaptivePostSolveReplay {
    SourceAdaptiveSpinResult result;
    std::map<std::string,std::uint64_t> timingsNanoseconds;
    // traceRecordStreamSha256 hashes encoded iteration records, not trace-file
    // container bytes.
    std::string checkpointSha256,traceRecordStreamSha256,encodedCheckpointSha256,diagnosticBytes;
    Eigen::VectorXd terminalState;
    KktResidual terminalKkt;
    bool optimizerInvoked=false;
    std::string binaryValidation="NOT_RUN",scientificAcceptance="NOT_GRANTED";
    std::string scope="SUPPLEMENTAL_DIAGNOSTIC_ONLY";
};
SourceAdaptivePostSolveReplay replaySourceAdaptivePostSolve(
    const SourceAdaptiveSpinRequest&,const std::string& materialCheckpointBytes,
    const std::string& expectedMaterialSettingsHash,const std::string& tracePath,int expectedPrePassCount);
}
#endif
