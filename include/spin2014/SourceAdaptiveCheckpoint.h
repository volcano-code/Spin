#pragma once
#include "spin2014/SourceAdaptiveSpinSolver.h"
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
namespace detail { struct SourceCheckpointAccess; }
struct SourceAdaptivePostSolveReplay;
SourceAdaptivePostSolveReplay replaySourceAdaptivePostSolve(
    const SourceAdaptiveSpinRequest&,const std::string&,const std::string&,const std::string&,int);
/** Validated boundary token. Created by the solver or source-aware codec.
 * Durable publication and resource supervision belong to an outer runner.
 * Full completed pass metadata is retained; sealed trace files stay external.
 * Copying this token does not consume it or reset global pass/trace limits.
 */
class SourceAdaptiveCheckpoint {
public:
    std::size_t completedPasses() const {return history_.passes.size();}
    std::uint64_t storedTraceBytes() const {return storedBytes_;}
    const std::string& requestHash() const {return history_.inputHash;}
private:
    friend struct detail::SourceCheckpointAccess;
    friend std::string encodeSourceAdaptiveCheckpoint(const SourceAdaptiveCheckpoint&);
    friend std::shared_ptr<const SourceAdaptiveCheckpoint> decodeSourceAdaptiveCheckpoint(
        const SourceAdaptiveSpinRequest&,const std::string&);
    // Read-only diagnostic access is deliberately narrower than the existing
    // solver checkpoint handoff: it cannot construct a resume control.
    friend SourceAdaptivePostSolveReplay replaySourceAdaptivePostSolve(
        const SourceAdaptiveSpinRequest&,const std::string&,const std::string&,const std::string&,int);
    SourceAdaptiveCheckpoint(const SourceAdaptiveSpinResult& history,std::uint64_t bytes,
        std::uint64_t maximumBytes,bool disk,std::vector<std::string> digests):history_(history),storedBytes_(bytes),
        maximumBytes_(maximumBytes),disk_(disk),traceDigests_(std::move(digests)){}
    const SourceAdaptiveSpinResult history_;
    const std::uint64_t storedBytes_,maximumBytes_;
    const bool disk_;
    const std::vector<std::string> traceDigests_;
};
struct SourceAdaptiveStageControl {
    // Execution slice only: not a new scientific stopping rule. Global
    // maximumAdaptivePasses is unchanged and bound by the request hash.
    int maximumNewPasses=1;
    std::shared_ptr<const SourceAdaptiveCheckpoint> resume;
};
struct SourceAdaptiveStageRun {
    SourceAdaptiveSpinResult result;
    std::shared_ptr<const SourceAdaptiveCheckpoint> checkpoint;
};
// Streaming resumes need a fresh exclusive attempt directory; older files
// remain immutable and referenced by the full retained history. Wall time,
// failed-attempt bytes and total filesystem usage need an external supervisor.
SourceAdaptiveStageRun solveSourceAdaptiveSpinStages(
    const SourceAdaptiveSpinRequest&,const SourceAdaptiveStageControl&);
}
#endif
