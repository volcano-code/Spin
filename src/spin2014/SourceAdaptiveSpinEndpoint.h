#pragma once
#include "spin2014/SourceAdaptiveSpinSolver.h"
#include "spin2014/SourceAdaptiveCheckpoint.h"

namespace spin2014::detail {
void finishSourceAdaptiveResult(SourceAdaptiveSpinResult&);
// Private trusted handoff, populated only within the actual solve invocation.
// Never reconstructed from caller-supplied leaves, hashes, status or cuboids.
struct SourceAdaptiveSpinEndpoint {
    std::optional<SourceAdaptiveDomain> domain;
    Eigen::VectorXd binary;
};
SourceAdaptiveSpinResult solveSourceAdaptiveSpinEndpoint(
    const SourceAdaptiveSpinRequest&,SourceAdaptiveSpinEndpoint*,
    const SourceAdaptiveStageControl* control=nullptr,
    std::shared_ptr<const SourceAdaptiveCheckpoint>* checkpoint=nullptr);
}
