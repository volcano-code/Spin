#include "spin/CanonicalFullRankHollowing.h"
#include "manufacturing/SlicerValidation.h"

namespace spin {
// Defined by compiling the frozen companion sources with entrypoint-only
// renames. CMake binds the exact sources and remaps the single transport call.
FullRankMomentEvaluation canonicalFullRankMaterialImplementation(
    const SpinDesignRequest&, const QuadricCoefficients10&,
    const QuadricMomentEvaluatorSettings&);
SpinMomentEvaluator canonicalFullRankAdapterImplementation(
    const SpinDesignRequest&, const QuadricMomentEvaluatorSettings&);
FullRankHollowingResult canonicalFullRankSolverImplementation(
    const SpinDesignRequest&, const QuadricMomentEvaluatorSettings&,
    const std::vector<QuadricCoefficients10>&, const SpinNewtonEventSink&);
FullRankHollowingResult canonicalFullRankSolverImplementation(
    const SpinDesignRequest&, const MassMomentProblem&, const std::string&,
    const QuadricMomentEvaluatorSettings&, const std::vector<QuadricCoefficients10>&,
    const SpinNewtonEventSink&);

namespace {
std::string bindCanonicalPolicy(const std::string& role,const std::string& value) {
    return manufacturing::computeContentSha256(
        "CanonicalMaterialOrdering-v1-permutation-only;"+role+';'+value);
}
}
FullRankMomentEvaluation evaluateCanonicalFullRankMaterialMoments(
    const SpinDesignRequest& r,const QuadricCoefficients10& q,
    const QuadricMomentEvaluatorSettings& settings) {
    auto result=canonicalFullRankMaterialImplementation(r,q,settings);
    result.scientificHash=bindCanonicalPolicy("evaluation",result.scientificHash);
    return result;
}
SpinMomentEvaluator makeCanonicalFullRankMomentEvaluator(
    const SpinDesignRequest& r,const QuadricMomentEvaluatorSettings& settings) {
    return canonicalFullRankAdapterImplementation(r,settings);
}
FullRankHollowingResult solveCanonicalFullRankQuadricHollowing(
    const SpinDesignRequest& r,const QuadricMomentEvaluatorSettings& settings,
    const std::vector<QuadricCoefficients10>& guesses, const SpinNewtonEventSink& eventSink) {
    auto result=canonicalFullRankSolverImplementation(r,settings,guesses,eventSink);
    result.inputHash=bindCanonicalPolicy("input",result.inputHash);
    result.outputHash=bindCanonicalPolicy("output",result.inputHash+';'+result.outputHash);
    return result;
}

FullRankHollowingResult solveCanonicalFullRankQuadricHollowing(
    const SpinDesignRequest& r,const MassMomentProblem& problem,
    const std::string& problemIdentity,const QuadricMomentEvaluatorSettings& settings,
    const std::vector<QuadricCoefficients10>& guesses, const SpinNewtonEventSink& eventSink) {
    auto result=canonicalFullRankSolverImplementation(r,problem,problemIdentity,settings,guesses,eventSink);
    result.inputHash=bindCanonicalPolicy("input",result.inputHash);
    result.outputHash=bindCanonicalPolicy("output",result.inputHash+';'+result.outputHash);
    return result;
}
}
