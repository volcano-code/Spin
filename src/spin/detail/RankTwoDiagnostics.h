#pragma once

namespace spin::detail {
// Observation only. No clocks, result fields, branches or numerical settings
// are supplied to the algorithm. Each thread defaults to a null sink.
using RankTwoDiagnosticSink=void(*)(const char*,void*) noexcept;
struct RankTwoDiagnosticTarget {RankTwoDiagnosticSink sink=nullptr;void* context=nullptr;};
inline thread_local RankTwoDiagnosticTarget rankTwoDiagnosticTarget;
inline void rankTwoDiagnostic(const char* phase) noexcept {
    const auto target=rankTwoDiagnosticTarget;
    if(target.sink)target.sink(phase,target.context);
}
class ScopedRankTwoDiagnostics {
public:
    ScopedRankTwoDiagnostics(RankTwoDiagnosticSink sink,void* context) noexcept
        :previous_(rankTwoDiagnosticTarget){rankTwoDiagnosticTarget={sink,context};}
    ~ScopedRankTwoDiagnostics(){rankTwoDiagnosticTarget=previous_;}
    ScopedRankTwoDiagnostics(const ScopedRankTwoDiagnostics&)=delete;
    ScopedRankTwoDiagnostics& operator=(const ScopedRankTwoDiagnostics&)=delete;
private:
    RankTwoDiagnosticTarget previous_;
};
} // namespace spin::detail
