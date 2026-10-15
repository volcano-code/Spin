#pragma once

#if defined(MAKE_IT_STAND_DIAGNOSTIC_TIMING)
#include <chrono>
#include <cstdio>

namespace spin2014 {
class DiagnosticStageTimer {
public:
    DiagnosticStageTimer(const char* phase,int pass)
        :phase_(phase),pass_(pass),start_(std::chrono::steady_clock::now()){}
    DiagnosticStageTimer(const DiagnosticStageTimer&)=delete;
    DiagnosticStageTimer& operator=(const DiagnosticStageTimer&)=delete;
    ~DiagnosticStageTimer() noexcept {
        const auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now()-start_).count();
        std::fprintf(stderr,"SPIN2014_DIAG_TIMING phase=%s pass=%d elapsed_ns=%lld\n",
            phase_,pass_,static_cast<long long>(elapsed));
    }
private:
    const char* phase_;
    int pass_;
    std::chrono::steady_clock::time_point start_;
};
}

#define SPIN2014_DIAG_JOIN_INNER(a,b) a##b
#define SPIN2014_DIAG_JOIN(a,b) SPIN2014_DIAG_JOIN_INNER(a,b)
#define SPIN2014_DIAG_TIMER(phase,pass) \
    ::spin2014::DiagnosticStageTimer SPIN2014_DIAG_JOIN(spin2014_diag_timer_,__LINE__){phase,pass}
#else
#define SPIN2014_DIAG_TIMER(phase,pass) do {} while(false)
#endif
