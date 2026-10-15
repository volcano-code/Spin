#pragma once
// Test harness only. Observations never change geometry/solver decisions.
#include <algorithm>
#include <cstdint>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <ctime>
namespace rank_two_test {
struct ClockPair {std::uint64_t wall,cpu;};
inline ClockPair clocks(){
    timespec wall{},cpu{};
    if(clock_gettime(CLOCK_MONOTONIC,&wall)||clock_gettime(CLOCK_THREAD_CPUTIME_ID,&cpu))
        throw std::runtime_error("phase profiler clock unavailable");
    return {std::uint64_t(wall.tv_sec)*1000000000+std::uint64_t(wall.tv_nsec),
        std::uint64_t(cpu.tv_sec)*1000000000+std::uint64_t(cpu.tv_nsec)};
}
struct PhaseTotal {std::uint64_t spans=0,wall=0,cpu=0,maxWall=0,maxCpu=0;};
class PhaseAccumulator {
public:
    explicit PhaseAccumulator(ClockPair begin):begin_(begin),last_(begin){}
    void advance(const std::string& next,ClockPair now){
        if(now.wall<last_.wall||now.cpu<last_.cpu)throw std::runtime_error("nonmonotone profile clock");
        if(next.empty()||next.size()>96)throw std::runtime_error("invalid bounded profile label");
        if(!phase_.empty()){
            auto& t=totals_[phase_];const auto w=now.wall-last_.wall,c=now.cpu-last_.cpu;
            ++t.spans;t.wall+=w;t.cpu+=c;t.maxWall=std::max(t.maxWall,w);t.maxCpu=std::max(t.maxCpu,c);
        }
        if(totals_.size()>128)throw std::runtime_error("phase profiler label budget exceeded");
        phase_=next;last_=now;
    }
    const auto& totals()const{return totals_;}
    ClockPair begin()const{return begin_;}
    ClockPair last()const{return last_;}
private:
    ClockPair begin_,last_;std::string phase_;std::map<std::string,PhaseTotal> totals_;
};
inline std::string phaseLabel(const char* raw){
    std::string token(raw);token=token.substr(0,token.find(' '));
    // These completion/fallback markers remain DISTINCT in raw diagnostics.
    // The bounded AFTER-event histogram groups them with the existing end
    // and begin of that same operation; it never claims exclusive function
    // timing or changes scientific work or the original 128-label budget.
    if(token=="native-fiber:raw-resultant-full-product-verified")
        return "native-fiber:raw-resultant-complete";
    if(token=="native-fiber:raw-resultant-full-product-direct-fallback")
        return "native-fiber:raw-resultant-begin";
    if(token=="algebraic-root:resultant-component-product-verified")
        return "algebraic-root:native-factor-product-verified";
    if(token=="native-fiber:component-square-free-certified")
        return "native-fiber:square-free-certified";
    if(token=="native-fiber:component-square-free-original-fallback")
        return "native-fiber:square-free-begin";
    if(token=="regular-strip:memo-hit" || token=="regular-strip:memo-miss")
        return "regular-edge:attempt";
    if(token=="algebraic-root:reload-complete")return "unattributed:after-root-complete";
    if(token=="mixed-intervals:edge-ready")return "unattributed:after-edge-complete";
    if(token.rfind("registry:subdivision:",0)==0){
        const auto tail=token.substr(token.rfind(':')+1);
        return tail=="begin"?"registry:subdivision-begin":"unattributed:after-subdivision-"+tail;
    }
    for(const auto* prefix:{"algebraic-root:frozen-degree-","algebraic-root:isolate-factor-degree-",
            "algebraic-root:factor-root-set-","algebraic-root:expression-degree-",
            "native-fiber:factor-root-set-"}){
        if(token.rfind(prefix,0)==0)return std::string(prefix).substr(0,std::string(prefix).size()-1);
    }
    const auto first=token.find(':');
    if(first!=std::string::npos)token=token.substr(0,token.find(':',first+1));
    return token;
}
class PhaseProfile {
public:
    explicit PhaseProfile(std::ostream& out):out_(out),data_(clocks()){
        data_.advance("harness:begin",data_.begin());
    }
    void stage(const char* name){
        data_.advance("observer:accounting",clocks());snapshot("stage");
        // These spans are harness validation/replays AFTER native construction,
        // not exclusive algorithm timing. Whole-output validation and the two
        // diagnostic replays use the existing import/validation category,
        // preserving the original cap even after all four owners complete.
        // No other stage identity or raw source diagnostic is changed.
        const auto stage=std::string(name);
        const char* label=(stage=="harness:validate-coverage-output"
            ||stage=="harness:source-projection-replay"
            ||stage=="harness:source-subdivision-replay")
            ?"harness:import-and-validate":name;
        data_.advance(label,clocks());
    }
    void event(const char* phase) noexcept {
        try {
            data_.advance("observer:logging",clocks());
            out_<<"WHOLE_COVERAGE_DIAGNOSTIC "<<phase<<'\n';
            data_.advance("observer:accounting",clocks());++events_;
            if(events_%4096==0||std::string(phase).rfind("mixed-intervals:round-counts:",0)==0)
                snapshot("event");
            data_.advance(phaseLabel(phase),clocks());
        }catch(...){valid_=false;out_<<"RANK_TWO_PHASE_ERROR observation unavailable\n";}
    }
    void finish(){
        data_.advance("observer:accounting",clocks());snapshot("final");
        if(!valid_)throw std::runtime_error("invalid actual phase observation");
    }
private:
    void snapshot(const char* reason){
        // Copy the counters, then close ONLY the snapshot's open interval.
        // The live accumulator accounts checkpoint writing at the next event.
        auto saved=data_;saved.advance("snapshot:boundary",clocks());++snapshots_;
        out_<<"RANK_TWO_PHASE_SNAPSHOT_BEGIN 1 "<<snapshots_<<' '<<events_<<' '
            <<saved.begin().wall<<' '<<saved.begin().cpu<<' '<<saved.last().wall<<' '
            <<saved.last().cpu<<' '<<reason<<'\n';
        for(const auto& item:saved.totals()){
            const auto& t=item.second;
            out_<<"RANK_TWO_PHASE_BUCKET "<<snapshots_<<' '<<item.first<<' '<<t.spans<<' '
                <<t.wall<<' '<<t.cpu<<' '<<t.maxWall<<' '<<t.maxCpu<<'\n';
        }
        out_<<"RANK_TWO_PHASE_SNAPSHOT_END "<<snapshots_<<'\n';out_.flush();
    }
    std::ostream& out_;PhaseAccumulator data_;std::uint64_t events_=0,snapshots_=0;bool valid_=true;
};
}
