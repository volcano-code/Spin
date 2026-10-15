// Replay one recorded terminal QP, not the nonlinear optimization. All records
// and L-BFGS history counts are verified before invoking the current QP kernel.
#include "spin2014/SlqpIterationStore.h"
#include "spin2014/SlqpAuditTrace.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>

void require(bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}
int main(int argc,char** argv) {
    try {
        require(argc==6 || argc==7,
                "usage: native_qp_replay TRACE BOUNDED HISTORY_CAPACITY FACES ITERATIONS [SELECTED_ITERATION]");
        const int bounded=std::stoi(argv[2]),capacity=std::stoi(argv[3]);
        const int selectedIteration=argc==7?std::stoi(argv[6]):-1;
        spin2014::TangentQpSettings settings;
        settings.maximumFaces=std::stoi(argv[4]);settings.maximumIterations=std::stoi(argv[5]);
        require(bounded>0 && capacity>0 && selectedIteration>=-1,"invalid replay controls");
        spin2014::LbfgsHistory history;history.capacity=capacity;
        spin2014::SlqpIterationReader reader(argv[1]);
        spin2014::SlqpIteration current,previous,selected;
        spin2014::LbfgsHistory selectedHistory;
        bool selectedFound=selectedIteration<0;
        std::size_t count=0;
        while(reader.next(current)) {
            require(!count || !previous.finalState,"records after terminal state");
            require(current.iteration==static_cast<int>(count),"iteration sequence mismatch");
            if(count) {
                if(previous.restoration.accepted)history=spin2014::LbfgsHistory{capacity};
                else if(previous.curvaturePairAccepted || previous.curvaturePairRejected) {
                    require(previous.acceptedStep>0 && previous.kkt.valid,"invalid curvature source");
                    const Eigen::VectorXd s=current.state-previous.state;
                    // Preserve native arithmetic order and the previous witness
                    // multipliers for both endpoints of the secant pair.
                    const Eigen::VectorXd y=(current.gradient+current.constraintJacobian.transpose()*previous.kkt.lambda)-
                        (previous.gradient+previous.constraintJacobian.transpose()*previous.kkt.lambda);
                    const bool accepted=history.update(s,y);
                    require(accepted==previous.curvaturePairAccepted && accepted!=previous.curvaturePairRejected,
                            "curvature acceptance mismatch");
                }
                require(history.s.size()==static_cast<std::size_t>(previous.historyAfter),"history after mismatch");
            }
            require(history.s.size()==static_cast<std::size_t>(current.historyBefore),"history before mismatch");
            const std::vector<unsigned char> empty;
            history.resetForMembership(count?previous.activeMembership:empty,current.activeMembership);
            require(history.s.size()==static_cast<std::size_t>(current.lbfgsHistorySize),"metric history mismatch");
            if(current.iteration==selectedIteration) {
                selected=current;
                selectedHistory=history;
                selectedFound=true;
            }
            previous=std::move(current);++count;
        }
        require(count>0 && previous.finalState,"missing terminal state");
        require(selectedFound,"selected iteration is not present in trace");
        if(selectedIteration>=0) previous=std::move(selected);
        const auto& replayHistory=selectedIteration>=0?selectedHistory:history;
        require(previous.state.size()==bounded+1 && previous.metricTangent.attempted &&
                !previous.metricTangent.trace.empty(),"selected iteration has no recorded metric QP");
        if(selectedIteration<0)
            require(!previous.metricTangent.checked,"terminal QP is unexpectedly already checked");
        else
            require(previous.metricTangent.checked,"selected iteration is not a completed checked QP");
        std::size_t recordedPeakFaces=0;
        for(const auto& step:previous.metricTangent.trace)
            if(step.workingSet.size()>recordedPeakFaces)recordedPeakFaces=step.workingSet.size();
        const auto start=std::chrono::steady_clock::now();
        const auto result=spin2014::solvePrimalBoxQp(previous.state,previous.gradient,
            previous.constraintJacobian,previous.constraints,bounded,replayHistory,settings);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        require(!result.trace.empty(),"missing replay trace");
        std::cout<<std::setprecision(17)<<"{\"schema\":\"spin2014-native-qp-replay-v1\",\"records\":"<<count
                 <<",\"selected_iteration\":"<<previous.iteration
                 <<",\"recorded_peak_working_faces\":"<<recordedPeakFaces
                 <<",\"history_pairs\":"<<replayHistory.s.size()<<",\"maximum_faces\":"<<settings.maximumFaces
                 <<",\"maximum_iterations\":"<<settings.maximumIterations<<",\"seconds\":"<<seconds
                 <<",\"last_action\":"<<std::quoted(result.trace.back().action)
                 <<",\"working_faces\":"<<result.workingSet.size()<<",\"steps\":"<<result.trace.size()
                 <<",\"original_terminal_action\":"<<std::quoted(previous.metricTangent.trace.back().action)
                 <<",\"original_working_set_equal\":"<<(result.workingSet==previous.metricTangent.workingSet?"true":"false")
                 <<",\"history\":[";
        auto vector=[](const Eigen::VectorXd& v){std::cout<<'[';for(Eigen::Index i=0;i<v.size();++i){if(i)std::cout<<',';std::cout<<v[i];}std::cout<<']';};
        for(std::size_t i=0;i<replayHistory.s.size();++i){if(i)std::cout<<',';std::cout<<"{\"s\":";vector(replayHistory.s[i]);std::cout<<",\"y\":";vector(replayHistory.y[i]);std::cout<<'}';}
        std::cout<<"],\"terminal\":{\"iteration\":"<<previous.iteration;
        previous.metricTangent=result;
        spin2014::writeSlqpAuditJson(std::cout,previous);
        std::cout<<"},\"scope\":\"FIXED_QP_REPLAY_NOT_NONLINEAR_OR_BINARY_VALIDATION\"}\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
