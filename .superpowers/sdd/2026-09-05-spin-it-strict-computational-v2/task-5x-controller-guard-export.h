#pragma once
// Controller-only serialization of actual callbacks/geometry/KKT, not a solver.
#include "spin2014/SourceCageGuardedSlqp.h"
#include "spin2014/SlqpAuditTrace.h"
#include "task-5w-controller-export-helpers.h"
#include <stdexcept>

void dumpGuarded(const spin2014::SourceCagePaperFunctional& f,
                 const spin2014::SourceCageGuardedSlqpAudit& audit) {
    using namespace spin2014;
    const auto& a=audit.facts();const auto& s=a.source;
    std::cout<<"{\"hash\":"<<std::quoted(a.scientificSha256)
        <<",\"final_status\":"<<int(a.finalStatus)<<",\"outcome\":"<<int(s.outcome)
        <<",\"numerical_available\":"<<s.numericalResultAvailable
        <<",\"final_evaluation_available\":"<<s.finalEvaluationAvailable
        <<",\"status\":"<<std::quoted(s.numericalResultAvailable?slqpStatusToken(s.numerical.status):"UNAVAILABLE")
        <<",\"checker_calls\":"<<a.checkerInvocationCount<<",\"final_rechecks\":"<<a.finalRecheckCount
        <<",\"callback_count\":"<<s.callbackCount<<",\"source_calls\":"<<s.sourceEvaluationCount
        <<",\"initial\":";arrayJson(s.initialState);
    std::cout<<",\"final\":";arrayJson(s.numerical.finalState);
    std::cout<<",\"controls\":";writeSlqpAuditControlsJson(std::cout,a.settings.source.solver);
    std::cout<<",\"calls\":[";
    for(std::size_t i=0;i<s.recordCount;++i) {
        if(i)std::cout<<',';const auto& r=s.records[i];
        const Eigen::VectorXd x=s.callbackStates.col(i);
        if(r.status==SourceCageSlqpEvaluationStatus::Evaluated && f.evaluate(x).outputSha256!=r.outputSha256.data())
            throw std::runtime_error("callback does not bind actual material evaluation");
        std::cout<<"{\"sequence\":"<<r.sequence<<",\"purpose\":"<<int(r.purpose)
            <<",\"status\":"<<int(r.status)<<",\"invoked\":"<<r.sourceInvoked
            <<",\"state_hash\":"<<std::quoted(r.stateSha256.data())
            <<",\"output_hash\":"<<std::quoted(r.outputSha256.data())
            <<",\"reason\":"<<std::quoted(r.reason.data())<<",\"state\":";
        arrayJson(x);std::cout<<'}';
    }
    std::cout<<"],\"guards\":[";
    for(std::size_t i=0;i<a.attempts.size();++i) {
        if(i)std::cout<<',';const auto& g=a.attempts[i];
        if(g.sourceOutputSha256[0] && f.evaluateValues(s.callbackStates.col(g.callbackRecord)).outputSha256!=g.sourceOutputSha256.data())
            throw std::runtime_error("geometry does not bind actual material state");
        std::cout<<"{\"callback\":"<<g.callbackRecord<<",\"first_stage\":"<<g.firstStage
            <<",\"stage_count\":"<<g.stageCount<<",\"invoked\":"<<g.checkerInvoked
            <<",\"completed\":"<<g.checkCompleted<<",\"outcome\":"<<int(g.outcome)<<'}';
    }
    std::cout<<"],\"stages\":[";
    for(std::size_t i=0;i<a.stages.size();++i) {
        if(i)std::cout<<',';const auto& g=a.stages[i];
        std::cout<<"{\"stage\":"<<int(g.stage)<<",\"outcome\":"<<int(g.outcome)
            <<",\"first\":"<<g.first<<",\"second\":"<<g.second
            <<",\"reason\":"<<std::quoted(g.reason.data())
            <<",\"distance\":"<<std::quoted(g.squaredDistanceExact)<<'}';
    }
    std::cout<<"],\"trace\":[";
    for(std::size_t i=0;i<s.numerical.trace.size();++i) {
        if(i)std::cout<<',';const auto& t=s.numerical.trace[i];const auto e=f.evaluate(t.state);
        if(e.totalObjective!=t.objective || !e.jointGradient.isApprox(t.gradient,0.) ||
           !e.constraints.isApprox(t.constraints,0.) || !e.jointConstraintJacobian.isApprox(t.constraintJacobian,0.))
            throw std::runtime_error("KKT trace does not bind actual material derivatives");
        std::cout<<"{\"iteration\":"<<t.iteration<<",\"objective\":"<<t.objective
            <<",\"merit\":"<<t.meritValue<<",\"accepted_step\":"<<t.acceptedStep
            <<",\"line_search_evaluations\":"<<t.lineSearchEvaluations;
        writeSlqpAuditJson(std::cout,t);std::cout<<'}';
    }
    std::cout<<']';
    if(s.finalEvaluationAvailable) {
        const auto e=f.evaluate(s.numerical.finalState);
        if(e.outputSha256!=s.finalEvaluation.outputSha256)
            throw std::runtime_error("final evaluation differs from actual endpoint");
        std::cout<<",\"material\":";dumpSource(f,s.numerical.finalState);
    }
    std::cout<<'}';
}
