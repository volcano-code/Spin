#include "spin2014/SlqpIterationStore.h"
#include "spin2014/SlqpAuditTrace.h"
#include <filesystem>
#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include <csignal>
using namespace spin2014;
void check(bool b,const char* m){if(!b)throw std::runtime_error(m);}
template<class F> void fails(F f){bool caught=false;try{f();}catch(const SlqpStorageError&){caught=true;}check(caught,"storage error required");}
int main(){try {
    char tmp[]="/tmp/slqp-store-test-XXXXXX";const auto dir=std::string(mkdtemp(tmp));
    // Synthetic native problem: a maximum-iteration accepted step needs a new final record.
    EqualityBoxProblem p{2,2,[](const Eigen::VectorXd& x){AdaptiveFillEvaluation e;
        e.objective=.5*x.squaredNorm();e.gradient=x;e.constraints.resize(0);
        e.constraintJacobian.resize(0,2);return e;}};
    SlqpLbfgsSettings s;s.algorithm=SlqpAlgorithm::AuditedKktRestorationV1;s.maximumIterations=1;
    const Eigen::VectorXd x=Eigen::VectorXd::Constant(2,.7);
    auto memory=solveSlqpLbfgs(p,x,s);s.iterationStorePath=dir+"/solve";
    auto disk=solveSlqpLbfgs(p,x,s);
    check(disk.trace.empty(),"disk history must not accumulate in vector");
    check(memory.traceHash==disk.traceHash,"scientific hash mismatch");
    check(disk.iterationStore->retainedRecords()==0&&disk.iterationStore->peakRetainedRecords()==1,"retention bound");
    std::size_t i=0;forEachSlqpIteration(disk,[&](const auto& t){check(encodeSlqpIteration(t)==encodeSlqpIteration(memory.trace.at(i++)),"complete record mismatch");});
    check(i==memory.trace.size()&&memory.trace.back().finalState,"final state lost");
    fails([&]{SlqpIterationStore collision(s.iterationStorePath,1000000);});
    auto impossible=p;impossible.evaluate=[](const Eigen::VectorXd& x){AdaptiveFillEvaluation e;
        e.objective=0;e.gradient=Eigen::VectorXd::Zero(2);e.constraints=Eigen::VectorXd::Ones(1);
        e.constraintJacobian=Eigen::MatrixXd::Zero(1,2);return e;};
    s.iterationStorePath.clear();const auto failureMemory=solveSlqpLbfgs(impossible,x,s);
    s.iterationStorePath=dir+"/failure";const auto failureDisk=solveSlqpLbfgs(impossible,x,s);
    check(failureDisk.status!=SlqpStatus::Converged&&failureDisk.traceHash==failureMemory.traceHash,"failure result mismatch");
    // Nondefault nested fields, NaN and infinity must survive the binary codec.
    SlqpIteration t=memory.trace.back();t.metricTangent.attempted=true;t.metricTangent.checked=true;
    t.metricTangent.status="quoted\"\nstatus";t.metricTangent.direction=x;t.metricTangent.multipliers=x;
    t.metricTangent.workingSet={1,2};t.metricTangent.trace.push_back({{3,4},"remove",7,.8,.9});
    t.metricTangentSelected=t.metricTangentAccepted=t.metricTangentFallbackAttempted=true;t.metricTangentTrialCount=9;
    t.restoration.attempted=true;t.restoration.accepted=true;t.restoration.direction=x;
    t.restoration.trials.push_back({.25,INFINITY,false,true,"trial"});t.objectiveTrials=t.restoration.trials;
    t.restoration.status="status";t.restoration.reason="reason";t.objective=NAN;
    t.historyBefore=2;t.historyAfter=3;t.curvaturePairAccepted=true;t.curvaturePairRejected=true;
    t.historyResetReason="reset";t.usedFirstOrderFallback=true;
    {SlqpIterationStore store(dir+"/fields",10000000);for(int n=0;n<1000;++n)store.append(t);
     check(store.retainedRecords()==1,"online retention");store.seal();}
    SlqpIterationReader reader(dir+"/fields");SlqpIteration decoded;int n=0;
    while(reader.next(decoded)){check(encodeSlqpIteration(decoded)==encodeSlqpIteration(t),"nested field lost");++n;}
    check(n==1000,"count mismatch");
    {SlqpIterationStore partial(dir+"/partial",1000000);partial.append(t);}
    fails([&]{SlqpIterationReader partial(dir+"/partial");});
    fails([&]{SlqpIterationStore limited(dir+"/limited",64);limited.append(t);limited.seal();});
    // Real kernel write failure; no mock, and the incomplete file remains.
    struct rlimit oldLimit;getrlimit(RLIMIT_FSIZE,&oldLimit);auto small=oldLimit;small.rlim_cur=128;
    auto oldSignal=std::signal(SIGXFSZ,SIG_IGN);setrlimit(RLIMIT_FSIZE,&small);
    fails([&]{SlqpIterationStore failed(dir+"/write-failure",1000000);failed.append(t);failed.seal();});
    setrlimit(RLIMIT_FSIZE,&oldLimit);std::signal(SIGXFSZ,oldSignal);
    fails([&]{SlqpIterationReader failed(dir+"/write-failure");});
    std::cout<<"PASS native online history; artifacts preserved at "<<dir<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
