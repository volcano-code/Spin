// Development characterization only; not production, not an author-case benchmark.
#include "spin/CanonicalFullRankHollowing.h"
#include <nlohmann/json.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
using Json=nlohmann::json;
namespace {
Json number(double x) { return std::isfinite(x)?Json(x):Json("NONFINITE"); }
template<class T> Json optional(const std::optional<T>& x) {return x?Json(*x):Json(nullptr);}
Json optional(const std::optional<double>& x) {return x?number(*x):Json(nullptr);}
template<class T> Json vector(const T& x) {
    Json result=Json::array(); for(Eigen::Index i=0;i<x.size();++i)result.push_back(number(x(i)));
    return result;
}
void require(bool b,const char* diagnostic) {if(!b)throw std::runtime_error(diagnostic);}
double arg(const char* value) {
    std::size_t end=0; const std::string text(value); const double x=std::stod(text,&end);
    require(end==text.size() && std::isfinite(x),"invalid numeric argument"); return x;
}
int budgetArg(const char* value) {
    const std::string text(value);
    require(text=="20" || text=="40" || text=="60","budget outside preregistered matrix");
    return std::stoi(text);
}
spin::TriangleMesh box(double radius,double z) {
    spin::TriangleMesh m; m.vertices.resize(8,3);m.faces.resize(12,3);
    const double a=-radius,b=radius;
    m.vertices<<a,a,a,b,a,a,b,b,a,a,b,a,a,a,b,b,a,b,b,b,b,a,b,b;
    m.vertices.col(2).array()+=z;
    m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
    return m;
}
Json event(const spin::SpinSolverEvent& e) {
    return {{"type",e.type},{"restart",e.restart},{"iteration",e.iteration},
        {"evaluation_index",e.evaluationIndex},{"restart_evaluation_index",e.restartEvaluationIndex},
        {"step",optional(e.step)},{"outcome",e.outcome},{"reason",e.reason},{"detail",e.detail},
        {"quadric_coefficients",optional(e.quadricCoefficients)},{"quadric_hash",e.quadricHash},
        {"normalized_residual",optional(e.normalizedEquationResidual)},
        {"physical_residual",optional(e.physicalRelativeResidual)},{"objective",optional(e.objective)},
        {"mu",optional(e.multiplierMu)},{"geometry_hash",optional(e.materialGeometryHash)},
        {"topology_hash",optional(e.materialTopologyHash)},
        {"profile",spin::spinEquationProfileToken(e.equationProfile)},
        {"raw_residual",optional(e.rawEquationResidual)}};
}
Json trace(const spin::SpinIterationTrace& t) {
    return {{"restart",t.restart},{"iteration",t.iteration},{"evaluations",t.functionEvaluations},
        {"normalized_residual",number(t.normalizedEquationResidual)},
        {"physical_residual",number(t.physicalRelativeResidual)},{"objective",number(t.objective)},
        {"mu",number(t.multiplierMu)},{"accepted_geodesic_step",number(t.acceptedGeodesicStep)},
        {"quadric_hash",t.quadricHash},{"geometry_hash",t.materialGeometryHash},
        {"topology_hash",t.materialTopologyHash},{"profile",spin::spinEquationProfileToken(t.equationProfile)},
        {"raw_residual",optional(t.rawEquationResidual)},{"quadric_coefficients",optional(t.quadricCoefficients)}};
}
}
int main(int argc,char** argv) {
    try {
        const auto started=std::chrono::steady_clock::now();
        require(argc==4,"usage: canonical_budget_probe <z:0|2> <weight:0.1|1|10> <budget:20|40|60>");
        const double z=arg(argv[1]),wc=arg(argv[2]); const int budget=budgetArg(argv[3]);
        require((z==0 || z==2) && (wc==.1 || wc==1 || wc==10),"case outside preregistered matrix");
        spin::SpinDesignRequest r;
        r.designDomain=box(.8,z);r.outerSurface=box(1,z);r.minimumWallThickness=.15;r.materialDensity=1;
        r.mode=spin::SpinMode::SpinningTop;r.centerOfMassWeight=wc;
        r.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.outerSurface.vertices,r.outerSurface.faces,1))-
            physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.designDomain.vertices,r.designDomain.faces,1));
        r.settings.equationProfile=spin::SpinEquationProfile::DualMonomialRawV1;
        r.settings.maximumRestarts=1;r.settings.maximumFunctionEvaluationsPerRestart=budget;
        spin::QuadricMomentEvaluatorSettings es;es.triangulation.icosphereSubdivisions=1;es.triangulation.truncationRadius=8;
        spin::QuadricCoefficients10::Vector p;p<<-.31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03;
        const auto seed=spin::QuadricCoefficients10::fromUnnormalized(p);const Eigen::Vector3d c(0,0,z);
        p(0)=seed.constant()-seed.linearVector().dot(c)+c.dot(seed.quadraticMatrix()*c);
        p.segment<3>(1)=seed.linearVector()-2*seed.quadraticMatrix()*c;p.tail<6>()=seed.values().tail<6>();
        const auto initial=spin::QuadricCoefficients10::fromUnnormalized(p);
        const auto out=spin::solveCanonicalFullRankQuadricHollowing(r,es,{initial});
        const auto& n=out.numerical;
        require(n.status!=spin::SpinSolverStatus::Infeasible,"search labeled mathematically infeasible");
        require(n.functionEvaluations<=budget && n.evaluationTimings.size()==static_cast<std::size_t>(n.functionEvaluations),"evaluation count mismatch");
        if(out.finalGeometryVerified)require(out.independentlyRecomputedRawResidual && *out.independentlyRecomputedRawResidual<=1e-12,"false final certification");
        if(n.status!=spin::SpinSolverStatus::Converged)require(!out.finalGeometryVerified && !out.finalEvaluation,"failed search exposed final material");
        const double solverOther=out.totalWallSeconds-n.evaluatorWallSeconds-out.finalVerificationSeconds;
        require(solverOther>=0,"nested timing double count");
        Json j={{"schema","m1-budget-canonical-v1"},{"z",z},{"weight",wc},{"evaluation_budget",budget},
            {"numerical_status",spin::spinSolverStatusToken(n.status)},{"numerical_diagnostic",n.diagnostic},
            {"final_verified",out.finalGeometryVerified},{"final_diagnostic",out.verificationDiagnostic},
            {"raw_residual",optional(n.rawEquationResidual)},{"independent_raw_residual",optional(out.independentlyRecomputedRawResidual)},
            {"physical_residual",number(n.physicalRelativeResidual)},
            {"normalized_residual",number(n.normalizedEquationResidual)},{"mu",number(n.multiplierMu)},
            {"multipliers",vector(n.multipliers)},{"initial_quadric",vector(initial.values())},
            {"quadric",n.quadric?vector(n.quadric->values()):Json(nullptr)},
            {"input_hash",out.inputHash},{"output_hash",out.outputHash},{"trace_hash",n.traceHash},{"event_hash",n.eventHash},
            {"final_geometry_hash",out.finalEvaluation?Json(out.finalEvaluation->geometryHash):Json(nullptr)},
            {"final_moments",out.independentTotalMoments?vector(out.independentTotalMoments->values()):Json(nullptr)},
            {"evaluations",n.functionEvaluations},{"restarts",n.restartCount},
            {"restart_evaluations",n.restartFunctionEvaluations},{"event_stream_version",n.eventStreamVersion},
            {"timing",{{"solve_wrapper_seconds",out.totalWallSeconds},
                {"coupled_evaluator_seconds",n.evaluatorWallSeconds},{"final_audit_seconds",out.finalVerificationSeconds},
                {"solver_other_seconds",solverOther},{"evaluator_reported_residual_seconds",n.residualEvaluationSeconds},
                {"evaluator_reported_derivative_seconds",n.derivativeEvaluationSeconds}}},
            {"events",Json::array()},{"trace",Json::array()},{"evaluation_timings",Json::array()}};
        for(const auto& e:n.events)j["events"].push_back(event(e));
        for(const auto& t:n.trace)j["trace"].push_back(trace(t));
        for(const auto& t:n.evaluationTimings)j["evaluation_timings"].push_back({
            {"evaluation_index",t.evaluationIndex},{"restart",t.restart},{"restart_evaluation_index",t.restartEvaluationIndex},
            {"returned",t.evaluatorReturned},{"wall_seconds",t.wallSeconds},
            {"reported_residual_seconds",optional(t.reportedResidualSeconds)},
            {"reported_derivative_seconds",optional(t.reportedDerivativeSeconds)}});
        j["timing"]["entry_through_record_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        std::cout<<"M1_RESULT "<<j.dump()<<'\n';return 0;
    } catch(const std::exception& e) {
        std::cerr<<"M1_HARNESS_ERROR "<<e.what()<<'\n';return 2;
    }
}
