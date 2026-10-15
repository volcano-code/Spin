// Diagnostic only. Include the byte-unchanged frozen solver, not a new solver.
#include MIS_FROZEN_SOLVER_CPP
#include "spin/CanonicalFullRankHollowing.h"
#include <Eigen/SVD>
#include <nlohmann/json.hpp>
#include <iostream>

using DiagJson = nlohmann::json;
using DiagVector = spin::QuadricCoefficients10::Vector;
using DiagMatrix = Eigen::Matrix<double,10,10>;

namespace {
void diagRequire(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
template<class T> DiagJson diagVector(const T& v) {
    DiagJson a=DiagJson::array();
    for(Eigen::Index i=0;i<v.size();++i) {
        diagRequire(std::isfinite(v(i)),"nonfinite diagnostic vector");a.push_back(v(i));
    }
    return a;
}
template<class T> DiagJson diagMatrix(const T& m) {
    DiagJson a=DiagJson::array();
    for(Eigen::Index i=0;i<m.rows();++i)a.push_back(diagVector(m.row(i)));
    return a;
}
double diagArg(const char* text) {
    std::size_t consumed=0;std::string str(text);double value=std::stod(str,&consumed);
    diagRequire(consumed==str.size() && std::isfinite(value),"invalid argument");return value;
}
spin::TriangleMesh diagBox(double radius,double z) {
    spin::TriangleMesh m;m.vertices.resize(8,3);m.faces.resize(12,3);const double a=-radius,b=radius;
    m.vertices<<a,a,a,b,a,a,b,b,a,a,b,a,a,a,b,b,a,b,b,b,b,a,b,b;m.vertices.col(2).array()+=z;
    m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
    return m;
}
struct DiagLinearization {
    DiagMatrix hessian;
    Eigen::MatrixXd matrix;
    Eigen::VectorXd rhs,solution,singular;
    DiagVector direction;
    Eigen::Index rank;
};
DiagLinearization diagLinearization(const NewtonState& state,const spin::QuadricCoefficients10& q) {
    // Literal Eq14 assembly from frozen QuadricHollowingSolver.cpp:708-797.
    DiagLinearization out;const Eigen::Index constraints=state.problem.constraints.size();
    out.hessian=state.problem.objectiveHessian;
    for(Eigen::Index i=0;i<constraints;++i)out.hessian-=state.lambda(i)*state.problem.constraintHessians[static_cast<std::size_t>(i)];
    out.matrix=Eigen::MatrixXd::Zero(11+constraints,11+constraints);
    out.rhs=Eigen::VectorXd::Zero(11+constraints);
    out.matrix.topLeftCorner(10,10)=out.hessian*state.map.momentJacobian+state.mu*DiagMatrix::Identity();
    out.matrix.block(0,10,10,constraints)=-state.problem.constraintJacobian.transpose();
    out.matrix.block(10,0,constraints,10)=state.problem.constraintJacobian*state.map.momentJacobian;
    out.matrix.block(0,10+constraints,10,1)=state.solverCoordinates;
    out.matrix.block(10+constraints,0,1,10)=2.0*state.solverCoordinates.transpose();
    out.rhs.head(10)=-state.optimalityResidual;out.rhs.segment(10,constraints)=-state.problem.constraints;
    Eigen::FullPivLU<Eigen::MatrixXd> lu(out.matrix);out.rank=lu.rank();out.solution=lu.solve(out.rhs);
    out.direction=spin::DualMonomialCoordinates10::fromGeometry(q).projectToTangent(out.solution.head(10));
    out.singular=Eigen::JacobiSVD<Eigen::MatrixXd>(out.matrix).singularValues();
    return out;
}
DiagJson diagPoint(const spin::QuadricCoefficients10& q,const spin::SpinMomentProblem& problem,
                   const spin::SpinMomentEvaluator& evaluator) {
    try {
        auto map=evaluator(q);auto state=evaluateState(problem,q,map,spin::SpinEquationProfile::DualMonomialRawV1);
        return {{"status","OK"},{"q",diagVector(q.values())},{"q_hash",q.hash()},
            {"moments",diagVector(map.rawMoments)},{"raw_residual",state.rawResidual},
            {"physical_residual",state.physicalResidual},{"phi",.5*state.rawResidual*state.rawResidual},
            {"lambda",diagVector(state.lambda)},{"mu",state.mu},{"geometry_hash",map.geometryHash},
            {"connectivity_hash",map.topologyHash}};
    }catch(const std::exception& error) {
        return {{"status","EVALUATION_FAILED"},{"q",diagVector(q.values())},{"q_hash",q.hash()},{"error",error.what()}};
    }
}
DiagJson diagState(std::size_t index,const std::vector<spin::QuadricCoefficients10>& seen,
                  const std::vector<std::optional<spin::SpinMomentMapEvaluation>>& maps,
                  const spin::SpinMomentProblem& problem,const spin::SpinMomentEvaluator& evaluator) {
    const auto& q=seen.at(index);const auto& map=maps.at(index).value();
    const auto s=spin::DualMonomialCoordinates10::fromGeometry(q);
    const auto state=evaluateState(problem,q,map,spin::SpinEquationProfile::DualMonomialRawV1);
    const auto linear=diagLinearization(state,q);const auto norm=linear.direction.norm();
    double qError=0;
    if(index+1<seen.size() && norm>1e-15) {
        const auto next=s.geodesicUpdate(linear.direction,std::min(norm,.5)).geometry();
        qError=(next.values()-seen[index+1].values()).cwiseAbs().maxCoeff();
    }
    Eigen::VectorXd orthogonality(state.lambda.size()+1);
    orthogonality.head(state.lambda.size())=state.problem.constraintJacobian*state.optimalityResidual;
    orthogonality(state.lambda.size())=state.solverCoordinates.dot(state.optimalityResidual);
    DiagJson row={{"evaluation_index",index+1},{"q",diagVector(q.values())},{"q_hash",q.hash()},
        {"raw_residual",state.rawResidual},{"physical_residual",state.physicalResidual},
        {"phi",.5*state.rawResidual*state.rawResidual},{"moments",diagVector(map.rawMoments)},
        {"native_jacobian",diagMatrix(map.momentJacobian)},{"dual_jacobian",diagMatrix(state.map.momentJacobian)},
        {"mu",state.mu},{"lambda",diagVector(state.lambda)},
        {"optimality_residual",diagVector(state.optimalityResidual)},
        {"constraints",diagVector(state.problem.constraints)},
        {"projection_orthogonality",diagVector(orthogonality)},
        {"system_matrix",diagMatrix(linear.matrix)},{"rhs",diagVector(linear.rhs)},
        {"linear_solution",diagVector(linear.solution)},{"singular_values",diagVector(linear.singular)},
        {"lu_rank",linear.rank},{"linear_relative_residual",(linear.matrix*linear.solution-linear.rhs).norm()/std::max(1e-300,linear.rhs.norm())},
        {"newton_direction",diagVector(linear.direction)},{"newton_norm",norm},{"next_trial_q_max_error",qError},
        {"geometry_hash",map.geometryHash},{"connectivity_hash",map.topologyHash},{"directions",DiagJson::array()}};
    for(int axis=0;axis<11;++axis) {
        DiagVector d=axis<10?s.projectToTangent(DiagVector::Unit(axis)):linear.direction;
        const double dn=d.norm();DiagJson direction={{"axis",axis},{"unscaled_norm",dn}};
        if(dn<=(axis<10?1e-12:1e-15)) {direction["skipped"]=true;row["directions"].push_back(direction);continue;}
        d/=dn;const DiagVector dm=state.map.momentJacobian*d;
        const DiagVector dg=linear.hessian*dm+state.mu*d;
        const double predicted=state.optimalityResidual.dot(dg)+state.problem.constraints.dot(state.problem.constraintJacobian*dm);
        direction["direction"]=diagVector(d);direction["moment_derivative"]=diagVector(dm);
        direction["predicted_phi_derivative"]=predicted;direction["samples"]=DiagJson::array();
        for(double h:{1e-2,1e-3,1e-4,1e-5,1e-6}) {
            direction["samples"].push_back({{"h",h},
                {"plus",diagPoint(s.geodesicUpdate(d,h).geometry(),problem,evaluator)},
                {"minus",diagPoint(s.geodesicUpdate(-d,h).geometry(),problem,evaluator)}});
        }
        row["directions"].push_back(std::move(direction));
    }
    return row;
}
}
int main(int argc,char** argv) {
    try {
        diagRequire(argc==3,"usage: direction_probe z weight");const double z=diagArg(argv[1]),weight=diagArg(argv[2]);
        diagRequire((z==0||z==2)&&(weight==.1||weight==1||weight==10),"case outside frozen matrix");
        spin::SpinDesignRequest r;r.designDomain=diagBox(.8,z);r.outerSurface=diagBox(1,z);r.minimumWallThickness=.15;r.materialDensity=1;
        r.mode=spin::SpinMode::SpinningTop;r.centerOfMassWeight=weight;
        r.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.outerSurface.vertices,r.outerSurface.faces,1))-
            physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.designDomain.vertices,r.designDomain.faces,1));
        r.settings.equationProfile=spin::SpinEquationProfile::DualMonomialRawV1;
        r.settings.maximumRestarts=1;r.settings.maximumFunctionEvaluationsPerRestart=60;
        spin::QuadricMomentEvaluatorSettings settings;settings.triangulation.icosphereSubdivisions=1;settings.triangulation.truncationRadius=8;
        DiagVector p;p<<-.31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03;
        const auto seed=spin::QuadricCoefficients10::fromUnnormalized(p);const Eigen::Vector3d c(0,0,z);
        p(0)=seed.constant()-seed.linearVector().dot(c)+c.dot(seed.quadraticMatrix()*c);
        p.segment<3>(1)=seed.linearVector()-2*seed.quadraticMatrix()*c;p.tail<6>()=seed.values().tail<6>();
        const auto initial=spin::QuadricCoefficients10::fromUnnormalized(p);
        const auto evaluator=spin::makeCanonicalFullRankMomentEvaluator(r,settings);
        std::vector<spin::QuadricCoefficients10> seen;
        std::vector<std::optional<spin::SpinMomentMapEvaluation>> maps;
        const auto recording=[&](const spin::QuadricCoefficients10& q) {
            seen.push_back(q);maps.push_back(std::nullopt);auto m=evaluator(q);maps.back()=m;return m;
        };
        const auto problem=spin::SpinMomentProblem::spinningTop(weight);
        const auto replay=spin::solveSpinMomentSystem(problem,r.settings,recording,{initial});
        diagRequire(!replay.trace.empty(),"replay has no accepted state");
        std::size_t late=0;for(const auto& t:replay.trace)if(t.rawEquationResidual && *t.rawEquationResidual>1e-12)late=static_cast<std::size_t>(t.functionEvaluations-1);
        DiagJson out={{"schema","m1-direction-diagnostic-v1"},{"z",z},{"weight",weight},
            {"trace_hash",replay.traceHash},{"event_hash",replay.eventHash},{"replay_evaluations",replay.functionEvaluations},
            {"replay_status",spin::spinSolverStatusToken(replay.status)},{"states",DiagJson::array()}};
        std::vector<std::size_t> selected={0};if(late!=0)selected.push_back(late);
        for(auto index:selected){std::cerr<<"probe state evaluation="<<index+1<<'\n';out["states"].push_back(diagState(index,seen,maps,problem,evaluator));}
        std::cout<<"M1_DIRECTION "<<out.dump()<<'\n';return 0;
    }catch(const std::exception& error){std::cerr<<"M1_DIRECTION_ERROR "<<error.what()<<'\n';return 2;}
}
