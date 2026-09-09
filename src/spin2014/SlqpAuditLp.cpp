#include "spin2014/SlqpKkt.h"
#include <Highs.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace spin2014 {
void validateSlqpAuditSettings(const SlqpAuditSettings& s) {
    if(!std::isfinite(s.restorationTrustBound)||s.restorationTrustBound<=0 ||
       s.lpIterationLimit<=0 || !std::isfinite(s.lpTimeLimit)||s.lpTimeLimit<=0 ||
       !std::isfinite(s.lpFeasibilityTolerance)||s.lpFeasibilityTolerance<1e-10 || s.lpFeasibilityTolerance>1e-4)
        throw std::invalid_argument("invalid audited SLQP LP controls");
}
namespace {
// Row-wise sparse assembly: O(n*m) for the witness, O(n*m+m) restoration.
struct SmallLp {
    HighsModel model;
    explicit SmallLp(int n) {
        auto& p=model.lp_;p.num_col_=n;p.num_row_=0;
        p.col_cost_.assign(n,0.);p.col_lower_.assign(n,-kHighsInf);p.col_upper_.assign(n,kHighsInf);
        p.a_matrix_.format_=MatrixFormat::kRowwise;p.a_matrix_.start_={0};
    }
    void row(const std::vector<std::pair<int,double>>& values,double lower,double upper) {
        auto& p=model.lp_;p.num_row_++;p.row_lower_.push_back(lower);p.row_upper_.push_back(upper);
        for(auto [i,v]:values) if(v!=0.) {p.a_matrix_.index_.push_back(i);p.a_matrix_.value_.push_back(v);}
        p.a_matrix_.start_.push_back(p.a_matrix_.value_.size());
    }
    bool solve(const SlqpAuditSettings& s,Eigen::VectorXd& x,std::string& status) {
        Highs h;
        bool ok=h.setOptionValue("output_flag",false)==HighsStatus::kOk;
        ok&=h.setOptionValue("threads",1)==HighsStatus::kOk;
        ok&=h.setOptionValue("parallel","off")==HighsStatus::kOk;
        ok&=h.setOptionValue("random_seed",0)==HighsStatus::kOk;
        ok&=h.setOptionValue("solver","simplex")==HighsStatus::kOk;
        ok&=h.setOptionValue("presolve","on")==HighsStatus::kOk;
        ok&=h.setOptionValue("simplex_iteration_limit",s.lpIterationLimit)==HighsStatus::kOk;
        ok&=h.setOptionValue("time_limit",s.lpTimeLimit)==HighsStatus::kOk;
        ok&=h.setOptionValue("primal_feasibility_tolerance",s.lpFeasibilityTolerance)==HighsStatus::kOk;
        ok&=h.setOptionValue("dual_feasibility_tolerance",s.lpFeasibilityTolerance)==HighsStatus::kOk;
        if(!ok) {status="LP_OPTION_ERROR";return false;}
        if(h.passModel(model)==HighsStatus::kError) {status="LP_MODEL_ERROR";return false;}
        const auto run=h.run();status=h.modelStatusToString(h.getModelStatus());
        if(run==HighsStatus::kError || h.getModelStatus()!=HighsModelStatus::kOptimal) return false;
        const auto& solution=h.getSolution();
        if(!solution.value_valid || solution.col_value.size()!=static_cast<std::size_t>(model.lp_.num_col_)) return false;
        x=Eigen::Map<const Eigen::VectorXd>(solution.col_value.data(),solution.col_value.size());
        return x.allFinite();
    }
};
}
KktResidual computeKktWitness(const AdaptiveFillEvaluation& e,const Eigen::VectorXd& x,
    Eigen::Index b,double tol,const SlqpAuditSettings& s) {
    validateSlqpAuditSettings(s);KktResidual k;
    if(!validSlqpEvaluation(e,x.size())||!x.allFinite()||b<0||b>x.size()||!std::isfinite(tol)||tol<0||tol>=.5) return k;
    const int m=e.constraints.size();SmallLp lp(m+1);lp.model.lp_.col_cost_[m]=1;lp.model.lp_.col_lower_[m]=0;
    for(int i=0;i<x.size();++i) {
        const bool lower=i<b&&x(i)<=tol,upper=i<b&&x(i)>=1.-tol;
        auto row=[&](double sign,double scale) {
            std::vector<std::pair<int,double>> v;
            for(int j=0;j<m;++j)v.emplace_back(j,sign*scale*e.constraintJacobian(j,i));
            v.emplace_back(m,-1.);lp.row(v,-kHighsInf,-sign*scale*e.gradient(i));
        };
        if(!lower)row(1,1.);if(!upper)row(-1,1.);
        if(lower)row(1,x(i));if(upper)row(-1,1.-x(i));
    }
    Eigen::VectorXd solution;std::string status;
    if(!lp.solve(s,solution,status)) {k.status="WITNESS_LP_UNKNOWN: "+status;return k;}
    Eigen::VectorXd lambda=solution.head(m),alpha=Eigen::VectorXd::Zero(x.size()),eta=alpha;
    Eigen::VectorXd r=e.gradient+e.constraintJacobian.transpose()*lambda;
    if(!r.allFinite()) {k.status="DERIVED_OVERFLOW";return k;}
    for(int i=0;i<b;++i) {if(x(i)<=tol)alpha(i)=std::max(0.,r(i));else if(x(i)>=1.-tol)eta(i)=std::max(0.,-r(i));}
    k=evaluateKkt(e,x,b,lambda,alpha,eta);if(k.valid)k.status="WITNESS_LP_OPTIMAL_FIRST_ORDER";return k;
}
RestorationEvidence solveRestoration(const AdaptiveFillEvaluation& e,const Eigen::VectorXd& x,
    Eigen::Index b,const SlqpAuditSettings& s) {
    validateSlqpAuditSettings(s);RestorationEvidence r;r.attempted=true;
    if(!validSlqpEvaluation(e,x.size())||!x.allFinite()||b<0||b>x.size()) {r.status="INVALID_DATA";return r;}
    const int n=x.size(),m=e.constraints.size();SmallLp lp(n+2*m);
    for(int i=0;i<n;++i) {
        lp.model.lp_.col_lower_[i]=i<b?std::max(-s.restorationTrustBound,-x(i)):-s.restorationTrustBound;
        lp.model.lp_.col_upper_[i]=i<b?std::min(s.restorationTrustBound,1.-x(i)):s.restorationTrustBound;
    }
    for(int j=0;j<2*m;++j) {lp.model.lp_.col_lower_[n+j]=0;lp.model.lp_.col_cost_[n+j]=1;}
    for(int j=0;j<m;++j) {
        std::vector<std::pair<int,double>> v;
        for(int i=0;i<n;++i)v.emplace_back(i,e.constraintJacobian(j,i));
        v.emplace_back(n+j,-1);v.emplace_back(n+m+j,1);lp.row(v,-e.constraints(j),-e.constraints(j));
    }
    Eigen::VectorXd solution;
    if(!lp.solve(s,solution,r.status))return r;
    r.direction=solution.head(n);
    r.predictedResidual=(e.constraintJacobian*r.direction+e.constraints).lpNorm<1>();
    if(!std::isfinite(r.predictedResidual))r.status="DERIVED_OVERFLOW";
    return r;
}
} // namespace spin2014
