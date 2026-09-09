#include "spin2014/SlqpKkt.h"
#include <algorithm>
#include <cmath>

namespace spin2014 {
bool validSlqpEvaluation(const AdaptiveFillEvaluation& e,Eigen::Index n) {
    return std::isfinite(e.objective) && e.gradient.size()==n &&
        e.constraintJacobian.cols()==n && e.constraintJacobian.rows()==e.constraints.size() &&
        e.gradient.allFinite() && e.constraints.allFinite() && e.constraintJacobian.allFinite();
}
KktResidual evaluateKkt(const AdaptiveFillEvaluation& e,const Eigen::VectorXd& x,
    Eigen::Index b,const Eigen::VectorXd& lambda,const Eigen::VectorXd& alpha,const Eigen::VectorXd& eta) {
    KktResidual k;k.lambda=lambda;k.alpha=alpha;k.eta=eta;
    if(!validSlqpEvaluation(e,x.size()) || !x.allFinite() || b<0 || b>x.size() ||
       lambda.size()!=e.constraints.size() || alpha.size()!=x.size() || eta.size()!=x.size() ||
       !lambda.allFinite() || !alpha.allFinite() || !eta.allFinite()) return k;
    const Eigen::VectorXd r=e.gradient+e.constraintJacobian.transpose()*lambda-alpha+eta;
    if(!r.allFinite()) {k.status="DERIVED_OVERFLOW";return k;}
    k.equalityPrimal=e.constraints.size()?e.constraints.cwiseAbs().maxCoeff():0.;
    k.stationarity=r.size()?r.cwiseAbs().maxCoeff():0.;
    k.boundPrimal=k.dualSign=k.complementarity=0.;
    for(Eigen::Index i=0;i<x.size();++i) {
        k.dualSign=std::max({k.dualSign,-alpha(i),-eta(i)});
        if(i>=b) {
            if(alpha(i)!=0. || eta(i)!=0.) {k.status="UNBOUNDED_COORDINATE_MULTIPLIER";k.full=INFINITY;return k;}
        } else {
            k.boundPrimal=std::max({k.boundPrimal,-x(i),x(i)-1.});
            const double lower=alpha(i)*x(i),upper=eta(i)*(1.-x(i));
            if(!std::isfinite(lower)||!std::isfinite(upper)) {k.status="DERIVED_OVERFLOW";k.full=INFINITY;return k;}
            k.complementarity=std::max({k.complementarity,std::abs(lower),std::abs(upper)});
        }
    }
    k.full=std::max({k.equalityPrimal,k.boundPrimal,k.stationarity,k.dualSign,k.complementarity});
    k.valid=std::isfinite(k.full);k.status=k.valid?"VALID_FIRST_ORDER_RESIDUAL":"DERIVED_OVERFLOW";
    return k;
}
} // namespace spin2014
