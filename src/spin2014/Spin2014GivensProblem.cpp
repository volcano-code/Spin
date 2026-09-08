#include "spin2014/Spin2014GivensProblem.h"
#include <cmath>
#include <stdexcept>

namespace {
// A dedicated second-order scalar for the eleven paper coordinates. Chain
// rules are analytic; it is deliberately independent of the Faster evaluator.
struct Scalar {
    double v;
    Eigen::Matrix<double,11,1> d = Eigen::Matrix<double,11,1>::Zero();
    Eigen::Matrix<double,11,11> h = Eigen::Matrix<double,11,11>::Zero();
    Scalar(double value=0.) : v(value) {}
};
Scalar operator+(const Scalar& a,const Scalar& b) { Scalar r(a.v+b.v);r.d=a.d+b.d;r.h=a.h+b.h;return r; }
Scalar operator-(const Scalar& a,const Scalar& b) { Scalar r(a.v-b.v);r.d=a.d-b.d;r.h=a.h-b.h;return r; }
Scalar operator*(const Scalar& a,const Scalar& b) {
    Scalar r(a.v*b.v);r.d=a.d*b.v+b.d*a.v;
    r.h=a.h*b.v+b.h*a.v+a.d*b.d.transpose()+b.d*a.d.transpose();return r;
}
Scalar chain(const Scalar& a,double value,double first,double second) {
    Scalar r(value);r.d=first*a.d;r.h=first*a.h+second*a.d*a.d.transpose();return r;
}
Scalar inverse(const Scalar& a) { double u=1./a.v;return chain(a,u,-u*u,2*u*u*u); }
Scalar operator/(const Scalar& a,const Scalar& b) { return a*inverse(b); }
Scalar sine(const Scalar& a) {return chain(a,std::sin(a.v),std::cos(a.v),-std::sin(a.v));}
Scalar cosine(const Scalar& a) {return chain(a,std::cos(a.v),-std::sin(a.v),-std::cos(a.v));}
bool finite(const Scalar& a) {return std::isfinite(a.v)&&a.d.allFinite()&&a.h.allFinite();}
}

namespace spin2014 {
const char* spin2014FormulationToken(Spin2014Formulation f) noexcept {
    switch(f) {case Spin2014Formulation::LegacyAngleEliminated:return "LEGACY_ANGLE_ELIMINATED";case Spin2014Formulation::PaperGivens:return "PAPER_GIVENS";}return "INVALID";
}
Spin2014GivensProblem::Spin2014GivensProblem(Spin2014Objective objective,double inertiaWeight,double centerWeight)
    : objective_(objective),inertiaWeight_(inertiaWeight),centerWeight_(centerWeight) {
    (void)SpinIt2014MomentProblem::fromObjective(objective,inertiaWeight,centerWeight);
}
Spin2014GivensEvaluation Spin2014GivensProblem::evaluate(const Vector11& state) const {
    if(!state.allFinite()||state(0)<=0.|| (hasAngle() && state(4)+state(5)<=0.))
        throw std::invalid_argument("invalid paper Givens mass, axial inertia or state");
    Scalar x[11];for(int i=0;i<11;++i) {x[i]=Scalar(state(i));x[i].d(i)=1.;}
    Scalar objective;
    std::vector<Scalar> constraints={x[1],x[2]};
    if(hasAngle()) {
        const Scalar c=cosine(x[10]),s=sine(x[10]);
        const Scalar delta=objective_==Spin2014Objective::YoYo?Scalar(0.):x[3]*x[3]/x[0];
        const Scalar ia=c*c*(x[5]+x[6])+s*s*(x[4]+x[6])+Scalar(2.)*c*s*x[9]-delta;
        const Scalar ib=s*s*(x[5]+x[6])+c*c*(x[4]+x[6])-Scalar(2.)*c*s*x[9]-delta;
        const Scalar ic=x[4]+x[5];
        objective=Scalar(inertiaWeight_)*((ia/ic)*(ia/ic)+(ib/ic)*(ib/ic));
        if(objective_==Spin2014Objective::YoYo) constraints.push_back(x[3]);
        constraints.push_back(x[8]);constraints.push_back(x[7]);
        constraints.push_back(c*s*(x[4]-x[5])+(c*c-s*s)*x[9]);
    }
    if(objective_==Spin2014Objective::TopMassAndLowering)
        objective=objective+Scalar(centerWeight_)*x[3]*x[3];
    else if(objective_==Spin2014Objective::TopLoweringOnly||!hasAngle())
        objective=objective+Scalar(centerWeight_)*(x[3]/x[0])*(x[3]/x[0]);
    if(!finite(objective)) throw std::invalid_argument("nonfinite paper Givens objective derivatives");
    Spin2014GivensEvaluation r;r.objective=objective.v;r.objectiveGradient=objective.d;r.objectiveHessian=objective.h;
    r.constraints.resize(constraints.size());r.constraintJacobian.resize(constraints.size(),11);
    for(std::size_t i=0;i<constraints.size();++i) {
        if(!finite(constraints[i]))throw std::invalid_argument("nonfinite paper Givens constraint derivatives");
        r.constraints(i)=constraints[i].v;r.constraintJacobian.row(i)=constraints[i].d.transpose();r.constraintHessians.push_back(constraints[i].h);
    }return r;
}
}
