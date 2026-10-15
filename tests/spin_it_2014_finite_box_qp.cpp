#include "spin2014/SlqpTangentQp.h"
#include <iostream>
#include <stdexcept>
#include <limits>
#ifdef TEST_PRIMAL
#define solveFiniteBoxQp solvePrimalBoxQp
#endif
#ifdef BASELINE_TANGENT
#define solveFiniteBoxQp solveTangentQp
#endif
void require(bool v, const char* m) { if(!v) throw std::runtime_error(m); }
void near(double v,double w) { require(std::isfinite(v)&&std::abs(v-w)<1e-9,"wrong analytic box solution"); }
int main() { try {
    using namespace spin2014;
    Eigen::MatrixXd a(1,2); a<<1.,1.;
    Eigen::VectorXd c=Eigen::VectorXd::Zero(1);
    const Eigen::MatrixXd empty(0,2); const Eigen::VectorXd zero(0);
    // Removing finite interior bounds gives d=(1,-1), not (.2,-.2).
    auto r=solveFiniteBoxQp(Eigen::Vector2d(.8,.2),Eigen::Vector2d(0.,2.),a,c,2,{});
    require(r.checked,"box QP not checked");near(r.direction[0],.2);near(r.direction[1],-.2);
    auto repeat=solveFiniteBoxQp(Eigen::Vector2d(.8,.2),Eigen::Vector2d(0.,2.),a,c,2,{});
    require(r.direction==repeat.direction && r.workingSet==repeat.workingSet,"nondeterministic solve");
    // Nonzero equality RHS and unbounded phi-like coordinate.
    c[0]=.1;
    r=solveFiniteBoxQp(Eigen::Vector2d(.8,0.),Eigen::Vector2d(-3.,1.),a,c,1,{});
    require(r.checked,"nonzero equality failed");near(r.direction[0],.2);near(r.direction[1],-.3);
    // Independent inverse metric [[1,2],[2,9]], B=[[1.8,-.4],[-.4,.2]].
    LbfgsHistory h;h.capacity=2;require(h.update(Eigen::Vector2d(1.,2.),Eigen::Vector2d(1.,0.)),"history rejected");
    r=solveFiniteBoxQp(Eigen::Vector2d(.5,.5),Eigen::Vector2d(0.,1.),empty,zero,1,h);
    require(r.checked,"coupled box failed");near(r.direction[0],-.5);near(r.direction[1],-6.);
    require(r.workingSet==std::vector<int>{0},"lower face encoding");near(r.multipliers[0],-1.5);
    r=solveFiniteBoxQp(Eigen::Vector2d(.5,.5),Eigen::Vector2d(0.,-1.),empty,zero,1,h);
    require(r.checked,"upper coupled box failed");near(r.direction[0],.5);near(r.direction[1],6.);
    require(r.workingSet==std::vector<int>{1},"upper face encoding");near(r.multipliers[0],1.5);
    // A tolerance-sized outward drift still gives EXACTLY zero feasible step.
    // Add the face and re-solve; do not clamp a free coordinate after checking.
    r=solveFiniteBoxQp(Eigen::Vector2d(0.,.5),Eigen::Vector2d(1e-12,1.),empty,zero,2,{});
    require(r.checked && r.direction[0]==0.,"tiny outward boundary drift gives zero step");
    near(r.direction[1],-.5);
    r=solveFiniteBoxQp(Eigen::Vector2d(1.,.5),Eigen::Vector2d(-1e-12,-1.),empty,zero,2,{});
    require(r.checked && r.direction[0]==0.,"tiny upper drift gives zero step");
    near(r.direction[1],.5);
    r=solveFiniteBoxQp(Eigen::Vector2d(1.,.5),Eigen::Vector2d(-1e-18,1.),empty,zero,2,{});
    require(r.checked && r.direction[0]==0.,"upper sub-ULP outward drift gives zero step");
    // Contradictory bounds/equality cannot be promoted to a certificate.
    a<<1.,0.;c[0]=2.;
    r=solveFiniteBoxQp(Eigen::Vector2d(.5,.5),Eigen::Vector2d::Ones(),a,c,2,{});
    require(!r.checked && r.status=="BOX_QP_UNKNOWN","false infeasibility/success claim");
    TangentQpSettings s;s.maximumFaces=0;
    require(!solveFiniteBoxQp(Eigen::Vector2d(.5,.5),Eigen::Vector2d(2.,0.),empty,zero,2,{},s).checked,"face cap ignored");
    s={};s.maximumIterations=1;
    require(!solveFiniteBoxQp(Eigen::Vector2d(.5,.5),Eigen::Vector2d(2.,0.),empty,zero,2,{},s).checked,"iteration cap ignored");
    s={};s.residualTolerance=std::numeric_limits<double>::quiet_NaN();
    require(!solveFiniteBoxQp(Eigen::Vector2d(.5,.5),Eigen::Vector2d(2.,0.),empty,zero,2,{},s).checked,"NaN accepted");
    require(!solveFiniteBoxQp(Eigen::Vector2d(-.1,.5),Eigen::Vector2d(2.,0.),empty,zero,2,{}).checked,"invalid box accepted");
    std::cout<<"finite box analytic tests passed\n";return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;} }
