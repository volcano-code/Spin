#include "spin2014/SlqpTangentQp.h"
#include "spin2014/SlqpLbfgsSolver.h"
#include <iostream>
#include <stdexcept>
using namespace spin2014;
static void require(bool v,const char* s){if(!v)throw std::runtime_error(s);}
int main(){try{
    TangentQpSettings automatic;automatic.maximumFaces=-1;
    Eigen::VectorXd x=Eigen::VectorXd::Constant(5,.5),g(5);g<<2.,-3.,4.,-5.,.25;
    Eigen::MatrixXd a(0,5);Eigen::VectorXd c(0);
    auto fixed=automatic;fixed.maximumFaces=2;
    auto capped=solvePrimalBoxQp(x,g,a,c,4,{},fixed);
    require(!capped.checked&&capped.trace.back().action=="FACE_LIMIT","explicit cap was bypassed");
    auto full=solvePrimalBoxQp(x,g,a,c,4,{},automatic);
    require(full.checked,"bounded-variable capacity not supported");
    Eigen::VectorXd expected(5);expected<<-.5,.5,-.5,.5,-.25;
    require((full.direction-expected).lpNorm<Eigen::Infinity>()<1e-12,"wrong analytic/free-variable solution");
    require(full.workingSet.size()==4,"bounded prefix not fully supported");
    fixed.maximumFaces=0;
    require(!solvePrimalBoxQp(x,g,a,c,4,{},fixed).checked,"legacy zero budget bypassed");
    fixed=automatic;fixed.maximumIterations=1;
    auto stopped=solvePrimalBoxQp(x,g,a,c,4,{},fixed);
    require(!stopped.checked&&stopped.trace.back().action=="ITERATION_LIMIT","iteration gate bypassed");
    auto free=solvePrimalBoxQp(x,g,a,c,0,{},automatic);
    require(free.checked&&(free.direction+g).norm()<1e-12,"zero bounded variables unsupported");
    Eigen::MatrixXd eq(1,2);eq<<1.,1.;Eigen::VectorXd residual(1);residual<<.1;
    auto phase=solvePrimalBoxQp(Eigen::Vector2d(.8,0.),Eigen::Vector2d(-3.,1.),eq,residual,1,{},automatic);
    require(phase.checked&&(phase.direction-Eigen::Vector2d(.2,-.3)).norm()<1e-9,"auto capacity not passed to Phase I");
    require(phase.trace.front().action.find("PHASE_I_")==0,"phase fixture did not exercise Phase I");
    fixed.maximumFaces=-2;
    require(!solvePrimalBoxQp(x,g,a,c,4,{},fixed).checked,"invalid negative policy accepted");
    require(!solveTangentQp(x,g,a,c,4,{},automatic).checked,"legacy tangent semantics changed");
    require(!solveFiniteBoxQp(x,g,a,c,4,{},automatic).checked,"legacy finite-box semantics changed");
    SlqpLbfgsSettings controls;controls.algorithm=SlqpAlgorithm::AuditedKktRestorationV1;
    controls.directionPolicy=SlqpDirectionPolicy::FiniteBoxPrimalV2;controls.metricTangent=automatic;
    validateSlqpDirectionSettings(controls);
    controls.directionPolicy=SlqpDirectionPolicy::FiniteBoxQpV1;
    bool rejected=false;try{validateSlqpDirectionSettings(controls);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"automatic capacity allowed in unsupported policy");
    std::cout<<"bounded-capacity analytic, phase I, compatibility and resource gates passed\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
