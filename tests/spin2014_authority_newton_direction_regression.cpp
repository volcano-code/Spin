#include "spin2014/QuadricAuthorityNewtonState.h"
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>

static void require(bool c,const char* message){if(!c)throw std::runtime_error(message);}
int main() {
 try {
    using namespace spin2014;
    using A=spin::QuadricCoefficientAuthority10;
    A::RationalVector q;q.fill("0/1");q[0]="1/1";
    const auto authority=A::fromNormalizedGeometry(q);
    const AuthorityMomentVector10 moments{"2/1","0/1","0/1","1/1","3/1","4/1","5/1","0/1","0/1","0/1"};
    AuthorityMomentMatrix10 jacobian;
    for(unsigned i=0;i<10;++i)for(unsigned j=0;j<10;++j)jacobian[i][j]=i==j?"1/1":"0/1";
    QuadricObjectiveSettings settings;settings.objective=Spin2014Objective::StaticLowering;
    const auto state=linearizeSpin2014AuthorityNewtonState(settings,authority,moments,jacobian,"hand-identity-map");
    const auto direction=solveSpin2014AuthorityNewtonDirection(settings,state,"1/1000000000000");
    require(direction.defined && direction.linearSolve.solved,"protected authority equation did not reach direction");
    require(direction.initialArcLength=="1/2","existing maximum geodesic step lost");
    require(mpq_class(direction.dualDirection[0])==0,"radial component not projected away");
    require(mpq_class(direction.dualDirection[3])==mpq_class(-2.0/3.0),"hand-derived represented direction changed");
    for(unsigned i:{1U,2U,4U,5U,6U,7U,8U,9U})
        require(mpq_class(direction.dualDirection[i])==0,"unexpected transverse direction");
    const auto updated=authority.geodesicUpdate(direction.dualDirection,direction.initialArcLength);
    require(mpq_class(updated.dual()[3])<0 && updated.geometry()!=authority.geometry(),"actual authority update absent");
    require(direction.stateIdentity==state.scientificHash && !direction.scientificHash.empty(),"direction not bound to actual state");
    // Original equation has a nonzero binary64 candidate error. A tighter
    // caller limit must reject it, not regularize or issue a direction.
    const auto tooStrict=solveSpin2014AuthorityNewtonDirection(settings,state,"1/1000000000000000000000000000000");
    require(!tooStrict.defined && !tooStrict.linearSolve.solved,"failed original-equation check issued direction");
    auto modified=state;modified.rightHandSide[3]="-1/1";
    bool rejected=false;
    try{(void)solveSpin2014AuthorityNewtonDirection(settings,modified,"1/1000000000000");}
    catch(const std::exception&){rejected=true;}
    require(rejected,"modified equation reused old state identity");
    auto wrongObjective=settings;wrongObjective.centerWeight=2;
    rejected=false;
    try{(void)solveSpin2014AuthorityNewtonDirection(wrongObjective,state,"1/1000000000000");}
    catch(const std::exception&){rejected=true;}
    require(rejected,"different objective reused old Newton state");
    std::cout<<"protected authority equation, tangent projection, geodesic update and stale-state rejection passed\n";
    std::cout<<"direction hash="<<direction.scientificHash<<'\n';
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
