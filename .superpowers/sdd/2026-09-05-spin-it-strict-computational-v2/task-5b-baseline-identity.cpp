#include "spin2014/SlqpLbfgsSolver.h"
#include <iostream>
namespace spin2014 {
SlqpLbfgsResult baselineSolveSlqpLbfgsPaper(const AdaptiveFillFunctional&,const Eigen::VectorXd&,double,const SlqpLbfgsSettings&);
}
int main() {
    voxel::GridSpec grid;grid.dimensions=Eigen::Vector3i(4,4,4);grid.minimum=Eigen::Vector3d(-1.5,-3,-1.5);grid.maximum=Eigen::Vector3d(1.5,3,1.5);
    std::vector<unsigned char> material(64,1),carveable(64,0);
    for(int z=1;z<3;++z)for(int y=1;y<3;++y)for(int x=1;x<3;++x)carveable[x+4*(y+4*z)]=1;
    auto source=voxel::VoxelDomain::create(grid,material,{},carveable);
    spin2014::AdaptiveVoxelDomain domain(source,2,2);
    spin2014::AdaptiveFillFunctional f(domain,spin2014::SpinIt2014MomentProblem::fromObjective(spin2014::Spin2014Objective::TopInertiaOnly),.03);
    Eigen::VectorXd x=Eigen::VectorXd::Constant(domain.leaves().size(),.2);
    spin2014::SlqpLbfgsSettings s;s.maximumIterations=180;s.equalityTolerance=1e-9;
    const auto base=spin2014::baselineSolveSlqpLbfgsPaper(f,x,1.7,s),current=spin2014::solveSlqpLbfgsPaper(f,x,1.7,s);
    std::cout<<"BASE080b1f1 paper trace hash="<<base.traceHash<<"\ncurrent paper trace hash="<<current.traceHash<<'\n';
    return base.traceHash==current.traceHash && base.beta.isApprox(current.beta,0.) && base.finalPhi==current.finalPhi ? 0:1;
}
