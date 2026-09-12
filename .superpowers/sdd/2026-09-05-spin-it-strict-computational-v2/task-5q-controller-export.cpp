#include "spin2014/SourceCageMomentFunctional.h"
#include "source_surface_test_utils.h"
#include <iostream>
#include <iomanip>

template<class T> void row(const T& xs) {
    for (Eigen::Index i=0;i<xs.size();++i) std::cout << ' ' << xs(i);
}
template<class T> void seq(const T& xs) { for (auto x:xs) std::cout << ' ' << x; }

int main() {
    using namespace spin2014;
    std::cout << std::setprecision(17);
    for (int fixture=0;fixture<2;++fixture) {
        const double scale=fixture ? 1.7 : 1.;
        const Eigen::Vector3d offset=fixture ? Eigen::Vector3d(-.2,.15,.6) : Eigen::Vector3d::Zero();
        auto outer=fixture ? surface_test::box({.31,.32,.3},{.69,.65,.7}) : surface_test::cube(.3,.7);
        auto design=fixture ? surface_test::box({.35,.36,.38},{.61,.6,.65}) : surface_test::cube(.34,.66);
        outer.v=(scale*outer.v).rowwise()+offset.transpose();
        design.v=(scale*design.v).rowwise()+offset.transpose();
        SourceSurfaceOctreeConfig sc;sc.maximumDepth=5;sc.rootOrigin=offset;sc.rootSide=scale;
        auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
        SourceAdaptiveDomain domain(source,2,5);
        SourceSurfaceCageConfig cc;cc.cageLevel=2;cc.rootOrigin=offset;cc.rootSide=scale;
        SourceCageGridConfig gc;gc.fineLevel=3;
        SourceCageFemConfig fc;SourceCageBbwConfig bc;
        const SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),fc),bc);
        if (!bbw.success()) throw std::runtime_error(bbw.terminalReason());
        const SourceCageBinding binding(domain,bbw,SourceCageBindingConfig{});
        const SourceConformingBoundary b(binding,SourceConformingBoundaryConfig{});
        const double density=fixture ? 2.3 : .8;
        const SourceCageMomentFunctional functional(b,density,SourceCageMomentFunctionalConfig{});
        Eigen::VectorXd beta(functional.betaCount());
        for (Eigen::Index i=0;i<beta.size();++i) beta(i)=.13+.017*i;
        auto parameters=functional.initialHandleParameters();
        for (std::size_t h=0;h<functional.handleCount();++h) {
            parameters(4*h)=1+.01*(int(h%5)-2);
            parameters(4*h+1)=scale*.006*(int(h%3)-1);
            parameters(4*h+2)=scale*.004*(int(h%5)-2);
            parameters(4*h+3)=scale*.005*(int(h%7)-3);
        }
        const auto r=functional.evaluate(beta,parameters);
        const auto v=functional.evaluateValues(beta,parameters);
        const auto& grid=bbw.fem().grid();const auto& cage=grid.cage();
        std::cout << "CASE " << fixture << ' ' << density << ' ' << beta.size() << ' ' << functional.handleCount() << '\n';
        std::cout << "IDENT " << functional.configSha256() << ' ' << r.outputSha256 << ' ' << v.outputSha256 << '\n';
        std::cout << "FLAGS " << r.hasJacobians << ' ' << v.hasJacobians << ' ' << v.betaJacobian.size() << ' ' << v.handleJacobian.size() << '\n';
        std::cout << "BETA";row(beta);std::cout << '\n';
        std::cout << "PARAM";row(parameters);std::cout << '\n';
        for (Eigen::Index i=0;i<grid.nodePositions().rows();++i) {
            std::cout << "NODE " << i;row(grid.nodePositions().row(i));row(bbw.finalWeights().row(i));std::cout << '\n';
        }
        for (std::size_t i=0;i<grid.cells().size();++i) {
            std::cout << "CELL " << i;seq(grid.cells()[i]);seq(grid.cellNodes()[i]);std::cout << '\n';
        }
        for (std::size_t h=0;h<functional.handleCount();++h) {
            std::cout << "HANDLE " << h;row(cage.boundaryVertices().row(cage.handleVertexIndices()[h]));std::cout << '\n';
        }
        for (Eigen::Index i=0;i<b.points().rows();++i) {
            std::cout << "POINT " << i;row(b.points().row(i));std::cout << '\n';
            std::cout << "MOVED " << i;row(r.deformedPoints.row(i));std::cout << '\n';
            std::cout << "VALUEPOINT " << i;row(v.deformedPoints.row(i));std::cout << '\n';
        }
        for (Eigen::Index i=0;i<outer.f.rows();++i) {
            std::cout << "OUTERFACE " << i;row(outer.f.row(i));std::cout << '\n';
        }
        for (std::size_t i=0;i<b.triangles().size();++i) {
            std::cout << "TRI " << i;seq(b.triangles()[i]);std::cout << '\n';
        }
        for (const auto& p:b.patches()) {
            const auto owner=[&](std::size_t id)->long long {return id==SourceConformingBoundary::noOwner?-1:binding.primitives()[id].betaIndex;};
            std::cout << "PATCH " << owner(p.negativeOwner) << ' ' << owner(p.positiveOwner) << ' ' << p.triangleBegin << ' ' << p.triangleCount << '\n';
        }
        std::cout << "MOMENTS";row(r.finalMoments.values());std::cout << '\n';
        std::cout << "VALUEMOMENTS";row(v.finalMoments.values());std::cout << '\n';
        for (Eigen::Index i=0;i<r.betaJacobian.cols();++i) {std::cout << "JB " << i;row(r.betaJacobian.col(i));std::cout << '\n';}
        for (Eigen::Index i=0;i<r.handleJacobian.cols();++i) {std::cout << "JH " << i;row(r.handleJacobian.col(i));std::cout << '\n';}
        std::cout << "END\n";
    }
}
