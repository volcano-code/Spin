#include "spin2014/SourceSurfaceOctree.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sys/resource.h>
using namespace surface_test;
using namespace spin2014;
SourceSurfaceOctree build(const Mesh& outer,const Mesh& design,int depth=5,std::size_t budget=5000000) {
    SourceSurfaceOctreeConfig c; c.maximumDepth=depth; c.maximumNodes=budget;
    return SourceSurfaceOctree(outer.v,outer.f,design.v,design.f,c);
}
void enclosure(const SourceSurfaceOctree& tree,const Ten& exact) {
    close(tree.designSurfaceMoments().values(),exact,"independent original surface moments");
    const double delta=exact(0)-tree.carveableMoments().mass();
    require(delta>=-1e-14 && delta<=tree.unresolvedBoundaryVolume()+1e-14,"conservative volume enclosure");
}
int main(int argc,char**) {
    try {
        auto outer=cube(0,1), design=cube(.25,.75);
        if(argc>1) {
            const auto start=std::chrono::steady_clock::now();
            auto tree=build(outer,design,9);
            const double h=1./512;
            close(tree.carveableMoments().values(),boxMoments(Eigen::Vector3d::Constant(.25+h),Eigen::Vector3d::Constant(.75-h)),"level9 all ten closed-touch moments");
            enclosure(tree,boxMoments(Eigen::Vector3d::Constant(.25),Eigen::Vector3d::Constant(.75)));
            require(tree.stats().maximumVisitedLevel==9 && tree.stats().levels.at(9)>0,"genuine depth9 traversal");
            require(tree.cells().size()==2097993,"fixed sparse full-child hierarchy count");
            require(std::abs(tree.unresolvedBoundaryVolume()-49153./8388608)<1e-14,"fixed closed boundary uncertainty");
            rusage usage{}; getrusage(RUSAGE_SELF,&usage);
            std::cout<<std::setprecision(17)<<"depth9 nodes="<<tree.cells().size()<<" internal="<<tree.stats().internalNodes<<" boundary="<<tree.stats().boundaryNodes<<" full="<<tree.stats().fullNodes<<" empty="<<tree.stats().emptyNodes<<" max_level="<<tree.stats().maximumVisitedLevel<<" stored_cell_bytes="<<sizeof(SourceSurfaceCell)*tree.cells().size()<<" peak_rss_KiB="<<usage.ru_maxrss<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" volume="<<tree.carveableMoments().mass()<<" uncertainty="<<tree.unresolvedBoundaryVolume()<<" partition="<<tree.partitionSha256()<<'\n';
            std::cout<<"depth9 moments="<<tree.carveableMoments().values().transpose()<<"\nlevel_counts=";
            for(auto count:tree.stats().levels) std::cout<<count<<',';
            std::cout<<'\n';
            return 0;
        }
        auto tree=build(outer,design);
        close(tree.outerSurfaceMoments().values(),boxMoments({0,0,0},{1,1,1}),"original outer moments");
        const double h=1./32;
        close(tree.carveableMoments().values(),boxMoments(Eigen::Vector3d::Constant(.25+h),Eigen::Vector3d::Constant(.75-h)),"depth5 all moments");
        auto same=build(outer,design);
        require(tree.sourceSha256()==same.sourceSha256() && tree.configSha256()==same.configSha256() && tree.partitionSha256()==same.partitionSha256(),"deterministic scientific identities");
        for(const auto& cell:tree.cells()) {
            auto query=tree.cell(cell.level,cell.coordinates);
            close(query.moments.values(),cell.moments.values(),"stored query");
            if(cell.kind==SourceSurfaceCellKind::Internal) {
                Ten sum=Ten::Zero(); double uncertainty=0;
                for(int child=0;child<8;++child) {
                    const auto& c=tree.cells().at(cell.firstChild+child);
                    require(c.parent==cell.index && c.level==cell.level+1,"lineage");
                    for(int axis=0;axis<3;++axis) require(c.coordinates[axis]==2*cell.coordinates[axis]+((child>>axis)&1),"fixed child order");
                    sum+=c.moments.values(); uncertainty+=c.unresolvedVolume;
                }
                close(cell.moments.values(),sum,"parent moment sum");
                require(cell.unresolvedVolume==uncertainty,"parent uncertainty sum");
            } else if(cell.level<5) {
                Ten sum=Ten::Zero();
                for(int child=0;child<8;++child) {
                    auto xyz=cell.coordinates; for(int a=0;a<3;++a) xyz[a]=xyz[a]*2+((child>>a)&1);
                    auto q=tree.cell(cell.level+1,xyz);
                    require(q.index==SourceSurfaceCell::noIndex && q.parent==cell.index,"virtual terminal descendant");
                    sum+=q.moments.values();
                }
                close(cell.moments.values(),sum,"terminal descendant partition");
            }
        }
        const auto source=tree.sourceSha256(), partition=tree.partitionSha256();
        design.v.setZero(); design.f.setZero(); outer.v.setZero(); outer.f.setZero();
        require(tree.sourceSha256()==source && tree.partitionSha256()==partition && tree.designVertices()(0,0)==.25 && tree.outerVertices()(7,0)==1,"owned source immutability");
        require(tree.designFaces()(0,1)==2 && tree.outerFaces()(0,1)==2,"owned face copy");
        close(tree.cell(0,{0,0,0}).moments.values(),boxMoments(Eigen::Vector3d::Constant(.25+h),Eigen::Vector3d::Constant(.75-h)),"caller mutation cannot alter saved query");
        outer=cube(0,1); design=cube(.25,.75);
        require(build(outer,cube(.25,.7)).sourceSha256()!=source,"changed original source identity");
        require(build(outer,design,4).configSha256()!=tree.configSha256(),"depth identity");
        const auto budgetVariant=build(outer,design,5,5000001);
        require(budgetVariant.sourceSha256()==tree.sourceSha256() && budgetVariant.configSha256()!=tree.configSha256() && budgetVariant.partitionSha256()!=tree.partitionSha256(),"configured node guard must bind scientific identity");
        auto shell=join(design,cube(7./16,9./16),true);
        auto cavity=build(outer,shell,5);
        enclosure(cavity,boxMoments({.25,.25,.25},{.75,.75,.75})-boxMoments({7./16,7./16,7./16},{9./16,9./16,9./16}));
        require(cavity.cell(5,{16,16,16}).kind==SourceSurfaceCellKind::Empty,"nested cavity actual material classification");
        auto tinyCavity=join(design,cube(.39,.40),true);
        auto enclosed=build(outer,tinyCavity,3);
        require(enclosed.cell(3,{3,3,3}).kind==SourceSurfaceCellKind::Boundary,"triangles entirely inside cell cannot be missed");
        auto islands=join(cube(.11,.19),cube(.81,.89));
        enclosure(build(outer,islands,5),boxMoments({.11,.11,.11},{.19,.19,.19})+boxMoments({.81,.81,.81},{.89,.89,.89}));
        auto tinyIsland=build(outer,cube(.39,.40),3);
        require(tinyIsland.cell(3,{3,3,3}).kind==SourceSurfaceCellKind::Boundary,"enclosed tiny island boundary");
        Mesh t; t.v.resize(4,3); t.v<<.1,.2,.15,.8,.25,.2,.2,.8,.3,.3,.35,.9;
        t.f.resize(4,3); t.f<<1,2,3,0,2,1,0,1,3,0,3,2;
        enclosure(build(outer,t,5),tetraMoments(t.v));
        // Valid material islands nested inside a reversed void remain legal.
        auto nested=join(shell,cube(.48,.52)); enclosure(build(outer,nested,4),boxMoments({.25,.25,.25},{.75,.75,.75})-boxMoments({7./16,7./16,7./16},{9./16,9./16,9./16})+boxMoments({.48,.48,.48},{.52,.52,.52}));
        require(build(outer,design,0).cells().size()==1,"depth zero boundary terminal");
        const double tiny=std::ldexp(1.,-21);
        auto depth20=build(outer,cube(tiny,3*tiny),20);
        require(depth20.stats().maximumVisitedLevel==20 && depth20.cells().size()<10000,"supported depth20 sparse traversal");
        enclosure(depth20,boxMoments(Eigen::Vector3d::Constant(tiny),Eigen::Vector3d::Constant(3*tiny)));
        // Query arithmetic must also preserve partitioning away from dyadic
        // origin/side defaults, using shared represented boundary coordinates.
        SourceSurfaceOctreeConfig shiftedConfig;
        shiftedConfig.rootOrigin={-.3,-.2,-.1}; shiftedConfig.rootSide=1.4; shiftedConfig.maximumDepth=4;
        auto shiftedOuter=box({-.2,-.1,0},{1,1,1});
        SourceSurfaceOctree shifted(shiftedOuter.v,shiftedOuter.f,design.v,design.f,shiftedConfig);
        enclosure(shifted,boxMoments({.25,.25,.25},{.75,.75,.75}));
        const auto savedOrigin=shifted.config().rootOrigin;
        shiftedConfig.rootOrigin.setZero(); shiftedConfig.maximumDepth=0;
        require(shifted.config().rootOrigin==savedOrigin && shifted.config().maximumDepth==4,"owned immutable config");
        rejects([&]{tree.cell(6,{0,0,0});},"invalid query depth");
        rejects([&]{tree.cell(-1,{0,0,0});},"negative query depth");
        rejects([&]{tree.cell(3,{8,0,0});},"invalid query coordinate");
        rejects([&]{build(outer,design,-1);},"negative depth");
        rejects([&]{build(outer,design,21);},"excess depth");
        bool resource=false; try { build(outer,design,5,8); } catch(const std::exception& e) {resource=std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos;}
        require(resource,"explicit node-budget fail closed");
        rejects([&]{build(outer,design,5,0);},"zero budget");
        auto config=SourceSurfaceOctreeConfig{};
        auto construct=[&]{SourceSurfaceOctree x(outer.v,outer.f,design.v,design.f,config);};
        config.rootSide=0; rejects(construct,"zero root side");
        config.rootSide=std::numeric_limits<double>::infinity(); rejects(construct,"nonfinite root");
        config.rootSide=.5; rejects(construct,"root clipping");
        config.rootSide=1; config.rootOrigin.x()=std::numeric_limits<double>::quiet_NaN(); rejects(construct,"NaN root");
        config.rootOrigin=Eigen::Vector3d::Constant(1e16); config.rootSide=4; config.maximumDepth=20; rejects(construct,"collapsed dyadic coordinates");
        config.rootOrigin.setZero(); config.rootSide=1e110; rejects(construct,"volume overflow");
        auto bad=design; bad.f.conservativeResize(11,3); rejects([&]{build(outer,bad);},"open mesh");
        bad=design; bad.f=design.f.rowwise().reverse(); rejects([&]{build(outer,bad);},"flipped root");
        rejects([&]{build(outer,join(design,cube(.4,.6)));},"misnested signs");
        rejects([&]{build(outer,join(cube(.2,.6),cube(.4,.8)));},"self intersection");
        rejects([&]{build(outer,cube(.8,1.1));},"outside root design");
        rejects([&]{build(cube(.1,.9),cube(.05,.2));},"outside outer within root");
        rejects([&]{build(join(outer,cube(.4,.6),true),cube(.45,.55));},"design inside outer void");
        bad=design; bad.v.conservativeResize(9,3); bad.v.row(8)<<.5,.5,.5; rejects([&]{build(outer,bad);},"unreferenced vertex");
        bad=design; bad.v(0,0)=std::numeric_limits<double>::quiet_NaN(); rejects([&]{build(outer,bad);},"NaN geometry");
        bad=outer; bad.f.conservativeResize(11,3); rejects([&]{build(bad,design);},"invalid original outer boundary");
        std::cout<<"source surface hierarchy, signed material, queries, hashes, validation passed; faster="<<MAKE_IT_STAND_HAS_SPIN_IT_FASTER<<" exact="<<MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS<<'\n';
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n'; return 1;}
}
