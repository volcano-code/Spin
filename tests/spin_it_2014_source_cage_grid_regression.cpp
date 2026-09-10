#include "spin2014/SourceCageGrid.h"
#include "source_surface_test_utils.h"
#include <algorithm>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <set>
#include <type_traits>
using namespace surface_test;
using namespace spin2014;
using Key=SourceCageGrid::Key;

SourceSurfaceCage cage(const Mesh& m,int level=3) {
    SourceSurfaceCageConfig c; c.cageLevel=level; return {m.v,m.f,c};
}
template<class F> void resource(F f) {
    bool ok=false;
    try {f();} catch(const std::length_error& e) {ok=std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos;}
    require(ok,"resource failure must be explicit RESOURCE_LIMIT");
}
void near(double a,double b) {require(std::abs(a-b)<=1e-12*std::max(1.,std::abs(b)),"analytic reconstruction to 1e-12");}
Eigen::MatrixXd fields(const Eigen::MatrixXd& p) {
    Eigen::MatrixXd v(p.rows(),8);
    for(Eigen::Index i=0;i<p.rows();++i) {
        const double x=p(i,0),y=p(i,1),z=p(i,2);
        v.row(i)<<1,x,y,z,x*y,x*z,y*z,x*y*z;
    }
    return v;
}
void reconstruct(const SourceCageGrid& g,const Eigen::MatrixXd& queries) {
    const auto got=g.interpolate(fields(g.nodePositions()),queries), expected=fields(queries);
    require(got.rows()==queries.rows() && got.cols()==8,"query by field shape preserved");
    for(Eigen::Index i=0;i<got.rows();++i) for(int k=0;k<8;++k) near(got(i,k),expected(i,k));
}
// Same independent analytic annulus extrusion fixture used by source-cage tests.
Mesh tunnel() {
    Mesh m; m.v.resize(16,3); m.f.resize(32,3);
    const double rings[2][4][2]={{{.1,.1},{.9,.1},{.9,.9},{.1,.9}},{{.3,.3},{.7,.3},{.7,.7},{.3,.7}}};
    for(int z=0;z<2;++z) for(int r=0;r<2;++r) for(int j=0;j<4;++j)
        m.v.row(8*z+4*r+j)<<rings[r][j][0],rings[r][j][1],z?.9:.1;
    int row=0;
    auto quad=[&](int a,int b,int c,int d) {m.f.row(row++)<<a,b,c; m.f.row(row++)<<a,c,d;};
    for(int j=0;j<4;++j) {
        int k=(j+1)%4;
        quad(j,k,k+8,j+8); quad(j+4,j+12,k+12,k+4);
        quad(j,j+4,k+4,k); quad(j+8,k+8,k+12,j+12);
    }
    return m;
}
void connectivity(const SourceCageGrid& g) {
    require(g.cells().size()==g.cellNodes().size(),"one connectivity row per fine cell");
    require(std::is_sorted(g.cells().begin(),g.cells().end()),"lexicographic fine cells");
    require(std::adjacent_find(g.cells().begin(),g.cells().end())==g.cells().end(),"unique cells");
    require(std::is_sorted(g.nodeKeys().begin(),g.nodeKeys().end()),"lexicographic nodes");
    require(std::adjacent_find(g.nodeKeys().begin(),g.nodeKeys().end())==g.nodeKeys().end(),"unique nodes");
    std::set<SourceCageGrid::Edge> fromCells;
    std::set<int> usedNodes;
    for(std::size_t i=0;i<g.cells().size();++i) for(int b=0;b<8;++b) {
        int node=g.cellNodes()[i][b]; require(node>=0 && std::size_t(node)<g.nodeKeys().size(),"valid node indices");
        usedNodes.insert(node);
        for(int a=0;a<3;++a) require(g.nodeKeys()[node][a]==g.cells()[i][a]+((b>>a)&1),"x-bit-first integer corners");
        for(int a=0;a<3;++a) if(!(b&(1<<a))) fromCells.insert({node,g.cellNodes()[i][b|(1<<a)]});
    }
    require(usedNodes.size()==g.nodeKeys().size(),"no orphan nodes");
    require(std::vector<SourceCageGrid::Edge>(fromCells.begin(),fromCells.end())==g.edges(),"unique sorted edges all belong to fine cells");
    for(const auto& e:g.edges()) {
        int distance=0; for(int a=0;a<3;++a) distance+=std::abs(int(g.nodeKeys()[e[1]][a])-int(g.nodeKeys()[e[0]][a]));
        require(distance==1,"every edge is one axis-neighbor step");
    }
    std::set<int> handles(g.handleNodeIndices().begin(),g.handleNodeIndices().end());
    require(handles.size()==g.cage().handleVertexIndices().size(),"all handles retained exactly once");
    for(std::size_t i=0;i<g.handleNodeIndices().size();++i)
        require(g.nodePositions().row(g.handleNodeIndices()[i])==g.cage().boundaryVertices().row(g.cage().handleVertexIndices()[i]),"coarse handles unchanged bit-exactly");
}
int main() {
    int failures=0;
    auto test=[&](const char* name,const std::function<void()>& run) {
        try {run(); std::cout<<"PASS "<<name<<'\n';}
        catch(const std::exception& e) {++failures; std::cerr<<"FAIL "<<name<<": "<<e.what()<<'\n';}
    };
    const auto source=cube(.25,.75); const auto coarse=cage(source);
    SourceCageGridConfig config; const SourceCageGrid g(coarse,config);
    test("refine_selected_volume_without_losing_full_cells_or_handles",[&] {
        require(coarse.selectedCells().size()==216 && coarse.handleVertexIndices().size()==218,"validated coarse analytic fixture");
        require(g.cells().size()==1728 && g.nodeKeys().size()==2197 && g.edges().size()==6084,"12^3 cells, 13^3 nodes, 3*12*13^2 edges");
        std::vector<Key> expected;
        for(unsigned x=2;x<14;++x) for(unsigned y=2;y<14;++y) for(unsigned z=2;z<14;++z) expected.push_back({x,y,z});
        require(g.cells()==expected,"exact literal fine cell union");
        require(g.nodePositions().minCoeff()==.125 && g.nodePositions().maxCoeff()==.875,"unchanged selected-volume envelope");
        connectivity(g);
    });
    test("multilinear_transfer_original_and_new_adaptive_vertices_without_resolve",[&] {
        Eigen::MatrixXd q(13,3); q.topRows(8)=cube(.125,.875).v;
        q.row(8)<<.3,.45,.6; q.row(9)<<.5,.45,.6; q.row(10)<<.5,.5,.6;
        q.row(11)<<.28125,.34375,.40625; q.row(12)<<.25,.375,.5;
        reconstruct(g,q); reconstruct(g,source.v); reconstruct(g,g.nodePositions());
        for(Eigen::Index i=0;i<g.nodePositions().rows();++i) {
            auto e=g.embed(g.nodePositions().row(i).transpose()); int ones=0;
            for(int b=0;b<8;++b) {require(e.coefficients[b]==0 || e.coefficients[b]==1,"nodal identity coefficients"); if(e.coefficients[b]==1) {++ones; require(e.nodeIndices[b]==i,"identity selects same node");}}
            require(ones==1,"one unit coefficient at every node");
        }
        Eigen::MatrixXd f(g.nodePositions().rows(),3);
        for(Eigen::Index i=0;i<f.rows();++i) {double a=.1+.01*(i%7),b=.2+.01*(i%11); f.row(i)<<a,b,1-a-b;}
        const auto out=g.interpolate(f,q);
        for(Eigen::Index i=0;i<out.rows();++i) {near(out.row(i).sum(),1); require(out.row(i).minCoeff()>0,"positive partition fixture stays positive");}
        f.setConstant(2); const auto nonpartition=g.interpolate(f,q);
        for(Eigen::Index i=0;i<nonpartition.rows();++i) for(int k=0;k<3;++k) near(nonpartition(i,k),2);
        const auto empty=g.interpolate(f,Eigen::MatrixXd(0,3)); require(empty.rows()==0 && empty.cols()==3,"empty query preserves columns");
    });
    test("closed_represented_ties_and_nextafter_do_not_snap_outside",[&] {
        const auto interior=g.embed({.3,.45,.6});
        require(interior.cellKey==Key{4,7,9},"interior selects literal fine cell");
        const double expected[8]={.064,.256,.016,.064,.096,.384,.024,.096};
        double sum=0;
        for(int b=0;b<8;++b) {near(interior.coefficients[b],expected[b]); require(std::isfinite(interior.coefficients[b]) && interior.coefficients[b]>=0,"finite nonnegative coefficients"); sum+=interior.coefficients[b];}
        near(sum,1);
        require(g.cells().at(interior.cellIndex)==interior.cellKey && g.cellNodes().at(interior.cellIndex)==interior.nodeIndices,"embedding index and connectivity refer to selected cell");
        const auto split=g.embed({.28125,.34375,.40625});
        require(split.cellKey==Key{4,5,6},"newly split analytic fine cell");
        for(double coefficient:split.coefficients) require(coefficient==.125,"newly split midpoint gives eight exact one-eighth coefficients");
        require(g.embed({.5,.45,.6}).cellKey==Key{7,7,9},"face chooses lexicographic first active incident cell");
        require(g.embed({.5,.5,.6}).cellKey==Key{7,7,9},"edge chooses first incident cell");
        require(g.embed({.5,.5,.5}).cellKey==Key{7,7,7},"corner chooses first incident cell");
        for(int a=0;a<3;++a) for(double boundary: {.125,.875}) {
            Eigen::Vector3d p=Eigen::Vector3d::Constant(.5); p(a)=boundary;
            g.embed(p); p(a)=std::nextafter(boundary,boundary==.125?0.:1.);
            rejects([&]{g.embed(p);},"one ULP outside must reject");
            p(a)=std::nextafter(boundary,.5); g.embed(p);
        }
        const Eigen::MatrixXd f=fields(g.nodePositions());
        for(Eigen::Vector3d p: {Eigen::Vector3d(.5,.45,.6),Eigen::Vector3d(.5,.5,.6),Eigen::Vector3d(.5,.5,.5)}) {
            Eigen::MatrixXd q(3,3); q.row(0)=p; p.x()=std::nextafter(.5,0.); q.row(1)=p; p.x()=std::nextafter(.5,1.); q.row(2)=p;
            const auto out=g.interpolate(f,q); for(int i=1;i<3;++i) require((out.row(i)-out.row(0)).cwiseAbs().maxCoeff()<1e-12,"continuous across shared faces/edges/corners");
        }
    });
    test("sparse_refinement_preserves_cavity_tunnel_and_disconnected_components",[&] {
        const SourceCageGrid hollow(cage(join(cube(.1,.9),cube(.3,.7),true)),config);
        const SourceCageGrid ring(cage(tunnel()),config);
        const SourceCageGrid islands(cage(join(cube(.03,.09),cube(.91,.97))),config);
        require(hollow.cells().size()==4032 && ring.cells().size()==3840 && islands.cells().size()==16,"exact coarse union refined eightfold");
        connectivity(hollow); connectivity(ring); connectivity(islands);
        rejects([&]{hollow.embed({.5,.5,.5});},"bounded cavity stays empty");
        for(double z:{0.,.5,1.}) rejects([&]{ring.embed({.5,.5,z});},"through-hole stays empty including open ends");
        rejects([&]{islands.embed({.5,.5,.5});},"component gap stays empty");
        islands.embed({.05,.05,.05}); islands.embed({.95,.95,.95});
        require(hollow.embed({.375,.5,.5}).cellKey==Key{5,7,7},"hole boundary chooses available incident cell");
        require(hollow.embed({.625,.5,.5}).cellKey==Key{10,7,7},"hole boundary skips absent lexicographically earlier cell");
    });
    test("tiny_sparse_level10_and_level20_avoid_root_lattice_expansion",[&] {
        auto small=config; small.fineLevel=10; small.maximumCells=8; small.maximumNodes=27; small.maximumEdges=54;
        SourceCageGrid s(cage(cube(100.25/512,100.75/512),9),small);
        require(s.cells().size()==8 && s.nodeKeys().size()==27 && s.edges().size()==54 && s.handleNodeIndices().size()==8,"sparse level10 literal counts");
        require(s.cells().front()==Key{200,200,200} && s.cells().back()==Key{201,201,201},"sparse level10 keys"); connectivity(s);
        small.fineLevel=20;
        resource([&]{SourceCageGrid huge(coarse,small);});
        SourceCageGrid max(cage(cube(100.25/524288,100.75/524288),19),small);
        require(max.cells().size()==8 && max.nodeKeys().size()==27 && max.edges().size()==54,"supported genuine sparse level20"); connectivity(max);
        reconstruct(max,max.nodePositions());
        std::cout<<"SYNTHETIC_REGRESSION sparse fineLevel20: 8 cells / 27 nodes / 54 edges; no optimization claim\n";
    });
    test("nondyadic_root_uses_FMA_shared_planes_and_actual_interval_lengths",[&] {
        SourceSurfaceCageConfig c; c.rootOrigin={-.3,-.2,-.1}; c.rootSide=1.7;
        Eigen::Vector3d lo,hi; for(int a=0;a<3;++a) {lo(a)=std::fma(.25,1.7,c.rootOrigin(a)); hi(a)=std::fma(.75,1.7,c.rootOrigin(a));}
        const auto mesh=box(lo,hi); const SourceCageGrid a(SourceSurfaceCage(mesh.v,mesh.f,c),config);
        require(a.nodeKeys().size()==2197,"nonempty affine grid"); connectivity(a);
        for(std::size_t i=0;i<a.nodeKeys().size();++i) for(int axis=0;axis<3;++axis)
            require(a.nodePositions()(i,axis)==std::fma(double(a.nodeKeys()[i][axis])/16,1.7,c.rootOrigin(axis)),"canonical FMA positions");
        Eigen::MatrixXd q(4,3); q<<.3,.45,.6,.5,.5,.5,.28125,.34375,.40625,.125,.875,.5;
        for(int i=0;i<q.rows();++i) for(int axis=0;axis<3;++axis) q(i,axis)=std::fma(q(i,axis),1.7,c.rootOrigin(axis));
        reconstruct(a,q); reconstruct(a,a.nodePositions()); reconstruct(a,mesh.v);
        require(a.embed(q.row(1).transpose()).cellKey==Key{7,7,7},"represented non-dyadic corner tie");
        for(int axis=0;axis<3;++axis) for(double t:{.125,.875}) {
            Eigen::Vector3d p=q.row(1).transpose(); double b=std::fma(t,1.7,c.rootOrigin(axis));
            p(axis)=std::nextafter(b,t==.125?-10.:10.); rejects([&]{a.embed(p);},"nondyadic nextafter outside rejects");
            p(axis)=std::nextafter(b,q(1,axis)); a.embed(p);
        }
    });
    for(int crossing=0;crossing<3;++crossing) test(
        crossing==0?"positive_binade_crossing_must_not_create_zero_width_fine_cells":
        crossing==1?"negative_binade_crossing_equal_ULP_requires_explicit_rejection":
                    "aligned_equal_ULP_root_documents_conservative_rejection",[&] {
        const double binade=std::ldexp(1.,48);
        SourceSurfaceCageConfig c;
        c.rootOrigin.x()=crossing==0?std::nextafter(binade,0.):crossing==1?-binade-.5:binade;
        Eigen::Vector3d lo,hi;
        for(int a=0;a<3;++a) {lo(a)=std::fma(.25,c.rootSide,c.rootOrigin(a)); hi(a)=std::fma(.75,c.rootSide,c.rootOrigin(a));}
        const auto mesh=box(lo,hi);
        // Construct outside the expected-rejection catch: an upstream failure
        // is a broken fixture, never evidence of the grid's behavioral RED.
        const SourceSurfaceCage upstream(mesh.v,mesh.f,c);
        require(upstream.boundaryAudit().facts().boundaryValid && upstream.selectedCells().size()==216,"binade fixture is an actual valid coarse source cage");
        for(int i=0;i<8;++i) require(std::fma(double(i)/8,1.,c.rootOrigin.x())<std::fma(double(i+1)/8,1.,c.rootOrigin.x()),"coarse planes are strictly increasing");
        int collapsed=0;
        for(int i=0;i<16;++i) if(std::fma(double(i)/16,1.,c.rootOrigin.x())==std::fma(double(i+1)/16,1.,c.rootOrigin.x())) ++collapsed;
        require(crossing==0?collapsed==7:collapsed==0,"independent FMA evidence distinguishes collapsed and otherwise valid equal-ULP grids");
        std::cout<<"binade fixture "<<crossing<<": upstream valid 216 cells; independently collapsed fine intervals="<<collapsed<<'\n';
        bool rejected=false;
        try {const SourceCageGrid invalid(upstream,config);}
        catch(const std::invalid_argument& e) {rejected=std::string(e.what()).find("numerically")!=std::string::npos;}
        require(rejected,"equal-ULP finer root must reject explicitly before publishing grid");
    });
    test("invalid_dimensions_nonfinite_inputs_and_each_resource_budget_reject",[&] {
        for(int level:{-1,0,3,21}) {auto c=config; c.fineLevel=level; rejects([&]{SourceCageGrid bad(coarse,c);},"fine level must be greater and at most20");}
        for(int which=0;which<4;++which) for(std::size_t value:{std::size_t(0),std::numeric_limits<std::size_t>::max()}) {
            auto c=config; if(which==0)c.maximumCells=value; if(which==1)c.maximumNodes=value; if(which==2)c.maximumEdges=value; if(which==3)c.maximumInterpolationValues=value;
            resource([&]{SourceCageGrid bad(coarse,c);});
        }
        for(int which=0;which<3;++which) {auto c=config; if(which==0)c.maximumCells=1727; if(which==1)c.maximumNodes=2196; if(which==2)c.maximumEdges=6083; resource([&]{SourceCageGrid bad(coarse,c);});}
        auto c=config; c.maximumInterpolationValues=7; const SourceCageGrid limited(coarse,c);
        Eigen::MatrixXd q(1,3); q<<.3,.45,.6; const auto f=fields(g.nodePositions());
        resource([&]{limited.interpolate(f,q);});
        rejects([&]{g.interpolate(Eigen::MatrixXd(1,8),q);},"mismatched node count");
        rejects([&]{g.interpolate(Eigen::MatrixXd(g.nodePositions().rows(),0),q);},"zero field columns");
        rejects([&]{g.interpolate(f,Eigen::MatrixXd(1,2));},"invalid query columns");
        auto bad=f; if(bad.size()) bad(0,0)=std::numeric_limits<double>::quiet_NaN();
        rejects([&]{g.interpolate(bad,q);},"NaN field");
        for(double v:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
            auto nonfinite=f; nonfinite(0,0)=v;
            rejects([&]{g.interpolate(nonfinite,Eigen::MatrixXd(0,3));},"nonfinite fields rejected even for empty queries");
            q(0,0)=v; rejects([&]{g.interpolate(f,q);},"nonfinite query"); rejects([&]{g.embed(q.row(0).transpose());},"nonfinite embedding");
        }
    });
    test("temporary_container_capacity_budget_rejects_before_construction",[&] {
        auto c=config; c.maximumEdges=std::set<SourceCageGrid::Edge>().max_size()+1;
        resource([&]{SourceCageGrid bad(coarse,c);});
    });
    test("copy_ownership_and_source_cage_config_output_identity_are_complete",[&] {
        auto mutableCage=coarse; auto mutableConfig=config; const SourceCageGrid owned(mutableCage,mutableConfig),again(coarse,config);
        mutableCage=cage(cube(.26,.74)); mutableConfig.fineLevel=5;
        require(owned.cage().sourceVertices()==source.v && owned.config().fineLevel==4,"copy-owned immutable cage and configuration");
        require(g.cells()==again.cells() && g.nodeKeys()==again.nodeKeys() && g.cellNodes()==again.cellNodes() && g.edges()==again.edges() && g.handleNodeIndices()==again.handleNodeIndices(),"byte-identical integer arrays and mapping");
        require(g.nodePositions().size()>0 && std::memcmp(g.nodePositions().data(),again.nodePositions().data(),g.nodePositions().size()*sizeof(double))==0,"byte-identical positions");
        require(g.sourceSha256()==coarse.sourceSha256() && g.cageConfigSha256()==coarse.configSha256() && g.cageOutputSha256()==coarse.outputSha256(),"upstream identities unchanged and exposed");
        require(g.configSha256().size()==64 && g.outputSha256().size()==64 && g.configSha256()==again.configSha256() && g.outputSha256()==again.outputSha256(),"stable complete hashes");
        SourceCageGrid different(mutableCage,config); require(g.sourceSha256()!=different.sourceSha256() && g.outputSha256()!=different.outputSha256(),"different source binds output despite identical cell domain");
        for(int which=0;which<5;++which) {
            auto c=config; if(which==0)++c.maximumCells; if(which==1)++c.maximumNodes; if(which==2)++c.maximumEdges; if(which==3)++c.maximumInterpolationValues; if(which==4)++c.fineLevel;
            SourceCageGrid d(coarse,c); require(d.configSha256()!=g.configSha256() && d.outputSha256()!=g.outputSha256(),"each config field bound to identity");
        }
        SourceSurfaceCageConfig cc; ++cc.maximumSelectedCells; SourceCageGrid d(SourceSurfaceCage(source.v,source.f,cc),config);
        require(d.cageConfigSha256()!=g.cageConfigSha256() && d.outputSha256()!=g.outputSha256(),"upstream cage configuration binds output");
        static_assert(std::is_same_v<decltype(g.nodePositions()),const Eigen::MatrixXd&>);
        static_assert(std::is_same_v<decltype(g.cage()),const SourceSurfaceCage&>);
    });
    std::cout<<"source cage grid SYNTHETIC_REGRESSION: "<<failures<<" failed groups; no BBW or physics authority\n";
    return failures?1:0;
}
