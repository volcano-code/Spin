#include "spin2014/SourceSurfaceCage.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <type_traits>
using namespace surface_test;
using namespace spin2014;
using Key=std::array<std::uint32_t,3>;

// Independent signed origin-tetra integration, not the cage/oracle integrator.
Ten signedMoments(const SourceSurfaceCage& cage) {
    Ten result=Ten::Zero();
    for(Eigen::Index i=0;i<cage.boundaryFaces().rows();++i) {
        Eigen::Matrix<double,3,3> p;
        for(int j=0;j<3;++j) p.row(j)=cage.boundaryVertices().row(cage.boundaryFaces()(i,j));
        const double volume=p.determinant()/6;
        const Eigen::Vector3d sum=p.colwise().sum();
        const Eigen::Matrix3d q=volume*(sum*sum.transpose()+p.transpose()*p)/20;
        Ten t; t<<volume,volume*sum/4,q(0,0),q(1,1),q(2,2),q(1,2),q(0,2),q(0,1);
        result+=t;
    }
    return result;
}
SourceSurfaceCage build(const Mesh& m,int level=3) {
    SourceSurfaceCageConfig c; c.cageLevel=level; return SourceSurfaceCage(m.v,m.f,c);
}
template<class F> void resourceRejects(F f) {
    bool rejected=false;
    try {f();} catch(const std::length_error& e) {
        rejected=std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos;
    }
    require(rejected,"exhausted budget must produce RESOURCE_LIMIT, no partial cage");
}
void audit(const SourceSurfaceCage& c,std::size_t cells,int handles,int triangles,int components) {
    require(c.selectedCells().size()==cells,"complete selected cell count");
    require(c.boundaryVertices().rows()==handles && c.handleVertexIndices().size()==std::size_t(handles),"every boundary corner is one handle, no sampling cap");
    require(c.boundaryFaces().rows()==triangles,"all boundary triangles retained");
    require(c.boundaryAudit().facts().boundaryValid,"independent exact boundary audit");
    require(c.boundaryAudit().facts().componentVolumeSigns.size()==std::size_t(components),"all boundary components retained");
    require(std::is_sorted(c.selectedCells().begin(),c.selectedCells().end()),"canonical selected keys");
    require(std::adjacent_find(c.selectedCells().begin(),c.selectedCells().end())==c.selectedCells().end(),"unique selected keys");
    require(std::is_sorted(c.boundaryVertexKeys().begin(),c.boundaryVertexKeys().end()),"canonical vertex keys");
    std::set<Key> uniqueKeys(c.boundaryVertexKeys().begin(),c.boundaryVertexKeys().end());
    require(uniqueKeys.size()==std::size_t(handles),"shared corners cannot be duplicated");
    std::map<std::array<int,2>,int> edges;
    std::set<std::array<int,3>> faces;
    for(int i=0;i<triangles;++i) {
        std::array<int,3> f{};
        for(int j=0;j<3;++j) {
            f[j]=c.boundaryFaces()(i,j);
            std::array<int,2> e{f[j],c.boundaryFaces()(i,(j+1)%3)};
            std::sort(e.begin(),e.end()); ++edges[e];
        }
        std::sort(f.begin(),f.end()); require(faces.insert(f).second,"no duplicate/internal face triangles");
    }
    for(auto e:edges) require(e.second==2,"two triangles per boundary edge");
    for(int i=0;i<handles;++i) require(c.handleVertexIndices()[i]==i,"complete canonical handle mapping");
}
// Analytic square annulus extrusion. Four ring corners counterclockwise;
// outward outer wall, inward tunnel wall, and oppositely oriented end annuli.
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
int main() {
    try {
        // Catches boundary-only selection, lost Full descendants, wrong closed
        // contact ties, internal faces, smoothing and the legacy 16-handle cap.
        auto source=cube(.25,.75); const Mesh saved=source;
        SourceSurfaceCageConfig config; config.cageLevel=3;
        SourceSurfaceCage c(source.v,source.f,config); audit(c,216,218,432,1);
        std::set<Key> expected;
        for(unsigned x=1;x<=6;++x) for(unsigned y=1;y<=6;++y) for(unsigned z=1;z<=6;++z) expected.insert({x,y,z});
        require(std::set<Key>(c.selectedCells().begin(),c.selectedCells().end())==expected,"literal closed-contact indices 1..6 on each axis");
        require(c.boundaryVertices().minCoeff()==.125 && c.boundaryVertices().maxCoeff()==.875,"exact level3 envelope extrema");
        close(signedMoments(c),boxMoments({.125,.125,.125},{.875,.875,.875}),"all ten independent cuboid moments, volume 27/64");
        require(source.v==saved.v && source.f==saved.f,"construction must not mutate source");
        auto two=build(source,2); audit(two,64,98,192,1);
        close(signedMoments(two),boxMoments({0,0,0},{1,1,1}),"level2 root envelope");
        auto zero=build(source,0); audit(zero,1,8,12,1);
        // Catches nondeterministic iteration, partial hashes and borrowed input.
        auto same=build(source);
        require(c.selectedCells()==same.selectedCells() && c.boundaryVertexKeys()==same.boundaryVertexKeys() && c.handleVertexIndices()==same.handleVertexIndices(),"deterministic keys and handles");
        require(std::memcmp(c.boundaryVertices().data(),same.boundaryVertices().data(),sizeof(double)*c.boundaryVertices().size())==0 && std::memcmp(c.boundaryFaces().data(),same.boundaryFaces().data(),sizeof(int)*c.boundaryFaces().size())==0,"byte-identical geometry");
        require(c.sourceSha256()==same.sourceSha256() && c.configSha256()==same.configSha256() && c.outputSha256()==same.outputSha256(),"deterministic identity");
        require(c.sourceSha256()==two.sourceSha256() && c.configSha256()!=two.configSha256() && c.outputSha256()!=two.outputSha256(),"level identity");
        auto changed=build(cube(.26,.74));
        require(c.sourceSha256()!=changed.sourceSha256() && c.outputSha256()!=changed.outputSha256(),"source bound into output identity even with same selected cells");
        for(int b=0;b<4;++b) {
            auto v=config;
            if(b==0) ++v.maximumSourceNodes; if(b==1) ++v.maximumSelectedCells;
            if(b==2) ++v.maximumBoundaryTriangles; if(b==3) ++v.maximumBoundaryVertices;
            SourceSurfaceCage budget(source.v,source.f,v);
            require(c.sourceSha256()==budget.sourceSha256() && c.configSha256()!=budget.configSha256() && c.outputSha256()!=budget.outputSha256(),"each resource budget in identity");
        }
        source.v.setZero(); source.f.setZero(); config.rootSide=9;
        require(c.sourceVertices()==saved.v && c.sourceFaces()==saved.f && c.config().rootSide==1,"copy-owned source/config");
        static_assert(std::is_same_v<decltype(c.boundaryVertices()),const Eigen::MatrixXd&>);
        static_assert(std::is_same_v<decltype(c.selectedCells()),const std::vector<Key>&>);
        // Catches dense root traversal: only 9 branches, 73 stored oracle nodes.
        const auto tiny=cube(100.25/512,100.75/512);
        SourceSurfaceCageConfig sparse; sparse.cageLevel=9; sparse.maximumSourceNodes=73;
        sparse.maximumSelectedCells=1; sparse.maximumBoundaryVertices=8; sparse.maximumBoundaryTriangles=12;
        SourceSurfaceCage level9(tiny.v,tiny.f,sparse); audit(level9,1,8,12,1);
        require(level9.selectedCells().front()==Key{100,100,100} && level9.sourceNodeCount()==73,"genuine sparse level9 geometry");
        close(signedMoments(level9),boxMoments(Eigen::Vector3d::Constant(100./512),Eigen::Vector3d::Constant(101./512)),"sparse cube moments");
        std::cout<<"level9 GEOMETRY: nodes="<<level9.sourceNodeCount()<<" selected=1 handles=8 triangles=12\n";
        // Catches cavity filling, wrong cavity winding, largest-shell culling.
        auto hollow=build(join(cube(.1,.9),cube(.3,.7),true)); audit(hollow,504,412,816,2);
        const auto& signs=hollow.boundaryAudit().facts().componentVolumeSigns;
        require(std::count(signs.begin(),signs.end(),-1)==1,"bounded cavity must have inward boundary");
        close(signedMoments(hollow),boxMoments({0,0,0},{1,1,1})-boxMoments({.375,.375,.375},{.625,.625,.625}),"cavity envelope volume 63/64");
        // 768 outer triangles -16 open endpatch triangles +128 tunnel walls.
        // 386 outer corners -2 hole interior endcorners +8*7 inner wall corners.
        auto torus=build(tunnel()); audit(torus,480,440,880,1);
        close(signedMoments(torus),boxMoments({0,0,0},{1,1,1})-boxMoments({.375,.375,0},{.625,.625,1}),"through-hole volume 15/16");
        require(torus.boundaryVertices().rows()-torus.boundaryFaces().rows()/2==0,"through-hole genus one");
        auto islands=build(join(cube(.03,.09),cube(.91,.97))); audit(islands,2,16,24,2);
        require(std::abs(signedMoments(islands)(0)-1./256)<1e-15,"disconnected union volume");
        // Catches silently repaired or usable singular unions. Original boxes
        // are disjoint, but their selected cells touch only at edge or vertex.
        for(const Key key: {Key{2,2,1},Key{2,2,2}}) {
            auto a=cube(1.25/8,1.75/8);
            Eigen::Vector3d lo,hi; for(int j=0;j<3;++j) {lo(j)=(key[j]+.25)/8; hi(j)=(key[j]+.75)/8;}
            auto m=join(a,box(lo,hi)); bool rejected=false;
            try {build(m);} catch(const SourceSurfaceCageGeometryError& e) {
                rejected=!e.boundaryAudit().facts().boundaryValid && std::string(e.what()).find("INVALID_GEOMETRY")!=std::string::npos;
            }
            require(rejected,"singular edge/vertex union must fail with exact audit diagnostics");
        }
        // Face, edge and corner ties in [.125,.25]^3 select all 3^3 cells.
        audit(build(cube(.125,.25)),27,56,108,1);
        // Catches inconsistent floating planes under affine non-dyadic roots.
        SourceSurfaceCageConfig affine; affine.cageLevel=3;
        affine.rootOrigin={-.3,-.2,-.1}; affine.rootSide=1.7;
        Eigen::Vector3d lo,hi;
        for(int j=0;j<3;++j) {lo(j)=std::fma(.25,1.7,affine.rootOrigin(j)); hi(j)=std::fma(.75,1.7,affine.rootOrigin(j));}
        auto m=box(lo,hi); SourceSurfaceCage shifted(m.v,m.f,affine); audit(shifted,216,218,432,1);
        for(int i=0;i<shifted.boundaryVertices().rows();++i) for(int j=0;j<3;++j)
            require(shifted.boundaryVertices()(i,j)==std::fma(double(shifted.boundaryVertexKeys()[i][j])/8,1.7,affine.rootOrigin(j)),"canonical FMA corner coordinates");
        auto rootChanged=SourceSurfaceCageConfig{}; rootChanged.cageLevel=3; rootChanged.rootSide=2;
        SourceSurfaceCage otherRoot(saved.v,saved.f,rootChanged);
        require(otherRoot.sourceSha256()==c.sourceSha256() && otherRoot.configSha256()!=c.configSha256(),"root side identity");
        rootChanged.rootSide=1; rootChanged.rootOrigin.setConstant(-.1);
        SourceSurfaceCage otherOrigin(saved.v,saved.f,rootChanged);
        require(otherOrigin.configSha256()!=c.configSha256(),"root origin identity");
        // Catches budget truncation/fallback and unchecked Full expansion.
        config=SourceSurfaceCageConfig{}; config.cageLevel=3;
        for(int b=0;b<4;++b) for(bool zeroBudget: {false,true}) {
            auto v=config;
            if(b==0) v.maximumSourceNodes=zeroBudget?0:1;
            if(b==1) v.maximumSelectedCells=zeroBudget?0:215;
            if(b==2) v.maximumBoundaryVertices=zeroBudget?0:217;
            if(b==3) v.maximumBoundaryTriangles=zeroBudget?0:431;
            resourceRejects([&]{SourceSurfaceCage fail(saved.v,saved.f,v);});
        }
        auto full=cube(0,1); auto low=config; low.cageLevel=5; low.maximumSelectedCells=1;
        resourceRejects([&]{SourceSurfaceCage fail(full.v,full.f,low);});
        // The exact audit uses signed-int face ids and at most INT_MAX/2
        // triangles. An unrepresentable budget must fail before extraction.
        auto oversized=config; oversized.maximumBoundaryTriangles=std::numeric_limits<int>::max();
        resourceRejects([&]{SourceSurfaceCage fail(saved.v,saved.f,oversized);});
        for(int b=0;b<3;++b) {
            auto v=config;
            if(b==0) v.maximumSelectedCells=std::numeric_limits<std::size_t>::max();
            if(b==1) v.maximumBoundaryVertices=std::numeric_limits<std::size_t>::max();
            if(b==2) v.maximumSourceNodes=std::numeric_limits<std::size_t>::max();
            resourceRejects([&]{SourceSurfaceCage fail(saved.v,saved.f,v);});
        }
        for(int depth: {-1,21}) {auto v=config; v.cageLevel=depth; rejects([&]{SourceSurfaceCage fail(saved.v,saved.f,v);},"invalid depth");}
        for(double side: {0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
            auto v=config; v.rootSide=side; rejects([&]{SourceSurfaceCage fail(saved.v,saved.f,v);},"invalid root side");
        }
        auto bad=config; bad.rootOrigin.x()=std::numeric_limits<double>::quiet_NaN();
        rejects([&]{SourceSurfaceCage fail(saved.v,saved.f,bad);},"invalid root origin");
        auto nan=saved; nan.v(0,0)=std::numeric_limits<double>::quiet_NaN(); rejects([&]{build(nan);},"nonfinite source");
        auto invalid=saved; invalid.f(0,0)=99; rejects([&]{build(invalid);},"invalid source index");
        auto open=saved; open.f.conservativeResize(11,3); rejects([&]{build(open);},"open source");
        rejects([&]{build(Mesh{});},"empty source");
        rejects([&]{build(cube(-.1,.5));},"source outside root");
        const auto microscopic=cube(100.25/1048576,100.75/1048576);
        auto maxDepth=sparse; maxDepth.cageLevel=20; maxDepth.maximumSourceNodes=161;
        audit(SourceSurfaceCage(microscopic.v,microscopic.f,maxDepth),1,8,12,1);
        std::cout<<"source surface cage regression passed; complete geometry only, no BBW/physics authority\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n'; return 1;}
}
