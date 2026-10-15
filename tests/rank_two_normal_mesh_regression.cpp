#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#if __has_include("spin/detail/RankTwoNormalMesh.h")
#include "spin/detail/RankTwoNormalMesh.h"
#include "verification/Sha256.h"
#if __has_include("spin/detail/FixedNormalSphere.h")
#include "spin/detail/FixedNormalSphere.h"
#endif

namespace {
void require(bool v,const char* message){if(!v)throw std::runtime_error(message);}
using Q=CGAL::Gmpq;
using Rays=std::vector<std::array<Q,3>>;
using Faces=std::vector<std::array<std::size_t,3>>;
void tests(const char* input,const char* output){
    using namespace spin;
    const Rays octahedron={{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}};
    const Faces faces={{{0,2,4}},{{2,1,4}},{{1,3,4}},{{3,0,4}},
                       {{2,0,5}},{{1,2,5}},{{3,1,5}},{{0,3,5}}};
    const auto valid=PreparedRankTwoNormalMesh::fromExact(octahedron,faces);
    require(valid.edges().size()==12 && valid.checkedFacePairs()==28,"incomplete octahedron admission");
    require(valid.face(3).vertexIds==faces[3],"normal indices changed");
    auto reject=[&](const Rays& rays,const Faces& triangles,const char* kind){
        bool refused=false;
        try {PreparedRankTwoNormalMesh::fromExact(rays,triangles);}
        catch(const RankTwoSourceError& e){
            refused=e.code()=="SOURCE_TOPOLOGY_INVALID";
            if(kind)refused=refused && e.witness().find(kind)!=std::string::npos;
        }
        require(refused,"invalid normal mesh not rejected with witness");
    };
    auto open=faces;open.pop_back();reject(octahedron,open,"edge");
    auto flipped=faces;std::swap(flipped[0][1],flipped[0][2]);reject(octahedron,flipped,"orientation");
    auto duplicate=faces;duplicate.push_back(faces[0]);reject(octahedron,duplicate,"duplicate face");
    auto zero=octahedron;zero[0]={0,0,0};reject(zero,faces,"zero ray");
    auto unused=octahedron;unused.push_back({1,1,1});reject(unused,faces,"unused vertex");
    reject(octahedron,Faces{{{0,2,4}},{{0,4,2}}},nullptr);
    // A combinatorially valid pentagonal bipyramid, but the equator is a
    // pentagram. All ten triangle determinants are positive and all links
    // are simple cycles. It covers spherical directions TWICE. Edge counts,
    // determinants and Euler characteristic alone must not admit it.
    const Rays twice={{{0,0,1}},{{0,0,-1}},{{3,0,0}},{{-3,2,0}},{{1,-3,0}},{{1,3,0}},{{-3,-2,0}}};
    Faces winding;
    for(std::size_t i=0;i<5;++i){const auto a=2+i,b=2+(i+1)%5;winding.push_back({a,b,0});winding.push_back({b,a,1});}
    reject(twice,winding,"cone intersection");
    auto scaled=octahedron;for(std::size_t i=0;i<scaled.size();++i)for(auto& q:scaled[i])q*=Q(static_cast<int>(i+2));
    const auto rescaled=PreparedRankTwoNormalMesh::fromExact(scaled,faces);
    require(rescaled.inputHash()!=valid.inputHash(),"raw normal authority not bound");
    require(rescaled.edges()==valid.edges(),"positive ray scaling changed topology");

    std::ifstream stream(input,std::ios::binary);require(bool(stream),"frozen normal input missing");
    std::string bytes((std::istreambuf_iterator<char>(stream)),std::istreambuf_iterator<char>());
    verification::Sha256 digest;digest.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    require(digest.finish()=="5ad29232268675e93c7f8be0de0709fefbb242310ed7ca2cfc7797b2ada21140",
        "frozen normal artifact identity mismatch");
    std::istringstream rows(bytes);rows.imbue(std::locale::classic());std::string line;Rays rays;Faces triangles;
    while(std::getline(rows,line)){
        std::istringstream row(line);row.imbue(std::locale::classic());char kind;
        require(bool(row>>kind),"empty normal input record");
        if(kind=='v'){
            double x,y,z;require(bool(row>>x>>y>>z) && std::isfinite(x) && std::isfinite(y) && std::isfinite(z),"bad saved ray");
            rays.push_back({Q(x),Q(y),Q(z)}); // stored binary64 authority, NOT decimal rationals
        }else if(kind=='f'){
            std::size_t a,b,c;require(bool(row>>a>>b>>c),"bad normal face");triangles.push_back({a,b,c});
        }else throw std::runtime_error("unknown frozen normal record");
        std::string extra;require(!(row>>extra),"extra frozen normal tokens");
    }
    const auto frozen=PreparedRankTwoNormalMesh::fromExact(rays,triangles);
#if __has_include("spin/detail/FixedNormalSphere.h")
    const auto generated=detail::fixedNormalSphere(2);
    require(generated.vertices.size()==rays.size() && generated.faces.size()==triangles.size(),"shared generator changed grid size");
    for(std::size_t i=0;i<rays.size();++i)for(int d=0;d<3;++d)
        require(Q(generated.vertices[i][d])==rays[i][d],"shared generator changed stored binary64 normal");
    for(std::size_t i=0;i<triangles.size();++i)for(int d=0;d<3;++d)
        require(std::size_t(generated.faces[i][d])==triangles[i][d],"shared generator changed face indices");
#else
    throw std::runtime_error("MISSING_FEATURE: shared original fixed-normal generator");
#endif
    require(frozen.rays().size()==162 && frozen.faces().size()==320 && frozen.edges().size()==480,
        "wrong frozen level-2 mesh");
    require(frozen.checkedFacePairs()==51040,"not every actual normal cone pair checked");
    const auto reloaded=PreparedRankTwoNormalMesh::reload(rays,triangles,frozen.canonicalRecord());
    require(reloaded.inputHash()==frozen.inputHash(),"frozen admission reload mismatch");
    bool failed=false;
    try {PreparedRankTwoNormalMesh::reload(rays,triangles,frozen.canonicalRecord()+"passed=true\n");}
    catch(const RankTwoSourceError& e){failed=e.code()=="INPUT_IDENTITY_MISMATCH";}
    require(failed,"forged whole-normal admission accepted");
    std::ifstream previous(output);require(!previous.good(),"refuse overwrite normal evidence");
    std::ofstream saved(output);require(bool(saved),"cannot save normal evidence");
    saved<<valid.canonicalRecord()<<frozen.canonicalRecord();saved.close();require(bool(saved),"normal save failed");
    std::cout<<"PASS actual normal sphere: 162 rays,320 faces,480 edges,51040 exact cone-pair checks\n"
        <<frozen.inputHash()<<'\n';
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoNormalMesh.h")
    if(argc!=3)throw std::invalid_argument("frozen normals and new evidence path required");
    tests(argv[1],argv[2]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: exact whole normal sphere admission");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
