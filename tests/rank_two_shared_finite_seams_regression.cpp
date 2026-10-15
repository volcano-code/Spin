#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoSharedFiniteSeams.h")
#include "spin/detail/RankTwoSharedFiniteSeams.h"
#include "fixtures/rank_two_source_cases.h"
#include <algorithm>
#include <sstream>
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){bool caught=false;
    try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}
    require(caught,"invalid seam request was not rejected with required classification");}
spin::PreparedRankTwoNormalMesh octa(bool reverse=false){
    spin::PreparedRankTwoNormalMesh::Rays rays={{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}};
    spin::PreparedRankTwoNormalMesh::Faces faces={{0,2,4},{2,1,4},{1,3,4},{3,0,4},{2,0,5},{1,2,5},{3,1,5},{0,3,5}};
    if(reverse)std::reverse(faces.begin(),faces.end());
    return spin::PreparedRankTwoNormalMesh::fromExact(rays,faces);
}
void run(const char* path){using namespace spin;
    const auto s=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"shared-seams-annulus-v1");
    const auto a=analyzeRankTwoArrangement(s,Q(3,2),octa());
    const auto m=buildRankTwoSharedFiniteSeams(s,a);
    require(m.points().size()==33&&m.intervals().size()==20,"annulus finite seams lost roots/empty intervals or duplicated north pole");
    std::size_t finite=0,moving=0;
    for(const auto& p:m.points())if(!p.endpoint){++moving;require(p.fraction==Q(1,2),"moving fraction altered");
        require(p.parentIdentities[0]!=p.parentIdentities[1],"moving point has no actual endpoint pair");}
    require(moving==16,"not one actual moving sample per finite interval");
    for(const auto& i:m.intervals()){
        if(!i.finite){require(i.points.empty(),"outside projection was given fake finite points");continue;}
        ++finite;require(i.points.size()==3,"shared interval does not own both ends plus actual interior sample");
        for(int side=0;side<2;++side){auto indices=m.consumerPoints(i.identity,side);
            auto expected=i.points;if(i.sourceReversed[side])std::reverse(expected.begin(),expected.end());
            require(indices==expected,"incident consumer did not use SAME master point indices in native order");}
        const auto& p=m.points()[i.points[1]];
        const auto x=p.coordinates[0].expression(),y=p.coordinates[1].expression();
        require(x.sign()==0||y.sign()==0,"octa seam midpoint left actual supporting normal line");
        const auto r=x*x+y*y,h=r-CORE::Expr(2);
        require((r+h*h-CORE::Expr(CORE::BigRat(9,4))).sign()==i.originalSigns[0]
            &&h.sign()==i.originalSigns[1]&&(CORE::Expr(CORE::BigRat(9,4))-r).sign()==i.originalSigns[2],
            "master sample does not preserve actual F/h/D signs");
    }
    require(finite==16,"annulus finite normal interval cover incomplete");
    const auto other=analyzeRankTwoArrangement(s,Q(3,2),octa(true));
    const auto reordered=buildRankTwoSharedFiniteSeams(s,other);
    require(m.points().size()==reordered.points().size(),"face-order duplicated shared points");
    for(std::size_t k=0;k<m.points().size();++k){const auto& p=m.points()[k];const auto& q=reordered.points()[k];
        require(p.identity==q.identity,"shared node identity depends on incident face order");
        for(int d=0;d<2;++d)require(p.coordinates[d].expression()==q.coordinates[d].expression(),"shared root changed on face reorder");}
    const auto saddle=PreparedRankTwoSource::fromExact({0,0,0,1,0,0,0,0,0,Q(-1,2)},"shared-seams-saddle-v1");
    const auto b=analyzeRankTwoArrangement(saddle,Q(1),octa());
    const auto n=buildRankTwoSharedFiniteSeams(saddle,b);
    require(n.points().size()==9,"saddle coincident original curves lost shared point");
    std::size_t hOwner=0;
    for(const auto& i:n.intervals())if(i.finite){require(i.originalSigns[1]==0,"saddle h coincident line sign lost");
        for(const auto& o:i.originalOwners)if(o.polynomial=="h"&&o.multiplicity==1)++hOwner;}
    require(hOwner==4,"coincident h multiplicity dropped");
    rejects([&]{buildRankTwoSharedFiniteSeams(s,b);},"INPUT_IDENTITY_MISMATCH");
    for(const auto fractions:std::vector<std::vector<Q>>{{Q(0)},{Q(1)},{Q(1,2),Q(1,2)},{Q(3,4),Q(1,4)}}){
        RankTwoSharedFiniteSeamSettings bad;bad.fractions=fractions;
        rejects([&]{buildRankTwoSharedFiniteSeams(s,a,bad);},"INPUT_IDENTITY_MISMATCH");}
    RankTwoSharedFiniteSeamSettings tiny;tiny.maximumPoints=1;
    rejects([&]{buildRankTwoSharedFiniteSeams(s,a,tiny);},"RESOURCE_LIMIT");
    tiny={};tiny.maximumChecks=0;rejects([&]{buildRankTwoSharedFiniteSeams(s,a,tiny);},"RESOURCE_LIMIT");
    std::ofstream out(path);require(bool(out),"cannot save actual shared seam evidence");
    out<<"rank-two-shared-finite-seams-test-v1\nCASE annulus\n"<<m.canonicalRecord()
        <<"CASE_END\nCASE reordered\n"<<reordered.canonicalRecord()<<"CASE_END\nCASE saddle\n"<<n.canonicalRecord()<<"CASE_END\n";
    out.flush();require(bool(out),"incomplete shared seam evidence");
    std::cout<<"PASS annulus33points16finite20total; saddle9points4coincident; both native consumers; negative gates\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoSharedFiniteSeams.h")
    if(argc!=2)throw std::invalid_argument("seam evidence path required");run(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: shared finite normal-edge master point table and native consumers");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
