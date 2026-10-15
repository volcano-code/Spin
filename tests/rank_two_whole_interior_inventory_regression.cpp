#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoWholeInteriorInventory.h")
#include "spin/detail/RankTwoWholeInteriorInventory.h"
#include "fixtures/rank_two_source_cases.h"
namespace {
using Q=CGAL::Gmpq;using E=CORE::Expr;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code){try{f();}catch(const spin::RankTwoSourceError& e){require(e.code()==code,"wrong rejection class");return;}throw std::runtime_error("incomplete owner inventory accepted");}
void tests(const char* path){using namespace spin;
    const auto normal=PreparedRankTwoNormalMesh::fromExact(
        {{{1,0,0}},{{-1,0,0}},{{0,1,0}},{{0,-1,0}},{{0,0,1}},{{0,0,-1}}},
        {{{0,2,4}},{{2,1,4}},{{1,3,4}},{{3,0,4}},{{2,0,5}},{{1,2,5}},{{3,1,5}},{{0,3,5}}});
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"shared-seams-annulus-v1");
    const auto whole=analyzeRankTwoArrangement(source,Q(3,2),normal);
    const auto inventory=enumerateRankTwoWholeMovingInteriors(source,whole);
    const auto& native=wholeArrangementRecords(whole);
    require(inventory.sourceIdentity()==source.inputHash()&&inventory.arrangementIdentity()==whole.scientificHash()
        &&inventory.normalIdentity()==normal.inputHash(),"whole inventory has foreign source/grid authority");
    require(inventory.owners().size()==8&&inventory.pointCount()==24,"active owner omitted or inactive owner dropped");
    const int xsign[4]={-1,1,1,-1},ysign[4]={-1,-1,1,1};
    for(std::size_t owner=0;owner<8;++owner){const auto& local=inventory.owners()[owner];
        require(local.arrangementIdentity()==native.normalFaces[owner].scientificHash(),"owner order or native identity lost");
        // Annulus has x-events +/-3/2, +/-sqrt((3+sqrt(2))/2),
        // +/-sqrt(2), +/-sqrt((3-sqrt(2))/2),0. Each quadrant has
        // 3,2,1,0 gaps on its four occupied x-slabs; not the global13.
        require(local.slabCount()==8&&local.points().size()==(owner<4?6u:0u),"wrong complete normal-overlay gap census");
        std::map<std::pair<int,int>,std::size_t> counts;std::set<std::string> regions;
        for(const auto& p:local.points()){
            const auto x=p.coordinates()[0].expression(),y=p.coordinates()[1].expression();
            const E h=x*x+y*y-2,D=E(9)/4-x*x-y*y,F=x*x+y*y+h*h-E(9)/4;
            require(owner<4&&x.sign()==xsign[owner]&&y.sign()==ysign[owner],"point escaped its original normal owner");
            require(D.sign()>0&&(p.regionKind()==RankTwoRegionKind::Full?(F.sign()>0&&h.sign()<0):F.sign()<0),"wrong active native material signs");
            const E upper=CORE::sqrt(D),lower=p.regionKind()==RankTwoRegionKind::Full?-upper:h;
            require((upper-lower).sign()>0,"native owner point lacks positive thickness");
            require(p.definition().abscissaFraction==Q(1,2)&&p.definition().ordinateFraction==Q(1,2),"answer-selected interior fraction");
            counts[{static_cast<int>(p.regionKind()),h.sign()}]++;regions.insert(p.definition().regionIdentity);
        }
        if(owner<4)require(counts[std::make_pair(3,-1)]==1&&counts[std::make_pair(2,-1)]==2
            &&counts[std::make_pair(2,1)]==3&&regions.size()==3,"Full/Cut split or native region missing");
        else require(regions.empty(),"negative normal owner invented active material");
    }
    RankTwoInteriorInventorySettings limited;limited.maximumPoints=23;
    rejects([&]{enumerateRankTwoWholeMovingInteriors(source,whole,limited);},"RESOURCE_LIMIT");
    limited.maximumCandidateChecks=0;
    rejects([&]{enumerateRankTwoWholeMovingInteriors(source,whole,limited);},"RESOURCE_LIMIT");
    const auto foreign=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"foreign-frame");
    rejects([&]{enumerateRankTwoWholeMovingInteriors(foreign,whole);},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{enumerateRankTwoWholeMovingInteriors(source,native.global);},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{reloadRankTwoWholeMovingInteriors(source,whole,inventory.canonicalRecord()+"caller_passed true\n");},"INPUT_IDENTITY_MISMATCH");
    std::ofstream out(path);require(bool(out),"cannot create complete owner evidence");
    out<<"rank-two-whole-interior-inventory-test-v1\nINVENTORY "<<inventory.scientificHash()<<'\n'<<inventory.canonicalRecord();
    out.flush();require(bool(out),"incomplete owner evidence");
    std::cout<<"PASS original normal owners=8 active=4 inactive=4 points=24 total-point-cap/foreign-source/forged-reload rejected\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoWholeInteriorInventory.h")
    if(argc!=2)throw std::invalid_argument("whole inventory evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: complete original normal-owner native interior composition");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
