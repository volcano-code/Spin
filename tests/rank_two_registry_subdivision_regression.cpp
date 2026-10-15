#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoArcRegistrySubdivisions.h")
#include "spin/detail/RankTwoArcRegistrySubdivisions.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include <algorithm>
#include <locale>
#include <sstream>
#include <set>
namespace {
using Q=CGAL::Gmpq;using C=spin::detail::CertifiedAlgebraicReal;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
template<class F>void rejects(F f,const char* code,const char* why){
    bool caught=false;try{f();}catch(const spin::RankTwoSourceError& e){caught=e.code()==code;}require(caught,why);
}
const spin::RankTwoRegistryPoint& point(const spin::RankTwoArcRegistry& r,const std::string& id){
    for(const auto& p:r.points())if(p.identity==id)return p;throw std::runtime_error("point absent");
}
void tests(const char* path){
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact({1,0,-1,1,-1,-1,0,0,0,0},"endpoint-neighborhood-method-v1");
    const auto global=analyzeRankTwoGlobalArrangement(source,Q(1));
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,-1,1},{-3,-1,1},{0,-3,1}}};
    const auto a=analyzeRankTwoNormalFace(source,global,cone);
    std::string event;RankTwoArcSamplingPlan plan;
    for(const auto& v:arrangementRecords(a).vertices)
        if(C::reload(v.coordinates[0]).expression()==1 && C::reload(v.coordinates[1]).expression()==0)event=v.identity;
    for(const auto& arc:arrangementRecords(a).arcs)
        if(arc.ends[0].location==RankTwoArcEndLocation::Finite && arc.ends[1].location==RankTwoArcEndLocation::Finite)
            plan.emplace(arc.identity,std::vector<Q>{Q(1,2)});
    const auto n=certifyRankTwoEndpointNeighborhood(a,event);
    const auto base=buildRankTwoArcRegistryWithNeighborhoods(a,plan,{n});
    // Frozen positive-y h graph: x^2+y^2+y-1=0. Select its actual regular
    // gate-to-midpoint edge, never an endpoint-fan edge or a success-tuned input.
    const RankTwoRegistryEdge* selected=nullptr;
    for(const auto& e:base.edges()){
        if(!e.endpointNeighborhoodIdentity.empty())continue;
        bool h=false;for(const auto& owner:e.owners)h|=owner.polynomial=="h";
        const auto& l=base.points()[e.points[0]],&r=base.points()[e.points[1]];
        if(h && (l.kind==RankTwoRegistryPointKind::SharedGate || r.kind==RankTwoRegistryPointKind::SharedGate)
            && l.coordinates[1].expression()>0 && r.coordinates[1].expression()>0){
            require(!selected,"frozen h edge selection is ambiguous");selected=&e;
        }
    }
    require(selected,"frozen h gate-to-midpoint edge absent");
    const auto& left=base.points()[selected->points[0]],&right=base.points()[selected->points[1]];
    RankTwoArcSubdivision first{selected->arcIdentity,{left.identity,right.identity},Q(1,2)};
    const auto one=buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{first});
    const RankTwoRegistryPoint* middle=nullptr;
    for(const auto& p:one.points())if(p.kind==RankTwoRegistryPointKind::ArcSubdivision){require(!middle,"extra subdivision");middle=&p;}
    require(middle && middle->endpointIdentities==first.pointIdentities && middle->fraction==Q(1,2),
        "actual parent DAG or fraction lost");
    const auto x=middle->coordinates[0].expression(),y=middle->coordinates[1].expression();
    require(x==(left.coordinates[0].expression()+right.coordinates[0].expression())/2,
        "subdivision uses interval midpoint rather than actual parent expressions");
    require((x*x+y*y+y-1).sign()==0 && y>0,"subdivision left actual h graph or selected another fiber branch");
    require(y!=(left.coordinates[1].expression()+right.coordinates[1].expression())/2,
        "new curve point was substituted by linear chord interpolation");
    RankTwoArcSubdivision second{first.arcIdentity,{middle->identity,right.identity},Q(1,2)};
    const auto result=buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{first,second});
    const auto reordered=buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{second,first});
    require(result.scientificHash()==reordered.scientificHash() && result.canonicalRecord()==reordered.canonicalRecord(),
        "subdivision request order changed scientific DAG or indices");
    require(result.points().size()==53 && result.edges().size()==68 && result.subdivisions().size()==2,
        "nested shared subdivision was omitted or duplicated");
    std::size_t bound=0;for(const auto& e:result.edges())bound+=!e.endpointNeighborhoodIdentity.empty();
    require(bound==8,"refinement lost complete endpoint evidence");
    const RankTwoRegistryPoint* nested=nullptr;
    for(const auto& p:result.points())if(p.kind==RankTwoRegistryPointKind::ArcSubdivision && p.identity!=middle->identity)nested=&p;
    require(nested && nested->coordinates[0].expression()==(x+right.coordinates[0].expression())/2,
        "nested subdivision reused original arc endpoints");
    for(const auto* p:{middle,nested}){
        const auto u=p->coordinates[0].expression(),v=p->coordinates[1].expression();
        require((u*u+v*v+v-1).sign()==0 && v>0,"nested graph membership failed");
    }
    // Original directed rings must consume exactly the same newly split chain.
    std::size_t consumers=0;
    for(const auto& region:result.regions())for(bool inner:{false,true}){
        const auto& loops=inner?region.innerLoops:region.outerLoops;
        for(std::size_t k=0;k<loops.size();++k){
            if(!loops[k].finite)continue;
            const auto indices=result.loopVertices(region.identity,inner,k);
            consumers+=std::count(indices.begin(),indices.end(),std::size_t(&point(result,nested->identity)-result.points().data()));
        }
    }
    require(consumers==2,"adjacent regions did not consume one shared nested point");
    auto bad=first;bad.pointIdentities[0]="missing";
    rejects([&]{buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{bad});},"INPUT_IDENTITY_MISMATCH","unknown parent accepted");
    bad=first;std::swap(bad.pointIdentities[0],bad.pointIdentities[1]);
    rejects([&]{buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{bad});},"INPUT_IDENTITY_MISMATCH","reversed monotone parents accepted");
    bad=first;bad.fraction=0;
    rejects([&]{buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{bad});},"INPUT_IDENTITY_MISMATCH","endpoint alias accepted");
    bad=first;for(const auto& arc:base.arcs())if(arc.identity!=first.arcIdentity && arc.finite){bad.arcIdentity=arc.identity;break;}
    rejects([&]{buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{bad});},"INPUT_IDENTITY_MISMATCH","foreign-arc parents accepted");
    rejects([&]{buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{first,first});},"INPUT_IDENTITY_MISMATCH","duplicate dependency accepted");
    rejects([&]{buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{second});},"INPUT_IDENTITY_MISMATCH","missing DAG parent accepted");
    for(const auto& e:base.edges())if(!e.endpointNeighborhoodIdentity.empty()){
        RankTwoArcSubdivision splitFan{e.arcIdentity,{base.points()[e.points[0]].identity,base.points()[e.points[1]].identity},Q(1,2)};
        rejects([&]{buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{splitFan});},"ARRANGEMENT_UNRESOLVED",
            "split invalidated endpoint chord but inherited old proof");break;
    }
    RankTwoSourceLimits tiny;tiny.maximumVertices=52;
    rejects([&]{buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{first,second},tiny);},"RESOURCE_LIMIT","DAG vertices ignored budget");
    const auto empty=buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},{});
    require(empty.scientificHash()==base.scientificHash() && empty.canonicalRecord()==base.canonicalRecord(),
        "empty subdivision changed historical registry-v2");
    require(result.canonicalRecord().find("rank-two-shared-arc-registry-v3\n")==0,
        "DAG semantics lack explicit evidence version");
    // All these coordinates are rational affine combinations of TWO frozen
    // source projections. Treating sibling nodes as independent algebraic
    // roots spuriously multiplies eliminant degrees and rejects this DAG.
    auto deep=one;std::vector<RankTwoArcSubdivision> dag{first};std::string center=middle->identity;
    Q centerFraction(1,2);
    const auto add=[&](const RankTwoArcSubdivision& request){
        std::set<std::string> old;for(const auto& p:deep.points())old.insert(p.identity);
        dag.push_back(request);deep=buildRankTwoArcRegistryWithSubdivisions(a,plan,{n},dag);
        std::string added;for(const auto& p:deep.points())if(!old.count(p.identity)){
            require(added.empty(),"DAG request published multiple points");added=p.identity;
        }
        require(!added.empty(),"DAG request lost its point");return added;
    };
    for(int depth=0;depth<6;++depth){
        const auto low=add({first.arcIdentity,{left.identity,center},Q(1,2)});
        const auto high=add({first.arcIdentity,{center,right.identity},Q(1,2)});
        center=add({first.arcIdentity,{low,high},Q(1,3)});
        centerFraction=centerFraction/2+Q(1,6);
        const auto& p=point(deep,center);
        std::ostringstream fractionText;fractionText.imbue(std::locale::classic());fractionText<<centerFraction;
        const CORE::Expr weight(CORE::BigRat(fractionText.str().c_str()));
        require(p.coordinates[0].expression()==(1-weight)
            *left.coordinates[0].expression()+weight
            *right.coordinates[0].expression(),"deep sibling dependency changed exact coordinate");
        require(p.coordinates[0].certificate().coefficients.size()<=5,
            "two-anchor affine dependency expanded into fictitious independent conjugates");
    }
    require(deep.subdivisions().size()==19,"deep DAG requests omitted");
    std::ofstream out(path);require(bool(out),"cannot save subdivision evidence");
    out<<"rank-two-registry-subdivision-test-v1\nREGISTRY "<<result.scientificHash()<<'\n'<<result.canonicalRecord();
    out.flush();require(bool(out),"incomplete subdivision evidence");
    std::cout<<"PASS actual gate-to-sample nested curve subdivision, original shared consumers, preserved endpoint proofs, DAG/order/budget rejection\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoArcRegistrySubdivisions.h")
    if(argc!=2)throw std::invalid_argument("subdivision evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: actual shared-point curve subdivision, including endpoint gates");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
