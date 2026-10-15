#include "spin/detail/RankTwoArcRegistrySubdivisions.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "fixtures/rank_two_source_cases.h"
#include "verification/Sha256.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
namespace {
using Q=CGAL::Gmpq;
void require(bool x,const char* why){if(!x)throw std::runtime_error(why);}
std::string id(const spin::RankTwoArrangement& a,const spin::RankTwoArcSubdivision& s){
    std::ostringstream row;row<<"rank-two-moving-arc-subdivision-v1\n"<<a.scientificHash()<<'\n'
        <<s.arcIdentity<<'\n'<<s.pointIdentities[0]<<'\n'<<s.pointIdentities[1]<<'\n'
        <<s.fraction.numerator()<<'/'<<s.fraction.denominator()<<'\n';
    const auto bytes=row.str();verification::Sha256 h;
    h.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());return h.finish();
}
void coord(std::ostream& out,const spin::detail::CertifiedAlgebraicReal& v){
    const auto& c=v.certificate();out<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
    for(const auto& x:c.coefficients)out<<x<<',';out<<'\n';
}
void tests(const char* path){
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"affine-leaf-dependency-method-v1");
    const auto a=analyzeRankTwoGlobalArrangement(source,Q(3,2));RankTwoArcSamplingPlan plan;
    for(const auto& arc:arrangementRecords(a).arcs)
        if(arc.ends[0].location==RankTwoArcEndLocation::Finite && arc.ends[1].location==RankTwoArcEndLocation::Finite)
            plan.emplace(arc.identity,std::vector<Q>{Q(1,4),Q(1,2),Q(3,4)});
    const auto initial=buildRankTwoArcRegistry(a,plan);
    std::ofstream out(path);require(bool(out),"cannot save affine dependency evidence");
    out<<"rank-two-native-affine-dependency-test-v1\nAUTHORITY "<<a.scientificHash()<<'\n';
    std::size_t count=0;
    // ALL actual quartic F arcs, not a selected well-conditioned branch.
    for(const auto& arc:initial.arcs()){
        if(arc.supportPolynomial.find("4,0:")==std::string::npos)continue;
        require(arc.finite && !arc.vertical,"frozen quartic arc bounds absent");
        const RankTwoRegistryPoint *quarter=nullptr,*threeQuarter=nullptr;
        for(const auto& p:initial.points())if(p.kind==RankTwoRegistryPointKind::ArcSample && p.arcIdentity==arc.identity){
            if(p.fraction==Q(1,4))quarter=&p;if(p.fraction==Q(3,4))threeQuarter=&p;
        }
        require(quarter && threeQuarter,"frozen primitive moving samples missing");
        RankTwoArcSubdivision first{arc.identity,{*arc.ends[0].vertexIdentity,quarter->identity},Q(1,2)};
        RankTwoArcSubdivision second{arc.identity,{id(a,first),threeQuarter->identity},Q(1,2)};
        // Composition of actual parent definitions: 0 ->1/4 gives1/8;
        // (1/8+3/4)/2=7/16. This is NOT a chord-coordinate interpolation.
        const auto expected=sampleRankTwoArc(a,arc.identity,Q(7,16));
        out<<"CASE_BEGIN "<<arc.identity<<" EXPECTED 7/16\n";
        coord(out,expected.coordinates[0]);coord(out,expected.coordinates[1]);out.flush();
        std::cerr<<"AFFINE_CASE_BEGIN "<<arc.identity<<'\n';
        const auto actual=buildRankTwoArcRegistryWithSubdivisions(a,plan,{}, {first,second});
        const RankTwoRegistryPoint* point=nullptr;
        for(const auto& p:actual.points())if(p.identity==id(a,second))point=&p;
        require(point && point->endpointIdentities==second.pointIdentities && point->fraction==Q(1,2),
            "flattening erased original moving parent identities");
        for(int axis=0;axis<2;++axis){
            coord(out,point->coordinates[axis]);
            require(point->coordinates[axis].certificate().coefficients.size()
                <=expected.coordinates[axis].certificate().coefficients.size(),
                "primitive affine leaves treated as independent conjugates; unnecessary degree growth");
            require(point->coordinates[axis].expression()==expected.coordinates[axis].expression(),
                "composed native projection or selected fiber branch changed");
        }
        out<<"REGISTRY_BEGIN\n"<<actual.canonicalRecord()<<"REGISTRY_END\nCASE_END\n";out.flush();++count;
    }
    // The frozen arrangement has inner AND outer quartic branches: eight
    // actual arcs. Count against the authoritative inventory, not four
    // quadrants mistaken for the complete native curve inventory.
    std::size_t nativeCount=0;
    for(const auto& arc:arrangementRecords(a).arcs)
        nativeCount+=arc.supportPolynomial.find("4,0:")!=std::string::npos;
    require(count==nativeCount && count==8,"not all eight original quartic native arcs checked");
    out<<"CASES 8\nprimitive_affine_dependency PASSED\nfull_carrier_completion NOT_GRANTED\n";out.flush();
    require(bool(out),"incomplete dependency evidence");
}
}
int main(int argc,char** argv){try{
    if(argc!=2)throw std::invalid_argument("affine dependency evidence path required");tests(argv[1]);return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
