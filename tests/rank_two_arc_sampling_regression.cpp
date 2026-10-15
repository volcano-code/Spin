#include <iostream>
#include <stdexcept>
#include <sstream>
#include <fstream>
#if __has_include("spin/detail/RankTwoArcSampling.h")
#include "spin/detail/RankTwoArcSampling.h"
#include "spin/detail/RankTwoArcStrip.h"
#include "spin/detail/RankTwoNormalOverlay.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "fixtures/rank_two_source_cases.h"
#include <CGAL/Arr_algebraic_segment_traits_2.h>
namespace {
using E=CORE::Expr;
using Q=CGAL::Gmpq;
void require(bool v,const char* reason){if(!v)throw std::runtime_error(reason);}
E expr(const Q& q){std::ostringstream s;s<<q;return E(CORE::BigRat(s.str().c_str()));}
bool onOriginal(const std::array<Q,25>& coefficients,
    const std::array<spin::detail::CertifiedAlgebraicReal,2>& coordinates){
    using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
    using P=Traits::Polynomial_2;
    CGAL::Gmpz denominator(1);
    for(const auto& c:coefficients)mpz_lcm(denominator.mpz(),denominator.mpz(),c.denominator().mpz());
    const P x=CGAL::shift(P(1),1,0),y=CGAL::shift(P(1),1,1);
    P p(0);
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)if(coefficients[5*i+j]!=0){
        const auto& c=coefficients[5*i+j];std::ostringstream integer;
        integer<<c.numerator()*(denominator/c.denominator());
        p+=P(CORE::BigInt(integer.str()))*CGAL::ipower(x,i)*CGAL::ipower(y,j);
    }
    // Independent checker: reconstruct the point from the returned coordinate
    // certificates, NOT from the sampler's supporting curve/branch. Query the
    // original F/h/D with the unbounded exact zero predicate. CORE expansion
    // of two dependent root expressions is exact but impractically costly.
    auto& kernel=Traits::CKvA_2::instance().kernel();
    const auto point=kernel.construct_algebraic_real_2_object()(
        spin::detail::rankTwoAlgebraicHandle(coordinates[0]),
        spin::detail::rankTwoAlgebraicHandle(coordinates[1]));
    return kernel.is_zero_at_2_object()(CGAL::make_square_free(p),point);
}
void tests(bool trace,const char* output,bool stripOnly=false){
    using namespace spin;
    const auto source=PreparedRankTwoSource::fromExact(rank_two_cases::annulus(),"arc-sample-method");
    const auto global=analyzeRankTwoGlobalArrangement(source,Q(3,2));
    RankTwoNormalFace cone;cone.vertexIds={0,1,2};cone.rays={{{0,0,1},{1,0,0},{0,1,0}}};
    const auto local=analyzeRankTwoNormalFace(source,global,cone);
    std::ofstream evidence(output);
    if(!evidence)throw std::runtime_error("cannot open independent arc evidence");
    evidence<<"rank-two-actual-arc-sampling-v1\n";
    std::size_t count=0,vertical=0,nonlinear=0,infiniteRejected=0;
    if(!stripOnly)for(const auto* a:{&global,&local}){
        const auto& records=arrangementRecords(*a);
        for(const auto& arc:records.arcs){
            if(arc.ends[0].location!=RankTwoArcEndLocation::Finite
                ||arc.ends[1].location!=RankTwoArcEndLocation::Finite){
                bool rejected=false;
                try{sampleRankTwoArc(*a,arc.identity,Q(1,2));}
                catch(const RankTwoSourceError& e){rejected=e.code()=="ARRANGEMENT_UNRESOLVED";}
                require(rejected,"unbounded arc was sampled without certified clipping");
                ++infiniteRejected;continue;
            }
            const auto endpoint=[&](int i)->const RankTwoVertexRecord&{
                for(const auto& v:records.vertices)if(v.identity==*arc.ends[i].vertexIdentity)return v;
                throw std::runtime_error("missing test endpoint");
            };
            for(const Q t:{Q(1,3),Q(1,2),Q(2,3)}){
                if(trace)std::cerr<<"SAMPLE_BEGIN "<<count<<' '<<arc.identity<<" t="<<t<<" vertical="<<arc.vertical<<'\n';
                const auto sample=sampleRankTwoArc(*a,arc.identity,t);
                if(trace)std::cerr<<"SAMPLE_RETURNED "<<count<<'\n';
                const auto& xy=sample.coordinates;
                const E x=xy[0].expression(),y=xy[1].expression();
                const int axis=arc.vertical?1:0;
                const E left=detail::CertifiedAlgebraicReal::reload(endpoint(0).coordinates[axis]).expression();
                const E right=detail::CertifiedAlgebraicReal::reload(endpoint(1).coordinates[axis]).expression();
                if(trace)std::cerr<<"CHECK_WEIGHTED "<<count<<'\n';
                require((xy[axis].expression()-((1-expr(t))*left+expr(t)*right)).sign()==0,
                    "arc sample substituted isolator midpoint for moving coordinate");
                require(sample.arcIdentity==arc.identity && sample.endpointIdentities[0]==endpoint(0).identity
                    && sample.endpointIdentities[1]==endpoint(1).identity && sample.fraction==t,
                    "sample lost parent dependency identity");
                for(const auto& owner:arc.owners){
                    if(trace)std::cerr<<"CHECK_OWNER "<<count<<' '<<owner.polynomial<<'\n';
                    const int i=owner.polynomial=="F"?0:owner.polynomial=="h"?1:owner.polynomial=="D"?2:-1;
                    if(i>=0)require(onOriginal(records.originalFhD[i],xy),"point is not on its exact owner curve");
                    if(i==0){
                        const E wx=detail::CertifiedAlgebraicReal::reload(arc.interiorWitness[0]).expression();
                        const E wy=detail::CertifiedAlgebraicReal::reload(arc.interiorWitness[1]).expression();
                        // Annulus F has two negative-y (and two positive-y)
                        // branches. F==0 and y sign alone would miss swapping
                        // the inner and outer circle of the same polynomial.
                        require((x*x+y*y-E(3)/2).sign()==(wx*wx+wy*wy-E(3)/2).sign(),
                            "sample switched inner/outer annulus component");
                    }
                }
                if(!arc.vertical){
                    const E witness=detail::CertifiedAlgebraicReal::reload(arc.interiorWitness[1]).expression();
                    require(y.sign()==witness.sign(),"sample switched upper/lower circle branch");
                    ++nonlinear;
                }else{require(x==detail::CertifiedAlgebraicReal::reload(endpoint(0).coordinates[0]).expression(),"vertical sample moved off its line");++vertical;}
                const auto again=sampleRankTwoArc(*a,arc.identity,t);
                require(again.identity==sample.identity && again.coordinates[0].expression()==x
                    && again.coordinates[1].expression()==y,"same arc dependency changed on replay");
                evidence<<"SAMPLE "<<count<<' '<<sample.identity<<' '<<sample.arrangementIdentity<<' '
                    <<sample.arcIdentity<<' '<<sample.endpointIdentities[0]<<' '<<sample.endpointIdentities[1]
                    <<' '<<t.numerator()<<'/'<<t.denominator()<<' '<<sample.vertical<<' '<<sample.interiorBranch<<'\n';
                for(const auto& value:xy){
                    const auto& c=value.certificate();
                    evidence<<"COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
                    for(const auto& coefficient:c.coefficients)evidence<<coefficient<<',';
                    evidence<<'\n';
                }
                ++count;
                if(trace)std::cerr<<"SAMPLE_END "<<count<<'\n';
            }
        }
    }
    if(!stripOnly){
        require(vertical>0 && nonlinear>0,"missing vertical or nonlinear actual-arrangement coverage");
        require(infiniteRejected>0,"missing unbounded arc rejection coverage");
    }
    detail::AlgebraicKernel k;
    const auto zero=detail::CertifiedAlgebraicReal::fromRoot(k.construct_algebraic_real_1_object()(CORE::BigRat(0)));
    require(!onOriginal(arrangementRecords(global).originalFhD[0],{zero,zero}),
        "independent checker accepted coordinates off the original F");
    // Consume the actual quartic F from the immutable annulus arrangement,
    // not an independently substituted circle. Bind two exact moving points
    // on its upper inner component to a rational, regular isolation strip.
    const auto& g=arrangementRecords(global);
    std::size_t certifiedStrips=0;
    evidence.flush();
    for(const auto& arc:g.arcs){
        if(arc.vertical || arc.ends[0].location!=RankTwoArcEndLocation::Finite
            || arc.ends[1].location!=RankTwoArcEndLocation::Finite)continue;
        bool isF=false;for(const auto& owner:arc.owners)isF|=owner.polynomial=="F";
        if(!isF)continue;
        if(trace)std::cerr<<"STRIP_SELECT "<<arc.identity<<'\n';
        const E wy=detail::CertifiedAlgebraicReal::reload(arc.interiorWitness[1]).expression();
        if(wy<=0 || wy>=1)continue;
        if(trace)std::cerr<<"STRIP_LEFT_SAMPLE\n";
        const auto l=sampleRankTwoArc(global,arc.identity,Q(1,2)-Q(1,65536));
        if(trace)std::cerr<<"STRIP_RIGHT_SAMPLE\n";
        const auto r=sampleRankTwoArc(global,arc.identity,Q(1,2)+Q(1,65536));
        if(trace)std::cerr<<"STRIP_SAMPLES_RETURNED\n";
        const detail::RankTwoStripBox box{Q(-1,4096),Q(1,4096),Q(89,100),Q(9,10)};
        for(const auto* point:{&l,&r}){
            if(trace)std::cerr<<"STRIP_BOX_CHECK "<<point->fraction<<'\n';
            const E px=point->coordinates[0].expression(),py=point->coordinates[1].expression();
            require(px>expr(box[0]) && px<expr(box[1]) && py>expr(box[2]) && py<expr(box[3]),
                "actual arc endpoints do not belong to locked regular strip");
            if(trace)std::cerr<<"STRIP_OWNER_CHECK "<<point->fraction<<'\n';
            require(onOriginal(g.originalFhD[0],point->coordinates),"strip endpoint not on original quartic");
        }
        if(trace)std::cerr<<"STRIP_CERTIFICATE\n";
        const auto proof=detail::certifyRankTwoRegularStrip(g.originalFhD[0],box,
            {g.originalFhD[1],g.originalFhD[2]},Q(1,1048576));
        // Input is z+2-x^2-y^2=0, hence h=x^2+y^2-2 is NEGATIVE here.
        // F uses h^2, but its sign does not change the stored h owner.
        require(proof.yDerivativeSign==-1 && proof.obstacleSigns==std::vector<int>{-1,1}
            && proof.chordDistanceUpperBound<=Q(1,1048576),"actual quartic strip not isolated from h/D");
        evidence<<"REGULAR_STRIP "<<global.scientificHash()<<' '<<arc.identity<<' '<<proof.inputIdentity;
        for(const auto& b:box)evidence<<' '<<b;
        evidence<<' '<<proof.derivativeLowerBound<<' '<<proof.slopeUpperBound<<' '
            <<proof.curvatureUpperBound<<' '<<proof.chordDistanceUpperBound<<'\n';
        for(const auto* point:{&l,&r}){
            evidence<<"STRIP_ENDPOINT "<<point->identity<<' '<<point->fraction<<'\n';
            for(const auto& value:point->coordinates){
                const auto& c=value.certificate();
                evidence<<"STRIP_COORD "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
                for(const auto& coefficient:c.coefficients)evidence<<coefficient<<',';
                evidence<<'\n';
            }
        }
        ++certifiedStrips;
    }
    require(certifiedStrips==1,"missing/duplicated actual regular quartic strip");
    bool rejected=false;
    try{sampleRankTwoArc(global,"not-an-arc",Q(1,2));}catch(const RankTwoSourceError& e){rejected=e.code()=="INPUT_IDENTITY_MISMATCH";}
    require(rejected,"unknown curve identity accepted");
    evidence<<"CHECKED "<<count<<" VERTICAL "<<vertical<<" OTHER "<<nonlinear
        <<" UNBOUNDED_REJECTED "<<infiniteRejected<<" REGULAR_STRIPS "<<certifiedStrips
        <<"\nfull_carrier_completion NOT_GRANTED\n";
    evidence.flush();require(bool(evidence),"arc evidence write failed");
    std::cout<<"PASS actual arrangement arc samples "<<count<<" vertical "<<vertical
        <<" nonlinear "<<nonlinear<<"; one local 2D strip; endpoint-sector/topology/mesh NOT_EXECUTED\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoArcSampling.h")
    const bool stripOnly=argc==3 && std::string(argv[2])=="--strip-profile";
    const bool trace=stripOnly || (argc==3 && std::string(argv[2])=="--trace");
    if(argc!=2 && !trace)throw std::invalid_argument("arc evidence path and optional --trace required");
    if(trace){
        const auto sink=[](const char* phase,void*)noexcept{std::cerr<<"PHASE "<<phase<<'\n';};
        spin::detail::ScopedRankTwoDiagnostics observer(sink,nullptr);tests(true,argv[1],stripOnly);
    }else tests(false,argv[1]);
    return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: certified moving samples from actual arrangement arcs");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
