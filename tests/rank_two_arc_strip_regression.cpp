#include <fstream>
#include <iostream>
#include <stdexcept>
#if __has_include("spin/detail/RankTwoArcStrip.h")
#include "spin/detail/RankTwoArcStrip.h"
#include "spin/RankTwoSource.h"
namespace {
using Q=CGAL::Gmpq;
using P=spin::detail::RankTwoStripPolynomial;
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
template<class F> void rejects(F f,const char* code,const char* why){
    bool rejected=false;
    try{f();}catch(const spin::RankTwoSourceError& e){rejected=e.code()==code;}
    require(rejected,why);
}
void tests(const char* path){
    using namespace spin::detail;
    // Hand-derived: y=x^2, x in [-1,1]; the horizontal chord y=1 has
    // exact maximum graph/chord distance 1. Losing Pxx makes this bound zero.
    P parabola{};parabola[1]=1;parabola[10]=-1;
    const RankTwoStripBox wide{Q(-1),Q(1),Q(-1),Q(2)};
    auto p=certifyRankTwoRegularStrip(parabola,wide,{},Q(1));
    require(p.yDerivativeSign==1 && p.derivativeLowerBound==1
        && p.slopeUpperBound==2 && p.curvatureUpperBound==2
        && p.chordDistanceUpperBound==1,"incorrect parabola derivative/chord proof");
    rejects([&]{certifyRankTwoRegularStrip(parabola,wide,{},Q(999,1000));},
        "ARRANGEMENT_UNRESOLVED","coarse chord received requested accuracy");
    // Offset parabola never touches the true graph, but this coarse chord
    // crosses it twice. Endpoint-only or finite sampling checks miss it.
    P obstruction=parabola;obstruction[0]=Q(-1,2);
    rejects([&]{certifyRankTwoRegularStrip(parabola,wide,{obstruction},Q(1));},
        "ARRANGEMENT_UNRESOLVED","other curve inside strip was ignored");
    P above{};above[1]=1;above[0]=-3;
    auto separated=certifyRankTwoRegularStrip(parabola,wide,{above},Q(1));
    require(separated.obstacleSigns==std::vector<int>{-1},"obstacle exclusion sign not certified");
    // x^2+y^2=1, upper regular branch. Analytic interval bounds:
    // |Py|>=3/2, |f'|<=1/3, |f''|<=40/27, error<=5/108.
    P circle{};circle[0]=-1;circle[10]=1;circle[2]=1;
    const RankTwoStripBox upper{Q(-1,4),Q(1,4),Q(3,4),Q(5,4)};
    const auto c=certifyRankTwoRegularStrip(circle,upper,{},Q(1,16));
    require(c.derivativeLowerBound==Q(3,2) && c.slopeUpperBound==Q(1,3)
        && c.curvatureUpperBound==Q(40,27) && c.chordDistanceUpperBound==Q(5,108),
        "circle certificate lacks independently derived interval bound");
    P reversed=circle;for(auto& q:reversed)q*=Q(-7);
    const auto r=certifyRankTwoRegularStrip(reversed,upper,{},Q(1,16));
    require(r.yDerivativeSign==-1 && r.slopeUpperBound==c.slopeUpperBound
        && r.chordDistanceUpperBound==c.chordDistanceUpperBound,
        "polynomial rescaling changed geometric certificate");
    rejects([&]{certifyRankTwoRegularStrip(circle,{Q(3,4),Q(5,4),Q(-1,4),Q(1,4)},{},Q(1));},
        "ARRANGEMENT_UNRESOLVED","critical endpoint accepted as regular strip");
    rejects([&]{certifyRankTwoRegularStrip(circle,{Q(-1,4),Q(1,4),Q(-2),Q(2)},{},Q(1));},
        "ARRANGEMENT_UNRESOLVED","two branches accepted as unique graph");
    rejects([&]{certifyRankTwoRegularStrip(circle,{Q(-1,4),Q(1,4),Q(2),Q(3)},{},Q(1));},
        "ARRANGEMENT_UNRESOLVED","empty strip accepted as curve graph");
    // Four corner signs alone bracket a graph, but at x=1/2 its height is
    // one, above this rectangle. The complete upper edge must be checked.
    P bump{};bump[1]=1;bump[5]=-4;bump[10]=4;
    rejects([&]{certifyRankTwoRegularStrip(bump,{Q(0),Q(1),Q(-1,4),Q(3,4)},{},Q(1));},
        "ARRANGEMENT_UNRESOLVED","corner checks substituted for whole-edge bracketing");
    // Nonzero mixed derivative. Omitting the 2*Pxy*f' term makes this
    // rational graph falsely have zero curvature.
    P mixed{};mixed[6]=1;mixed[1]=2;mixed[0]=-1;
    const auto m=certifyRankTwoRegularStrip(mixed,{Q(0),Q(1),Q(1,4),Q(3,4)},{},Q(1,16));
    require(m.curvatureUpperBound==Q(3,8) && m.chordDistanceUpperBound==Q(3,64),
        "mixed derivative contribution lost");
    P root{};root[2]=1;root[5]=-1;
    const auto close=certifyRankTwoRegularStrip(root,{Q(1,16),Q(5,64),Q(3,16),Q(5,16)},{},Q(1,800));
    require(close.derivativeLowerBound==Q(3,8) && close.curvatureUpperBound==Q(1024,27)
        && close.chordDistanceUpperBound==Q(1,864),"near-critical regular graph incorrectly bounded");
    rejects([&]{certifyRankTwoRegularStrip(parabola,{Q(1),Q(-1),Q(-1),Q(2)},{},Q(1));},
        "INPUT_IDENTITY_MISMATCH","reversed box accepted");
    rejects([&]{certifyRankTwoRegularStrip(P{},wide,{},Q(1));},
        "INPUT_IDENTITY_MISMATCH","zero polynomial accepted");
    rejects([&]{certifyRankTwoRegularStrip(parabola,wide,{},Q(0));},
        "INPUT_IDENTITY_MISMATCH","zero requested error accepted");
    // Positive separation smaller than double roundoff remains strict.
    P flat{};flat[1]=1;Q tiny(1);for(int i=0;i<100;++i)tiny/=2;
    P near=flat;near[0]=-tiny;
    const auto gap=certifyRankTwoRegularStrip(flat,{Q(-1),Q(1),-tiny/4,tiny/4},{near},tiny);
    require(gap.chordDistanceUpperBound==0 && gap.obstacleSigns[0]==-1,
        "positive exact tiny gap lost to floating conversion");
    // Exact translation of y=x^2: y=(x-1)^2+1 on [1,1.01].
    // Monomial interval cancellation falsely makes both horizontal-edge
    // signs UNKNOWN. Whole-edge Bernstein bounds must prove this regular
    // graph without changing its derivative/error bound or accepting only
    // its corners. The hand-derived chord error is 2*(.01)^2/8=1/40000.
    P shifted{};shifted[1]=1;shifted[0]=-2;shifted[5]=2;shifted[10]=-1;
    const auto translated=certifyRankTwoRegularStrip(shifted,
        {Q(1),Q(101,100),Q(99,100),Q(51,50)},{},Q(1,10000));
    require(translated.yDerivativeSign==1 && translated.derivativeLowerBound==1
        && translated.slopeUpperBound==Q(1,50) && translated.curvatureUpperBound==2
        && translated.chordDistanceUpperBound==Q(1,40000),
        "translation cancellation blocked exact whole-edge graph bracketing");
    // Original annulus at the negative-x projection extreme. The F quartic
    // is strictly positive throughout this WHOLE chart box: squared radius
    // >=(149/100)^2 > (3+sqrt(2))/2. Expanded monomial intervals lose the
    // cancellation between 4*r^4 and -12*r^2, although no curve is present.
    P projection{};projection[0]=Q(-9,4);projection[10]=1;projection[2]=1;
    P originalF{};originalF[0]=7;originalF[10]=-12;originalF[2]=-12;
    originalF[20]=4;originalF[12]=8;originalF[4]=4;
    const RankTwoStripBox extreme{Q(-713,409600),Q(0),Q(-151,100),Q(-149,100)};
    const auto excluded=certifyRankTwoRegularStrip(projection,extreme,{originalF},Q(1,1048576));
    require(excluded.obstacleSigns==std::vector<int>{1},"actual F obstacle not excluded at D extreme");
    P actualCrossing{};actualCrossing[0]=Q(3,2);actualCrossing[1]=1;
    rejects([&]{certifyRankTwoRegularStrip(projection,extreme,{actualCrossing},Q(1,1048576));},
        "ARRANGEMENT_UNRESOLVED","whole-box cancellation repair accepted an actual crossing");
    std::ofstream out(path);require(bool(out),"cannot open strip evidence");
    out<<"rank-two-regular-strip-kernel-v1\n"
        <<"PARABOLA "<<p.derivativeLowerBound<<' '<<p.slopeUpperBound<<' '
        <<p.curvatureUpperBound<<' '<<p.chordDistanceUpperBound<<'\n'
        <<"CIRCLE "<<c.derivativeLowerBound<<' '<<c.slopeUpperBound<<' '
        <<c.curvatureUpperBound<<' '<<c.chordDistanceUpperBound<<'\n'
        <<"TINY_GAP "<<tiny<<'\n'
        <<"TRANSLATED_PARABOLA "<<translated.derivativeLowerBound<<' '
        <<translated.slopeUpperBound<<' '<<translated.curvatureUpperBound<<' '
        <<translated.chordDistanceUpperBound<<"\nEXTREME_OBSTACLE "<<excluded.obstacleSigns.front()
        <<"\nfull_carrier_completion NOT_GRANTED\n";
    out.flush();require(bool(out),"strip evidence write failed");
    std::cout<<"PASS regular strip: exact graph, derivative bounds, obstacles and chord error; endpoint sectors/3D NOT_EXECUTED\n";
}
}
#endif
int main(int argc,char** argv){try{
#if __has_include("spin/detail/RankTwoArcStrip.h")
    if(argc!=2)throw std::invalid_argument("evidence path required");tests(argv[1]);return 0;
#else
    throw std::runtime_error("MISSING_FEATURE: exact regular-arc isolation and chord bound");
#endif
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
