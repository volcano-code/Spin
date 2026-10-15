#include <iostream>
#include <stdexcept>

#if __has_include("spin/RankTwoSource.h")
#include "spin/RankTwoSource.h"
#include "fixtures/rank_two_source_cases.h"
#include <array>
#include <cmath>
#include <limits>
#include <string>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F&& operation, const char* code) {
    try { operation(); }
    catch (const spin::RankTwoSourceError& e) {
        require(e.code() == code, "wrong rejection classification");
        require(!e.witness().empty(), "missing rejection witness");
        return;
    }
    throw std::runtime_error("invalid source was accepted");
}
using Q=CGAL::Gmpq;
using Polynomial=std::array<Q,25>;
Polynomial add(Polynomial a,const Polynomial& b) {
    for(int i=0;i<25;++i)a[i]+=b[i];return a;
}
Polynomial scale(Polynomial a,const Q& b) {
    for(auto& v:a)v*=b;return a;
}
Polynomial multiply(const Polynomial& a,const Polynomial& b) {
    Polynomial out{};
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)
        for(int k=0;k+i<5;++k)for(int l=0;l+j<5;++l)
            out[(i+k)*5+j+l]+=a[i*5+j]*b[k*5+l];
    return out;
}
void graphIdentities(const spin::PreparedRankTwoSource& source) {
    const auto& q=source.exactCoefficients();
    const auto& graph=source.graph();
    const int powers[6]={0,5,1,10,6,2};
    Polynomial h{};
    for(int j=0;j<6;++j)h[powers[j]]=graph.height[j];
    std::array<Polynomial,3> x{};
    for(int i=0;i<3;++i){
        x[i]=scale(h,graph.kernel[i]);
        x[i][5]+=graph.plane[i][0];x[i][1]+=graph.plane[i][1];
    }
    // Integrate no shared production polynomial helper into this oracle:
    // directly substitute each complete coordinate polynomial into q.
    Polynomial surface{};surface[0]=q[0];
    const Q A[3][3]={{q[4],q[9],q[8]}, {q[9],q[5],q[7]}, {q[8],q[7],q[6]}};
    Polynomial sphere{};sphere[0]=-Q(9,4);
    Polynomial projectedGradient{};
    for(int i=0;i<3;++i){
        surface=add(surface,scale(x[i],q[i+1]));
        sphere=add(sphere,multiply(x[i],x[i]));
        Polynomial g{};g[0]=q[i+1];
        for(int j=0;j<3;++j){
            surface=add(surface,scale(multiply(x[i],x[j]),A[i][j]));
            g=add(g,scale(x[j],2*A[i][j]));
        }
        Polynomial declared{};
        declared[0]=graph.gradient[i][0];declared[5]=graph.gradient[i][1];
        declared[1]=graph.gradient[i][2];
        require(g==declared,"gradient graph polynomial differs");
        projectedGradient=add(projectedGradient,scale(g,graph.kernel[i]));
        Q orthogonal=0;
        for(int j=0;j<3;++j)orthogonal+=A[i][j]*graph.kernel[j];
        require(orthogonal==0,"graph kernel not exact");
    }
    require(surface==Polynomial{},"q(X) not identically zero");
    Polynomial beta{};beta[0]=graph.beta;
    require(projectedGradient==beta && graph.beta>0,"oriented gradient identity failed");
    require(sphere==graph.spherePolynomial(Q(3,2)),"sphere polynomial coefficient error");
    for(int j=0;j<2;++j){
        Q perpendicular=0;
        for(int i=0;i<3;++i)perpendicular+=graph.kernel[i]*graph.plane[i][j];
        require(perpendicular==0,"plane not perpendicular to kernel");
    }
}
void tests() {
    using Q=CGAL::Gmpq;
    using C=spin::RankTwoExactCoefficients;
    using P=spin::PreparedRankTwoSource;
    // Removing exact-rank admission, sign preservation or identity binding
    // must fail these hand-derived fixtures; no eigenvalue threshold oracle.
    const C disk=rank_two_cases::disk();
    const auto source=P::fromExact(disk,"target-frame-v1");
    auto twice=disk; for(auto& x:twice)x*=2;
    const auto scaled=P::fromExact(twice,"target-frame-v1");
    require(source.geometryHash()==scaled.geometryHash(),"positive scale changed geometry");
    require(source.inputHash()!=scaled.inputHash(),"raw source identity was discarded");
    require(source.numericViewHash().size()==64,"numeric view lacks binding digest");
    require(source.numericViewHash()==P::fromExact(disk,"target-frame-v1").numericViewHash(),
            "numeric view binding not deterministic");
    require(source.numericViewHash()!=scaled.numericViewHash(),
            "numeric view digest not bound to raw authority");
    auto negative=disk; for(auto& x:negative)x=-x;
    require(source.geometryHash()!=P::fromExact(negative,"target-frame-v1").geometryHash(),
            "negative scale lost the material side");
    require(source.geometryHash()!=P::fromExact(disk,"another-frame").geometryHash(),
            "coordinate identity omitted from geometry");
    require(source.exactCoefficients()==disk,"authority changed on construction");
    require(source.inputHash().size()==64 && source.geometryHash().size()==64,
            "missing SHA-256 identity");
    require(std::abs(source.normalizedView().values().norm()-1)<1e-14,"invalid numeric view");
    rejects([&]{P::fromExact(disk,"");},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{P::fromExact(disk,std::string("x\0y",3));},"INPUT_IDENTITY_MISMATCH");
    auto full=disk; full[6]=Q(1)/Q(CGAL::Gmpz("10000000000000000000000000000000000000000"));
    rejects([&]{P::fromExact(full,"f");},"EXACT_STRATUM_MISMATCH");
    auto cylinder=disk; cylinder[3]=0;
    rejects([&]{P::fromExact(cylinder,"f");},"EXACT_STRATUM_MISMATCH");
    auto rankOne=disk; rankOne[5]=0;
    rejects([&]{P::fromExact(rankOne,"f");},"EXACT_STRATUM_MISMATCH");
    rejects([&]{P::fromExact(C{},"f");},"EXACT_STRATUM_MISMATCH");
    spin::QuadricCoefficients10::Vector raw;
    raw << 0,0,0,1,-1,-1,0,0,0,0;
    const auto represented=spin::QuadricCoefficients10::fromUnnormalized(raw);
    const auto exact=P::fromRepresented(represented,"f");
    for(int i=0;i<10;++i)
        require(exact.exactCoefficients()[i]==Q(represented.values()[i]),
                "represented values were renormalized or decimalized");
    raw[6]=std::ldexp(1.0,-500);
    rejects([&]{P::fromRepresented(spin::QuadricCoefficients10::fromUnnormalized(raw),"f");},
            "EXACT_STRATUM_MISMATCH");
    // Huge/small positive scales must neither overflow nor become zero.
    Q huge(CGAL::Gmpz("1"+std::string(400,'0')));
    auto large=disk,small=disk;
    for(int i=0;i<10;++i){large[i]*=huge;small[i]/=huge;}
    require(P::fromExact(large,"target-frame-v1").geometryHash()==source.geometryHash(),
            "large scale changed authority");
    require(P::fromExact(small,"target-frame-v1").geometryHash()==source.geometryHash(),
            "small scale changed authority");
    // R=1/7*[[-3,-2,6],[6,-3,2],[2,6,3]] is rational orthogonal;
    // q=49*(z'-x'^2-2*y'^2), with (x',y',z')=R*x.
    const C oblique{0,14,42,21,-81,-22,-44,24,-6,30};
    const auto tilted=P::fromExact(oblique,"f");
    require(tilted.normalizedView().values().allFinite(),"oblique numeric view invalid");
    require(tilted.graph().kernelCoordinate==1,"wrong maximum-principal-minor chart");
    require(tilted.graph().kernel==std::array<Q,3>{Q(1,3),1,Q(1,2)},"rotated kernel differs");
    require(tilted.graph().beta==Q(343,6),"rotated orientation differs");
    graphIdentities(source);graphIdentities(tilted);
    graphIdentities(P::fromExact(negative,"f"));
    graphIdentities(P::fromExact(C{0,1,1,1,-3,-3,-8,4,4,-1},"f"));
    // For the disk, nonzero unit coefficients have exact squared value 1/3.
    const auto& bounds=source.normalizedViewErrorBounds();
    for(int i=0;i<10;++i){
        require(bounds[i]>=0 && bounds[i]<Q("1/1000000000000"),"unbounded numeric view error");
        const Q value(std::abs(source.normalizedView().values()[i]));
        if(disk[i]!=0){
            require((value-bounds[i])*(value-bounds[i])<=Q(1,3)
                    && (value+bounds[i])*(value+bounds[i])>=Q(1,3),
                    "numeric-view interval excludes exact normalization");
        }
    }
    std::cout << "PASS exact source admission, identity, scaling, represented rank, graph identities\n"
              << source.geometryHash() << '\n' << tilted.geometryHash() << '\n';
}
}
#endif

int main() {
    try {
#if __has_include("spin/RankTwoSource.h")
        tests();
        return 0;
#else
        throw std::runtime_error("MISSING_FEATURE: exact rank-two source factory");
#endif
    } catch(const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
