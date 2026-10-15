#include <iostream>
#include <stdexcept>
#include <fstream>
#if __has_include("spin/detail/RankTwoCarrierKernels.h")
#include "spin/detail/RankTwoCarrierKernels.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "spin/RankTwoSource.h"
#include <array>
#include <vector>
namespace {
using namespace spin::detail;
using E=CORE::Expr;
using Q=CGAL::Gmpq;
void require(bool v,const char* reason){if(!v)throw std::runtime_error(reason);}
template<class F> void rejects(F f,const char* code){
    try{f();}catch(const spin::RankTwoSourceError& e){require(e.code()==code,"wrong failure class");return;}
    throw std::runtime_error("invalid carrier kernel input accepted");
}
CertifiedAlgebraicReal root(const AlgebraicPolynomial& p,std::size_t i){
    AlgebraicKernel k;std::vector<AlgebraicReal> values;
    k.solve_1_object()(p,false,std::back_inserter(values));
    return CertifiedAlgebraicReal::fromRoot(values.at(i));
}
void thickness(){
    const std::array<bool,3> boundary{true,true,true},internal{false,false,false};
    // Disk vertices all on F=0: upper and lower coincide everywhere on its
    // coarse interpolant, even though the analytic disk contains material.
    const E disk=(CORE::sqrt(E(10))-1)/2;
    rejects([&]{certifyRankTwoTriangleThickness({disk,disk,disk},{disk,disk,disk},boundary);},"SOURCE_TOPOLOGY_INVALID");
    // Equator-only projected ball has the same defect with zero heights.
    rejects([&]{certifyRankTwoTriangleThickness({0,0,0},{0,0,0},boundary);},"SOURCE_TOPOLOGY_INVALID");
    // Positive total prism volume is not enough: edge 0--1 is pinched.
    rejects([&]{certifyRankTwoTriangleThickness({0,0,0},{0,0,1},internal);},"SOURCE_TOPOLOGY_INVALID");
    certifyRankTwoTriangleThickness({0,0,0},{0,0,1},{true,false,false});
    certifyRankTwoTriangleThickness({-1,-1,-1},{1,1,1},internal);
    rejects([&]{certifyRankTwoTriangleThickness({0,0,0},{-1,2,2},internal);},"SOURCE_TOPOLOGY_INVALID");
    // No tolerance may erase a genuine strictly positive algebraic gap.
    const E tiny=E(1)/CORE::BigInt("1267650600228229401496703205376");
    certifyRankTwoTriangleThickness({0,0,0},{tiny,tiny,tiny},internal);
    // Same high-degree expression representation is exactly equal, whereas
    // a distinct arbitrarily small gap must use the original exact path.
    E shared(0);for(const int p:{2,3,5,7,11,13,17,19,23,29,31,37})shared+=CORE::sqrt(E(p));
    require(rankTwoExactLiftDifference(shared,shared).sign()==0,"shared seam expression not exactly zero");
    certifyRankTwoTriangleThickness({shared,shared,0},{shared,shared,1},{true,false,false});
    require(rankTwoExactLiftDifference(E(1)+tiny,E(1)).sign()>0,"distinct positive gap erased by identity shortcut");
    require(rankTwoExactLiftDifference(E(1),E(1)+tiny).sign()<0,"distinct negative gap erased by identity shortcut");
}
void interpolation(){
    // Explicit coefficients avoid CGAL's two-iterator template being selected
    // for a pair of ints (which would become vector(count,value)).
    const auto a=root(AlgebraicPolynomial(CORE::BigInt(-1),CORE::BigInt(3)),0);
    const auto b=root(AlgebraicPolynomial(CORE::BigInt(-5),CORE::BigInt(7)),0);
    const auto c=sampleRankTwoAlgebraicCoordinate(a,b,Q(2,5));
    require((35*c.expression()-17).sign()==0,"rational weighted coordinate rounded or wrong weight");
    const auto p=root(AlgebraicPolynomial(-2,0,1),1),n=root(AlgebraicPolynomial(-2,0,1),0);
    const auto zero=sampleRankTwoAlgebraicCoordinate(p,n,Q(1,2));
    require(zero.expression().sign()==0,"dependent endpoints lost resultant zero root");
    for(const Q t:{Q(1,3),Q(2,3)}){
        const auto weighted=sampleRankTwoAlgebraicCoordinate(p,n,t);
        const E value=weighted.expression();
        require((9*value*value-2).sign()==0 && value.sign()==(t<Q(1,2)?1:-1),
            "antipodal scaling lost sign or rational scale");
    }
    const std::vector<CORE::BigInt> coefficients{7,0,-12,0,4};
    const AlgebraicPolynomial quartic(coefficients.begin(),coefficients.end());
    const auto innerNegative=root(quartic,1),innerPositive=root(quartic,2);
    const auto nearCenter=sampleRankTwoAlgebraicCoordinate(innerNegative,innerPositive,Q(1,2)-Q(1,65536));
    const E scaled=32768*nearCenter.expression();
    require((4*scaled*scaled*scaled*scaled-12*scaled*scaled+7).sign()==0
        && scaled<0 && scaled*scaled<E(3)/2,"dependent quartic sample changed branch or fraction");
    const auto q=root(AlgebraicPolynomial(-3,0,1),1);
    const auto midpoint=sampleRankTwoAlgebraicCoordinate(p,q,Q(1,2));
    const E z=midpoint.expression();
    require((16*z*z*z*z-40*z*z+1).sign()==0 && z>E(3)/2 && z<2,
        "wrong conjugate of algebraic midpoint selected");
    const auto reflected=sampleRankTwoAlgebraicCoordinate(q,p,Q(1,2));
    require((reflected.expression()-z).sign()==0,"endpoint reversal changed geometric point");
    require(rankTwoAlgebraicHandle(midpoint).compare(rankTwoAlgebraicHandle(reflected))==CGAL::EQUAL,
        "arithmetic coordinate cannot be passed back to curve analysis exactly");
    auto refined=p.certificate();refined.lower="11/8";refined.upper="23/16";
    const auto replay=sampleRankTwoAlgebraicCoordinate(CertifiedAlgebraicReal::reload(refined),q,Q(1,2));
    require(replay.certificate().rootIdentity==midpoint.certificate().rootIdentity,
        "isolator refinement changed weighted sample identity");
    require(sampleRankTwoAlgebraicCoordinate(p,p,Q(1,3)).expression()==p.expression(),
        "equal-coordinate endpoints lost their value");
    const auto outside=affineRankTwoAlgebraicCoordinate(p,q,Q(-1,2));
    require(outside.expression()==(3*p.expression()-q.expression())/2,
        "signed affine resultant selected a wrong conjugate or inward interval");
    const auto opposite=affineRankTwoAlgebraicCoordinate(p,n,Q(-1,2));
    require(opposite.expression()==2*p.expression(),"signed antipodal ruler lost its branch or scale");
    const auto rationalOutside=affineRankTwoAlgebraicCoordinate(a,b,Q(-1,2));
    require(7*rationalOutside.expression()==1,"signed rational ruler uses wrong parent weight");
    rejects([&]{affineRankTwoAlgebraicCoordinate(a,b,Q(-1));},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{affineRankTwoAlgebraicCoordinate(a,b,Q(0));},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{sampleRankTwoAlgebraicCoordinate(a,b,Q(-1,3));},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{sampleRankTwoAlgebraicCoordinate(a,b,Q(1));},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{sampleRankTwoAlgebraicCoordinate(a,b,Q(0));},"INPUT_IDENTITY_MISMATCH");
    rejects([&]{sampleRankTwoAlgebraicCoordinate(a,b,Q(1,2),0);},"INPUT_IDENTITY_MISMATCH");
}
void selectedConjugateFactors(){
    // Same primitive quartic that supplies the real annulus endpoint
    // abscissae. Its four roots are +/-sqrt((3+/-sqrt(2))/2).
    // Mixing the two DISTINCT negative roots at1/16 and1/8 has degree8:
    // the degree16 all-conjugate eliminant also includes the equal-parent
    // quartic and the antipodal-parent scaled quartic. Neither contains our
    // selected root. Keeping those factors consumes the downstream cap.
    const std::vector<CORE::BigInt> coefficients{7,0,-12,0,4};
    const AlgebraicPolynomial p(coefficients.begin(),coefficients.end());
    const auto outer=root(p,0),inner=root(p,1);
    const auto a=sampleRankTwoAlgebraicCoordinate(outer,inner,Q(1,16));
    const auto b=sampleRankTwoAlgebraicCoordinate(outer,inner,Q(1,8));
    require(a.certificate().coefficients.size()<=9 && b.certificate().coefficients.size()<=9,
        "nonselected equal/antipodal conjugate factors consume downstream coordinate budget");
    const auto c=sampleRankTwoAlgebraicCoordinate(a,b,Q(1,2));
    // Independently derived rational weights: (15/16+7/8)/2=29/32.
    const E expected=(29*outer.expression()+3*inner.expression())/32;
    require((c.expression()-expected).sign()==0,
        "selected-factor reduction changed chained actual affine coordinate");
    require(c.expression()<inner.expression() && c.expression()>outer.expression(),
        "selected-factor reduction chose a conjugate outside the actual endpoints");
    std::cout<<"SELECTED_FACTOR_CHAIN PASS distinct quartic roots; exact29/32 and3/32 weights\n";
}
void affineDependencyReuse(){
    // Observable cost regression: rebuilding a new CORE polynomial root
    // discards the ALREADY PROVED exact affine dependency. Source roots are
    // created before observing the sample, not mocked or approximated.
    const auto a=root(AlgebraicPolynomial(-2,0,1),1);
    const auto b=root(AlgebraicPolynomial(-3,0,1),1);
    int expressionRebuilds=0;
    ScopedRankTwoDiagnostics observer([](const char* phase,void* context) noexcept {
        if(std::string(phase)=="algebraic-root:reload-expression")
            ++*static_cast<int*>(context);
    },&expressionRebuilds);
    const auto value=sampleRankTwoAlgebraicCoordinate(a,b,Q(2,3));
    require(expressionRebuilds==0,"proved affine sample redundantly reconstructs CORE polynomial root");
    const E z=value.expression();
    // Independently derived: 9z^2=14+4sqrt(6), hence
    // 81z^4-252z^2+100=0. All-conjugate eliminant identity remains intact.
    require(value.certificate().coefficients==std::vector<std::string>{"100","0","-252","0","81"},
        "affine expression reuse changed frozen selected-root certificate");
    require((81*z*z*z*z-252*z*z+100).sign()==0 && z>E(3)/2 && z<2,
        "affine expression reuse lost independent polynomial or conjugate");
    const auto restored=CertifiedAlgebraicReal::reload(value.certificate());
    require(expressionRebuilds==1,"external affine certificate reload bypassed CORE reconstruction");
    require((restored.expression()-z).sign()==0,"saved root differs from actual affine expression");
    std::cout<<"AFFINE_DEPENDENCY_REUSE PASS literal81z4-252z2+100 and full public reload\n";
}
void nativeHandleReuse(){
    const auto value=root(AlgebraicPolynomial(-2,0,1),1);
    const auto reloaded=CertifiedAlgebraicReal::reload(value.certificate());
    std::size_t solves=0;
    ScopedRankTwoDiagnostics sink([](const char* phase,void* raw)noexcept{
        if(std::string(phase)=="algebraic-handle:solve-frozen")++*static_cast<std::size_t*>(raw);
    },&solves);
    for(const auto& coordinate:{value,reloaded})for(int i=0;i<3;++i){
        auto handle=rankTwoAlgebraicHandle(coordinate);
        require(handle.compare(CORE::BigRat(1))==CGAL::LARGER&&handle.compare(CORE::BigRat(2))==CGAL::SMALLER,
            "native handle reuse selected wrong conjugate or interval");
        require(handle.is_root_of(AlgebraicPolynomial(-2,0,1)),"native handle reuse lost exact polynomial substitution");
        require(coordinate.certificate().rootIdentity==value.certificate().rootIdentity,
            "native query mutated saved certificate identity");
    }
    std::cerr<<"NATIVE_HANDLE_EXTRA_SOLVES "<<solves<<'\n';
    require(solves==0,"MISSING_NATIVE_HANDLE_REUSE: repeated native queries re-isolate an already authenticated root");
    const auto oldEscape=CORE::get_static_EscapePrec();CORE::get_static_EscapePrec()=CORE::extLong(64);
    try{rejects([&]{rankTwoAlgebraicHandle(value);},"ARRANGEMENT_UNRESOLVED");}
    catch(...){CORE::get_static_EscapePrec()=oldEscape;throw;}
    CORE::get_static_EscapePrec()=oldEscape;
    requireExactAlgebraicEnvironment();
    std::cout<<"NATIVE_HANDLE_REUSE PASS native and fully validated public reload; exact conjugate preserved\n";
}
#ifdef RANK_TWO_ACTUAL_PARENT_AFFINE
void actualParentAffine(const char* path){
    std::ifstream in(path);std::string schema;in>>schema;
    require(schema=="RANK_TWO_R188_FAILED_AFFINE_V1","wrong actual capture schema");
    std::vector<CertifiedAlgebraicReal> parents;
    for(int j=0;j<2;++j){
        std::size_t count=0;in>>count;require(count==(j==0?73u:17u),"actual parent degree changed");
        AlgebraicRootCertificate c;c.coefficients.resize(count);
        for(auto& s:c.coefficients)in>>s;
        in>>c.rootIndex>>c.lower>>c.upper>>c.rootIdentity;
        require(bool(in),"truncated original parent capture");
        parents.push_back(CertifiedAlgebraicReal::reload(c));
    }
    std::string extra;require(!(in>>extra),"unexpected caller-selected factor or answer");
    const auto beforeA=parents[0].certificate(),beforeB=parents[1].certificate();
    ScopedRankTwoDiagnostics diagnostics([](const char* phase,void*) noexcept {
        if(std::string(phase).find("weighted-coordinate:authenticated-parent ")==0)
            std::cout<<"PARENT_COMPONENT "<<phase<<'\n';
    },nullptr);
    // The actual failed edge5--8; the expected value is the original exact
    // affine expression, not a supplied selected factor, isolator midpoint,
    // or a certificate produced by the sampler under test.
    const auto midpoint=sampleRankTwoAlgebraicCoordinate(parents[0],parents[1],Q(1,2));
    require((2*midpoint.expression()-parents[0].expression()-parents[1].expression()).sign()==0,
        "authenticated parent construction changed actual edge midpoint");
    const auto restored=CertifiedAlgebraicReal::reload(midpoint.certificate());
    require(restored.expression()==midpoint.expression(),"actual edge midpoint does not survive public reload");
    const auto reverse=sampleRankTwoAlgebraicCoordinate(parents[1],parents[0],Q(1,2));
    require(reverse.expression()==midpoint.expression(),"actual parent reversal chose another conjugate");
    require(parents[0].certificate().coefficients==beforeA.coefficients
        &&parents[0].certificate().rootIndex==beforeA.rootIndex
        &&parents[0].certificate().rootIdentity==beforeA.rootIdentity
        &&parents[1].certificate().coefficients==beforeB.coefficients
        &&parents[1].certificate().rootIndex==beforeB.rootIndex
        &&parents[1].certificate().rootIdentity==beforeB.rootIdentity,
        "parent factor membership rewrote original certificate authority");
    require(rankTwoAlgebraicHandle(parents[0]).polynomial().degree()==72
        &&rankTwoAlgebraicHandle(parents[1]).polynomial().degree()==16,
        "parent membership mutated caller native handle polynomial");
    const auto outside=affineRankTwoAlgebraicCoordinate(parents[0],parents[1],Q(-1,2));
    require(2*outside.expression()==3*parents[0].expression()-parents[1].expression(),
        "signed actual parent factor construction selected wrong affine root");
    std::cout<<"ACTUAL_PARENT_AFFINE PASS edge5--8; original72/16 authority; exact1/2; publicreload; reversal; signed\n";
    const auto& c=midpoint.certificate();
    std::cout<<"ACTUAL_MIDPOINT "<<c.rootIdentity<<' '<<c.rootIndex<<' '<<c.lower<<' '<<c.upper<<' ';
    for(const auto& s:c.coefficients)std::cout<<s<<',';
    std::cout<<'\n';
    // Eisenstein at2 and3: these irreducible parents truly need289 paired
    // conjugates. Exact factor authentication must not bypass the256 cap.
    std::vector<CORE::BigInt> irreducible(18,0);irreducible[17]=1;irreducible[0]=-2;
    const auto u=root(AlgebraicPolynomial(irreducible.begin(),irreducible.end()),0);
    irreducible[0]=-3;
    const auto v=root(AlgebraicPolynomial(irreducible.begin(),irreducible.end()),0);
    rejects([&]{sampleRankTwoAlgebraicCoordinate(u,v,Q(1,2));},"RESOURCE_LIMIT");
    std::cout<<"GENUINE_PARENT_CAP PASS irreducible17x17 rejected without relaxation\n";
}
#endif
}
#endif
int main(int argc,char** argv){
    try{
#if __has_include("spin/detail/RankTwoCarrierKernels.h")
        nativeHandleReuse();affineDependencyReuse();thickness();interpolation();selectedConjugateFactors();
#ifdef RANK_TWO_ACTUAL_PARENT_AFFINE
        require(argc==2,"actual capture required");actualParentAffine(argv[1]);
#endif
        std::cout<<"PASS carrier thickness and exact moving-coordinate kernels; full carrier NOT_EXECUTED\n";
        return 0;
#else
        throw std::runtime_error("MISSING_FEATURE: exact carrier thickness and moving-coordinate kernels");
#endif
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
