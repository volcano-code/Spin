#pragma once
// TEST-ONLY arithmetic feasibility probe, not a production provider.
// Full factors, including complex-only factors, content and multiplicities
// contribute to the resultant. No selected/minimal factor substitutes p.
#include <flint.h>
#include <fmpz.h>
#include <fmpz_poly.h>
#include <fmpz_poly_factor.h>
#include <memory>
#include <sstream>
namespace rank_two_test {
using Z=CORE::BigInt;
using A=spin::detail::AlgebraicPolynomial;
using B=CGAL::Polynomial<A>;
using BT=CGAL::Polynomial_traits_d<B>;
struct NativePoly {
    fmpz_poly_t value;
    NativePoly(){fmpz_poly_init(value);}
    ~NativePoly(){fmpz_poly_clear(value);}
};
struct NativeFactors {
    fmpz_poly_factor_t value;
    NativeFactors(){fmpz_poly_factor_init(value);}
    ~NativeFactors(){fmpz_poly_factor_clear(value);}
};
inline Z integer(const fmpz* value){
    char* data=fmpz_get_str(nullptr,10,value);
    if(!data)throw std::runtime_error("probe integer allocation failed");
    std::unique_ptr<char,decltype(&flint_free)> owned(data,&flint_free);return Z(data);
}
inline A squareFreeProbe(const A& input,const char* label){
    if(std::string(FLINT_VERSION)!="3.3.1"||flint_get_num_threads()!=1||input.degree()<1)
        throw std::runtime_error("probe square-free input/backend contract mismatch");
    NativePoly p;
    for(int i=0;i<=input.degree();++i){
        std::ostringstream text;text<<input[i];
        fmpz_t coefficient;fmpz_init(coefficient);
        const int bad=fmpz_set_str(coefficient,text.str().c_str(),10);
        if(bad){fmpz_clear(coefficient);throw std::runtime_error("probe square-free integer parse failed");}
        fmpz_poly_set_coeff_fmpz(p.value,i,coefficient);fmpz_clear(coefficient);
    }
    NativeFactors factors;fmpz_poly_factor_squarefree(factors.value,p.value);
    const Z content=integer(&factors.value->c);
    A product(content),candidate(1);
    std::cout<<"SQUAREFREE_PROBE_CONTENT "<<label<<' '<<content<<'\n';
    for(slong i=0;i<factors.value->num;++i){
        const int degree=fmpz_poly_degree(factors.value->p+i),multiplicity=factors.value->exp[i];
        if(degree<1||multiplicity<1||degree>input.degree()||multiplicity>input.degree())
            throw std::runtime_error("probe invalid square-free degree/multiplicity");
        std::vector<Z> coefficients;
        for(int j=0;j<=degree;++j)coefficients.push_back(integer(factors.value->p[i].coeffs+j));
        const A factor(coefficients.begin(),coefficients.end());
        product*=CGAL::ipower(factor,multiplicity);candidate*=factor;
        std::cout<<"SQUAREFREE_PROBE_FACTOR "<<label<<' '<<multiplicity<<' '<<coefficients.size();
        for(const auto& c:coefficients)std::cout<<' '<<c;std::cout<<'\n';
    }
    if(product!=input||candidate.degree()<1)
        throw std::runtime_error("probe square-free decomposition does not equal original input");
    // TEST ONLY: independent checker and exact original-CGAL comparison must
    // additionally establish square-freeness. Product identity alone is NOT
    // a generic square-free certificate, nor authority to change production.
    std::cout<<"SQUAREFREE_PROBE_PRODUCT_VERIFIED "<<label<<' '<<input.degree()<<'\n';
    return candidate;
}
inline void squareFreeLiteralChecks(){
    // Literal primitive square-free parts, derived before implementation.
    // Missing multiplicity groups, signed content or complex-only factors
    // must not change the complete root set.
    if(squareFreeProbe(A(-12,0,0,0,9,0,-3),"signed-repeated-complex")!=A(-2,0,-1,0,1))
        throw std::runtime_error("square-free candidate lost a complete factor");
    if(squareFreeProbe(A(-2,6,-8,8,-6,2),"unequal-multiplicities")!=A(-1,1,-1,1))
        throw std::runtime_error("square-free candidate lost unequal multiplicities");
    if(squareFreeProbe(A(4,-4,-8,8,4,-4),"quadratic-resultant")!=A(-1,0,1))
        throw std::runtime_error("square-free candidate changed conjugate resultant roots");
    std::cout<<"PASS literal square-free: signed content, complete factors, unequal multiplicities\n";
}
inline A factorResultantProbe(const B& swapped,const A& input,const char* label){
    if(std::string(FLINT_VERSION)!="3.3.1"||flint_get_num_threads()!=1)
        throw std::runtime_error("probe factor backend contract mismatch");
    NativePoly p;
    for(int i=0;i<=input.degree();++i){
        std::ostringstream text;text<<input[i];
        fmpz_t coefficient;fmpz_init(coefficient);
        const int bad=fmpz_set_str(coefficient,text.str().c_str(),10);
        if(bad){fmpz_clear(coefficient);throw std::runtime_error("probe integer parse failed");}
        fmpz_poly_set_coeff_fmpz(p.value,i,coefficient);fmpz_clear(coefficient);
    }
    NativeFactors factors;fmpz_poly_factor(factors.value,p.value);
    const Z content=integer(&factors.value->c);
    A product(content),result(CGAL::ipower(content,swapped.degree()));
    std::cout<<"RESULTANT_PROBE_CONTENT "<<label<<' '<<content<<' '<<swapped.degree()<<'\n';
    for(slong i=0;i<factors.value->num;++i){
        const int degree=fmpz_poly_degree(factors.value->p+i),multiplicity=factors.value->exp[i];
        if(degree<1||multiplicity<1||degree>input.degree()||multiplicity>input.degree())
            throw std::runtime_error("probe invalid factor degree/multiplicity");
        std::vector<Z> coefficients;
        for(int j=0;j<=degree;++j)coefficients.push_back(integer(factors.value->p[i].coeffs+j));
        const A factor(coefficients.begin(),coefficients.end());
        const A part=CGAL::resultant(swapped,BT::Swap()(B(factor),0,1));
        product*=CGAL::ipower(factor,multiplicity);result*=CGAL::ipower(part,multiplicity);
        std::cout<<"RESULTANT_PROBE_FACTOR "<<label<<' '<<multiplicity<<' '<<coefficients.size();
        for(const auto& c:coefficients)std::cout<<' '<<c;std::cout<<'\n';
    }
    if(product!=input)throw std::runtime_error("probe factor product does not equal original input");
    std::cout<<"RESULTANT_PROBE_PRODUCT_VERIFIED "<<label<<' '<<input.degree()<<'\n';
    return result;
}
inline void factorResultantLiteralChecks(){
    const B x=CGAL::shift(B(1),1,0),y=CGAL::shift(B(1),1,1);
    // Hand-derived -3*(x^2-2)^2*(x^2+1), includes complex-only factor.
    const A first(-12,0,0,0,9,0,-3);
    const auto one=factorResultantProbe(BT::Swap()(x-y,0,1),first,"linear-content-complex");
    if(one!=first)throw std::runtime_error("linear resultant lost signed content, complex factor or multiplicity");
    // 2*(x-1)^3*(x^2+1). Res_x(x^2-y,p)=4*(1-y)^3*(1+y)^2.
    const A second(-2,6,-8,8,-6,2),expected(4,-4,-8,8,4,-4);
    const auto two=factorResultantProbe(BT::Swap()(x*x-y,0,1),second,"quadratic-content-power");
    if(two!=expected)throw std::runtime_error("resultant content exponent or repeated factor wrong");
    std::cout<<"PASS literal resultant: signed content, complex-only factor and multiplicities\n";
}
} // namespace rank_two_test
