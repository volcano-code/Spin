#pragma once
#include <mpfr.h>
#include <gmpxx.h>
#include <utility>
#include <stdexcept>

// Outward rounded intervals; every operation encloses the exact rational
// expression. No ambient floating-point rounding mode is changed.
struct CertifiedInterval {
    mpfr_t lo, hi;
    CertifiedInterval(long value=0) {
        mpfr_init2(lo,256); mpfr_init2(hi,256);
        mpfr_set_si(lo,value,MPFR_RNDD); mpfr_set_si(hi,value,MPFR_RNDU);
    }
    explicit CertifiedInterval(const mpq_class &value) : CertifiedInterval() {
        mpfr_set_q(lo,value.get_mpq_t(),MPFR_RNDD);
        mpfr_set_q(hi,value.get_mpq_t(),MPFR_RNDU);
    }
    CertifiedInterval(const CertifiedInterval &v) : CertifiedInterval() {
        mpfr_set(lo,v.lo,MPFR_RNDD); mpfr_set(hi,v.hi,MPFR_RNDU);
    }
    CertifiedInterval(CertifiedInterval &&v) noexcept : CertifiedInterval() {
        mpfr_swap(lo,v.lo); mpfr_swap(hi,v.hi);
    }
    CertifiedInterval &operator=(CertifiedInterval v) {
        mpfr_swap(lo,v.lo); mpfr_swap(hi,v.hi); return *this;
    }
    ~CertifiedInterval() { mpfr_clear(lo); mpfr_clear(hi); }
};
inline CertifiedInterval operator+(const CertifiedInterval &a,const CertifiedInterval &b) {
    CertifiedInterval r; mpfr_add(r.lo,a.lo,b.lo,MPFR_RNDD);
    mpfr_add(r.hi,a.hi,b.hi,MPFR_RNDU); return r;
}
inline CertifiedInterval operator*(const CertifiedInterval &a,const CertifiedInterval &b) {
    CertifiedInterval r; mpfr_t lower,upper;
    mpfr_init2(lower,256);mpfr_init2(upper,256);
    mpfr_set_inf(r.lo,1);mpfr_set_inf(r.hi,-1);
    for(int i=0;i<2;++i)for(int j=0;j<2;++j) {
        mpfr_mul(lower,i?a.hi:a.lo,j?b.hi:b.lo,MPFR_RNDD);
        mpfr_mul(upper,i?a.hi:a.lo,j?b.hi:b.lo,MPFR_RNDU);
        if(mpfr_less_p(lower,r.lo))mpfr_set(r.lo,lower,MPFR_RNDD);
        if(mpfr_greater_p(upper,r.hi))mpfr_set(r.hi,upper,MPFR_RNDU);
    }
    mpfr_clear(lower);mpfr_clear(upper);return r;
}
inline CertifiedInterval operator/(const CertifiedInterval &a,unsigned long b) {
    if(!b)throw std::invalid_argument("ZERO_INTERVAL_DIVISOR");
    CertifiedInterval r;mpfr_div_ui(r.lo,a.lo,b,MPFR_RNDD);
    mpfr_div_ui(r.hi,a.hi,b,MPFR_RNDU);return r;
}
