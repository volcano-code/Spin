#pragma once
#include <CGAL/CORE_Expr.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Gmpq.h>
#include <memory>
#include <string>
#include <vector>

namespace spin::detail {
using AlgebraicKernel=CGAL::Algebraic_kernel_d_1<CORE::BigInt>;
using AlgebraicReal=AlgebraicKernel::Algebraic_real_1;
using AlgebraicPolynomial=AlgebraicKernel::Polynomial_1;
using AlgebraicBivariatePolynomial=CGAL::Polynomial<AlgebraicPolynomial>;
class RankTwoNativeFiberWitness;
struct RankTwoNativeFiberAccess;

struct AlgebraicRootCertificate {
    std::string schema="rank-two-algebraic-root-v1";
    std::vector<std::string> coefficients; // primitive, ascending, leading > 0
    std::size_t rootIndex=0; // zero-based among distinct real roots
    std::string lower,upper; // exact canonical dyadics; strict open isolator
    std::string rootIdentity; // independent of isolator refinement
};

// Scoped, opaque proof minted ONLY by complete exact construction. No writable
// proof fields, caller serialization or process-global cache. Native handles
// can refine their exact isolators; their mathematical values remain fixed.
class RankTwoNativeRootSet {
public:
    static RankTwoNativeRootSet fromPolynomial(const AlgebraicPolynomial&);
    // Mint a scoped proof from ACTUAL complete elimination construction.
    // original receives the same full square-free polynomial as the old
    // resultant→square-free→fromPolynomial chain, never a selected factor.
    // No caller-authored component/proof input or process-global cache.
    static RankTwoNativeRootSet fromResultant(const AlgebraicBivariatePolynomial&,
        const AlgebraicBivariatePolynomial&,AlgebraicPolynomial& original);
    const std::vector<AlgebraicReal>& roots() const;
private:
    friend class CertifiedAlgebraicReal;
    friend struct RankTwoNativeFiberAccess;
    RankTwoNativeRootSet withNativeFiberWitness(
        std::shared_ptr<const RankTwoNativeFiberWitness>) const;
    struct Impl;
    explicit RankTwoNativeRootSet(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
};

// Arithmetic-bearing adapter. CGAL root handles themselves only support
// comparison; no approximate midpoint is substituted for a root.
class CertifiedAlgebraicReal {
public:
    static CertifiedAlgebraicReal fromRoot(const AlgebraicReal&,int maximumIsolationBits=512);
    static CertifiedAlgebraicReal fromRoot(const AlgebraicReal&,const RankTwoNativeRootSet&,
        int maximumIsolationBits=512);
    static CertifiedAlgebraicReal reload(const AlgebraicRootCertificate&,int maximumIsolationBits=512);
    const CORE::Expr& expression() const;
    const AlgebraicRootCertificate& certificate() const;
private:
    friend struct RankTwoNativeFiberAccess;
    std::shared_ptr<const RankTwoNativeFiberWitness> nativeFiberWitness() const;
    friend AlgebraicReal rankTwoAlgebraicHandle(const CertifiedAlgebraicReal&);
    // Only the sampler may reuse its proved affine dependency. No caller can
    // attach an arbitrary expression to a valid root certificate.
    friend CertifiedAlgebraicReal sampleRankTwoAlgebraicCoordinate(
        const CertifiedAlgebraicReal&,const CertifiedAlgebraicReal&,const CGAL::Gmpq&,int);
    friend CertifiedAlgebraicReal affineRankTwoAlgebraicCoordinate(
        const CertifiedAlgebraicReal&,const CertifiedAlgebraicReal&,const CGAL::Gmpq&,int);
    static CertifiedAlgebraicReal fromProvenAffineRoot(const AlgebraicReal&,
        const CertifiedAlgebraicReal&,const CertifiedAlgebraicReal&,const CGAL::Gmpq&,int,
        const RankTwoNativeRootSet* =nullptr);
    static CertifiedAlgebraicReal fromRootWithExpression(const AlgebraicReal&,int,const CORE::Expr*,
        const RankTwoNativeRootSet* =nullptr);
    struct Impl;
    explicit CertifiedAlgebraicReal(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
};

// A finite CORE escape threshold can cache UNKNOWN as zero. Do not silently
// change process-global configuration or reuse such a polluted expression.
void requireExactAlgebraicEnvironment();
// Complete exact distinct-real set, authenticated against the WHOLE input.
// Internal arithmetic helper only: no caller-provided proof, persisted cache
// or substituted output polynomial. Legacy builds retain native Solve_1.
std::vector<AlgebraicReal> rankTwoCompleteNativeRealRoots(const AlgebraicPolynomial&);
// Complete primitive square-free part of the ORIGINAL polynomial. Native
// candidates require exact decomposition/product AND an independently
// verified degree-preserving modular Bezout certificate. Otherwise use the
// original exact CGAL operation; no approximate or selected-factor result.
AlgebraicPolynomial rankTwoSquareFreeNativePolynomial(const AlgebraicPolynomial&);
// Private arithmetic operand only. Authenticates the FULL frozen parent
// factor product, then exact membership of this already certified value on a
// detached handle. Never rewrites its original certificate/root ordinal,
// accepts caller-selected factors, or supplies a reduced whole-input proof.
// Non-native builds retain the original polynomial and resource behavior.
AlgebraicPolynomial rankTwoAuthenticatedParentComponent(const CertifiedAlgebraicReal&);
// Raw resultant eliminating the OUTER power of the nested CGAL polynomial.
// Preserve ALL integer coefficients, content and sign. Native builds use the
// fixed exact backend only in the bounded ordinary-fiber degree scope; failed
// or unsupported calculations retain the original exact CGAL operation.
AlgebraicPolynomial rankTwoCompleteNativeResultant(
    const AlgebraicBivariatePolynomial&,const AlgebraicBivariatePolynomial&);
} // namespace spin::detail
