#pragma once
#include "spin/QuadricCoefficients.h"
#include "spin/QuadricInverseGaussMap.h"
#include <CGAL/Gmpq.h>
#include <memory>
#include <string>
#include <vector>
namespace spin {
/** RegularParaboloid is an explicit rank-two analytic normal-edge extension,
 * not the paper's full-rank projective mesh or a complete face construction. */
enum class TruncationRootProfile { FullRank, RegularParaboloid };
struct TruncationRootInterval {
    CGAL::Gmpq lower, upper;
};
/** Complete regular positive-branch radius-root isolation on n(t)=
 * (1-t)*firstNormal+t*secondNormal, 0<t<1. Each disjoint interval
 * contains exactly one polynomial root, certified admissible. Empty with
 * certified=true proves no crossings; certified=false is UNKNOWN, never
 * an empty-edge certificate. Tangencies and radius-root endpoints remain
 * rejected; regular-paraboloid ideal normal endpoints may be supported. */
struct TruncationRootIsolation {
    bool certified = false;
    std::vector<TruncationRootInterval> intervals;
    std::string diagnostic;
};
/** Complete one-dimensional finite-radius membership partition. roots are
 * ordered along the supplied n(t); inside has roots.size()+1 entries for the
 * open regions between the exact algebraic roots, not just the bracket gaps.
 * Root brackets are open, have nonroot rational endpoints, and each contain
 * one simple admissible root. For persistent root IDs use canonical endpoint
 * order; reversing the edge reverses both roots and inside.
 *
 * This does NOT determine connectivity inside a normal-sphere triangle.
 * certified=false is UNKNOWN and carries no partial partition. */
struct TruncationEdgePartition {
    bool certified = false;
    std::vector<TruncationRootInterval> roots;
    std::vector<bool> inside;
    std::string diagnostic;
};
/** Local derivative of the unique positive-branch analytic radius crossing.
 * The source point remains the original 384-bit enclosure/representative;
 * a separate, narrower proof enclosure certifies its implicit differential.
 * This is not a derivative of rational rounding or an entire source mesh. */
struct TruncationRationalDifferentialResult {
    using Matrix=QuadricRationalFiniteDifferentialResult::Matrix;
    bool defined=false;
    QuadricRationalPointResult point;
    Matrix coefficientJacobian,lower,upper;
    std::array<std::string,3> derivativePointLower,derivativePointUpper;
    std::string diagnostic;
};
/** Reusable exact certificate context for many normal edges on one fixed
 * quadric and truncation radius. Quadric-only rational algebra is prepared
 * once; every edge still receives its own exact Sturm-Tarski proof. */
class TruncationRootCertifier {
public:
    TruncationRootCertifier(const QuadricCoefficients10&, double radius);
    /** Opt-in full-rank coefficient authority. Radius remains the frozen
     * represented input; coefficient preparation never uses its double view. */
    TruncationRootCertifier(const QuadricCoefficientAuthority10&, double radius);
    TruncationRootCertifier(const QuadricCoefficients10&, double radius,
        TruncationRootProfile profile);
    bool certify(const Eigen::Vector3d& firstNormal,
                 const Eigen::Vector3d& secondNormal,
                 std::string* diagnostic=nullptr) const;
    TruncationRootIsolation isolate(const Eigen::Vector3d& firstNormal,
                                   const Eigen::Vector3d& secondNormal,
                                   int maximumDepth=64) const;
    TruncationEdgePartition partition(const Eigen::Vector3d& firstNormal,
                                     const Eigen::Vector3d& secondNormal,
                                     int maximumDepth=64) const;
    /** One regular full-rank crossing, using the existing exact isolator's
     * Sturm-Tarski admissibility oracle and bounded coordinate refinement.
     * Nearby inadmissible squared-equation roots need not first be separated
     * by the legacy all-polynomial-root isolation interface.
     * Encloses each analytic coordinate to absolute width <=2^-384 and
     * returns a rational approximation. Not an exact analytic root; no
     * binary64 root seed, parameterization change or rank-two fallback. */
    QuadricRationalPointResult rationalPoint(const Eigen::Vector3d& firstNormal,
                                          const Eigen::Vector3d& secondNormal,
                                          bool compactRationalCoordinates = false) const;
    TruncationRationalDifferentialResult rationalDifferential(
        const Eigen::Vector3d& firstNormal,const Eigen::Vector3d& secondNormal,
        bool compactRationalCoordinates = false) const;
private:
    QuadricRationalPointResult rationalPointImpl(
        const Eigen::Vector3d&,const Eigen::Vector3d&,bool,unsigned enclosureBits) const;
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};

/** Sufficient exact-rational certificate: exactly one admissible positive
 * inverse-Gauss radius root on the open normal edge. Conservative rejection
 * of unresolved endpoint, multiplicity, interval or arithmetic budgets.
 */
bool certifySingleTruncationRoot(const QuadricCoefficients10&,
    const Eigen::Vector3d& firstNormal,const Eigen::Vector3d& secondNormal,
    double radius,std::string* diagnostic=nullptr);
}
