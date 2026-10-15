#pragma once

#include <array>
#include <string>

namespace spin {

/** Fixed-precision coefficient/dual state. This primitive does not certify
 * source geometry, derivatives, material or a Newton solution. */
class QuadricCoefficientAuthority10 {
public:
    using RationalVector = std::array<std::string, 10>;
    static constexpr unsigned precisionBits = 1024;
    static constexpr const char* policy = "COEFFICIENT_MPFR1024_V1";

    /** Preserve supplied normalized rational geometry; do not round-trip a
     * saved binary64 endpoint through a new normalization. */
    static QuadricCoefficientAuthority10 fromNormalizedGeometry(const RationalVector& geometry);
    static QuadricCoefficientAuthority10 fromUnnormalizedGeometry(const RationalVector& geometry);
    /** Reload both authoritative states without recomputing either. Reject a
     * policy/precision mismatch or an inconsistent geometry/dual pair. */
    static QuadricCoefficientAuthority10 fromSavedState(
        const RationalVector& geometry, const RationalVector& dual,
        const std::string& savedPolicy, unsigned savedPrecisionBits);
    QuadricCoefficientAuthority10 geodesicUpdate(
        const RationalVector& dualDirection, const std::string& rationalArcLength) const;

    const RationalVector& geometry() const noexcept { return geometry_; }
    const RationalVector& dual() const noexcept { return dual_; }
    std::array<double, 10> geometryDisplay() const;

private:
    QuadricCoefficientAuthority10(RationalVector geometry, RationalVector dual);
    RationalVector geometry_;
    RationalVector dual_;
};

} // namespace spin
