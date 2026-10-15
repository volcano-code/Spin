#include "spin/QuadricCoefficientAuthority.h"

#include <gmp.h>
#include <mpfr.h>
#include <stdexcept>
#include <utility>
#include <vector>

namespace spin {
namespace {

constexpr mpfr_rnd_t rounding = MPFR_RNDN;

struct Real {
    mpfr_t value;
    Real() { mpfr_init2(value, QuadricCoefficientAuthority10::precisionBits); mpfr_set_zero(value, 0); }
    Real(const Real& other) : Real() { mpfr_set(value, other.value, rounding); }
    Real& operator=(const Real& other) { mpfr_set(value, other.value, rounding); return *this; }
    ~Real() { mpfr_clear(value); }
};

using Vector = std::array<Real, 10>;
using RationalVector = QuadricCoefficientAuthority10::RationalVector;

std::string integerText(mpz_srcptr value) {
    std::vector<char> buffer(mpz_sizeinbase(value, 10) + 3);
    mpz_get_str(buffer.data(), 10, value);
    return std::string(buffer.data());
}

std::string rationalText(mpq_srcptr value) {
    return integerText(mpq_numref(value)) + "/" + integerText(mpq_denref(value));
}

std::string parse(const std::string& text, Real* result) {
    mpq_t value;
    mpq_init(value);
    if (text.empty() || mpq_set_str(value, text.c_str(), 10) != 0 ||
        mpz_sgn(mpq_denref(value)) == 0) {
        mpq_clear(value);
        throw std::invalid_argument("coefficient or arc must be a finite rational with nonzero denominator");
    }
    mpq_canonicalize(value);
    const auto canonical = rationalText(value);
    mpfr_set_q(result->value, value, rounding);
    mpq_clear(value);
    if (!mpfr_number_p(result->value))
        throw std::invalid_argument("coefficient or arc exceeds the MPFR number range");
    return canonical;
}

Vector read(const RationalVector& rationals, RationalVector* canonical = nullptr) {
    Vector result;
    for (unsigned i = 0; i < 10; ++i) {
        const auto text = parse(rationals[i], &result[i]);
        if (canonical) canonical->at(i) = text;
    }
    return result;
}

Real norm(const Vector& input) {
    Real sum, product;
    for (const auto& x : input) {
        mpfr_mul(product.value, x.value, x.value, rounding);
        mpfr_add(sum.value, sum.value, product.value, rounding);
    }
    mpfr_sqrt(sum.value, sum.value, rounding);
    return sum;
}

void normalize(Vector* input) {
    const auto length = norm(*input);
    if (!mpfr_number_p(length.value) || mpfr_sgn(length.value) <= 0)
        throw std::invalid_argument("coefficient vector must be finite and nonzero");
    for (auto& x : *input) mpfr_div(x.value, x.value, length.value, rounding);
}

RationalVector write(const Vector& input) {
    RationalVector result;
    mpq_t value;
    mpq_init(value);
    for (unsigned i = 0; i < 10; ++i) {
        // mpfr_get_q captures the actual represented 1024-bit value exactly.
        mpfr_get_q(value, input[i].value);
        result[i] = rationalText(value);
    }
    mpq_clear(value);
    return result;
}

Vector dualFromGeometry(const Vector& geometry) {
    Vector dual = geometry;
    for (unsigned i = 7; i < 10; ++i) mpfr_mul_2ui(dual[i].value, dual[i].value, 1, rounding);
    normalize(&dual);
    return dual;
}

} // namespace

QuadricCoefficientAuthority10::QuadricCoefficientAuthority10(RationalVector geometry, RationalVector dual)
    : geometry_(std::move(geometry)), dual_(std::move(dual)) {}

QuadricCoefficientAuthority10 QuadricCoefficientAuthority10::fromNormalizedGeometry(
    const RationalVector& geometry) {
    RationalVector canonical;
    const auto q = read(geometry, &canonical);
    Real error, tolerance;
    const auto length = norm(q);
    mpfr_sub_ui(error.value, length.value, 1, rounding);
    mpfr_abs(error.value, error.value, rounding);
    parse("1/1000000000000", &tolerance);
    if (!mpfr_number_p(length.value) || mpfr_cmp(error.value, tolerance.value) > 0)
        throw std::invalid_argument("saved geometry coefficients must already have unit norm");
    return QuadricCoefficientAuthority10(std::move(canonical), write(dualFromGeometry(q)));
}

QuadricCoefficientAuthority10 QuadricCoefficientAuthority10::fromUnnormalizedGeometry(
    const RationalVector& geometry) {
    auto q = read(geometry);
    const auto length = norm(q);
    Real minimum;
    parse("1/1000000000000000", &minimum);
    if (!mpfr_number_p(length.value) || mpfr_cmp(length.value, minimum.value) <= 0)
        throw std::invalid_argument("raw geometry coefficient norm is below the existing minimum");
    normalize(&q);
    return QuadricCoefficientAuthority10(write(q), write(dualFromGeometry(q)));
}

QuadricCoefficientAuthority10 QuadricCoefficientAuthority10::fromSavedState(
    const RationalVector& geometry, const RationalVector& dual,
    const std::string& savedPolicy, unsigned savedPrecisionBits) {
    if (savedPolicy != policy || savedPrecisionBits != precisionBits)
        throw std::invalid_argument("saved coefficient policy or precision mismatch");
    // Validate the existing normalized-geometry contract without replacing
    // the geometry. A reload must preserve the exact represented dual, too.
    const auto validated = fromNormalizedGeometry(geometry);
    RationalVector canonicalDual;
    const auto s = read(dual, &canonicalDual);
    if (write(s) != canonicalDual)
        throw std::invalid_argument("saved dual is not exactly representable at the declared precision");
    const auto expected = read(validated.dual());
    Real tolerance, error;
    mpfr_set_ui_2exp(tolerance.value, 1, -1000, rounding);
    const auto length = norm(s);
    mpfr_sub_ui(error.value, length.value, 1, rounding);
    mpfr_abs(error.value, error.value, rounding);
    if (mpfr_cmp(error.value, tolerance.value) > 0)
        throw std::invalid_argument("saved dual must have unit norm at the declared precision");
    for (unsigned i = 0; i < 10; ++i) {
        mpfr_sub(error.value, s[i].value, expected[i].value, rounding);
        mpfr_abs(error.value, error.value, rounding);
        if (mpfr_cmp(error.value, tolerance.value) > 0)
            throw std::invalid_argument("saved geometry and dual are inconsistent");
    }
    return QuadricCoefficientAuthority10(validated.geometry(), std::move(canonicalDual));
}

QuadricCoefficientAuthority10 QuadricCoefficientAuthority10::geodesicUpdate(
    const RationalVector& dualDirection, const std::string& rationalArcLength) const {
    Real arc;
    parse(rationalArcLength, &arc);
    if (mpfr_sgn(arc.value) < 0) throw std::invalid_argument("geodesic arc must be nonnegative");
    const auto s = read(dual_);
    auto tangent = read(dualDirection);
    Real pairing, product;
    for (unsigned i = 0; i < 10; ++i) {
        mpfr_mul(product.value, s[i].value, tangent[i].value, rounding);
        mpfr_add(pairing.value, pairing.value, product.value, rounding);
    }
    for (unsigned i = 0; i < 10; ++i) {
        mpfr_mul(product.value, s[i].value, pairing.value, rounding);
        mpfr_sub(tangent[i].value, tangent[i].value, product.value, rounding);
    }
    normalize(&tangent);
    if (mpfr_zero_p(arc.value)) return *this;
    Real cosine, sine, first, second;
    mpfr_sin_cos(sine.value, cosine.value, arc.value, rounding);
    Vector updated;
    for (unsigned i = 0; i < 10; ++i) {
        mpfr_mul(first.value, cosine.value, s[i].value, rounding);
        mpfr_mul(second.value, sine.value, tangent[i].value, rounding);
        mpfr_add(updated[i].value, first.value, second.value, rounding);
    }
    normalize(&updated);
    Vector geometry = updated;
    for (unsigned i = 7; i < 10; ++i) mpfr_div_2ui(geometry[i].value, geometry[i].value, 1, rounding);
    normalize(&geometry);
    return QuadricCoefficientAuthority10(write(geometry), write(updated));
}

std::array<double, 10> QuadricCoefficientAuthority10::geometryDisplay() const {
    const auto q = read(geometry_);
    std::array<double, 10> result;
    for (unsigned i = 0; i < 10; ++i) result[i] = mpfr_get_d(q[i].value, rounding);
    return result;
}

} // namespace spin
