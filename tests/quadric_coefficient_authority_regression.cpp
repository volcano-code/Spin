#include "spin/QuadricCoefficientAuthority.h"
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>
#include <string>

using Authority = spin::QuadricCoefficientAuthority10;

static void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

static Authority::RationalVector axis(int index) {
    Authority::RationalVector out;
    out.fill("0/1");
    out.at(index) = "1/1";
    return out;
}

int main() {
    try {
        const auto initial = Authority::fromNormalizedGeometry(axis(0));
        const auto changed = initial.geodesicUpdate(axis(1), "1/10000000000000000000000000");
        // A binary64 cosine rounds to one at this step. Authority must retain
        // both the positive new coefficient and the quadratic norm correction.
        require(mpq_class(changed.geometry()[0]) < 1, "sub-double norm correction was lost");
        const mpq_class h("1/10000000000000000000000000");
        const mpq_class observed(changed.geometry()[1]);
        require(observed > h * mpq_class(999, 1000) &&
                observed < h * mpq_class(1001, 1000), "small geodesic update was lost");
        require(changed.geometryDisplay()[0] == 1., "display is not the binary64 derivative");
        require(changed.geometry() != initial.geometry(), "authority was replaced by its display");
        const auto restored = Authority::fromSavedState(
            changed.geometry(), changed.dual(), Authority::policy, Authority::precisionBits);
        require(restored.geometry() == changed.geometry() && restored.dual() == changed.dual(),
                "saved authority did not preserve geometry and dual bits");
        require(restored.geodesicUpdate(axis(1), "1/10000000000000000000000000").geometry() ==
                changed.geodesicUpdate(axis(1), "1/10000000000000000000000000").geometry(),
                "reloaded state changed the subsequent geodesic");
        auto wrongDual = changed.dual(); wrongDual[2] = "1/100";
        bool savedRejected = false;
        try { Authority::fromSavedState(changed.geometry(), wrongDual, Authority::policy, 1024); }
        catch (const std::invalid_argument&) { savedRejected = true; }
        require(savedRejected, "mismatched saved dual was admitted");
        for (const auto bits : {512u, 2048u}) {
            savedRejected = false;
            try { Authority::fromSavedState(changed.geometry(), changed.dual(), Authority::policy, bits); }
            catch (const std::invalid_argument&) { savedRejected = true; }
            require(savedRejected, "wrong saved precision was admitted");
        }
        savedRejected = false;
        try { Authority::fromSavedState(changed.geometry(), changed.dual(), "LEGACY_DOUBLE", 1024); }
        catch (const std::invalid_argument&) { savedRejected = true; }
        require(savedRejected, "wrong saved coefficient policy was admitted");
        auto unrepresentableDual = initial.dual();
        const mpq_class hidden = mpq_class(1) + mpq_class(1, mpz_class(1) << 1100);
        unrepresentableDual[0] = hidden.get_str();
        savedRejected = false;
        try { Authority::fromSavedState(initial.geometry(), unrepresentableDual, Authority::policy, 1024); }
        catch (const std::invalid_argument&) { savedRejected = true; }
        require(savedRejected, "saved dual was silently rounded to fit the claimed precision");
        require(initial.geodesicUpdate(axis(1), "0/1").geometry() == initial.geometry(),
                "zero arc changed saved authority");

        auto weighted = axis(0);
        weighted[0] = "3/5";
        weighted[9] = "4/5";
        const auto cross = Authority::fromNormalizedGeometry(weighted);
        const mpq_class pairing = mpq_class(cross.dual()[9]) / mpq_class(cross.dual()[0]);
        const mpq_class error = abs(pairing - mpq_class(8, 3));
        require(error < mpq_class(1, mpz_class(1) << 1000), "dual cross-term factor is wrong");
        require(cross.geometry() == weighted, "normalized input bits were not preserved");
        auto scaled = axis(0); scaled[0] = "2/1";
        require(Authority::fromUnnormalizedGeometry(scaled).geometry() == axis(0),
                "raw coefficient normalization is wrong");
        auto tiny = axis(0); tiny[0] = "1/10000000000000000";
        bool tinyRejected = false;
        try { Authority::fromUnnormalizedGeometry(tiny); }
        catch (const std::invalid_argument&) { tinyRejected = true; }
        require(tinyRejected, "legacy raw-coefficient minimum norm was relaxed");

        auto zero = axis(0); zero[0] = "0/1";
        for (const auto& invalid : {zero, axis(0)}) {
            auto q = invalid;
            if (q == axis(0)) q[0] = "2/1";
            bool rejected = false;
            try { Authority::fromNormalizedGeometry(q); }
            catch (const std::invalid_argument&) { rejected = true; }
            require(rejected, "invalid normalized input was admitted");
        }
        bool rejected = false;
        try { initial.geodesicUpdate(axis(1), "-1/2"); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "negative arc was admitted");
        rejected = false;
        try { initial.geodesicUpdate(axis(0), "1/2"); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "zero tangent was admitted");
        auto invalid = axis(0); invalid[0] = "1/0";
        rejected = false;
        try { Authority::fromNormalizedGeometry(invalid); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "zero denominator was admitted");
        std::cout << "PASS: MPFR1024 coefficients, dual pairing, geodesic, invalid inputs\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
