#include "spin/QuadricStratumTangent.h"

#include <Eigen/SVD>
#include <iostream>
#include <stdexcept>

namespace {
using V = spin::QuadricCoefficients10::Vector;
void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}
void check(const char* name, const V& raw,
           spin::QuadricSurfaceSubtype expected, int dimension) {
    const auto q = spin::QuadricCoefficients10::fromUnnormalized(raw);
    const auto tangent = spin::evaluateQuadricStratumTangent(q);
    require(tangent.defined, std::string(name)+": undefined");
    require(tangent.subtype == expected, std::string(name)+": wrong subtype");
    require(tangent.tangentDimension == dimension, std::string(name)+": wrong dimension");
    require(tangent.constraintRows.cols() == 10 &&
            tangent.constraintRows.allFinite(), std::string(name)+": malformed rows");
    require(tangent.basis.rows() == 10 && tangent.basis.cols() == dimension &&
            tangent.basis.allFinite(), std::string(name)+": malformed tangent basis");
    require((tangent.basis.transpose() * tangent.basis -
             Eigen::MatrixXd::Identity(dimension, dimension)).norm() <= 1e-10,
            std::string(name)+": tangent basis is not orthonormal");
    require((tangent.constraintRows * tangent.basis).norm() <= 1e-10,
            std::string(name)+": tangent basis violates coefficient constraints");
    require(tangent.contains(V::Zero()), std::string(name)+": zero direction rejected");
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        tangent.constraintRows, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const int rank = static_cast<int>((svd.singularValues().array() > 1e-11).count());
    require(rank + tangent.tangentDimension == 10, std::string(name)+": rank/dimension mismatch");
    const double h = 1e-7;
    const auto translated = q.pullback(Eigen::Matrix3d::Identity(),
                                       Eigen::Vector3d(0.17, -0.11, 0.09) * h);
    const V translationDirection = (translated.values() - q.values()) / h;
    require(tangent.contains(translationDirection, 1e-6),
            std::string(name)+": translation left the fixed stratum tangent");
    std::cout << name << " subtype=" << quadricSurfaceSubtypeToken(tangent.subtype)
              << " constraints=" << rank << " dimension=" << tangent.tangentDimension
              << " diagnostic=" << tangent.diagnostic << '\n';
}
}

int main() {
    try {
        V plane; plane << 0, 1, 2, 3, 0, 0, 0, 0, 0, 0;
        V cylinder; cylinder << -1, 0, 0, 0, 1, 2, 0, 0, 0, 0;
        V paraboloid; paraboloid << 0, 0, 0, -1, 1, 2, 0, 0, 0, 0;
        V ellipsoid; ellipsoid << -1, 0, 0, 0, 1, 2, 3, 0, 0, 0;
        check("plane", plane, spin::QuadricSurfaceSubtype::Plane, 3);
        check("elliptic-cylinder", cylinder, spin::QuadricSurfaceSubtype::EllipticCylinder, 7);
        check("paraboloid", paraboloid, spin::QuadricSurfaceSubtype::Paraboloid, 8);
        check("ellipsoid", ellipsoid, spin::QuadricSurfaceSubtype::Ellipsoid, 9);
        V cone = V::Zero();
        cone(4) = 1.0;
        cone(5) = 1.0;
        cone(6) = -1.0;
        V hyperbolicCylinder = V::Zero();
        hyperbolicCylinder(0) = -1.0;
        hyperbolicCylinder(4) = 1.0;
        hyperbolicCylinder(5) = -1.0;
        check("hyperbolic-cylinder", hyperbolicCylinder,
              spin::QuadricSurfaceSubtype::HyperbolicCylinder, 7);
        V intersectingPlanes = hyperbolicCylinder;
        intersectingPlanes(0) = 0.0;
        for (const auto& unsupported : {cone, intersectingPlanes}) {
            bool tangentRejected = false;
            try {
                (void)spin::evaluateQuadricStratumTangent(
                    spin::QuadricCoefficients10::fromUnnormalized(unsupported));
            } catch (const std::invalid_argument&) {
                tangentRejected = true;
            }
            require(tangentRejected,
                    "unsupported degenerate type must not fabricate a fixed-stratum tangent basis");
        }
        std::cout << "[PASS] fixed quadric stratum tangent dimensions\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
