#include "spin/QuadricMomentEvaluator.h"

#include <array>
#include <iostream>
#include <stdexcept>

namespace {
spin::TriangleMesh tetrahedron(double scale) {
    spin::TriangleMesh mesh;
    mesh.vertices.resize(4,3);
    mesh.vertices << 0,0,0, scale,0,0, 0,scale,0, 0,0,scale;
    mesh.faces.resize(4,3);
    mesh.faces << 0,2,1, 0,1,3, 0,3,2, 1,2,3;
    return mesh;
}
}

int main() {
    try {
        spin::SpinDesignRequest request;
        request.outerSurface=tetrahedron(2);
        request.designDomain=tetrahedron(.5);
        request.designDomain.vertices.array()+=.25;
        spin::QuadricCoefficients10::Vector coefficients;
        coefficients << 0,0,0,-1,1,1,0,0,0,0;
        const auto quadric=spin::QuadricCoefficients10::fromUnnormalized(coefficients);
        std::cout << "fixed_quadric=" << quadric.hash() << '\n';
        std::array<spin::QuadricMomentEvaluatorSettings,4> settings;
        settings[0].triangulation.maximumRadiusDoublings=8;
        settings[2].triangulation.truncationRadius=16;
        settings[3].derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteReducedStratum;
        // Frozen public-API outputs before the source-policy hash repair.
        const std::array<const char*,4> historicalHashes={{
            "003259bfe7e738f42e3025619eef315b4080bfa7bcee368393445e5c2bf2668a",
            "219a73cbb48c7b8dab46bba67ac043a96c8ae7aaac902e11b1796548ff4d5e2c",
            "964086ca5c58e5f0a939d72182b86ee72835a2993953ddd3ac8e41b4a9f32e4f",
            "be8aa9875dd9686494280bdce61397c192b3bb7913df3975f99a70e597ca766b"
        }};
        bool collision=false;
        for(std::size_t i=0;i<settings.size();++i) {
            auto& policy=settings[i];
            const auto baseline=spin::computeQuadricEvaluatorRequestHash(request,policy);
            if(baseline!=historicalHashes[i])
                throw std::runtime_error("disabled source policy changed a historical hash");
            if(i==0 && baseline!=spin::computeSpinDesignRequestHash(request))
                throw std::runtime_error("historical request compatibility changed");
            policy.triangulation.fixedNormalParaboloid=true;
            const auto fixed=spin::computeQuadricEvaluatorRequestHash(request,policy);
            std::cout << "variant=" << i << " off=" << baseline << " on=" << fixed << '\n';
            collision|=baseline==fixed;
            if(fixed!=spin::computeQuadricEvaluatorRequestHash(request,policy))
                throw std::runtime_error("fixed-normal request identity is nondeterministic");
            policy.triangulation.fixedNormalParaboloid=false;
            if(baseline!=spin::computeQuadricEvaluatorRequestHash(request,policy))
                throw std::runtime_error("disabled source policy does not restore identity");
        }
        if(collision)
            throw std::runtime_error("fixedNormalParaboloid omitted from public evaluator request hash");
        std::cout << "[PASS] fixed-normal source policy has a distinct public request identity\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
