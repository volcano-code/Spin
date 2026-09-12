#pragma once
#include "spin2014/SourceCageBinding.h"
#include <limits>

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceConformingBoundaryConfig {
    std::size_t maximumFaceTreeNodes=10000000;
    std::size_t maximumPatches=2000000;
    std::size_t maximumPoints=4000000;
    std::size_t maximumTriangles=16000000;
};
/** Immutable conforming REST surfaces of the binding's actual Full cuboids.
 * IMPLEMENTATION_DEFAULT: sparse dyadic face overlay and shared center fans.
 * Each triangle is stored once in +axis orientation, referenced by primitive
 * IDs in binding order. This guarantees indexed primitive closure; it does not
 * certify global manifoldness, deformed legality, wall thickness or a solver.
 * Points retain original outer rows first; ALL internal indices are reassigned
 * as outerPointCount()+index into sorted unique materialPointKeys(), at D+1.
 * The original outer triangular topology remains through binding().source().
 * Positive checked caps include outer points. Resource failures throw
 * length_error with RESOURCE_LIMIT; invalid geometry/frame throws
 * invalid_argument; nonfinite arithmetic throws overflow_error.
 */
class SourceConformingBoundary {
public:
    using Key=SourceCageBinding::Key;
    using Triangle=std::array<std::size_t,3>;
    static constexpr std::size_t noOwner=std::numeric_limits<std::size_t>::max();
    struct Patch {
        int axis=0;
        // Plane is on SOURCE depth-D lattice; multiply by 2 for point keys.
        std::uint32_t plane=0;
        int level=0;
        // Tile coordinates at level in u=(axis+1)%3,v=(axis+2)%3.
        std::array<std::uint32_t,2> coordinates{};
        std::size_t negativeOwner=noOwner,positiveOwner=noOwner;
        std::size_t triangleBegin=0,triangleCount=0;
    };
    struct OrientedTriangle {
        std::size_t triangle=0;
        int orientation=1; // +1 as stored, -1 swaps triangle vertices 1 and 2.
    };
    SourceConformingBoundary(const SourceCageBinding&,const SourceConformingBoundaryConfig&);
    const SourceCageBinding& binding() const;
    const SourceConformingBoundaryConfig& config() const;
    const Eigen::MatrixXd& points() const;
    const std::vector<Key>& materialPointKeys() const;
    const std::vector<Patch>& patches() const;
    const std::vector<Triangle>& triangles() const;
    const std::vector<std::vector<OrientedTriangle>>& primitiveBoundaries() const;
    std::size_t outerPointCount() const;
    int latticeLevel() const;
    const std::string& configSha256() const;
    const std::string& outputSha256() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
};
}
#endif
