#pragma once
#include "verification/ExactMaterialRegion.h"

namespace verification {
struct ExactCuboid { Eigen::Vector3d lower,upper; };
/** IMPLEMENTATION_DEFAULT resource guards, not author or allocator limits. */
struct ExactCuboidUnionSettings {
    std::size_t maximumInputCuboids=10000;
    std::size_t maximumOutputVertices=1000000;
    std::size_t maximumOutputTriangles=2000000;
    std::size_t maximumIntermediateComplexity=4000000;
};
struct ExactNefComplexity {
    std::size_t vertices=0,halfedges=0,halffacets=0,volumes=0;
    std::size_t total() const {return vertices+halfedges+halffacets+volumes;}
};
struct ExactCuboidUnionFacts {
    std::string algorithm="exact-balanced-cuboid-union-binary64-roundtrip-v1";
    std::string backend=ExactMaterialRegionFacts{}.backend;
    std::string complexityDefinition="vertices+halfedges+halffacets+volumes; live explicit union/roundtrip Nef operands/results including imported region; excludes existing auditor internal temporaries; not allocator cap";
    std::string status="INCOMPLETE_UNKNOWN",diagnostic,inputDigest,outputDigest,meshDigest;
    ExactCuboidUnionSettings settings;
    ExactNefComplexity finalComplexity;
    std::size_t inputCuboids=0,booleanOperations=0,peakComplexity=0;
    std::optional<bool> simple,roundtripEqual;
};
/** Immutable regular-closed algebra product. Empty union is explicit empty,
 * never physical/solver success. No lossy coordinate conversion or repairs.
 * Failed operation exposes diagnostics only, no partial mesh or region.
 */
class ExactCuboidUnion {
public:
    const ExactCuboidUnionFacts& facts() const {return facts_;}
    const Eigen::MatrixXd& vertices() const {return vertices_;}
    const Eigen::MatrixXi& faces() const {return faces_;}
    const std::shared_ptr<const ExactMaterialRegion>& region() const {return region_;}
private:
    ExactCuboidUnion()=default;
    ExactCuboidUnionFacts facts_;
    Eigen::MatrixXd vertices_{0,3};
    Eigen::MatrixXi faces_{0,3};
    std::shared_ptr<const ExactMaterialRegion> region_;
    friend ExactCuboidUnion makeExactCuboidUnion(const std::vector<ExactCuboid>&,const ExactCuboidUnionSettings&);
};
void validateExactCuboidUnionSettings(const ExactCuboidUnionSettings&);
std::string hashExactCuboidUnionSettings(const ExactCuboidUnionSettings&);
// Malformed settings/coordinates throw invalid_argument. Operational failures
// are explicit UNKNOWN/NONMANIFOLD/NONREPRESENTABLE diagnostics, not infeasible.
ExactCuboidUnion makeExactCuboidUnion(const std::vector<ExactCuboid>&,const ExactCuboidUnionSettings& = {});
}
