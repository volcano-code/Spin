#pragma once
#include "RankTwoNormalEvents.h"
#include <optional>

namespace spin {
enum class RankTwoRegionKind { Outside,Empty,Cut,Full };
enum class RankTwoArcEndLocation { Finite,LeftInfinity,RightInfinity,BottomInfinity,TopInfinity };
struct RankTwoVertexRecord {
    std::string identity;
    std::array<detail::AlgebraicRootCertificate,2> coordinates;
    bool isolated=false;
};
struct RankTwoArcEndRecord {
    RankTwoArcEndLocation location=RankTwoArcEndLocation::Finite;
    std::optional<std::string> vertexIdentity;
    std::optional<detail::AlgebraicRootCertificate> finiteAbscissa;
};
struct RankTwoArcRecord {
    std::string identity,supportPolynomial,leftRegion,rightRegion;
    bool vertical=false;
    int branch=-1;
    std::array<RankTwoArcEndRecord,2> ends; // increasing lexicographic direction
    std::array<detail::AlgebraicRootCertificate,2> interiorWitness;
    std::vector<RankTwoPolynomialOwner> owners;
};
struct RankTwoDirectedArcRecord {
    std::string identity;
    bool reversed=false,atInfinity=false;
};
struct RankTwoRegionRecord {
    std::string identity;
    RankTwoRegionKind kind=RankTwoRegionKind::Outside;
    std::array<CGAL::Gmpq,2> witness;
    std::array<int,3> originalSigns{};
    std::vector<int> normalSigns;
    bool unbounded=false;
    std::vector<std::vector<RankTwoDirectedArcRecord>> outerLoops,innerLoops;
    std::vector<std::string> isolatedVertices;
};
struct RankTwoArrangementRecords {
    std::array<std::array<CGAL::Gmpq,25>,3> originalFhD;
    std::vector<std::array<CGAL::Gmpq,25>> normalHalfplanes;
    std::vector<RankTwoVertexRecord> vertices;
    std::vector<RankTwoArcRecord> arcs;
    std::vector<RankTwoRegionRecord> regions;
};
const RankTwoArrangementRecords& arrangementRecords(const RankTwoArrangement&);

// Reconstructs against the supplied authority, cone and cached global DCEL.
// No caller-supplied counts, flags or IDs can grant validity. This is not yet
// a cheap deserializer: replay is deliberately required at this stage.
RankTwoArrangement reloadRankTwoNormalFace(const PreparedRankTwoSource&,
    const RankTwoArrangement& global,const RankTwoNormalFace&,const std::string& saved,
    const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
} // namespace spin
