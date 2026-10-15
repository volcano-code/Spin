#pragma once
#include "RankTwoNormalMesh.h"
#include "RankTwoNormalEvents.h"
#include "RankTwoNormalVertices.h"
#include "RankTwoSharedIntervals.h"

namespace spin {
struct RankTwoSharedEdgeRecord {
    std::array<std::size_t,2> edge{};
    std::array<std::size_t,2> incidentFaces{};
    std::vector<RankTwoNormalBoundaryEvent> events;
    std::vector<RankTwoSharedInterval> intervals;
};
struct RankTwoWholeArrangementRecords {
    RankTwoArrangement global;
    PreparedRankTwoNormalMesh normals;
    std::vector<RankTwoArrangement> normalFaces; // same order as normals.faces()
    std::vector<RankTwoSharedEdgeRecord> sharedEdges; // canonical edge order
    std::vector<RankTwoNormalVertexRecord> normalVertices; // input normal index order
    std::size_t totalArrangementCells=0;
};

// Exact stored-grid overload. No re-normalization or inferred radius. The
// settings-only adapter must use the same original fixed normal generator.
RankTwoArrangement analyzeRankTwoArrangement(const PreparedRankTwoSource&,
    const CGAL::Gmpq& radius,const PreparedRankTwoNormalMesh&,
    const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
const RankTwoWholeArrangementRecords& wholeArrangementRecords(const RankTwoArrangement&);
RankTwoArrangement reloadRankTwoArrangement(const PreparedRankTwoSource&,
    const CGAL::Gmpq& radius,const PreparedRankTwoNormalMesh&,const std::string& saved,
    const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
} // namespace spin
