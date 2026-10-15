#pragma once
#include "spin/RankTwoSource.h"
#include <array>
#include <memory>
#include <vector>

namespace spin {
namespace detail {class CertifiedAlgebraicReal;}
struct RankTwoNormalFace;
struct RankTwoNormalBoundaryEvent;
struct RankTwoNormalVertexRecord;
struct RankTwoArrangementRecords;
class PreparedRankTwoNormalMesh;
struct RankTwoWholeArrangementRecords;
struct RankTwoPolynomialOwner {
    std::string polynomial;
    int multiplicity=0;
};
struct RankTwoZeroFactor {
    std::string polynomial; // primitive integer polynomial, indexed by x/y powers
    std::vector<RankTwoPolynomialOwner> owners;
};

// Global F/h/D analysis is shared before per-normal-face overlays. This class
// is not a mesh and cannot be passed to Boolean or to the motion evaluator.
class RankTwoArrangement {
public:
    static RankTwoArrangement reload(const PreparedRankTwoSource&,const CGAL::Gmpq& radius,
        const std::string& saved,const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
    std::size_t faceCount() const;
    std::size_t certifiedFaceCount() const;
    std::size_t cutComponents() const;
    std::size_t cutBoundaryLoops() const;
    std::size_t fullComponents() const;
    std::size_t isolatedPoints() const;
    const std::vector<RankTwoZeroFactor>& factors() const;
    const std::string& canonicalRecord() const;
    const std::string& scientificHash() const;
    // Local arrangement authority; whole collections must use their members.
    const std::string& sourceIdentity() const;
    const CGAL::Gmpq& truncationRadius() const;
private:
    struct Impl;
    explicit RankTwoArrangement(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoArrangement analyzeRankTwoGlobalArrangement(
        const PreparedRankTwoSource&,const CGAL::Gmpq&,const RankTwoSourceLimits&);
    friend RankTwoArrangement analyzeRankTwoNormalFace(const PreparedRankTwoSource&,
        const RankTwoArrangement&,const RankTwoNormalFace&,const RankTwoSourceLimits&);
    friend const std::vector<RankTwoNormalBoundaryEvent>& normalBoundaryEvents(const RankTwoArrangement&);
    friend const std::vector<RankTwoNormalVertexRecord>& normalVertexRecords(const RankTwoArrangement&);
    friend const RankTwoArrangementRecords& arrangementRecords(const RankTwoArrangement&);
    friend const std::array<detail::CertifiedAlgebraicReal,2>& certifiedRankTwoVertex(
        const RankTwoArrangement&,const std::string&);
    friend RankTwoArrangement analyzeRankTwoArrangement(const PreparedRankTwoSource&,
        const CGAL::Gmpq&,const PreparedRankTwoNormalMesh&,const RankTwoSourceLimits&);
    friend const RankTwoWholeArrangementRecords& wholeArrangementRecords(const RankTwoArrangement&);
};

RankTwoArrangement analyzeRankTwoGlobalArrangement(
    const PreparedRankTwoSource&,const CGAL::Gmpq& radius,
    const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
} // namespace spin
