#pragma once
#include "RankTwoNormalOverlay.h"
#include <vector>

namespace spin {
// Admission of the stored rays as ONE embedded, oriented spherical
// triangulation. Positive ray scaling is allowed but never applied to the
// authority. This is not a source-surface or floating mesh certificate.
class PreparedRankTwoNormalMesh {
public:
    using Rays=std::vector<std::array<CGAL::Gmpq,3>>;
    using Faces=std::vector<std::array<std::size_t,3>>;
    static PreparedRankTwoNormalMesh fromExact(const Rays&,const Faces&,
        const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
    static PreparedRankTwoNormalMesh reload(const Rays&,const Faces&,const std::string& saved,
        const RankTwoSourceLimits& limits=RankTwoSourceLimits{});
    const Rays& rays() const;
    const Faces& faces() const;
    const std::vector<std::array<std::size_t,2>>& edges() const;
    RankTwoNormalFace face(std::size_t) const;
    std::size_t checkedFacePairs() const;
    const std::string& inputHash() const;
    const std::string& canonicalRecord() const;
private:
    struct Impl;
    explicit PreparedRankTwoNormalMesh(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
};
} // namespace spin
