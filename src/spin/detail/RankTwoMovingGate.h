#pragma once
#include "RankTwoArrangementRecords.h"
#include "RankTwoCarrierKernels.h"
#include <memory>

namespace spin {
// One shared moving projection, tied to two ACTUAL immutable event IDs.
// Private construction prevents substituting an unrelated algebraic value.
class RankTwoMovingGate {
public:
    static RankTwoMovingGate fromVertices(const RankTwoArrangement&,
        const std::string& endpointIdentity,const std::string& anchorIdentity,
        int projectionAxis,const CGAL::Gmpq& fraction,int maximumIsolationBits=512);
    // A ruler on the opposite side of an actual extreme event:
    // coordinate = endpoint + fraction*(endpoint-anchor), 0<fraction<1.
    // The stored fraction is NEGATIVE and carries a distinct gate identity.
    // This does not certify native membership; the endpoint arc proof does.
    static RankTwoMovingGate fromReflectedVertices(const RankTwoArrangement&,
        const std::string& endpointIdentity,const std::string& anchorIdentity,
        int projectionAxis,const CGAL::Gmpq& fraction,int maximumIsolationBits=512);
    const std::string& arrangementIdentity() const;
    const std::string& endpointIdentity() const;
    const std::string& anchorIdentity() const;
    int projectionAxis() const;
    const CGAL::Gmpq& fraction() const;
    const detail::CertifiedAlgebraicReal& coordinate() const;
    const std::string& identity() const;
private:
    static RankTwoMovingGate fromVertexWeight(const RankTwoArrangement&,
        const std::string&,const std::string&,int,const CGAL::Gmpq&,int);
    struct Impl;
    explicit RankTwoMovingGate(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
};
} // namespace spin
