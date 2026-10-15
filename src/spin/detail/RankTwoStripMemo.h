#pragma once
#include "RankTwoArcStrip.h"
#include <memory>

namespace spin::detail {
// Pure rational strip calculations only. Private live entries bind the full
// polynomial, ORDERED obstacles, exact box and accuracy. No caller-supplied
// witness, persistence, global/thread cache or registry/endpoint proof.
// A refiner owns one instance; every actual edge still runs its fresh guards.
class RankTwoRegularStripMemo {
public:
    RankTwoRegularStripMemo();
    ~RankTwoRegularStripMemo();
    RankTwoRegularStripMemo(const RankTwoRegularStripMemo&)=delete;
    RankTwoRegularStripMemo& operator=(const RankTwoRegularStripMemo&)=delete;
    // Successful pure parsing only; bind raw factor text AND chart. Returned
    // references remain valid for this memo's lifetime. No owner/edge proof.
    const RankTwoStripPolynomial& parsedChart(const std::string&,int projectionAxis);
    // Complete raw source + chart + ORDERED raw obstacles index the SAME
    // owned pure rational context/query store, never a registry proof. Names
    // and coefficients are owned by this memo; caller strings may change.
    RankTwoRegularStripCertificate certifyChart(const std::string&,int projectionAxis,
        const RankTwoStripBox&,const std::vector<std::string>& orderedObstacles,
        const CGAL::Gmpq& maximumChordDistance);
    // Uses the very same owned context/query store as the legacy APIs.
    RankTwoRegularStripAttempt tryCertifyChart(const std::string&,int projectionAxis,
        const RankTwoStripBox&,const std::vector<std::string>& orderedObstacles,
        const CGAL::Gmpq& maximumChordDistance);
    RankTwoRegularStripCertificate certify(const RankTwoStripPolynomial&,
        const RankTwoStripBox&,const std::vector<RankTwoStripPolynomial>&,
        const CGAL::Gmpq& maximumChordDistance);
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
