#pragma once

#include <CGAL/Gmpq.h>

#include <vector>
#include <map>

namespace spin::detail {

using ExactRationalPolynomial = std::vector<CGAL::Gmpq>;

/** Sum sign(value(root)) over the distinct real roots of polynomial in (lo,hi).
 *
 * Preconditions: polynomial is square-free and neither interval endpoint is
 * a root. Coefficients are in ascending order. The result is exact.
 */
int sturmTarskiSignSum(const ExactRationalPolynomial& polynomial,
                       const ExactRationalPolynomial& value,
                       const CGAL::Gmpq& lo,const CGAL::Gmpq& hi);

/** Same signed remainder chains, reusable only for one immutable polynomial.
 * Each query still evaluates the exact signs at its own nonroot endpoints;
 * no interval answer, source geometry or global cache is reused. */
class PreparedSturmTarskiCounter {
public:
    explicit PreparedSturmTarskiCounter(ExactRationalPolynomial polynomial);
    int signSum(const ExactRationalPolynomial& value,
                const CGAL::Gmpq& lo,const CGAL::Gmpq& hi);
private:
    ExactRationalPolynomial polynomial_;
    std::map<ExactRationalPolynomial,std::vector<ExactRationalPolynomial>> chains_;
};

} // namespace spin::detail
