#include "spin/QuadricAuthorityLinearSolve.h"
#include <gmpxx.h>
#include <mpfr.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace spin {
namespace {
std::string text(const mpq_class& q) {
    return q.get_num().get_str() + "/" + q.get_den().get_str();
}
mpq_class read(const std::string& s) {
    mpq_class q;
    if (s.empty() || mpq_set_str(q.get_mpq_t(), s.c_str(), 10) != 0 ||
        q.get_den() == 0) throw std::invalid_argument("nonfinite rational equation input");
    q.canonicalize();
    if (text(q) != s) throw std::invalid_argument("equation input must be a canonical rational");
    return q;
}
mpq_class magnitude(const mpq_class& q) { return q < 0 ? mpq_class(-q) : q; }
std::string outward(const mpq_class& q, mpfr_rnd_t rounding) {
    mpfr_t value;
    mpfr_init2(value, 1024);
    mpfr_set_q(value, q.get_mpq_t(), rounding);
    if (!mpfr_number_p(value)) {
        mpfr_clear(value);
        throw std::invalid_argument("authority residual exceeds MPFR range");
    }
    mpq_class bound;
    mpfr_get_q(bound.get_mpq_t(), value);
    mpfr_clear(value);
    return text(bound);
}
} // namespace

SpinNewtonAuthorityLinearSolveResult solveSpinNewtonAuthorityLinearSystem(
    const std::vector<std::vector<std::string>>& matrix,
    const std::vector<std::string>& rightHandSide,
    const std::string& absoluteResidualLimit) {
    SpinNewtonAuthorityLinearSolveResult result;
    try {
        const auto size = matrix.size();
        if (size == 0 || rightHandSide.size() != size)
            throw std::invalid_argument("authority equation must be nonempty and square");
        const mpq_class limit = read(absoluteResidualLimit);
        if (limit <= 0) throw std::invalid_argument("absolute residual limit must be positive");
        std::vector<std::vector<mpq_class>> a(size, std::vector<mpq_class>(size));
        std::vector<mpq_class> b(size);
        Eigen::MatrixXd display(size, size);
        Eigen::VectorXd rhs(size);
        mpq_class aNorm = 0, bNorm = 0;
        for (std::size_t row = 0; row < size; ++row) {
            if (matrix[row].size() != size)
                throw std::invalid_argument("authority equation must be nonempty and square");
            mpq_class rowNorm = 0;
            for (std::size_t col = 0; col < size; ++col) {
                a[row][col] = read(matrix[row][col]);
                display(row, col) = a[row][col].get_d();
                rowNorm += magnitude(a[row][col]);
            }
            aNorm = std::max(aNorm, rowNorm);
            b[row] = read(rightHandSide[row]);
            rhs(row) = b[row].get_d();
            bNorm = std::max(bNorm, magnitude(b[row]));
        }
        if (!display.allFinite() || !rhs.allFinite())
            throw std::invalid_argument("authority equation has no finite double compatibility view");
        result.matrixAuthority = matrix;
        result.rightHandSideAuthority = rightHandSide;
        result.absoluteResidualLimit = absoluteResidualLimit;
        result.candidate = solveSpinNewtonLinearSystem(display, rhs);
        if (!result.candidate.solved) {
            result.diagnostic = "protected candidate rejected: " + result.candidate.diagnostic;
            return result;
        }
        // Capture each represented double exactly; never use a rounded double
        // A*x or a compatibility matrix as the authority for the check.
        std::vector<mpq_class> x(size);
        mpq_class xNorm = 0, residual = 0;
        for (std::size_t col = 0; col < size; ++col) {
            x[col] = mpq_class(result.candidate.solution(col));
            xNorm = std::max(xNorm, magnitude(x[col]));
        }
        for (std::size_t row = 0; row < size; ++row) {
            mpq_class value = -b[row];
            for (std::size_t col = 0; col < size; ++col) value += a[row][col] * x[col];
            residual = std::max(residual, magnitude(value));
        }
        const mpq_class denominator = aNorm * xNorm + bNorm;
        if (denominator <= 0) throw std::invalid_argument("zero backward-error denominator");
        const mpq_class backwardError = residual / denominator;
        result.residualInfinityLower = outward(residual, MPFR_RNDD);
        result.residualInfinityUpper = outward(residual, MPFR_RNDU);
        result.backwardErrorUpper = outward(backwardError, MPFR_RNDU);
        const mpq_class protectedLimit(128.0 * std::numeric_limits<double>::epsilon() *
                                      static_cast<double>(size));
        result.solved = read(result.residualInfinityUpper) <= limit &&
                        read(result.backwardErrorUpper) <= protectedLimit;
        result.diagnostic = result.solved
            ? "protected candidate passed authority equation and MPFR1024 outward residual bounds"
            : "protected candidate failed authority equation residual bounds";
    } catch (const std::exception& error) {
        result.solved = false;
        result.diagnostic = error.what();
    }
    return result;
}
} // namespace spin
