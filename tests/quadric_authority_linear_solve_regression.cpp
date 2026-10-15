#include "spin/QuadricAuthorityLinearSolve.h"
#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    try {
        // Catch a double residual falsely becoming zero after cancellation.
        const auto cancellation = spin::solveSpinNewtonAuthorityLinearSystem(
            {{"3/1"}}, {"100000000000000000000/1"}, "1/1000000000000");
        require(cancellation.candidate.solved, "protected double candidate should exist");
        require(cancellation.candidate.residualInfinityNorm == 0.0,
                "fixture must expose double residual cancellation");
        require(!cancellation.solved && cancellation.residualInfinityUpper == "4096/1",
                "authority residual must reject the falsely exact double solve");

        const auto fraction = spin::solveSpinNewtonAuthorityLinearSystem(
            {{"3/1"}}, {"1/1"}, "1/1000000000000");
        require(fraction.solved, "accurate protected solve should survive authority check");
        require(fraction.residualInfinityUpper == "1/18014398509481984",
                "residual must use exact represented solution, not rounded A*x");
        require(fraction.precisionBits == 1024 &&
                fraction.matrixAuthority[0][0] == "3/1", "authority identity lost");

        const auto tiny = spin::solveSpinNewtonAuthorityLinearSystem(
            {{"10000000000000000000000001/10000000000000000000000000"}},
            {"1/1"}, "1/100000000000000000000000000");
        require(tiny.candidate.solved && tiny.candidate.residualInfinityNorm == 0.0 &&
                !tiny.solved, "sub-double authority coefficient must survive residual validation");
        const auto equilibrated = spin::solveSpinNewtonAuthorityLinearSystem(
            {{"1000000000000/1", "0/1"}, {"0/1", "1/1000000"}},
            {"2000000000000/1", "-3/1000000"}, "1/1000000000000");
        require(equilibrated.solved && equilibrated.candidate.usedEquilibration &&
                equilibrated.residualInfinityUpper == "0/1",
                "authority check must retain reversible protected equilibration");

        const auto singular = spin::solveSpinNewtonAuthorityLinearSystem(
            {{"1/1", "1/1"}, {"1/1", "1/1"}}, {"2/1", "2/1"}, "1/1000000000000");
        require(!singular.solved && !singular.candidate.solved,
                "singular protected solve must not be replaced by a pseudoinverse");
        const auto malformed = spin::solveSpinNewtonAuthorityLinearSystem(
            {{"NaN"}}, {"1/1"}, "1/1000000000000");
        require(!malformed.solved && malformed.candidate.solution.size() == 0,
                "invalid authority must be rejected before solving");
        const auto noncanonical = spin::solveSpinNewtonAuthorityLinearSystem(
            {{"6/2"}}, {"1/1"}, "1/1000000000000");
        require(!noncanonical.solved, "noncanonical authority must not be silently accepted");
        const auto invalidLimit = spin::solveSpinNewtonAuthorityLinearSystem(
            {{"3/1"}}, {"1/1"}, "0/1");
        require(!invalidLimit.solved, "zero residual limit must be rejected");
        const auto dimensions = spin::solveSpinNewtonAuthorityLinearSystem(
            {{"1/1", "0/1"}}, {"1/1"}, "1/1000000000000");
        require(!dimensions.solved, "nonsquare authority must be rejected");
        std::cout << "authority linear solve: cancellation rejected; ordinary solve accepted; protected failures retained\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
