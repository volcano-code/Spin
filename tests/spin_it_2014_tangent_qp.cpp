#include "spin2014/SlqpTangentQp.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double actual, double expected, const char* message) {
    require(std::isfinite(actual) && std::abs(actual - expected) < 1e-9, message);
}
spin2014::LbfgsHistory metric() {
    spin2014::LbfgsHistory h;
    h.capacity = 2;
    require(h.update(Eigen::Vector2d(1., 2.), Eigen::Vector2d(1., 0.)), "pair rejected");
    return h;
}
void checked(const spin2014::TangentQpResult& r) {
    require(r.attempted && r.checked && r.status == "QP_FIRST_ORDER_CHECKED", "QP not checked");
    require(r.direction.allFinite() && r.multipliers.allFinite(), "nonfinite solution");
    require(std::isfinite(r.equalityResidual) && r.equalityResidual <= 1e-9, "equality error");
    require(std::isfinite(r.stationarity) && r.stationarity <= 1e-9, "stationarity error");
    require(!r.trace.empty(), "missing trace");
}
void unknown(const spin2014::TangentQpResult& r) {
    require(r.attempted && !r.checked && r.status == "QP_UNKNOWN", "failure incorrectly checked");
}
}
int main() {
    try {
        const Eigen::MatrixXd empty(0, 2);
        const Eigen::VectorXd zero(0);
        const auto h = metric(); // Inverse metric [[1,2],[2,9]], Hessian [[1.8,-.4],[-.4,.2]].
        auto lower = spin2014::solveTangentQp(Eigen::Vector2d(0., .5), Eigen::Vector2d(0., 1.), empty, zero, 1, h);
        checked(lower);
        require(lower.direction[0] == 0., "fixed lower face not exactly zero");
        near(lower.direction[1], -5., "lower bound metric solution");
        require(lower.workingSet == std::vector<int>{0}, "missing lower face");
        near(lower.multipliers[0], -2., "lower multiplier sign");
        auto upper = spin2014::solveTangentQp(Eigen::Vector2d(1., .5), Eigen::Vector2d(0., -1.), empty, zero, 1, h);
        checked(upper);
        require(upper.direction[0] == 0., "fixed upper face not exactly zero");
        near(upper.direction[1], 5., "upper bound metric solution");
        near(upper.multipliers[0], 2., "upper multiplier sign");

        auto released = spin2014::solveTangentQp(Eigen::Vector2d(0., .5), Eigen::Vector2d(-1., 1.), empty, zero, 1, {});
        checked(released);
        near(released.direction[0], 1., "wrong-sign face not released");
        require(released.workingSet.empty() && released.trace.size() >= 2, "release trace missing");

        Eigen::MatrixXd a(1, 2); a << 1., 1.;
        auto equality = spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d(0., 2.), a, Eigen::VectorXd::Zero(1), 2, h);
        checked(equality);
        near(equality.direction[0], 5./7., "equality tangent x");
        near(equality.direction[1], -5./7., "equality tangent y");
        near(equality.multipliers[0], -11./7., "equality multiplier");

        a << 1., 0.;
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(0., .5), Eigen::Vector2d(0., 1.), a, Eigen::VectorXd::Ones(1), 1, h));
        const auto rankFailure = spin2014::solveTangentQp(Eigen::Vector2d(0., .5), Eigen::Vector2d(0., 1.), a, Eigen::VectorXd::Zero(1), 1, h);
        unknown(rankFailure);
        require(rankFailure.direction.size() == 0 && rankFailure.multipliers.size() == 0, "unsolved rank failure exposes solve data");
        // Free coordinates at numerical box endpoints are not box faces.
        auto unbounded = spin2014::solveTangentQp(Eigen::Vector2d(.5, 0.), Eigen::Vector2d(0., 1.), empty, zero, 0, {});
        checked(unbounded); near(unbounded.direction[1], -1., "unbounded coordinate became face");

        auto repeat = spin2014::solveTangentQp(Eigen::Vector2d(0., .5), Eigen::Vector2d(-1., 1.), empty, zero, 1, {});
        require(repeat.trace.size() == released.trace.size(), "nondeterministic trace length");
        for (std::size_t i=0; i<repeat.trace.size(); ++i) {
            require(repeat.trace[i].workingSet == released.trace[i].workingSet && repeat.trace[i].action == released.trace[i].action && repeat.trace[i].changedIndex == released.trace[i].changedIndex, "nondeterministic trace");
        }
        // Coupling makes a previously released face outward: it must be added
        // back, and its tiny numerical direction represented as exactly zero.
        spin2014::LbfgsHistory coupled; coupled.capacity = 2;
        require(coupled.update(Eigen::Vector3d(3., -3., 2.), Eigen::Vector3d(-1., -1., 1.)), "coupled pair rejected");
        const auto added = spin2014::solveTangentQp(Eigen::Vector3d::Zero(), Eigen::Vector3d(-3., -4., -1.), Eigen::MatrixXd(0, 3), zero, 3, coupled);
        checked(added);
        near(added.direction[0], 268./69., "re-added face solution x");
        require(added.direction[1] == 0., "roundoff on fixed face");
        near(added.direction[2], 181./69., "re-added face solution z");
        require(added.workingSet == std::vector<int>{1}, "wrong re-added face");
        bool sawAdd = false;
        for (const auto& step : added.trace)
            sawAdd = sawAdd || (step.action == "ADD_FACE" && step.changedIndex == 1);
        require(sawAdd, "outward inactive face was not added");
        // Retaining more than one pair exercises the complete BFGS recursion;
        // the final secant identity gives B*[0,1]=[1,2] independently.
        auto multiple = h;
        require(multiple.update(Eigen::Vector2d(0., 1.), Eigen::Vector2d(1., 2.)), "second pair rejected");
        const auto twoPairs = spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d(-1., -2.), empty, zero, 0, multiple);
        checked(twoPairs);
        near(twoPairs.direction[0], 0., "multiple-pair secant x");
        near(twoPairs.direction[1], 1., "multiple-pair secant y");

        spin2014::TangentQpSettings settings;
        settings.maximumFaces = 0;
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(0., .5), Eigen::Vector2d::Ones(), empty, zero, 1, h, settings));
        settings = {}; settings.maximumIterations = 1;
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(0., .5), Eigen::Vector2d(-1., 1.), empty, zero, 1, {}, settings));
        settings.maximumIterations = 0;
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 0, {}, settings));
        settings = {}; settings.residualTolerance = std::numeric_limits<double>::quiet_NaN();
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 0, {}, settings));
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(-.1, .5), Eigen::Vector2d::Ones(), empty, zero, 1, h));
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 3, h));
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::VectorXd::Ones(3), empty, zero, 1, h));
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d(std::numeric_limits<double>::infinity(), 0.), empty, zero, 1, h));
        auto invalid = h; invalid.inverseCurvature[0] = 2.;
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 1, invalid));
        invalid = h; invalid.y.clear();
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 1, invalid));
        invalid = h; invalid.s[0] = Eigen::VectorXd::Ones(3);
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 1, invalid));
        invalid = h; invalid.capacity = 0;
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 1, invalid));
        invalid = h; invalid.y[0] *= -1.;
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 1, invalid));
        invalid = h; invalid.s[0][0] = std::numeric_limits<double>::quiet_NaN();
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 1, invalid));
        for (double tolerance : {-1., .5, std::numeric_limits<double>::infinity()}) {
            settings = {}; settings.activeTolerance = tolerance;
            unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 1, h, settings));
        }
        settings = {}; settings.maximumFaces = -1;
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 1, h, settings));
        settings = {}; settings.maximumIterations = -1;
        unknown(spin2014::solveTangentQp(Eigen::Vector2d(.5, .5), Eigen::Vector2d::Ones(), empty, zero, 1, h, settings));
        // A strict tolerance must reject unverified floating-point cancellation.
        settings = {}; settings.residualTolerance = 1e-30;
        unknown(spin2014::solveTangentQp(Eigen::Vector3d::Zero(), Eigen::Vector3d(-3., -4., -1.), Eigen::MatrixXd(0, 3), zero, 3, coupled, settings));
        std::cout << "spin2014 tangent QP tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
