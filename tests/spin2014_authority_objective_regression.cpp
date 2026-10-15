#include "spin2014/QuadricAuthorityMomentProblem.h"
#include <gmpxx.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
static mpq_class rational(const std::string& s) {
    mpq_class q(s); q.canonicalize(); return q;
}
static std::string text(const mpq_class& q) {
    return q.get_num().get_str()+"/"+q.get_den().get_str();
}
static double relative(const mpq_class& actual,const mpq_class& expected) {
    const mpq_class difference=actual-expected;
    return std::abs(difference.get_d())/std::max(1.0,std::abs(expected.get_d()));
}
static void near(const std::string& actual, const mpq_class& expected) {
    mpq_class error = rational(actual) - expected;
    if (error < 0) error = -error;
    mpz_class denominator; mpz_ui_pow_ui(denominator.get_mpz_t(), 10, 290);
    require(error < mpq_class(1, denominator), "MPFR1024 objective disagrees with hand rational value");
}

int main() {
    try {
        using namespace spin2014;
        const std::array<std::string, 10> moments{
            "2/1", "0/1", "0/1", "1/1", "3/1", "4/1", "5/1", "0/1", "0/1", "0/1"};
        const std::array<Spin2014Objective, 5> modes{
            Spin2014Objective::YoYo, Spin2014Objective::TopInertiaOnly,
            Spin2014Objective::TopLoweringOnly, Spin2014Objective::TopMassAndLowering,
            Spin2014Objective::StaticLowering};
        const std::array<mpq_class, 5> expected{
            mpq_class(257,98), mpq_class(257,98), mpq_class(563,196),
            mpq_class(355,98), mpq_class(1,4)};
        const std::array<std::size_t, 5> count{5,4,4,4,2};
        for (std::size_t i=0; i<modes.size(); ++i) {
            QuadricObjectiveSettings settings; settings.objective=modes[i];
            const auto value=evaluateSpin2014AuthorityMomentProblem(settings,moments);
            require(value.precisionBits==1024 && value.rawMomentsAuthority==moments,
                    "objective must bind actual high-precision moments");
            near(value.objective,expected[i]);
            require(value.constraints.size()==count[i], "wrong mode constraints");
            near(value.constraints[0],0); near(value.constraints[1],0);
            if (i==0) near(value.constraints[2],1);
            if (i==4) {
                near(value.objectiveGradient[0],mpq_class(-1,4));
                near(value.objectiveGradient[3],mpq_class(1,2));
                near(value.objectiveHessian[0][0],mpq_class(3,8));
                near(value.objectiveHessian[0][3],mpq_class(-1,2));
                near(value.objectiveHessian[3][3],mpq_class(1,2));
            }
        }
        QuadricObjectiveSettings staticSettings;
        staticSettings.objective=Spin2014Objective::StaticLowering;
        auto tiny=moments;
        tiny[3]="10000000000000000000000001/10000000000000000000000000";
        const auto a=evaluateSpin2014AuthorityMomentProblem(staticSettings,moments);
        const auto b=evaluateSpin2014AuthorityMomentProblem(staticSettings,tiny);
        require(rational(b.objective)>rational(a.objective) &&
                rational(b.objectiveGradient[3])>rational(a.objectiveGradient[3]),
                "sub-double moment update must not disappear in objective or gradient");
        staticSettings.solverMassUnitInDeclaredUnits=3;
        staticSettings.solverLengthUnitInDeclaredUnits=2;
        const auto units=evaluateSpin2014AuthorityMomentProblem(staticSettings,moments);
        near(units.objective,1); near(units.objectiveGradient[3],2);
        staticSettings.objective=Spin2014Objective::TopMassAndLowering;
        const auto massUnits=evaluateSpin2014AuthorityMomentProblem(staticSettings,moments);
        near(massUnits.objective,mpq_class(257,98)+36);

        auto asymmetric=moments;
        asymmetric[1]="1/10"; asymmetric[2]="1/5";
        asymmetric[7]="1/4"; asymmetric[8]="1/5"; asymmetric[9]="1/6";
        double maximumDerivativeError=0;
        for (const auto mode:modes) {
            QuadricObjectiveSettings settings; settings.objective=mode;
            const auto base=evaluateSpin2014AuthorityMomentProblem(settings,asymmetric);
            if (mode!=Spin2014Objective::StaticLowering) {
                const std::size_t firstAxis=mode==Spin2014Objective::YoYo?3:2;
                near(base.constraints[firstAxis],mpq_class(-3,20));
                near(base.constraints[firstAxis+1],mpq_class(-3,20));
            }
            for (unsigned direction=0;direction<10;++direction) {
                for (unsigned scale=0;scale<3;++scale) {
                    const mpq_class step(1,10000UL*(1UL<<scale));
                    auto plus=asymmetric,minus=asymmetric;
                    plus[direction]=text(rational(plus[direction])+step);
                    minus[direction]=text(rational(minus[direction])-step);
                    const auto p=evaluateSpin2014AuthorityMomentProblem(settings,plus);
                    const auto m=evaluateSpin2014AuthorityMomentProblem(settings,minus);
                    const mpq_class gradientFD=(rational(p.objective)-rational(m.objective))/(2*step);
                    maximumDerivativeError=std::max(maximumDerivativeError,
                        relative(gradientFD,rational(base.objectiveGradient[direction])));
                    for (unsigned row=0;row<10;++row) {
                        const mpq_class hessianFD=(rational(p.objectiveGradient[row])-
                            rational(m.objectiveGradient[row]))/(2*step);
                        maximumDerivativeError=std::max(maximumDerivativeError,
                            relative(hessianFD,rational(base.objectiveHessian[row][direction])));
                    }
                    for (std::size_t row=0;row<base.constraints.size();++row) {
                        const mpq_class constraintFD=(rational(p.constraints[row])-
                            rational(m.constraints[row]))/(2*step);
                        maximumDerivativeError=std::max(maximumDerivativeError,
                            relative(constraintFD,rational(base.constraintJacobian[row][direction])));
                        for (unsigned col=0;col<10;++col) {
                            const mpq_class secondFD=(rational(p.constraintJacobian[row][col])-
                                rational(m.constraintJacobian[row][col]))/(2*step);
                            maximumDerivativeError=std::max(maximumDerivativeError,
                                relative(secondFD,rational(base.constraintHessians[row][col][direction])));
                        }
                    }
                }
            }
        }
        require(maximumDerivativeError<=1e-6,"three-scale objective/constraint derivatives exceed tolerance");
        auto invalid=moments; invalid[0]="0/1";
        bool rejected=false;
        try { (void)evaluateSpin2014AuthorityMomentProblem(staticSettings,invalid); }
        catch (const std::exception&) { rejected=true; }
        require(rejected,"nonpositive mass must be rejected");
        invalid=moments; invalid[0]="4/2";
        rejected=false;
        try { (void)evaluateSpin2014AuthorityMomentProblem(staticSettings,invalid); }
        catch (const std::exception&) { rejected=true; }
        require(rejected,"noncanonical authority must be rejected");
        std::cout << "asymmetric ten-direction three-scale maximum derivative error=" << maximumDerivativeError << '\n';
        std::cout << "five authority objectives: hand values, derivatives, units and sub-double updates passed\n";
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
