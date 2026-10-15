#include "spin2014/QuadricAuthorityNewtonState.h"
#include "spin/QuadricAuthorityLinearSolve.h"
#include <gmpxx.h>
#include <algorithm>
#include <iostream>
#include <locale>
#include <stdexcept>

static void require(bool c,const char* message) { if(!c) throw std::runtime_error(message); }
static mpq_class read(const std::string& s) { mpq_class q(s); q.canonicalize(); return q; }
static std::string text(const mpq_class& q) {return q.get_num().get_str()+"/"+q.get_den().get_str();}
static void near(const std::string& actual,const mpq_class& expected) {
    mpq_class e=read(actual)-expected; if(e<0)e=-e;
    mpz_class denominator;mpz_ui_pow_ui(denominator.get_mpz_t(),10,290);
    require(e<mpq_class(1,denominator),"authority Eq.14 disagrees with hand-derived system");
}
struct GroupedNumbers:std::numpunct<char> {
    char do_thousands_sep()const override {return '_';}
    std::string do_grouping()const override {return "\1";}
};
int main() {
 try {
    using namespace spin2014;
    using A=spin::QuadricCoefficientAuthority10;
    A::RationalVector q;q.fill("0/1");q[0]="1/1";
    const auto authority=A::fromNormalizedGeometry(q);
    AuthorityMomentVector10 moments{"2/1","0/1","0/1","1/1","3/1","4/1","5/1","0/1","0/1","0/1"};
    AuthorityMomentMatrix10 jacobian;
    for(unsigned i=0;i<10;++i)for(unsigned j=0;j<10;++j)jacobian[i][j]=i==j?"1/1":"0/1";
    QuadricObjectiveSettings settings;settings.objective=Spin2014Objective::StaticLowering;
    const auto state=linearizeSpin2014AuthorityNewtonState(settings,authority,moments,jacobian,"algebraic-identity-map");
    require(state.precisionBits==1024 && state.matrix.size()==13,"state layout wrong");
    near(state.multiplierMu,mpq_class(1,4));
    near(state.multipliers[0],0);near(state.multipliers[1],0);
    near(state.rawEquationResidual,mpq_class(1,2));
    for(unsigned row=0;row<13;++row) {
        near(state.rightHandSide[row],row==3?mpq_class(-1,2):mpq_class(0));
        for(unsigned col=0;col<13;++col) {
            mpq_class expected=0;
            if(row<10 && row==col)expected=mpq_class(1,4);
            if(row==0 && col==3)expected=mpq_class(-1,2);
            if(row==3 && col==3)expected=mpq_class(3,4);
            if((row==1 && col==10)||(row==2 && col==11))expected=-1;
            if((row==10 && col==1)||(row==11 && col==2)||(row==0 && col==12))expected=1;
            if(row==12 && col==0)expected=2;
            near(state.matrix[row][col],expected);
        }
    }
    const auto direction=spin::solveSpinNewtonAuthorityLinearSystem(state.matrix,state.rightHandSide,"1/1000000000000");
    require(direction.solved,"actual protected linear kernel rejected hand Eq.14");
    require(std::abs(direction.candidate.solution(3)+2.0/3.0)<1e-15,"wrong restored direction");
    auto changed=moments;changed[3]="10000000000000000000000001/10000000000000000000000000";
    const auto tiny=linearizeSpin2014AuthorityNewtonState(settings,authority,changed,jacobian,"algebraic-identity-map");
    require(tiny.multiplierMu!=state.multiplierMu && tiny.rawEquationResidual!=state.rawEquationResidual &&
            tiny.scientificHash!=state.scientificHash,"sub-double moments lost in multipliers/state identity");
    const auto originalLocale=std::locale();
    std::locale::global(std::locale(originalLocale,new GroupedNumbers));
    const auto localized=linearizeSpin2014AuthorityNewtonState(settings,authority,moments,jacobian,"algebraic-identity-map");
    std::locale::global(originalLocale);
    require(localized.scientificHash==state.scientificHash,"scientific state identity must not depend on host numeric locale");
    auto radial=q;radial[0]="0/1";radial[1]="1/1";
    bool rejected=false;
    try {(void)linearizeSpin2014AuthorityNewtonState(settings,A::fromNormalizedGeometry(radial),moments,jacobian,"rank-deficient-map");}
    catch(const std::exception&) {rejected=true;}
    require(rejected,"rank-deficient multiplier system must not be regularized");
    rejected=false;
    try {(void)linearizeSpin2014AuthorityNewtonState(settings,authority,moments,jacobian,"");}
    catch(const std::exception&) {rejected=true;}
    require(rejected,"missing material mapping identity accepted");

    A::RationalVector nonDiagonal{"-2/1","1/10","-1/5","3/10","3/1","4/1","5/1","1/5","-1/10","3/20"};
    const auto curved=A::fromUnnormalizedGeometry(nonDiagonal);
    double maximumError=0;
    for(const auto mode:{Spin2014Objective::YoYo,Spin2014Objective::TopInertiaOnly,
            Spin2014Objective::TopLoweringOnly,Spin2014Objective::TopMassAndLowering,
            Spin2014Objective::StaticLowering}) {
        settings.objective=mode;
        const auto base=linearizeSpin2014AuthorityNewtonState(settings,curved,moments,jacobian,"affine-moment-test-map");
        for(unsigned axis=0;axis<10;++axis)for(unsigned scale=0;scale<3;++scale) {
            A::RationalVector positive,negative;positive.fill("0/1");negative.fill("0/1");
            positive[axis]="1/1";negative[axis]="-1/1";
            const mpq_class step(1,10000UL*(1UL<<scale));
            const auto p=curved.geodesicUpdate(positive,text(step));
            const auto n=curved.geodesicUpdate(negative,text(step));
            auto plusMoments=moments,minusMoments=moments;
            for(unsigned k=0;k<10;++k) {
                plusMoments[k]=text(read(moments[k])+read(p.geometry()[k])-read(curved.geometry()[k]));
                minusMoments[k]=text(read(moments[k])+read(n.geometry()[k])-read(curved.geometry()[k]));
            }
            const auto plus=evaluateSpin2014AuthorityMomentProblem(settings,plusMoments);
            const auto minus=evaluateSpin2014AuthorityMomentProblem(settings,minusMoments);
            const auto count=base.multipliers.size();
            for(unsigned row=0;row<10+count;++row) {
                mpq_class fp,fm;
                if(row<10) {
                    fp=read(plus.objectiveGradient[row])+read(base.multiplierMu)*read(p.dual()[row]);
                    fm=read(minus.objectiveGradient[row])+read(base.multiplierMu)*read(n.dual()[row]);
                    for(std::size_t c=0;c<count;++c) {
                        fp-=read(base.multipliers[c])*read(plus.constraintJacobian[c][row]);
                        fm-=read(base.multipliers[c])*read(minus.constraintJacobian[c][row]);
                    }
                } else {fp=read(plus.constraints[row-10]);fm=read(minus.constraints[row-10]);}
                mpq_class predicted=0;
                for(unsigned k=0;k<10;++k)predicted+=read(base.matrix[row][k])*
                    (read(p.dual()[k])-read(n.dual()[k]))/(2*step);
                const mpq_class actual=(fp-fm)/(2*step),difference=actual-predicted;
                const double error=std::abs(difference.get_d())/std::max(1.0,std::abs(predicted.get_d()));
                maximumError=std::max(maximumError,error);
            }
        }
    }
    require(maximumError<=1e-6,"five-mode Eq.14 dual-coordinate derivatives exceed tolerance");
    std::cout<<"five-mode ten-direction three-scale assembly difference="<<maximumError<<'\n';
    std::cout<<"MPFR1024 fixed-state Eq.14: hand matrix, multipliers, residual and protected direction passed\n";
    std::cout<<"state hash="<<state.scientificHash<<'\n';
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
