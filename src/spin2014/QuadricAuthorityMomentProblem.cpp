#include "spin2014/QuadricAuthorityMomentProblem.h"
#include "spin2014/QuadricAuthorityNewtonState.h"
#include "spin/MomentInertiaAlgebra.h"
#include "manufacturing/SlicerValidation.h"
#include <gmpxx.h>
#include <mpfr.h>
#include <stdexcept>
#include <numeric>
#include <sstream>
#include <locale>
#include <limits>
#include <new>

namespace spin2014 { namespace {
struct Real {
    mpfr_t value;
    Real() { mpfr_init2(value,1024); mpfr_set_zero(value,0); }
    Real(const Real& other):Real() { mpfr_set(value,other.value,MPFR_RNDN); }
    Real& operator=(const Real& other) { mpfr_set(value,other.value,MPFR_RNDN); return *this; }
    ~Real() { mpfr_clear(value); }
    static Real constant(double d) { Real r; mpfr_set_d(r.value,d,MPFR_RNDN); return r; }
};
Real operator+(const Real& a,const Real& b) { Real r; mpfr_add(r.value,a.value,b.value,MPFR_RNDN); return r; }
Real operator-(const Real& a,const Real& b) { Real r; mpfr_sub(r.value,a.value,b.value,MPFR_RNDN); return r; }
Real operator-(const Real& a) { Real r; mpfr_neg(r.value,a.value,MPFR_RNDN); return r; }
Real operator*(const Real& a,const Real& b) { Real r; mpfr_mul(r.value,a.value,b.value,MPFR_RNDN); return r; }
Real operator/(const Real& a,const Real& b) {
    if (mpfr_zero_p(b.value)) throw std::domain_error("zero authority objective denominator");
    Real r; mpfr_div(r.value,a.value,b.value,MPFR_RNDN); return r;
}
std::string rationalText(const mpq_class& q) { return q.get_num().get_str()+"/"+q.get_den().get_str(); }
Real parse(const std::string& text) {
    mpq_class q;
    if (mpq_set_str(q.get_mpq_t(),text.c_str(),10)!=0 || q.get_den()==0)
        throw std::invalid_argument("authority moments must be finite canonical rationals");
    q.canonicalize();
    if (rationalText(q)!=text) throw std::invalid_argument("authority moments must be canonical rationals");
    Real result; mpfr_set_q(result.value,q.get_mpq_t(),MPFR_RNDN);
    if (!mpfr_number_p(result.value)) throw std::domain_error("authority moment outside MPFR range");
    return result;
}
std::string write(const Real& r) {
    if (!mpfr_number_p(r.value)) throw std::domain_error("nonfinite authority objective evaluation");
    mpq_class q; mpfr_get_q(q.get_mpq_t(),r.value); return rationalText(q);
}
struct Jet {
    Real v;
    std::array<Real,10> g;
    std::array<std::array<Real,10>,10> h;
    static Jet constant(const Real& v) { Jet result; result.v=v; return result; }
    static Jet variable(const Real& v,std::size_t index) {
        Jet result=constant(v); result.g[index]=Real::constant(1); return result;
    }
};
Jet operator+(const Jet& a,const Jet& b) {
    Jet r; r.v=a.v+b.v;
    for (unsigned i=0;i<10;++i) { r.g[i]=a.g[i]+b.g[i];
        for (unsigned j=0;j<10;++j) r.h[i][j]=a.h[i][j]+b.h[i][j]; }
    return r;
}
Jet operator-(const Jet& a) {
    Jet r; r.v=-a.v;
    for (unsigned i=0;i<10;++i) { r.g[i]=-a.g[i];
        for (unsigned j=0;j<10;++j) r.h[i][j]=-a.h[i][j]; }
    return r;
}
Jet operator-(const Jet& a,const Jet& b) { return a+(-b); }
Jet operator*(const Jet& a,const Jet& b) {
    Jet r; r.v=a.v*b.v;
    for (unsigned i=0;i<10;++i) { r.g[i]=a.g[i]*b.v+b.g[i]*a.v;
        for (unsigned j=0;j<10;++j) r.h[i][j]=a.h[i][j]*b.v+b.h[i][j]*a.v+
            a.g[i]*b.g[j]+b.g[i]*a.g[j]; }
    return r;
}
Jet operator*(double scale,const Jet& a) { return Jet::constant(Real::constant(scale))*a; }
Jet reciprocal(const Jet& a) {
    const Real inverse=Real::constant(1)/a.v;
    Jet r; r.v=inverse;
    for (unsigned i=0;i<10;++i) { r.g[i]=-(inverse*inverse)*a.g[i];
        for (unsigned j=0;j<10;++j) r.h[i][j]=Real::constant(2)*inverse*inverse*inverse*
            (a.g[i]*a.g[j])-inverse*inverse*a.h[i][j]; }
    return r;
}
Jet operator/(const Jet& a,const Jet& b) { return a*reciprocal(b); }
AuthorityMomentVector10 gradient(const Jet& value) {
    AuthorityMomentVector10 r; for (unsigned i=0;i<10;++i) r[i]=write(value.g[i]); return r;
}
AuthorityMomentMatrix10 hessian(const Jet& value) {
    AuthorityMomentMatrix10 r;
    for (unsigned i=0;i<10;++i) for (unsigned j=0;j<10;++j) r[i][j]=write(value.h[i][j]);
    return r;
}
using RealVector10=std::array<Real,10>;
using RealMatrix10=std::array<RealVector10,10>;
Real dot(const RealVector10& a,const RealVector10& b) {
    Real result;
    for(unsigned i=0;i<10;++i) result=result+a[i]*b[i];
    return result;
}
Real magnitude(const Real& value) { Real r;mpfr_abs(r.value,value.value,MPFR_RNDN);return r; }
Real norm(const RealVector10& value) {
    Real r=dot(value,value);mpfr_sqrt(r.value,r.value,MPFR_RNDN);return r;
}
Real binaryPower(long exponent) {
    Real r;mpfr_set_ui_2exp(r.value,1,exponent,MPFR_RNDN);return r;
}
/** Column-pivoted QR with reorthogonalization, fixed 1024-bit arithmetic.
 * No normal-equation solve, pseudoinverse, diagonal regularizer or precision
 * escalation. Rank and original projection residual guards fail closed.
 */
std::vector<Real> projectMultipliers(const std::vector<RealVector10>& columns,
                                    const RealVector10& gradient,Real* residual) {
    const auto count=columns.size();
    if(count==0 || count>10) throw std::invalid_argument("invalid authority multiplier dimension");
    auto work=columns;
    std::vector<RealVector10> basis(count);
    std::vector<std::vector<Real>> r(count,std::vector<Real>(count));
    std::vector<std::size_t> permutation(count);std::iota(permutation.begin(),permutation.end(),0);
    Real initialNorm;
    for(const auto& col:columns) { const auto length=norm(col);
        if(mpfr_cmp(length.value,initialNorm.value)>0)initialNorm=length; }
    if(mpfr_sgn(initialNorm.value)<=0)throw std::domain_error("zero authority multiplier projection");
    const Real rankFloor=initialNorm*binaryPower(-900);
    std::vector<Real> projected(count);
    for(std::size_t k=0;k<count;++k) {
        std::size_t pivot=k;Real best=norm(work[k]);
        for(std::size_t j=k+1;j<count;++j) { const auto length=norm(work[j]);
            if(mpfr_cmp(length.value,best.value)>0) {best=length;pivot=j;} }
        if(pivot!=k) {
            std::swap(work[k],work[pivot]);std::swap(permutation[k],permutation[pivot]);
            for(std::size_t i=0;i<k;++i)std::swap(r[i][k],r[i][pivot]);
        }
        for(unsigned pass=0;pass<2;++pass)for(std::size_t i=0;i<k;++i) {
            const auto correction=dot(basis[i],work[k]);r[i][k]=r[i][k]+correction;
            for(unsigned a=0;a<10;++a)work[k][a]=work[k][a]-basis[i][a]*correction;
        }
        const auto length=norm(work[k]);
        if(!mpfr_number_p(length.value) || mpfr_cmp(length.value,rankFloor.value)<=0)
            throw std::domain_error("authority multiplier projection is rank deficient or uncertifiable at 1024 bits");
        r[k][k]=length;
        for(unsigned a=0;a<10;++a)basis[k][a]=work[k][a]/length;
        projected[k]=-dot(basis[k],gradient);
        for(std::size_t j=k+1;j<count;++j) {
            for(unsigned pass=0;pass<2;++pass) {
                const auto correction=dot(basis[k],work[j]);r[k][j]=r[k][j]+correction;
                for(unsigned a=0;a<10;++a)work[j][a]=work[j][a]-basis[k][a]*correction;
            }
        }
    }
    std::vector<Real> pivoted(count),answer(count);
    for(std::size_t offset=0;offset<count;++offset) {
        const std::size_t i=count-1-offset;Real value=projected[i];
        for(std::size_t j=i+1;j<count;++j)value=value-r[i][j]*pivoted[j];
        pivoted[i]=value/r[i][i];
    }
    for(std::size_t i=0;i<count;++i)answer[permutation[i]]=pivoted[i];
    RealVector10 optimality=gradient;
    for(std::size_t j=0;j<count;++j)for(unsigned a=0;a<10;++a)
        optimality[a]=optimality[a]+columns[j][a]*answer[j];
    Real scale=Real::constant(1)+norm(gradient);
    for(std::size_t j=0;j<count;++j)scale=scale+norm(columns[j])*magnitude(answer[j]);
    *residual=Real();
    for(const auto& col:columns) {const auto value=magnitude(dot(col,optimality));
        if(mpfr_cmp(value.value,residual->value)>0)*residual=value;}
    const auto limit=binaryPower(-850)*initialNorm*scale;
    if(!mpfr_number_p(residual->value) || mpfr_cmp(residual->value,limit.value)>0)
        throw std::domain_error("authority multiplier solve failed original projection residual guard");
    return answer;
}
} // namespace

QuadricAuthorityMomentEvaluation evaluateSpin2014AuthorityMomentProblem(
    const QuadricObjectiveSettings& settings,const AuthorityMomentVector10& rawMoments) {
    validateSpin2014QuadricObjectiveSettings(settings);
    const Real massScale=Real::constant(settings.solverMassUnitInDeclaredUnits);
    const Real lengthScale=Real::constant(settings.solverLengthUnitInDeclaredUnits);
    const std::array<Real,3> scales{massScale,massScale*lengthScale,massScale*lengthScale*lengthScale};
    std::array<Jet,10> r;
    for (unsigned i=0;i<10;++i) {
        const Real value=parse(rawMoments[i]);
        r[i]=Jet::variable(value,i)*Jet::constant(scales[i==0?0:(i<4?1:2)]);
    }
    if (mpfr_sgn(r[0].v.value)<=0) throw std::domain_error("authority objective requires positive mass");
    const Jet height=r[3]/r[0];
    Jet objective;
    std::vector<Jet> constraints{r[1],r[2]};
    if (settings.objective==Spin2014Objective::StaticLowering) {
        objective=settings.centerWeight*(height*height);
    } else {
        const auto inertia=spin::detail::centralInertia10(r);
        Real denominator; mpfr_abs(denominator.value,inertia.i33.v.value,MPFR_RNDN);
        const Real cutoff=Real::constant(1e-14);
        if (mpfr_cmp(denominator.value,cutoff.value)<=0)
            throw std::domain_error("authority axial central inertia denominator is zero");
        objective=settings.inertiaWeight*spin::detail::transverseInertiaRatio(inertia);
        if (settings.objective==Spin2014Objective::YoYo) constraints.push_back(r[3]);
        constraints.push_back(inertia.i13); constraints.push_back(inertia.i23);
        if (settings.objective==Spin2014Objective::TopLoweringOnly)
            objective=objective+settings.centerWeight*(height*height);
        else if (settings.objective==Spin2014Objective::TopMassAndLowering)
            objective=objective+settings.centerWeight*(r[3]*r[3]);
    }
    QuadricAuthorityMomentEvaluation result;
    result.objectiveProfile=spin2014QuadricObjectiveProfile();
    result.objectiveIdentity=spin2014QuadricObjectiveHash(settings);
    result.rawMomentsAuthority=rawMoments;
    result.objective=write(objective.v);
    result.objectiveGradient=gradient(objective);
    result.objectiveHessian=hessian(objective);
    for (const auto& c:constraints) {
        result.constraints.push_back(write(c.v));
        result.constraintJacobian.push_back(gradient(c));
        result.constraintHessians.push_back(hessian(c));
    }
    return result;
}

QuadricAuthorityNewtonState linearizeSpin2014AuthorityNewtonState(
    const QuadricObjectiveSettings& settings,
    const spin::QuadricCoefficientAuthority10& coefficientAuthority,
    const AuthorityMomentVector10& totalRawMoments,
    const AuthorityMomentMatrix10& nativeMomentJacobian,
    const std::string& evaluationMapIdentity) {
    if(evaluationMapIdentity.empty()) throw std::invalid_argument("authority Newton state requires a material mapping identity");
    QuadricAuthorityNewtonState state;
    state.coefficientPolicy=spin::QuadricCoefficientAuthority10::policy;
    state.evaluationMapIdentity=evaluationMapIdentity;
    state.multiplierMethod="MPFR1024_PIVOTED_REORTHOGONALIZED_QR_V1";
    state.geometryAuthority=coefficientAuthority.geometry();state.dualAuthority=coefficientAuthority.dual();
    state.nativeMomentJacobian=nativeMomentJacobian;
    state.objective=evaluateSpin2014AuthorityMomentProblem(settings,totalRawMoments);
    RealVector10 p,s,inverseDual,gradient;
    RealMatrix10 nativeJacobian,dualJacobian,lagrangianHessian;
    for(unsigned i=0;i<10;++i) {
        p[i]=parse(state.geometryAuthority[i]);s[i]=parse(state.dualAuthority[i]);
        inverseDual[i]=s[i]*Real::constant(i<7?1:0.5);
        gradient[i]=parse(state.objective.objectiveGradient[i]);
        for(unsigned j=0;j<10;++j) {
            nativeJacobian[i][j]=parse(nativeMomentJacobian[i][j]);
            lagrangianHessian[i][j]=parse(state.objective.objectiveHessian[i][j]);
        }
    }
    const Real normalization=norm(inverseDual);
    for(unsigned row=0;row<10;++row)for(unsigned col=0;col<10;++col) {
        Real sum;
        for(unsigned k=0;k<10;++k) {
            const Real projection=(Real::constant(k==col?1:0)-p[k]*p[col])*
                                  Real::constant(col<7?1:0.5)/normalization;
            sum=sum+nativeJacobian[row][k]*projection;
        }
        dualJacobian[row][col]=sum;state.dualMomentJacobian[row][col]=write(sum);
    }
    const auto count=state.objective.constraints.size();
    std::vector<RealVector10> columns(count+1);
    for(std::size_t j=0;j<count;++j)for(unsigned a=0;a<10;++a)
        columns[j][a]=-parse(state.objective.constraintJacobian[j][a]);
    columns[count]=s;
    Real projectionResidual;
    const auto multipliers=projectMultipliers(columns,gradient,&projectionResidual);
    state.multiplierProjectionResidual=write(projectionResidual);
    for(std::size_t j=0;j<count;++j)state.multipliers.push_back(write(multipliers[j]));
    const auto mu=multipliers[count];state.multiplierMu=write(mu);
    RealVector10 optimality=gradient;
    Real squaredResidual;
    for(unsigned a=0;a<10;++a) {
        for(std::size_t j=0;j<count+1;++j)optimality[a]=optimality[a]+columns[j][a]*multipliers[j];
        state.optimalityResidual[a]=write(optimality[a]);squaredResidual=squaredResidual+optimality[a]*optimality[a];
    }
    for(std::size_t j=0;j<count;++j) {
        const auto constraint=parse(state.objective.constraints[j]);
        squaredResidual=squaredResidual+constraint*constraint;
        for(unsigned a=0;a<10;++a)for(unsigned b=0;b<10;++b)
            lagrangianHessian[a][b]=lagrangianHessian[a][b]-multipliers[j]*
                parse(state.objective.constraintHessians[j][a][b]);
    }
    mpfr_sqrt(squaredResidual.value,squaredResidual.value,MPFR_RNDN);
    state.rawEquationResidual=write(squaredResidual);
    const auto dimension=11+count;
    std::vector<std::vector<Real>> system(dimension,std::vector<Real>(dimension));
    std::vector<Real> rhs(dimension);
    for(unsigned a=0;a<10;++a) {
        rhs[a]=-optimality[a];
        for(unsigned b=0;b<10;++b) {
            for(unsigned k=0;k<10;++k)system[a][b]=system[a][b]+lagrangianHessian[a][k]*dualJacobian[k][b];
            if(a==b)system[a][b]=system[a][b]+mu;
        }
        for(std::size_t j=0;j<count;++j)system[a][10+j]=columns[j][a];
        system[a][10+count]=s[a];system[10+count][a]=Real::constant(2)*s[a];
    }
    for(std::size_t j=0;j<count;++j) {
        rhs[10+j]=-parse(state.objective.constraints[j]);
        for(unsigned b=0;b<10;++b)for(unsigned k=0;k<10;++k)
            system[10+j][b]=system[10+j][b]+parse(state.objective.constraintJacobian[j][k])*dualJacobian[k][b];
    }
    std::ostringstream identity;
    identity.imbue(std::locale::classic());
    const auto bind=[&](const std::string& value) {identity<<value.size()<<':'<<value<<';';};
    bind(state.coefficientPolicy);bind(state.multiplierMethod);bind(evaluationMapIdentity);bind(state.objective.objectiveIdentity);
    for(const auto& value:state.geometryAuthority)bind(value);
    for(const auto& value:state.dualAuthority)bind(value);
    for(const auto& value:totalRawMoments)bind(value);
    for(const auto& row:nativeMomentJacobian)for(const auto& value:row)bind(value);
    for(const auto& value:state.multipliers)bind(value);
    bind(state.multiplierMu);bind(state.rawEquationResidual);bind(state.multiplierProjectionResidual);
    state.matrix.resize(dimension,std::vector<std::string>(dimension));state.rightHandSide.resize(dimension);
    for(std::size_t a=0;a<dimension;++a) {
        state.rightHandSide[a]=write(rhs[a]);bind(state.rightHandSide[a]);
        for(std::size_t b=0;b<dimension;++b) {state.matrix[a][b]=write(system[a][b]);bind(state.matrix[a][b]);}
    }
    state.scientificHash=manufacturing::computeContentSha256(identity.str());
    return state;
}
QuadricAuthorityNewtonDirection solveSpin2014AuthorityNewtonDirection(
    const QuadricObjectiveSettings& settings,const QuadricAuthorityNewtonState& supplied,
    const std::string& absoluteLinearResidualLimit) {
    const auto authority=spin::QuadricCoefficientAuthority10::fromSavedState(
        supplied.geometryAuthority,supplied.dualAuthority,supplied.coefficientPolicy,supplied.precisionBits);
    const auto state=linearizeSpin2014AuthorityNewtonState(settings,authority,
        supplied.objective.rawMomentsAuthority,supplied.nativeMomentJacobian,supplied.evaluationMapIdentity);
    if(state.scientificHash!=supplied.scientificHash || state.matrix!=supplied.matrix ||
       state.rightHandSide!=supplied.rightHandSide || state.rawEquationResidual!=supplied.rawEquationResidual ||
       state.multiplierMu!=supplied.multiplierMu || state.multipliers!=supplied.multipliers)
        throw std::invalid_argument("authority Newton direction rejected stale or modified state");
    QuadricAuthorityNewtonDirection result;
    result.stateIdentity=state.scientificHash;
    result.dualDirection.fill("0/1");
    result.linearSolve=spin::solveSpinNewtonAuthorityLinearSystem(
        state.matrix,state.rightHandSide,absoluteLinearResidualLimit);
    if(!result.linearSolve.solved) {
        result.diagnostic=result.linearSolve.diagnostic;return result;
    }
    RealVector10 s,d;
    for(unsigned i=0;i<10;++i) {
        s[i]=parse(state.dualAuthority[i]);
        d[i]=Real::constant(result.linearSolve.candidate.solution(i));
    }
    const auto radial=dot(s,d);
    for(unsigned i=0;i<10;++i)d[i]=d[i]-s[i]*radial;
    const auto length=norm(d);
    result.directionNorm=write(length);
    if(mpfr_cmp(length.value,parse("1/1000000000000000").value)<=0) {
        result.diagnostic="Newton direction vanished before convergence";return result;
    }
    result.initialArcLength=write(mpfr_cmp(length.value,Real::constant(0.5).value)<0?length:Real::constant(0.5));
    std::ostringstream identity;identity.imbue(std::locale::classic());
    const auto bind=[&](const std::string& value){identity<<value.size()<<':'<<value<<';';};
    bind("MPFR1024_EQ14_PROTECTED_DUAL_DIRECTION_V1");bind(result.stateIdentity);
    bind(absoluteLinearResidualLimit);bind(result.linearSolve.residualInfinityLower);
    bind(result.linearSolve.residualInfinityUpper);bind(result.linearSolve.backwardErrorUpper);
    for(unsigned i=0;i<10;++i){result.dualDirection[i]=write(d[i]);bind(result.dualDirection[i]);}
    bind(result.directionNorm);bind(result.initialArcLength);
    result.scientificHash=manufacturing::computeContentSha256(identity.str());
    result.defined=true;
    result.diagnostic="protected authority equation passed; MPFR1024 tangent direction ready for actual trial";
    return result;
}
QuadricAuthorityNewtonPoint evaluateSpin2014AuthorityNewtonPoint(
    const spin::QuadricAuthorityMaterialContext& context,const QuadricObjectiveSettings& settings,
    const spin::QuadricCoefficientAuthority10& authority,const spin::QuadricMomentEvaluatorSettings& evaluation) {
    QuadricAuthorityNewtonPoint point;
    point.evaluatorRequestIdentity=spin::computeQuadricEvaluatorRequestHash(context.request(),evaluation);
    point.material=spin::evaluateQuadricMaterialMoments(context,authority,evaluation);
    point.state=linearizeSpin2014AuthorityNewtonState(settings,authority,point.material.totalMomentsTarget,
        point.material.derivative.differential.coefficientJacobian,point.material.scientificHash);
    return point;
}

QuadricAuthorityNewtonTrial evaluateSpin2014AuthorityNewtonTrial(
    const spin::QuadricAuthorityMaterialContext& context,const QuadricObjectiveSettings& settings,
    const QuadricAuthorityNewtonPoint& base,const spin::QuadricCoefficientAuthority10::RationalVector& direction,
    const std::string& arcLength,const spin::QuadricMomentEvaluatorSettings& evaluation) {
    const auto step=parse(arcLength);
    if(mpfr_sgn(step.value)<=0 || mpfr_cmp(step.value,Real::constant(0.5).value)>0)
        throw std::invalid_argument("authority trial requires positive arc at most the frozen 1/2");
    if(base.evaluatorRequestIdentity!=spin::computeQuadricEvaluatorRequestHash(context.request(),evaluation) ||
       base.material.contextIdentity!=context.identity() || !base.material.derivative.defined ||
       base.state.evaluationMapIdentity!=base.material.scientificHash)
        throw std::invalid_argument("authority trial rejected stale material mapping identity");
    const auto authority=spin::QuadricCoefficientAuthority10::fromSavedState(
        base.state.geometryAuthority,base.state.dualAuthority,base.state.coefficientPolicy,base.state.precisionBits);
    const auto rebuilt=linearizeSpin2014AuthorityNewtonState(settings,authority,base.material.totalMomentsTarget,
        base.material.derivative.differential.coefficientJacobian,base.material.scientificHash);
    if(rebuilt.scientificHash!=base.state.scientificHash ||
       rebuilt.rawEquationResidual!=base.state.rawEquationResidual || rebuilt.multiplierMu!=base.state.multiplierMu)
        throw std::invalid_argument("authority trial rejected stale numerical point");
    const auto candidate=authority.geodesicUpdate(direction,arcLength);
    QuadricAuthorityNewtonTrial result;
    result.arcLength=arcLength;result.baseStateIdentity=rebuilt.scientificHash;
    result.point=evaluateSpin2014AuthorityNewtonPoint(context,settings,candidate,evaluation);
    const auto raw=parse(result.point.state.rawEquationResidual);
    const bool badStationary=mpfr_cmp(raw.value,parse("1/1000000000000").value)<=0 &&
        mpfr_sgn(parse(result.point.state.multiplierMu).value)<=0;
    const bool decreased=mpfr_cmp(raw.value,parse(rebuilt.rawEquationResidual).value)<0;
    result.topologyChanged=base.material.material.topologyHash!=result.point.material.material.topologyHash;
    result.accepted=decreased && !badStationary;
    result.reason=badStationary?"stationary_multiplier_nonpositive":
        !decreased?"no_residual_decrease":result.topologyChanged?"topology_transition_relinearized":"residual_decreased";
    std::ostringstream identity;identity.imbue(std::locale::classic());
    const auto bind=[&](const std::string& value){identity<<value.size()<<':'<<value<<';';};
    bind("MPFR1024_ACTUAL_MATERIAL_NEWTON_TRIAL_V1");bind(result.baseStateIdentity);
    bind(result.point.state.scientificHash);bind(arcLength);bind(result.reason);
    for(const auto& v:direction)bind(v);
    result.scientificHash=manufacturing::computeContentSha256(identity.str());
    return result;
}
QuadricAuthorityNewtonAdvance advanceSpin2014AuthorityNewtonPoint(
    const spin::QuadricAuthorityMaterialContext& context,const QuadricObjectiveSettings& settings,
    const QuadricAuthorityNewtonPoint& base,const spin::QuadricMomentEvaluatorSettings& evaluation,
    const std::string& absoluteLinearResidualLimit,unsigned budget,const QuadricAuthorityEventSink& eventSink) {
    if(budget==0 || budget>60)throw std::invalid_argument("remaining authority evaluation budget must be 1..60");
    if(base.evaluatorRequestIdentity!=spin::computeQuadricEvaluatorRequestHash(context.request(),evaluation) ||
       base.material.contextIdentity!=context.identity() || base.state.evaluationMapIdentity!=base.material.scientificHash)
        throw std::invalid_argument("authority backtracking rejected stale material context");
    QuadricAuthorityNewtonAdvance result;
    result.direction=solveSpin2014AuthorityNewtonDirection(settings,base.state,absoluteLinearResidualLimit);
    if(!result.direction.defined)result.status="linear_solve_or_direction_failed";
    else {
        const auto authority=spin::QuadricCoefficientAuthority10::fromSavedState(
            base.state.geometryAuthority,base.state.dualAuthority,base.state.coefficientPolicy,base.state.precisionBits);
        mpq_class arc(result.direction.initialArcLength);
        const mpq_class minimum(std::numeric_limits<double>::epsilon());
        while(result.productionEvaluations<budget && arc>minimum) {
            QuadricAuthorityNewtonTrialRecord record;
            record.arcLength=rationalText(arc);
            const auto candidate=authority.geodesicUpdate(result.direction.dualDirection,record.arcLength);
            record.geometryAuthority=candidate.geometry();record.dualAuthority=candidate.dual();
            ++result.productionEvaluations;
            QuadricAuthorityNewtonEvent event;
            event.type="trial_evaluation_start";event.evaluationIndex=result.productionEvaluations;
            event.geometryAuthority=record.geometryAuthority;event.dualAuthority=record.dualAuthority;
            event.arcLength=record.arcLength;event.stateIdentity=base.state.scientificHash;
            if(eventSink)eventSink(event);
            try {
                auto trial=evaluateSpin2014AuthorityNewtonTrial(context,settings,base,
                    result.direction.dualDirection,record.arcLength,evaluation);
                record.reason=trial.reason;record.accepted=trial.accepted;
                record.rawEquationResidual=trial.point.state.rawEquationResidual;
                record.trialIdentity=trial.scientificHash;
                result.accepted=trial.accepted;
                if(trial.accepted)result.finalTrial=std::move(trial);
            } catch(const std::bad_alloc&) {throw;}
              catch(const std::exception& error) {
                record.reason="geometry_or_derivative_failure";record.diagnostic=error.what();
            }
            result.trials.push_back(std::move(record));
            const auto& recorded=result.trials.back();
            event.type="trial_evaluation_terminal";event.status=recorded.reason;
            event.diagnostic=recorded.diagnostic;event.rawResidual=recorded.rawEquationResidual;
            event.trialIdentity=recorded.trialIdentity;event.accepted=recorded.accepted;
            if(eventSink)eventSink(event);
            if(result.accepted){result.status="accepted";break;}
            if(result.trials.back().reason=="stationary_multiplier_nonpositive") {
                result.status="stationary_multiplier_nonpositive_restart";break;
            }
            arc/=2;
        }
        if(result.status.empty())result.status=result.productionEvaluations>=budget?
            "evaluation_budget_exhausted":"line_search_exhausted";
    }
    std::ostringstream identity;identity.imbue(std::locale::classic());
    const auto bind=[&](const std::string& value){identity<<value.size()<<':'<<value<<';';};
    bind("MPFR1024_ACTUAL_NEWTON_BACKTRACKING_V1");bind(context.identity());bind(base.state.scientificHash);
    bind(result.direction.scientificHash);bind(result.direction.diagnostic);bind(result.status);
    bind(std::to_string(budget));bind(std::to_string(result.productionEvaluations));
    for(const auto& record:result.trials) {
        bind(record.arcLength);bind(record.reason);bind(record.diagnostic);bind(record.rawEquationResidual);bind(record.trialIdentity);
        for(const auto& v:record.geometryAuthority)bind(v);
        for(const auto& v:record.dualAuthority)bind(v);
    }
    result.scientificHash=manufacturing::computeContentSha256(identity.str());
    return result;
}
std::string classifySpin2014AuthorityNumericalState(
    const QuadricAuthorityNewtonState& state,bool initialState) {
    const auto read=[](const std::string& text) {
        mpq_class q;
        if(mpq_set_str(q.get_mpq_t(),text.c_str(),10)!=0 || q.get_den()==0)
            throw std::invalid_argument("invalid authority numerical decision input");
        q.canonicalize();
        if(rationalText(q)!=text)
            throw std::invalid_argument("noncanonical authority numerical decision input");
        return q;
    };
    const auto raw=read(state.rawEquationResidual),mu=read(state.multiplierMu);
    if(raw<0)throw std::invalid_argument("authority residual norm cannot be negative");
    bool sensitive=false;
    for(const auto& row:state.nativeMomentJacobian)
        for(const auto& entry:row)if(read(entry)!=0)sensitive=true;
    if(initialState && mu<=0)return "initial_multiplier_nonpositive";
    if(!sensitive)return initialState?"initial_zero_sensitivity":"zero_sensitivity_restart";
    if(raw<=mpq_class(1,1000000000000L))
        return mu>0?"stationary_candidate":"stationary_multiplier_nonpositive_restart";
    return "continue";
}
std::vector<spin::QuadricCoefficientAuthority10> makeSpin2014AuthorityInitialSequence(
    spin::SpinMode mode,const spin::SpinSolverSettings& settings) {
    if(mode!=spin::SpinMode::SpinningTop && mode!=spin::SpinMode::YoYo)
        throw std::invalid_argument("authority full-rank initialization only supports top or yo-yo");
    if(settings.maximumRestarts<1 || settings.maximumRestarts>8 ||
       settings.maximumFunctionEvaluationsPerRestart!=60 || settings.equationTolerance!=1e-12 ||
       settings.lineSearchReduction!=0.5 || settings.maximumGeodesicStep!=0.5)
        throw std::invalid_argument("authority initialization requires the frozen numerical protocol");
    const auto starts=spin::deterministicInitialQuadrics(mode,settings);
    std::vector<spin::QuadricCoefficientAuthority10> result;
    result.reserve(starts.size());
    for(const auto& start:starts) {
        spin::QuadricCoefficientAuthority10::RationalVector rational;
        for(unsigned i=0;i<10;++i)rational[i]=rationalText(mpq_class(start.values()(i)));
        result.push_back(spin::QuadricCoefficientAuthority10::fromUnnormalizedGeometry(rational));
    }
    return result;
}
QuadricAuthorityNewtonAttempt runSpin2014AuthorityNewtonAttempt(
    const spin::QuadricAuthorityMaterialContext& context,const QuadricObjectiveSettings& settings,
    const spin::QuadricCoefficientAuthority10& initial,const spin::QuadricMomentEvaluatorSettings& evaluation,
    const std::string& linearLimit,unsigned budget,const QuadricAuthorityEventSink& eventSink) {
    if(budget==0 || budget>60)throw std::invalid_argument("authority attempt budget must be 1..60");
    validateSpin2014QuadricObjectiveSettings(settings);
    if(mpfr_sgn(parse(linearLimit).value)<=0)throw std::invalid_argument("linear residual limit must be positive");
    QuadricAuthorityNewtonAttempt result;
    result.initialGeometryAuthority=initial.geometry();result.initialDualAuthority=initial.dual();
    QuadricAuthorityNewtonEvent event;
    event.type="attempt_start";event.geometryAuthority=initial.geometry();event.dualAuthority=initial.dual();
    if(eventSink)eventSink(event);
    std::optional<QuadricAuthorityNewtonPoint> current;
    ++result.productionEvaluations;
    event.type="initial_evaluation_start";event.evaluationIndex=result.productionEvaluations;
    if(eventSink)eventSink(event);
    try {current=evaluateSpin2014AuthorityNewtonPoint(context,settings,initial,evaluation);}
    catch(const std::bad_alloc&) {throw;}
    catch(const std::exception& error) {result.status="initial_geometry_or_derivative_failure";result.diagnostic=error.what();}
    event.type="initial_evaluation_terminal";event.status=result.status;event.diagnostic=result.diagnostic;
    if(current) {
        event.stateIdentity=current->state.scientificHash;event.materialIdentity=current->material.scientificHash;
        event.rawResidual=current->state.rawEquationResidual;
    }
    if(eventSink)eventSink(event);
    if(current) {
        result.initialTotalMoments=current->material.totalMomentsTarget;
        result.initialObjective=current->state.objective.objective;
        result.initialMaterialIdentity=current->material.scientificHash;
        result.initialStateIdentity=current->state.scientificHash;
        result.initialRawResidual=current->state.rawEquationResidual;
        result.initialMu=current->state.multiplierMu;
        auto screen=classifySpin2014AuthorityNumericalState(current->state,true);
        while(true) {
            result.finalStateIdentity=current->state.scientificHash;
            result.finalRawResidual=current->state.rawEquationResidual;
            if(screen=="stationary_candidate") {
                result.status=screen;result.stationaryCandidate=std::move(current);break;
            }
            if(screen!="continue") {result.status=screen;break;}
            if(result.productionEvaluations>=budget) {result.status="evaluation_budget_exhausted";break;}
            const unsigned offset=result.productionEvaluations;
            const auto forward=[&](const QuadricAuthorityNewtonEvent& observed) {
                auto translated=observed;translated.evaluationIndex+=offset;
                if(eventSink)eventSink(translated);
            };
            auto step=advanceSpin2014AuthorityNewtonPoint(context,settings,*current,evaluation,
                linearLimit,budget-result.productionEvaluations,forward);
            result.productionEvaluations+=step.productionEvaluations;
            result.stepIdentities.push_back(step.scientificHash);
            for(auto& record:step.trials)result.trials.push_back(std::move(record));
            if(!step.accepted) {result.status=step.status;break;}
            if(!step.finalTrial || !step.finalTrial->accepted)
                throw std::logic_error("accepted authority step omitted committable material");
            current=std::move(step.finalTrial->point);
            ++result.acceptedUpdates;
            result.acceptedStateIdentities.push_back(current->state.scientificHash);
            event.type="state_accepted";event.evaluationIndex=result.productionEvaluations;
            event.geometryAuthority=current->state.geometryAuthority;event.dualAuthority=current->state.dualAuthority;
            event.stateIdentity=current->state.scientificHash;event.materialIdentity=current->material.scientificHash;
            event.rawResidual=current->state.rawEquationResidual;event.accepted=true;
            event.status="accepted";event.diagnostic.clear();
            if(eventSink)eventSink(event);
            screen=classifySpin2014AuthorityNumericalState(current->state,false);
        }
    }
    std::ostringstream identity;identity.imbue(std::locale::classic());
    const auto bind=[&](const std::string& value){identity<<value.size()<<':'<<value<<';';};
    bind("MPFR1024_ACTUAL_NEWTON_ATTEMPT_V1");bind(context.identity());
    bind(spin::computeQuadricEvaluatorRequestHash(context.request(),evaluation));
    bind(spin2014QuadricObjectiveHash(settings));bind(initial.policy);
    bind(linearLimit);bind(std::to_string(budget));bind(std::to_string(result.productionEvaluations));
    bind(std::to_string(result.acceptedUpdates));bind(result.status);bind(result.diagnostic);
    bind(result.initialStateIdentity);bind(result.finalStateIdentity);
    bind(result.initialObjective);bind(result.initialMaterialIdentity);
    bind(result.initialTotalMoments?"initial_material_present":"initial_material_absent");
    if(result.initialTotalMoments)for(const auto& value:*result.initialTotalMoments)bind(value);
    for(const auto& v:result.initialGeometryAuthority)bind(v);
    for(const auto& v:result.initialDualAuthority)bind(v);
    for(const auto& value:result.stepIdentities)bind(value);
    for(const auto& value:result.acceptedStateIdentities)bind(value);
    result.scientificHash=manufacturing::computeContentSha256(identity.str());
    event.type="attempt_terminal";event.evaluationIndex=result.productionEvaluations;
    event.status=result.status;event.diagnostic=result.diagnostic;
    event.stateIdentity=result.finalStateIdentity;event.rawResidual=result.finalRawResidual;
    event.accepted=false;
    if(eventSink)eventSink(event);
    return result;
}
QuadricAuthorityNewtonRestarts runSpin2014AuthorityNewtonRestarts(
    const spin::QuadricAuthorityMaterialContext& context,const QuadricObjectiveSettings& settings,
    const spin::QuadricMomentEvaluatorSettings& evaluation,const std::string& linearLimit,unsigned budget,
    const QuadricAuthorityEventSink& eventSink) {
    if(budget==0 || budget>60)throw std::invalid_argument("restart evaluation budget must be 1..60");
    validateSpin2014QuadricObjectiveSettings(settings);
    const auto& request=context.request();
    const bool yo=settings.objective==Spin2014Objective::YoYo;
    if(settings.objective==Spin2014Objective::StaticLowering ||
       request.mode!=(yo?spin::SpinMode::YoYo:spin::SpinMode::SpinningTop))
        throw std::invalid_argument("authority rotation objective/request mode mismatch");
    const auto starts=makeSpin2014AuthorityInitialSequence(request.mode,request.settings);
    QuadricAuthorityNewtonRestarts result;
    result.seed=request.settings.deterministicSeed;
    result.protocolScope=budget==60?"CONFIGURED_FULL_ATTEMPT_BUDGET":"DIAGNOSTIC_REDUCED_ATTEMPT_BUDGET";
    for(const auto& initial:starts) {
        const auto forward=[&](const QuadricAuthorityNewtonEvent& observed) {
            auto event=observed;event.restart=static_cast<unsigned>(result.attempts.size());
            event.sequence=static_cast<unsigned>(result.events.size());
            result.events.push_back(event);
            if(eventSink)eventSink(event);
        };
        auto attempt=runSpin2014AuthorityNewtonAttempt(context,settings,initial,evaluation,linearLimit,budget,forward);
        result.productionEvaluations+=attempt.productionEvaluations;
        const bool candidate=attempt.stationaryCandidate.has_value();
        result.attempts.push_back(std::move(attempt));
        if(candidate) {result.stationaryAttemptIndex=static_cast<unsigned>(result.attempts.size()-1);break;}
    }
    result.status=result.stationaryAttemptIndex?"stationary_candidate":"restart_budget_exhausted";
    std::ostringstream identity;identity.imbue(std::locale::classic());
    const auto bind=[&](const std::string& value){identity<<value.size()<<':'<<value<<';';};
    bind("MPFR1024_REQUEST_SEEDED_NEWTON_RESTARTS_V1");bind(context.identity());
    bind(spin::computeQuadricEvaluatorRequestHash(context.request(),evaluation));
    bind(spin2014QuadricObjectiveHash(settings));bind(std::to_string(result.seed));
    bind(std::to_string(starts.size()));bind(std::to_string(budget));bind(result.protocolScope);
    bind(std::to_string(result.productionEvaluations));bind(result.status);
    for(const auto& attempt:result.attempts)bind(attempt.scientificHash);
    result.scientificHash=manufacturing::computeContentSha256(identity.str());
    QuadricAuthorityNewtonEvent event;
    event.type="run_terminal";event.status=result.status;event.evaluationIndex=result.productionEvaluations;
    event.restart=static_cast<unsigned>(result.attempts.size());event.sequence=static_cast<unsigned>(result.events.size());
    result.events.push_back(event);
    if(eventSink)eventSink(event);
    return result;
}
} // namespace spin2014
