#include "spin/QuadricInverseGaussMap.h"
#include "spin/RationalSourceRepresentative.h"

#include <Eigen/Core>
#include <Eigen/LU>

#include <CGAL/Gmpq.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>


namespace
{


using Quadric = spin::QuadricCoefficients10;
using Vector10 = Quadric::Vector;
using Jacobian = spin::QuadricInverseGaussMap::CoefficientJacobian;


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=quadric-inverse-gauss-map "
              << message << '\n';
    return EXIT_FAILURE;
}

Jacobian sphereOracle(double coefficientScale, double radius,
                      const Eigen::Vector3d& unitNormal);

int rationalFinitePointRegression()
{
    // A binary64 conversion here would erase the derivative-scale motion.
    // On this literal unit sphere, every coordinate at (1,1,1) is sqrt(1/3).
    using Q=CGAL::Gmpq;
    Vector10 values=Vector10::Zero();
    values(0)=-.5;values(4)=values(5)=values(6)=.5;
    const spin::QuadricInverseGaussMap map(Quadric::fromNormalized(values));
    const auto point=map.rationalPoint(Eigen::Vector3d(1,1,1));
    if(!point.defined)return fail("rational finite source point rejected: "+point.diagnostic);
    CGAL::Gmpz denominator(1);
    mpz_mul_2exp(denominator.mpz(),denominator.mpz(),384);
    const Q bound(CGAL::Gmpz(1),denominator);
    for(int axis=0;axis<3;++axis) {
        const Q x(point.coordinates[axis].c_str());
        const Q lo(point.lower[axis].c_str()),hi(point.upper[axis].c_str());
        if(!(lo>0 && 3*lo*lo<=1 && 3*hi*hi>=1 && lo<=x && x<=hi && hi-lo<=bound))
            return fail("rational point lacks the independent sqrt(1/3) enclosure");
        if(x==Q(CGAL::to_double(x)))return fail("rational point was rounded through binary64");
    }
    // These coordinates are exactly rational despite the inverse-Gauss root.
    const auto rational=map.rationalPoint(Eigen::Vector3d(2,-3,6));
    if(!rational.defined)return fail("rational axis fixture rejected");
    for(int axis=0;axis<3;++axis) {
        const Q expected(axis==0?2:(axis==1?-3:6),7);
        const Q lo(rational.lower[axis].c_str()),hi(rational.upper[axis].c_str());
        if(!(lo<=expected && expected<=hi && hi-lo<=bound))
            return fail("signed rational root enclosure selects the wrong branch");
    }
    // Exact sign brackets independently check severe center/root cancellation
    // on BOTH sides of rank two, without using the implementation's formula.
    for(double epsilon:{-1e-14,1e-14}) {
        Vector10 almostRankTwo;
        almostRankTwo<<.125,0,0,1,1,2,epsilon,0,0,0;
        const auto q=Quadric::fromUnnormalized(almostRankTwo);
        const auto p=spin::QuadricInverseGaussMap(q).rationalPoint(Eigen::Vector3d::UnitZ());
        if(!p.defined)return fail("rational cancellation point rejected");
        const Q lo(p.lower[2].c_str()),hi(p.upper[2].c_str());
        const Q a(q.values()(6)),b(q.values()(3)),c(q.values()(0));
        if(!(lo<0 && hi<0 && b+2*a*lo>0 && b+2*a*hi>0 &&
             a*lo*lo+b*lo+c<=0 && a*hi*hi+b*hi+c>=0 && hi-lo<=bound))
            return fail("rational cancellation bracket is not the positive-gradient root");
        for(int axis=0;axis<2;++axis)
            if(Q(p.coordinates[axis].c_str())!=0)
                return fail("axis normal acquired an off-axis source coordinate");
    }
    // Three genuinely small source perturbations would be lost if point
    // construction returned binary64 first. Use the ACTUAL represented
    // coefficient chords, not ideal pre-normalization steps.
    const auto oracle=sphereOracle(.5,1.,Eigen::Vector3d::Constant(1/std::sqrt(3.)));
    for(double h:{1e-15,5e-16,2.5e-16}) {
        auto rawPositive=values,rawNegative=values;
        rawPositive<<-.5+h,0,0,0,.5,.5,.5,0,0,0;
        rawNegative<<-.5-h,0,0,0,.5,.5,.5,0,0,0;
        const auto qp=Quadric::fromUnnormalized(rawPositive);
        const auto qm=Quadric::fromUnnormalized(rawNegative);
        const auto p=spin::QuadricInverseGaussMap(qp).rationalPoint(Eigen::Vector3d(1,1,1));
        const auto m=spin::QuadricInverseGaussMap(qm).rationalPoint(Eigen::Vector3d(1,1,1));
        if(!p.defined || !m.defined)return fail("rational secant endpoint rejected");
        for(int axis=0;axis<3;++axis) {
            Q prediction(0);
            for(int column=0;column<10;++column)
                prediction+=Q(oracle(axis,column))*(Q(qp.values()(column))-Q(qm.values()(column)))/(2*Q(h));
            const Q observed=(Q(p.coordinates[axis].c_str())-Q(m.coordinates[axis].c_str()))/(2*Q(h));
            if(std::abs(CGAL::to_double(observed-prediction))>1e-6)
                return fail("three-step rational source secant loses small coefficient motion");
        }
    }
    if(map.rationalPoint(Eigen::Vector3d::Zero()).defined ||
       map.rationalPoint(Eigen::Vector3d(0,0,1e-16)).defined)
        return fail("rational source policy bypasses the unchanged normal gate");
    values<<.125,0,0,1,1,2,0,0,0,0;
    const auto rankTwo=spin::QuadricInverseGaussMap(Quadric::fromUnnormalized(values))
        .rationalPoint(Eigen::Vector3d::UnitZ());
    if(rankTwo.defined)return fail("full-rank rational path silently substituted rank two");
    // Explicit invocation must leave the historical binary64 path unchanged.
    const auto legacy=map.point(Eigen::Vector3d(1,1,1));
    if(!legacy.valid() || (*legacy.point-Eigen::Vector3d::Constant(1/std::sqrt(3.))).norm()>1e-15)
        return fail("explicit rational path changed the legacy map");
    return EXIT_SUCCESS;
}

int compactRationalPointRegression()
{
    // Catch interval escape, negative rounding, binary64-first selection and
    // replacing exact/narrow non-dyadic coordinates with a nearby dyadic.
    using Q=CGAL::Gmpq;
    CGAL::Gmpz denominator(1);
    mpz_mul_2exp(denominator.mpz(),denominator.mpz(),1024);
    const Q ulp(CGAL::Gmpz(1),denominator);
    for(const Q center:{Q(1,3),Q(-1,3)}) {
        const Q lo=center-ulp,hi=center+ulp;
        const Q x=spin::rationalSourceRepresentative(lo,hi,true);
        if(!(lo<=x && x<=hi) || denominator%x.denominator()!=0)
            return fail("compact representative escapes enclosure or dyadic grid");
        if(spin::rationalSourceRepresentative(lo,hi,false)!=center)
            return fail("legacy interval midpoint changed");
        if(spin::rationalSourceRepresentative(center,center,true)!=center ||
           spin::rationalSourceRepresentative(center-ulp/100,center+ulp/100,true)!=center)
            return fail("compact policy loses exact/narrow fallback");
    }
    try {
        spin::rationalSourceRepresentative(Q(1),Q(-1),true);
        return fail("reversed source interval accepted");
    } catch(const std::invalid_argument&) {}
    Vector10 values=Vector10::Zero();
    values(0)=-.5;values(4)=values(5)=values(6)=.5;
    const spin::QuadricInverseGaussMap map(Quadric::fromNormalized(values));
    const auto old=map.rationalPoint(Eigen::Vector3d(1,-1,1));
    const auto compact=map.rationalPoint(Eigen::Vector3d(1,-1,1),true);
    if(!compact.defined || compact.lower!=old.lower || compact.upper!=old.upper)
        return fail("compact policy changes analytic bounds or source admission");
    for(int axis=0;axis<3;++axis) {
        const Q x(compact.coordinates[axis].c_str());
        const Q lo(compact.lower[axis].c_str()),hi(compact.upper[axis].c_str());
        if(!(lo<=x && x<=hi) || denominator%x.denominator()!=0)
            return fail("actual compact finite point is not enclosed dyadic");
    }
    // This sphere's old sqrt enclosures already have dyadic midpoints.
    // Equal representatives are valid; the signed 1/3 fixtures above require
    // actual selection rather than asserting all inputs must change.
    return EXIT_SUCCESS;
}


bool close(
    const Eigen::Vector3d& actual,
    const Eigen::Vector3d& expected,
    double tolerance
)
{
    return (actual - expected).norm() <= tolerance;
}


Jacobian sphereOracle(double coefficientScale, double radius,
                      const Eigen::Vector3d& unitNormal)
{
    Jacobian expected = Jacobian::Zero();
    expected.col(0) = -unitNormal / (2.0 * coefficientScale * radius);
    for (int coordinate = 0; coordinate < 3; ++coordinate)
    {
        Eigen::Vector3d basis = Eigen::Vector3d::Zero();
        basis(coordinate) = 1.0;
        expected.col(1 + coordinate) = -basis / (2.0 * coefficientScale);

        const double pointCoordinate = radius * unitNormal(coordinate);
        expected.col(4 + coordinate) =
            unitNormal * pointCoordinate * pointCoordinate /
                (2.0 * coefficientScale * radius) -
            basis * pointCoordinate / coefficientScale;
    }
    const std::array<std::array<int, 2>, 3> crossCoordinates = {{
        {{1, 2}}, {{0, 2}}, {{0, 1}}
    }};
    for (int cross = 0; cross < 3; ++cross)
    {
        const int first = crossCoordinates[cross][0];
        const int second = crossCoordinates[cross][1];
        Eigen::Vector3d firstBasis = Eigen::Vector3d::Zero();
        Eigen::Vector3d secondBasis = Eigen::Vector3d::Zero();
        firstBasis(first) = 1.0;
        secondBasis(second) = 1.0;
        const double firstPoint = radius * unitNormal(first);
        const double secondPoint = radius * unitNormal(second);
        expected.col(7 + cross) =
            unitNormal * firstPoint * secondPoint /
                (coefficientScale * radius) -
            (firstBasis * secondPoint + secondBasis * firstPoint) /
                coefficientScale;
    }
    return expected;
}


Jacobian centralDifference(
    const Vector10& normalizedCoefficients,
    const Eigen::Vector3d& normal,
    double step
)
{
    Jacobian result;
    for (int coordinate = 0; coordinate < normalizedCoefficients.size();
         ++coordinate)
    {
        Vector10 positive = normalizedCoefficients;
        Vector10 negative = normalizedCoefficients;
        positive(coordinate) += step;
        negative(coordinate) -= step;
        const spin::QuadricInverseGaussPointResult positivePoint =
            spin::QuadricInverseGaussMap(
                Quadric::fromUnnormalized(positive)
            ).point(normal);
        const spin::QuadricInverseGaussPointResult negativePoint =
            spin::QuadricInverseGaussMap(
                Quadric::fromUnnormalized(negative)
            ).point(normal);
        if (!positivePoint.valid() || !negativePoint.valid())
        {
            throw std::runtime_error(
                "finite-difference point unexpectedly became invalid"
            );
        }
        result.col(coordinate) =
            (*positivePoint.point - *negativePoint.point) / (2.0 * step);
    }
    return result;
}


Eigen::Vector3d implicitPoint(const Quadric& q, const Eigen::Vector3d& normal,
                              const Eigen::Vector3d& initial)
{
    using V = Eigen::Matrix<long double,4,1>;
    using M = Eigen::Matrix<long double,4,4>;
    const auto a=q.quadraticMatrix().cast<long double>().eval();
    const auto b=q.linearVector().cast<long double>().eval();
    const auto n=normal.cast<long double>().eval();
    V state; state.head<3>()=initial.cast<long double>();
    state(3)=(2*a*state.head<3>()+b).dot(n)/n.squaredNorm();
    for(int iteration=0;iteration<30;++iteration) {
        const auto x=state.head<3>().eval();
        const auto gradient=(2*a*x+b).eval();
        V residual;residual.head<3>()=gradient-state(3)*n;
        residual(3)=x.dot(a*x)+b.dot(x)+static_cast<long double>(q.constant());
        if(residual.norm()<1e-17L)return x.cast<double>();
        M jac=M::Zero();jac.topLeftCorner<3,3>()=2*a;
        jac.topRightCorner<3,1>()=-n;jac.bottomLeftCorner<1,3>()=gradient.transpose();
        auto solve=jac.fullPivLu();
        if(solve.rank()!=4)throw std::runtime_error("independent point oracle singular");
        state-=solve.solve(residual);
    }
    throw std::runtime_error("independent point oracle did not converge");
}

int algebraicTransitionRegression()
{
    // Numerical admissibility must not jump at exactly rank two.
    for(double epsilon: {-1e-10,0.,1e-10}) {
        Vector10 values;values<<.125,0,0,1,1,2,epsilon,0,0,0;
        const spin::QuadricInverseGaussMap map(Quadric::fromUnnormalized(values));
        const auto below=map.algebraicCoefficientDifferential(Eigen::Vector3d(0,0,1e-8));
        if(below.valid() || below.point.status!=spin::QuadricInverseGaussPointStatus::BelowNumericalThreshold)
            return fail("algebraic rank transition bypasses normal-magnitude cutoff");
        const auto above=map.algebraicCoefficientDifferential(Eigen::Vector3d(0,0,1e-6));
        if(!above.valid())return fail("admissible scaled normal rejected at rank transition");
    }
    // Catches cancellation in A^-1(n/w-b/2), not merely a type classifier.
    for(double epsilon: {1e-6,1e-8,1e-10,1e-12,1e-14,0.,
                         -1e-6,-1e-8,-1e-10,-1e-12,-1e-14}) {
        Vector10 values;values<<.125,0,0,1,1,2,epsilon,0,0,0;
        const auto q=Quadric::fromUnnormalized(values);
        const auto result=spin::QuadricInverseGaussMap(q)
            .algebraicCoefficientDifferential(Eigen::Vector3d(0,0,1));
        if(!result.valid())return fail("regular algebraic transition rejected: "+result.diagnostic);
        const long double a=q.values()(6),b=q.values()(3),c=q.values()(0);
        const double z=static_cast<double>(-2*c/(b+sqrtl(b*b-4*a*c)));
        const double error=(*result.point.point-Eigen::Vector3d(0,0,z)).norm();
        std::cout<<"algebraic_transition="<<epsilon<<" point_error="<<error<<'\n';
        if(error>1e-14)return fail("near-paraboloid finite point cancellation");
    }
    // Actual new map, not only an implicit oracle, on both sides of rank two.
    Vector10 values;values<<.125,.5,-.25,1.,1.,2.,0.,0.,0.,0.;
    const auto q=Quadric::fromUnnormalized(values);
    const Eigen::Vector3d normal(1.,-.75,1.);
    const auto base=spin::QuadricInverseGaussMap(q).algebraicCoefficientDifferential(normal);
    if(!base.valid())return fail("algebraic rank-two base missing");
    double firstError=0.;
    for(double step: {1e-4,3e-5,1e-5}) {
        Jacobian difference;
        for(int col=0;col<10;++col) {
            auto plus=q.values(),minus=q.values();plus(col)+=step;minus(col)-=step;
            const auto p=spin::QuadricInverseGaussMap(Quadric::fromUnnormalized(plus))
                .algebraicCoefficientDifferential(normal);
            const auto m=spin::QuadricInverseGaussMap(Quadric::fromUnnormalized(minus))
                .algebraicCoefficientDifferential(normal);
            if(!p.valid() || !m.valid())return fail("actual rank-changing map unavailable");
            difference.col(col)=(*p.point.point-*m.point.point)/(2*step);
            if((difference.col(col)-base.coefficientJacobian->col(col)).norm()/
                std::max(1.,difference.col(col).norm())>1e-6)
                return fail("actual algebraic map coefficient derivative mismatch");
            const auto independent=implicitPoint(Quadric::fromUnnormalized(plus),normal,*base.point.point);
            if((*p.point.point-independent).norm()>1e-13)
                return fail("algebraic point disagrees with independent implicit solver");
        }
        const double error=(difference-*base.coefficientJacobian).norm()/std::max(1.,difference.norm());
        if(step==1e-4)firstError=error;
        if(step==1e-5 && error>firstError*.2)return fail("actual map derivative lacks step convergence");
        std::cout<<"actual_algebraic_step="<<step<<" full_direction_error="<<error<<'\n';
    }
    return EXIT_SUCCESS;
}

int paraboloidFiniteMapRegression()
{
    for(int kind=0;kind<3;++kind) {
        Vector10 values=Vector10::Zero();
        Eigen::Vector3d expected,normal;
        if(kind<2) {
            values<<.125,.5,-.25,1.,1.,kind==0?2.:-2.,0.,0.,0.,0.;
            expected=Eigen::Vector3d(.25,-.125,kind==0?-.375:-.3125);
            normal=Eigen::Vector3d(1.,kind==0?-.75:.25,1.);
        } else {
            values<<.125,.25,.5,1.,0.,0.,0.,.25,.5,0.;
            expected=Eigen::Vector3d(0.,0.,-.125);
            normal=Eigen::Vector3d(.125,.4375,1.);
        }
        const auto q=Quadric::fromUnnormalized(values);
        const spin::QuadricInverseGaussMap map(q);
        const auto result=map.paraboloidCoefficientDifferential(normal);
        if(!result.valid())return fail("regular paraboloid full point derivative missing: "+result.diagnostic);
        if(!close(*result.point.point,expected,1e-13))return fail("paraboloid analytic point mismatch");
        if(((*result.coefficientJacobian)*q.values()).norm()>1e-12)
            return fail("paraboloid normalization radial derivative nonzero");
        double firstError=0.,lastError=0.;
        for(double step: {1e-4,3e-5,1e-5}) {
            Jacobian difference;
            for(int column=0;column<10;++column) {
                auto plus=q.values(),minus=q.values();plus(column)+=step;minus(column)-=step;
                difference.col(column)=(implicitPoint(Quadric::fromUnnormalized(plus),normal,expected)-
                    implicitPoint(Quadric::fromUnnormalized(minus),normal,expected))/(2*step);
                if((difference.col(column)-result.coefficientJacobian->col(column)).norm()/
                    std::max(1.,difference.col(column).norm())>1e-6)
                    return fail("one paraboloid coefficient direction exceeds derivative gate");
            }
            const double error=(difference-*result.coefficientJacobian).norm()/std::max(1.,difference.norm());
            if(error>1e-6)return fail("paraboloid full-coefficient finite difference mismatch");
            if(step==1e-4)firstError=error;
            lastError=error;
            std::cout<<"paraboloid="<<kind<<" step="<<step<<" full_direction_error="<<error<<'\n';
        }
        if(lastError>firstError*.2)return fail("paraboloid three-step error lacks convergence");
        if(map.point(normal).valid())return fail("legacy full-rank entry silently changed");
        if(map.paraboloidCoefficientDifferential(-normal).valid())return fail("opposite normal accepted");
        if(map.paraboloidCoefficientDifferential(Eigen::Vector3d::Zero()).valid())return fail("zero normal accepted");
    }
    Vector10 cylinder=Vector10::Zero();cylinder(0)=-1.;cylinder(4)=1.;cylinder(5)=2.;
    if(spin::QuadricInverseGaussMap(Quadric::fromUnnormalized(cylinder))
        .paraboloidCoefficientDifferential(Eigen::Vector3d(1.,0.,1.)).valid())
        return fail("cylinder promoted to regular paraboloid");
    cylinder(3)=1.;cylinder(5)=0.;
    if(spin::QuadricInverseGaussMap(Quadric::fromUnnormalized(cylinder))
        .paraboloidCoefficientDifferential(Eigen::Vector3d(1.,0.,1.)).valid())
        return fail("rank-one matrix promoted to regular paraboloid");
    cylinder(5)=2.;
    if(spin::QuadricInverseGaussMap(Quadric::fromUnnormalized(cylinder))
        .paraboloidCoefficientDifferential(Eigen::Vector3d::UnitX()).valid())
        return fail("ideal normal promoted to finite point");
    cylinder(6)=3.;
    if(spin::QuadricInverseGaussMap(Quadric::fromUnnormalized(cylinder))
        .paraboloidCoefficientDifferential(Eigen::Vector3d(1.,0.,1.)).valid())
        return fail("full rank matrix promoted to paraboloid");
    return EXIT_SUCCESS;
}

} // namespace


int main()
{
    if(compactRationalPointRegression()!=EXIT_SUCCESS)return EXIT_FAILURE;
    if(rationalFinitePointRegression()!=EXIT_SUCCESS)return EXIT_FAILURE;
    if(algebraicTransitionRegression()!=EXIT_SUCCESS)return EXIT_FAILURE;
    if(paraboloidFiniteMapRegression()!=EXIT_SUCCESS)return EXIT_FAILURE;
    Vector10 coefficients = Vector10::Zero();
    coefficients(0) = -4.0;
    coefficients(4) = 1.0;
    coefficients(5) = 1.0;
    coefficients(6) = 1.0;
    const Quadric sphere = Quadric::fromUnnormalized(coefficients);
    const spin::QuadricInverseGaussMap sphereMap(sphere);
    const Eigen::Vector3d sphereNormal(2.0, -3.0, 6.0);
    const Eigen::Vector3d sphereUnitNormal = sphereNormal.normalized();
    const spin::QuadricInverseGaussPointResult spherePoint =
        sphereMap.point(sphereNormal);
    const spin::QuadricInverseGaussDifferentialResult sphereDifferential =
        sphereMap.coefficientDifferential(sphereNormal);
    const double sphereScale = 1.0 / std::sqrt(19.0);
    const Jacobian expectedSphereDerivative =
        sphereOracle(sphereScale, 2.0, sphereUnitNormal);
    if (
        !spherePoint.valid() ||
        !close(*spherePoint.point, 2.0 * sphereUnitNormal, 2e-14) ||
        !sphereDifferential.point.valid() ||
        !sphereDifferential.valid() ||
        (*sphereDifferential.coefficientJacobian -
            expectedSphereDerivative).norm() > 2e-12 ||
        (*sphereDifferential.coefficientJacobian * sphere.values()).norm() >
            2e-13
    )
    {
        return fail(
            "repeated-eigenvalue sphere point, literal derivatives, or "
            "normalized radial-null contract is incorrect"
        );
    }

    const Eigen::Vector3d axisNormal = Eigen::Vector3d::UnitX();
    const spin::QuadricInverseGaussDifferentialResult axisSphereDifferential =
        sphereMap.coefficientDifferential(axisNormal);
    if (!axisSphereDifferential.valid())
    {
        return fail("axis-aligned sphere differential is undefined");
    }
    const Jacobian axisSphereFiniteDifference = centralDifference(
        sphere.values(),
        axisNormal,
        1e-4
    );
    const Jacobian axisSphereOracle =
        sphereOracle(sphereScale, 2.0, axisNormal);
    int analyticZeroEntryCount = 0;
    double analyticZeroAbsoluteError = 0.0;
    for (int row = 0; row < axisSphereOracle.rows(); ++row)
    {
        for (int column = 0; column < axisSphereOracle.cols(); ++column)
        {
            if (axisSphereOracle(row, column) == 0.0)
            {
                ++analyticZeroEntryCount;
                analyticZeroAbsoluteError = std::max(
                    analyticZeroAbsoluteError,
                    std::abs(
                        (*axisSphereDifferential.coefficientJacobian)(
                            row,
                            column
                        ) - axisSphereFiniteDifference(row, column)
                    )
                );
            }
        }
    }
    if (
        analyticZeroEntryCount == 0 || analyticZeroAbsoluteError > 1e-8
    )
    {
        return fail("explicit analytic-zero sphere entries violate absolute rule");
    }

    const spin::QuadricInverseGaussPointResult scaledNormalPoint =
        sphereMap.point(3.25 * sphereNormal);
    const spin::QuadricInverseGaussDifferentialResult scaledNormalDerivative =
        sphereMap.coefficientDifferential(3.25 * sphereNormal);
    if (
        !scaledNormalPoint.valid() || !scaledNormalDerivative.valid() ||
        !close(*scaledNormalPoint.point, *spherePoint.point, 2e-14) ||
        (*scaledNormalDerivative.coefficientJacobian -
            *sphereDifferential.coefficientJacobian).norm() > 2e-12
    )
    {
        return fail("positive normal rescaling changed the finite map");
    }

    const spin::QuadricInverseGaussPointResult oppositeSpherePoint =
        sphereMap.point(-sphereNormal);
    const spin::QuadricInverseGaussPointResult complementedPoint =
        spin::QuadricInverseGaussMap(sphere.inverted()).point(sphereNormal);
    if (
        !oppositeSpherePoint.valid() || !complementedPoint.valid() ||
        !close(
            *complementedPoint.point,
            *oppositeSpherePoint.point,
            2e-14
        ) ||
        close(*complementedPoint.point, *spherePoint.point, 1e-6)
    )
    {
        return fail("normal sign or complemented-quadric semantics changed");
    }

    Eigen::Matrix3d asymmetricA;
    asymmetricA <<
        1.1, 0.17, -0.08,
        0.17, 0.8, 0.12,
        -0.08, 0.12, 1.3;
    const Eigen::Vector3d center(0.2, -0.1, 0.3);
    coefficients.setZero();
    coefficients(0) = center.dot(asymmetricA * center) - 1.0;
    coefficients.segment<3>(1) = -2.0 * asymmetricA * center;
    coefficients(4) = asymmetricA(0, 0);
    coefficients(5) = asymmetricA(1, 1);
    coefficients(6) = asymmetricA(2, 2);
    coefficients(7) = asymmetricA(1, 2);
    coefficients(8) = asymmetricA(0, 2);
    coefficients(9) = asymmetricA(0, 1);
    const Quadric asymmetric = Quadric::fromUnnormalized(coefficients);
    const Eigen::Vector3d asymmetricNormal(0.4, -0.2, 0.7);
    const spin::QuadricInverseGaussDifferentialResult asymmetricDerivative =
        spin::QuadricInverseGaussMap(asymmetric).coefficientDifferential(
            asymmetricNormal
        );
    if (!asymmetricDerivative.valid())
    {
        return fail("asymmetric rotated/translated ellipsoid is undefined");
    }

    const std::array<double, 3> steps = {{3e-4, 1e-4, 3e-5}};
    std::array<double, 3> relativeErrors{};
    for (std::size_t index = 0; index < steps.size(); ++index)
    {
        const Jacobian finiteDifference = centralDifference(
            asymmetric.values(),
            asymmetricNormal,
            steps[index]
        );
        relativeErrors[index] =
            (finiteDifference -
                *asymmetricDerivative.coefficientJacobian).norm() /
            std::max(1.0, finiteDifference.norm());
        if (relativeErrors[index] > 1e-6)
        {
            return fail("three-step normalized finite difference exceeded 1e-6");
        }
    }
    if (
        !(relativeErrors[1] < relativeErrors[0] &&
          relativeErrors[2] < relativeErrors[1])
    )
    {
        return fail("central finite differences do not show the expected trend");
    }

    coefficients.setZero();
    coefficients(0) = -1.0;
    coefficients(4) = 1.0;
    coefficients(5) = 1.0;
    coefficients(6) = -1.0;
    const Quadric hyperboloid = Quadric::fromUnnormalized(coefficients);
    const spin::QuadricInverseGaussMap hyperboloidMap(hyperboloid);
    const spin::QuadricInverseGaussDifferentialResult hyperboloidFinite =
        hyperboloidMap.coefficientDifferential(Eigen::Vector3d(1.0, 0.2, 0.1));
    if (!hyperboloidFinite.valid())
    {
        return fail("admissible hyperboloid differential is undefined");
    }
    std::array<double, 3> hyperboloidRelativeErrors{};
    for (std::size_t index = 0; index < steps.size(); ++index)
    {
        const Jacobian finiteDifference = centralDifference(
            hyperboloid.values(),
            Eigen::Vector3d(1.0, 0.2, 0.1),
            steps[index]
        );
        hyperboloidRelativeErrors[index] =
            (finiteDifference -
                *hyperboloidFinite.coefficientJacobian).norm() /
            std::max(1.0, finiteDifference.norm());
        if (hyperboloidRelativeErrors[index] > 1e-6)
        {
            return fail("hyperboloid finite difference exceeded 1e-6");
        }
    }
    if (
        !(hyperboloidRelativeErrors[1] < hyperboloidRelativeErrors[0] &&
          hyperboloidRelativeErrors[2] < hyperboloidRelativeErrors[1])
    )
    {
        return fail("hyperboloid central finite differences lack convergence");
    }
    const spin::QuadricInverseGaussDifferentialResult hyperboloidIdeal =
        hyperboloidMap.coefficientDifferential(Eigen::Vector3d(1.0, 0.0, 1.0));
    const spin::QuadricInverseGaussPointResult hyperboloidInadmissible =
        hyperboloidMap.point(Eigen::Vector3d(0.0, 0.0, 1.0));
    if (
        !hyperboloidFinite.valid() ||
        hyperboloidIdeal.point.status !=
            spin::QuadricInverseGaussPointStatus::IdealOrInadmissibleNormal ||
        hyperboloidIdeal.status !=
            spin::QuadricInverseGaussDifferentialStatus::PointUndefined ||
        hyperboloidIdeal.coefficientJacobian.has_value() ||
        hyperboloidInadmissible.status !=
            spin::QuadricInverseGaussPointStatus::IdealOrInadmissibleNormal
    )
    {
        return fail("hyperboloid finite/ideal normal states are incorrect");
    }

    const Eigen::Vector3d nearIdealNormal(
        std::sqrt(1.0 + 1e-8),
        0.0,
        1.0
    );
    const spin::QuadricInverseGaussDifferentialResult nearIdeal =
        hyperboloidMap.coefficientDifferential(nearIdealNormal);
    if (!nearIdeal.point.valid() || !nearIdeal.valid() ||
        !nearIdeal.coefficientJacobian->allFinite())
    {
        return fail("near-ideal differential was not certified after reversible scaling: " +
                    nearIdeal.diagnostic);
    }
    const Jacobian nearIdealFiniteDifference = centralDifference(
        hyperboloid.values(),
        nearIdealNormal,
        1e-12
    );
    if ((
            *nearIdeal.coefficientJacobian - nearIdealFiniteDifference
        ).norm() / nearIdealFiniteDifference.norm() > 1e-3)
    {
        return fail(
            "near-ideal equilibrated differential disagrees with its independent "
            "fixed-normal finite difference");
    }

    Vector10 capturedNearConeCoefficients;
    capturedNearConeCoefficients <<
        -1.1089457960913767e-10, 3.1374659509406226e-08,
        -1.2306930695891306e-07, -1.4403414898221336e-07,
        .40315348093351583, .40315365134925046, -.8215439149561166,
        -3.6688686058282268e-06, 1.8991362047213163e-07,
        -4.039908646624158e-09;
    const Quadric capturedNearCone =
        Quadric::fromUnnormalized(capturedNearConeCoefficients);
    const Eigen::Vector3d capturedNormal(
        -.50599266182565983, -.27144044986507432,
        .81871332489194881);
    const spin::QuadricInverseGaussDifferentialResult capturedDifferential =
        spin::QuadricInverseGaussMap(capturedNearCone)
            .coefficientDifferential(capturedNormal);
    if (!capturedDifferential.valid())
    {
        return fail(
            "captured near-cone inverse-Gauss derivative is rejected despite "
            "its independently verified backward-stable solution: " +
            capturedDifferential.diagnostic);
    }
    Jacobian capturedFiniteDifference;
    capturedFiniteDifference <<
        1354468.63144981909, -1.29292687567228502, .206735839521987428,
        -.118734812184552762, -.121495553475417915, -.0351793834629522717,
        -.0770670408749348725, -.104137902197929767, -.193529090635252121,
        -.130754911490996839,
        726616.296273806559, -.0282737500528720345, -1.12931671288042782,
        -.0636906143651330861, -.0655771327438884959, -.0184729093933767039,
        -.0413432767542367634, -.0552737841236514447, -.104137864415949269,
        -.0696135760797206767,
        1075462.31916096425, -.0418481499934221388, .164155976377306653,
        .514333731908959906, -.0970604919118753962, -.0279328017403069501,
        -.0614824948883612012, -.0828825734335558065, -.154499815724463886,
        -.104137882061271916;
    if (
        (*capturedDifferential.coefficientJacobian -
            capturedFiniteDifference).norm() /
            capturedFiniteDifference.norm() > 1e-6
    )
    {
        return fail(
            "captured near-cone derivative differs from the independent "
            "80-digit central-difference reference");
    }

    Vector10 cherchiNewtonNearConeCoefficients;
    cherchiNewtonNearConeCoefficients <<
        -1.127002572903594e-20, -2.300216466123981e-14,
        1.3354520904472402e-15, -1.8028129010386688e-15,
        0.40315369019808345, 0.40315352253990416,
        -0.8215438754828222, -6.524260075274312e-13,
        1.372198932213066e-14, -9.629360740025284e-17;
    const Quadric cherchiNewtonNearCone =
        Quadric::fromUnnormalized(cherchiNewtonNearConeCoefficients);
    const spin::QuadricInverseGaussMap cherchiNewtonNearConeMap(
        cherchiNewtonNearCone
    );
    const spin::QuadricInverseGaussPointResult cherchiNewtonNearConePoint =
        cherchiNewtonNearConeMap.point(Eigen::Vector3d::UnitX());
    if (
        !cherchiNewtonNearConePoint.valid() ||
        !cherchiNewtonNearConePoint.point->allFinite()
    )
    {
        return fail(
            "current-source Cherchi Newton hyperboloid was classified as "
            "degenerate before evaluating a finite inverse-Gauss point: " +
            cherchiNewtonNearConePoint.diagnostic);
    }

    Vector10 secondCherchiNewtonNearConeCoefficients;
    secondCherchiNewtonNearConeCoefficients <<
        -2.8364972803646595e-21, -8.11768576886572e-15,
        3.615646720485601e-16, -4.875767860908935e-16,
        0.4031536902205084, 0.4031535225168989,
        -0.821543875483107, -2.3102175102938004e-13,
        4.8400302740371675e-15, -2.4278465740028797e-17;
    const spin::QuadricInverseGaussPointResult secondCherchiNearConePoint =
        spin::QuadricInverseGaussMap(
            Quadric::fromUnnormalized(
                secondCherchiNewtonNearConeCoefficients
            )
        ).point(Eigen::Vector3d::UnitX());
    if (
        !secondCherchiNearConePoint.valid() ||
        !secondCherchiNearConePoint.point->allFinite()
    )
    {
        return fail(
            "second current-source Cherchi Newton hyperboloid was rejected "
            "by an absolute centered-level cutoff: " +
            secondCherchiNearConePoint.diagnostic);
    }

    Vector10 illConditionedCherchiCoefficients;
    illConditionedCherchiCoefficients <<
        -1.127002572903594e-20, -2.300216466123981e-14,
        1.3354520904472402e-15, -1.8028129010386688e-15,
        0.40315369019808345, 0.40315352253990416,
        -0.8215438754828222, -6.524260075274312e-13,
        1.372198932213066e-14, -9.629360740025284e-17;
    const Eigen::Vector3d illConditionedCherchiNormal(
        -0.52573111211913359, 0.85065080835203999, 0.0);
    const spin::QuadricInverseGaussDifferentialResult
        illConditionedCherchiDifferential =
            spin::QuadricInverseGaussMap(
                Quadric::fromUnnormalized(illConditionedCherchiCoefficients)
            ).coefficientDifferential(illConditionedCherchiNormal);
    if (
        !illConditionedCherchiDifferential.valid() ||
        !illConditionedCherchiDifferential.coefficientJacobian->allFinite()
    )
    {
        return fail(
            "captured Cherchi near-cone fixed-normal differential must pass "
            "the unchanged 1e-10 residual gate with a certified solve: " +
            illConditionedCherchiDifferential.diagnostic);
    }
    const Jacobian illConditionedCherchiFiniteDifference = centralDifference(
        Quadric::fromUnnormalized(illConditionedCherchiCoefficients).values(),
        illConditionedCherchiNormal,
        1e-22
    );
    if (
        (*illConditionedCherchiDifferential.coefficientJacobian -
            illConditionedCherchiFiniteDifference).norm() /
            illConditionedCherchiFiniteDifference.norm() > 5e-3
    )
    {
        return fail(
            "exact-rational captured near-cone differential disagrees with "
            "its independent fixed-normal finite difference");
    }

    Vector10 scaleSensitiveCherchiCoefficients;
    scaleSensitiveCherchiCoefficients <<
        6.105054942668466e-26, -1.3858848384049204e-17,
        -1.0991071633590135e-17, -2.9048607649655397e-18,
        0.39824375457412803, 0.3982438519335089,
        -0.8263193972913035, 5.167709185408485e-18,
        -1.9731758134344458e-17, 6.320903791626323e-19;
    const Eigen::Vector3d scaleSensitiveCherchiNormal(
        0.23722869, -0.50420944, -0.83035859);
    const spin::QuadricInverseGaussDifferentialResult
        scaleSensitiveCherchiDifferential = spin::QuadricInverseGaussMap(
            Quadric::fromUnnormalized(scaleSensitiveCherchiCoefficients)
        ).coefficientDifferential(scaleSensitiveCherchiNormal);
    if (!scaleSensitiveCherchiDifferential.valid() ||
        !scaleSensitiveCherchiDifferential.coefficientJacobian->allFinite() ||
        !(scaleSensitiveCherchiDifferential.rawReciprocalCondition < 1e-12) ||
        !(scaleSensitiveCherchiDifferential.equilibratedReciprocalCondition >
            1e-12))
    {
        return fail(
            "captured near-cone fixed-normal system that becomes well-conditioned "
            "under reversible row/column scaling must retain its exact derivative: " +
            scaleSensitiveCherchiDifferential.diagnostic +
            "; raw_rcond=" + std::to_string(
                scaleSensitiveCherchiDifferential.rawReciprocalCondition) +
            "; equilibrated_rcond=" + std::to_string(
            scaleSensitiveCherchiDifferential.equilibratedReciprocalCondition));
    }
    const Jacobian scaleSensitiveFiniteDifference = centralDifference(
        Quadric::fromUnnormalized(scaleSensitiveCherchiCoefficients).values(),
        scaleSensitiveCherchiNormal,
        1e-26
    );
    if ((
            *scaleSensitiveCherchiDifferential.coefficientJacobian -
            scaleSensitiveFiniteDifference
        ).norm() / std::max(1.0, scaleSensitiveFiniteDifference.norm()) > 5e-3)
    {
        return fail(
            "captured scale-sensitive differential disagrees with its independent "
            "fixed-normal finite difference");
    }

    Vector10 rankScaledCoefficients = Vector10::Zero();
    rankScaledCoefficients << -1e-30, 0, 0, 0, 1, 1, -1, 0, 0, 0;
    const spin::QuadricInverseGaussDifferentialResult rankScaledDifferential =
        spin::QuadricInverseGaussMap(
            Quadric::fromUnnormalized(rankScaledCoefficients)
        ).coefficientDifferential(Eigen::Vector3d::UnitX());
    if (!rankScaledDifferential.valid() ||
        !(rankScaledDifferential.rawReciprocalCondition < 1e-12) ||
        !(rankScaledDifferential.equilibratedReciprocalCondition > 1e-12))
    {
        return fail(
            "near-cone differential must distinguish an exactly full-rank scaled "
            "system from the raw-scale rank decision: " +
            rankScaledDifferential.diagnostic + "; raw_rcond=" +
            std::to_string(rankScaledDifferential.rawReciprocalCondition) +
            "; equilibrated_rcond=" + std::to_string(
                rankScaledDifferential.equilibratedReciprocalCondition) +
            "; valid=" + (rankScaledDifferential.valid() ? "true" : "false"));
    }

    const spin::QuadricInverseGaussPointResult zeroNormal =
        sphereMap.point(Eigen::Vector3d::Zero());
    const spin::QuadricInverseGaussPointResult nonfiniteNormal =
        sphereMap.point(Eigen::Vector3d(
            std::numeric_limits<double>::quiet_NaN(),
            0.0,
            1.0
        ));
    const spin::QuadricInverseGaussPointResult infiniteNormal =
        sphereMap.point(Eigen::Vector3d(
            std::numeric_limits<double>::infinity(),
            0.0,
            1.0
        ));
    const spin::QuadricInverseGaussPointResult belowThresholdNormal =
        sphereMap.point(1e-8 * sphereUnitNormal);
    if (
        zeroNormal.status != spin::QuadricInverseGaussPointStatus::InvalidNormal ||
        zeroNormal.point.has_value() ||
        nonfiniteNormal.status !=
            spin::QuadricInverseGaussPointStatus::InvalidNormal ||
        nonfiniteNormal.point.has_value() ||
        infiniteNormal.status !=
            spin::QuadricInverseGaussPointStatus::InvalidNormal ||
        infiniteNormal.point.has_value() ||
        belowThresholdNormal.status !=
            spin::QuadricInverseGaussPointStatus::BelowNumericalThreshold ||
        belowThresholdNormal.point.has_value()
    )
    {
        return fail("invalid or magnitude-threshold normal contract is implicit");
    }

    coefficients.setZero();
    coefficients(4) = 1.0;
    coefficients(5) = 1.0;
    coefficients(6) = 1.0;
    const spin::QuadricInverseGaussPointResult degenerateScale =
        spin::QuadricInverseGaussMap(
            Quadric::fromUnnormalized(coefficients)
        ).point(Eigen::Vector3d::UnitX());
    coefficients(6) = 0.0;
    const spin::QuadricInverseGaussMap singularMap(
        Quadric::fromUnnormalized(coefficients)
    );
    const spin::QuadricInverseGaussPointResult singularQuadratic =
        singularMap.point(Eigen::Vector3d::UnitX());
    const spin::QuadricInverseGaussDifferentialResult singularDifferential =
        singularMap.coefficientDifferential(Eigen::Vector3d::UnitX());
    if (
        degenerateScale.status !=
            spin::QuadricInverseGaussPointStatus::DegenerateQuadricScale ||
        degenerateScale.diagnostic.empty() ||
        singularQuadratic.status !=
            spin::QuadricInverseGaussPointStatus::SingularQuadraticMatrix ||
        singularQuadratic.diagnostic.empty() ||
        singularDifferential.point.valid() || singularDifferential.valid() ||
        singularDifferential.status !=
            spin::QuadricInverseGaussDifferentialStatus::PointUndefined ||
        singularDifferential.coefficientJacobian.has_value() ||
        singularDifferential.diagnostic.empty()
    )
    {
        return fail("degenerate finite coefficient states are not explicit");
    }

    bool rejectedNonfiniteCoefficients = false;
    bool rejectedZeroCoefficients = false;
    try
    {
        coefficients.setConstant(
            std::numeric_limits<double>::infinity()
        );
        (void)Quadric::fromUnnormalized(coefficients);
    }
    catch (const std::invalid_argument&)
    {
        rejectedNonfiniteCoefficients = true;
    }
    try
    {
        (void)Quadric::fromUnnormalized(Vector10::Zero());
    }
    catch (const std::invalid_argument&)
    {
        rejectedZeroCoefficients = true;
    }
    if (!rejectedNonfiniteCoefficients || !rejectedZeroCoefficients)
    {
        return fail("existing coefficient validation was bypassed");
    }

    std::cout << "[PASS] scenario=quadric-inverse-gauss-map"
              << " fd_h3e-4=" << relativeErrors[0]
              << " fd_h1e-4=" << relativeErrors[1]
              << " fd_h3e-5=" << relativeErrors[2]
              << " hyperboloid_fd_h3e-4=" << hyperboloidRelativeErrors[0]
              << " hyperboloid_fd_h1e-4=" << hyperboloidRelativeErrors[1]
              << " hyperboloid_fd_h3e-5=" << hyperboloidRelativeErrors[2]
              << " analytic_zero_entries=" << analyticZeroEntryCount
              << " analytic_zero_abs=" << analyticZeroAbsoluteError
              << " radial_null="
              << (*sphereDifferential.coefficientJacobian * sphere.values()).norm()
              << '\n';
    return EXIT_SUCCESS;
}
