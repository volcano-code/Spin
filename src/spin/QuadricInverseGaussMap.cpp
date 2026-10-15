#include "spin/QuadricInverseGaussMap.h"
#include "spin/RationalSourceRepresentative.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <CGAL/Gmpq.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/number_utils.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <locale>
#include <limits>
#include <sstream>


namespace
{


constexpr double minimumNormalNorm = 1e-15;
constexpr double minimumWSquared = 1e-14;
constexpr double minimumBorderedReciprocalCondition = 1e-12;

// Sqrt_extension::to_double evaluates a0+a1*sqrt(root) after rounding the
// pieces separately. Cancellation in ideal-edge roots can destroy the FD
// accuracy of the actual map. Round the exact algebraic value instead.
double representedAlgebraic(
    const CGAL::Sqrt_extension<CGAL::Gmpq,CGAL::Gmpq>& value)
{
    using Rational=CGAL::Gmpq;
    using Algebraic=CGAL::Sqrt_extension<Rational,Rational>;
    static_assert(sizeof(double)==sizeof(std::uint64_t) &&
                  std::numeric_limits<double>::is_iec559 &&
                  std::numeric_limits<double>::digits==53,
                  "ideal-edge representation requires IEEE binary64");
    if(value==Algebraic(0))return 0;
    const bool negative=value<Algebraic(0);
    const Algebraic magnitude=negative ? -value : value;
    const auto decode=[](std::uint64_t bits) {
        double result;std::memcpy(&result,&bits,sizeof(result));return result;
    };
    std::uint64_t low=0,high=UINT64_C(0x7fefffffffffffff);
    if(magnitude>Algebraic(Rational(decode(high))))
        return std::numeric_limits<double>::quiet_NaN();
    // Positive finite binary64 bit patterns are ordered numerically.
    while(low<high) {
        const auto middle=low+(high-low+1)/2;
        if(Algebraic(Rational(decode(middle)))<=magnitude)low=middle;
        else high=middle-1;
    }
    const Rational lower(decode(low));
    if(magnitude!=Algebraic(lower)) {
        const Rational midpoint=(lower+Rational(decode(low+1)))/2;
        if(magnitude>Algebraic(midpoint) ||
           (magnitude==Algebraic(midpoint) && (low&1)))++low;
    }
    return negative ? -decode(low) : decode(low);
}


bool exactInverseGaussDenominator(
    const spin::QuadricCoefficients10& quadric,
    double& denominator
)
{
    using Rational = CGAL::Gmpq;
    const auto& coefficients = quadric.values();
    const Rational a00(coefficients(4));
    const Rational a11(coefficients(5));
    const Rational a22(coefficients(6));
    const Rational a12(coefficients(7));
    const Rational a02(coefficients(8));
    const Rational a01(coefficients(9));
    const Rational determinant =
        a00 * (a11 * a22 - a12 * a12) -
        a01 * (a01 * a22 - a12 * a02) +
        a02 * (a01 * a12 - a11 * a02);
    if (determinant == 0)
    {
        return false;
    }

    const Rational b0(coefficients(1));
    const Rational b1(coefficients(2));
    const Rational b2(coefficients(3));
    const Rational adj00 = a11 * a22 - a12 * a12;
    const Rational adj01 = a02 * a12 - a01 * a22;
    const Rational adj02 = a01 * a12 - a02 * a11;
    const Rational adj11 = a00 * a22 - a02 * a02;
    const Rational adj12 = a01 * a02 - a00 * a12;
    const Rational adj22 = a00 * a11 - a01 * a01;
    const Rational quadraticForm =
        b0 * (adj00 * b0 + adj01 * b1 + adj02 * b2) +
        b1 * (adj01 * b0 + adj11 * b1 + adj12 * b2) +
        b2 * (adj02 * b0 + adj12 * b1 + adj22 * b2);
    const Rational centeredLevelNumerator =
        4 * Rational(coefficients(0)) * determinant - quadraticForm;
    denominator = CGAL::to_double(-centeredLevelNumerator / determinant);
    return true;
}


bool solveExactBorderedSystem(
    const Eigen::Matrix4d& system,
    const Eigen::Matrix<double, 4, 10>& rightHandSide,
    Eigen::Matrix<double, 4, 10>& solution
)
{
    using Rational = CGAL::Gmpq;
    std::array<std::array<Rational, 4>, 4> matrix{};
    std::array<std::array<Rational, 10>, 4> rhs{};
    for (int row = 0; row < 4; ++row)
    {
        for (int column = 0; column < 4; ++column)
        {
            matrix[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(column)] =
                Rational(system(row, column));
        }
        for (int column = 0; column < 10; ++column)
        {
            rhs[static_cast<std::size_t>(row)]
               [static_cast<std::size_t>(column)] =
                Rational(rightHandSide(row, column));
        }
    }
    const auto magnitude = [](const Rational& value) {
        return value < 0 ? -value : value;
    };

    for (int pivotColumn = 0; pivotColumn < 4; ++pivotColumn)
    {
        int pivotRow = pivotColumn;
        for (int row = pivotColumn + 1; row < 4; ++row)
        {
            if (magnitude(matrix[static_cast<std::size_t>(row)]
                                [static_cast<std::size_t>(pivotColumn)]) >
                magnitude(matrix[static_cast<std::size_t>(pivotRow)]
                                [static_cast<std::size_t>(pivotColumn)]))
            {
                pivotRow = row;
            }
        }
        if (matrix[static_cast<std::size_t>(pivotRow)]
                  [static_cast<std::size_t>(pivotColumn)] == 0)
        {
            return false;
        }
        if (pivotRow != pivotColumn)
        {
            std::swap(
                matrix[static_cast<std::size_t>(pivotRow)],
                matrix[static_cast<std::size_t>(pivotColumn)]
            );
            std::swap(
                rhs[static_cast<std::size_t>(pivotRow)],
                rhs[static_cast<std::size_t>(pivotColumn)]
            );
        }
        const Rational pivot = matrix[static_cast<std::size_t>(pivotColumn)]
                                     [static_cast<std::size_t>(pivotColumn)];
        for (int row = pivotColumn + 1; row < 4; ++row)
        {
            const Rational factor =
                matrix[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(pivotColumn)] / pivot;
            for (int column = pivotColumn + 1; column < 4; ++column)
            {
                matrix[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(column)] -=
                    factor * matrix[static_cast<std::size_t>(pivotColumn)]
                                [static_cast<std::size_t>(column)];
            }
            matrix[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(pivotColumn)] = Rational(0);
            for (int column = 0; column < 10; ++column)
            {
                rhs[static_cast<std::size_t>(row)]
                   [static_cast<std::size_t>(column)] -=
                    factor * rhs[static_cast<std::size_t>(pivotColumn)]
                                 [static_cast<std::size_t>(column)];
            }
        }
    }

    std::array<std::array<Rational, 10>, 4> exactSolution{};
    for (int rightColumn = 0; rightColumn < 10; ++rightColumn)
    {
        for (int row = 3; row >= 0; --row)
        {
            Rational value = rhs[static_cast<std::size_t>(row)]
                                [static_cast<std::size_t>(rightColumn)];
            for (int column = row + 1; column < 4; ++column)
            {
                value -= matrix[static_cast<std::size_t>(row)]
                                [static_cast<std::size_t>(column)] *
                    exactSolution[static_cast<std::size_t>(column)]
                                 [static_cast<std::size_t>(rightColumn)];
            }
            exactSolution[static_cast<std::size_t>(row)]
                         [static_cast<std::size_t>(rightColumn)] =
                value / matrix[static_cast<std::size_t>(row)]
                              [static_cast<std::size_t>(row)];
        }
    }

    for (int row = 0; row < 4; ++row)
    {
        for (int column = 0; column < 10; ++column)
        {
            Rational residual =
                -Rational(rightHandSide(row, column));
            for (int inner = 0; inner < 4; ++inner)
            {
                residual += Rational(system(row, inner)) *
                    exactSolution[static_cast<std::size_t>(inner)]
                                 [static_cast<std::size_t>(column)];
            }
            if (residual != 0)
            {
                return false;
            }
            solution(row, column) = CGAL::to_double(
                exactSolution[static_cast<std::size_t>(row)]
                            [static_cast<std::size_t>(column)]
            );
        }
    }
    return solution.allFinite();
}


spin::QuadricInverseGaussPointResult regularParaboloidPoint(
    const spin::QuadricCoefficients10& q, const Eigen::Vector3d& normal)
{
    using R=CGAL::Gmpq;
    using V=std::array<R,3>;
    using M=std::array<V,3>;
    using Status=spin::QuadricInverseGaussPointStatus;
    spin::QuadricInverseGaussPointResult result;
    const auto fail=[&](Status status,const char* message) {
        result.status=status;result.diagnostic=message;return result;
    };
    if(!normal.allFinite() || !std::isfinite(normal.norm()) || normal.norm()<=minimumNormalNorm)
        return fail(Status::InvalidNormal,"regular paraboloid requires a finite nonzero normal");
    const auto adj=[](const M& a) {
        M out{};
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            out[i][j]=a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]-
                      a[(j+1)%3][(i+2)%3]*a[(j+2)%3][(i+1)%3];
        return out;
    };
    const auto dot=[](const V& a,const V& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];};
    M a{};V b{},n{};
    const Eigen::Matrix3d matrix=q.quadraticMatrix();
    for(int i=0;i<3;++i) {
        b[i]=R(q.linearVector()(i));n[i]=R(normal(i));
        for(int j=0;j<3;++j)a[i][j]=R(matrix(i,j));
    }
    const auto cofactor=adj(a);
    R determinant=0;for(int i=0;i<3;++i)determinant+=a[0][i]*cofactor[i][0];
    if(determinant!=0)return fail(Status::SingularQuadraticMatrix,"explicit paraboloid path requires exact rank two");
    V u{};bool found=false;
    for(int column=0;column<3 && !found;++column) {
        for(int row=0;row<3;++row)u[row]=cofactor[row][column];
        found=dot(u,u)!=0;
    }
    if(!found)return fail(Status::SingularQuadraticMatrix,"quadric rank is less than two");
    const R ub=dot(u,b),un=dot(u,n);
    if(ub==0)return fail(Status::DegenerateQuadricScale,"rank-two quadric has no linear null-axis term");
    if(un==0)return fail(Status::IdealOrInadmissibleNormal,"paraboloid normal is ideal");
    const R t=ub/(2*un);
    if(t<=0)return fail(Status::IdealOrInadmissibleNormal,"paraboloid normal has opposite gradient orientation");
    M bordered=a;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)bordered[i][j]+=u[i]*u[j];
    const auto inverseNumerator=adj(bordered);
    R denominator=0;for(int i=0;i<3;++i)denominator+=bordered[0][i]*inverseNumerator[i][0];
    if(denominator==0)return fail(Status::SingularQuadraticMatrix,"regular paraboloid nullspace solve is singular");
    V rhs{},x{};
    for(int i=0;i<3;++i)rhs[i]=t*n[i]-b[i]/2;
    for(int i=0;i<3;++i)x[i]=dot(inverseNumerator[i],rhs)/denominator;
    const auto polynomial=[&](const V& point) {
        R value=R(q.constant())+dot(b,point);
        for(int i=0;i<3;++i)value+=point[i]*dot(a[i],point);
        return value;
    };
    const R shift=-polynomial(x)/ub;
    for(int i=0;i<3;++i)x[i]+=shift*u[i];
    if(polynomial(x)!=0)return fail(Status::NonFinitePoint,"exact paraboloid point residual is nonzero");
    for(int i=0;i<3;++i)
        if(2*dot(a[i],x)+b[i]!=2*t*n[i])
            return fail(Status::NonFinitePoint,"exact paraboloid normal residual is nonzero");
    const Eigen::Vector3d represented(CGAL::to_double(x[0]),CGAL::to_double(x[1]),CGAL::to_double(x[2]));
    if(!represented.allFinite())return fail(Status::NonFinitePoint,"paraboloid point overflows binary64");
    result.status=Status::Valid;result.point=represented;
    result.diagnostic="exact represented-coefficient regular paraboloid finite point; rounded to binary64";
    return result;
}

spin::QuadricInverseGaussPointResult algebraicRegularPoint(
    const spin::QuadricCoefficients10& q, const Eigen::Vector3d& normal,
    std::array<CGAL::Sqrt_extension<CGAL::Gmpq,CGAL::Gmpq>,3>* exactOutput=nullptr,
    const spin::QuadricCoefficientAuthority10* authority=nullptr)
{
    using R=CGAL::Gmpq;
    using E=CGAL::Sqrt_extension<R,R>;
    using V=std::array<R,3>;
    using M=std::array<V,3>;
    using Status=spin::QuadricInverseGaussPointStatus;
    spin::QuadricInverseGaussPointResult result;
    const auto fail=[&](Status status,const char* message) {
        result.status=status;result.diagnostic=message;return result;
    };
    if(!normal.allFinite() || !std::isfinite(normal.norm()) || normal.norm()<=minimumNormalNorm)
        return fail(Status::InvalidNormal,"algebraic finite map requires a finite normal with norm above 1e-15");
    M a{},inverse{};V b{},n{},inverseNormal{},center{};
    const auto matrix=q.quadraticMatrix();
    const auto coefficient=[&](int i) { return authority ? R(authority->geometry().at(i)) : R(q.values()(i)); };
    const std::array<std::array<int,3>,3> matrixIndices{{{{4,9,8}},{{9,5,7}},{{8,7,6}}}};
    for(int i=0;i<3;++i) {
        b[i]=coefficient(i+1);n[i]=R(normal(i));
        for(int j=0;j<3;++j)a[i][j]=authority ? coefficient(matrixIndices[i][j]) : R(matrix(i,j));
    }
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        inverse[i][j]=a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]-
                      a[(j+1)%3][(i+2)%3]*a[(j+2)%3][(i+1)%3];
    R determinant=0;for(int i=0;i<3;++i)determinant+=a[0][i]*inverse[i][0];
    if(determinant==0) {
        if(authority)return fail(Status::SingularQuadraticMatrix,"authority finite source requires exact full rank");
        const auto point=regularParaboloidPoint(q,normal);
        if(!point.valid())return point;
        // adj(A) has a nonzero nullspace column for regular rank two.
        // gradient=2*t*n, t=(u.b)/(2*u.n); the legacy w is 1/t.
        // Apply the same numerical gate on both sides of a rank transition.
        for(int column=0;column<3;++column) {
            R magnitude=0,ub=0,un=0;
            for(int row=0;row<3;++row) {
                const R u=inverse[row][column];
                magnitude+=u*u;ub+=u*b[row];un+=u*n[row];
            }
            if(magnitude==0)continue;
            if(4*un*un/(ub*ub)<=R(minimumWSquared))
                return fail(Status::BelowNumericalThreshold,"positive wSquared is within the historical 1e-14 cutoff");
            return point;
        }
        return fail(Status::SingularQuadraticMatrix,"regular rank-two nullspace column missing");
    }
    R denominator=-4*coefficient(0),normalForm=0;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
        inverse[i][j]/=determinant;
        inverseNormal[i]+=inverse[i][j]*n[j];
        center[i]-=inverse[i][j]*b[j]/2;
        denominator+=b[i]*inverse[i][j]*b[j];
    }
    if(denominator==0)return fail(Status::DegenerateQuadricScale,"algebraic finite map has exactly zero centered level");
    for(int i=0;i<3;++i)normalForm+=n[i]*inverseNormal[i];
    const R wSquared=4*normalForm/denominator;
    if(wSquared<=0)return fail(Status::IdealOrInadmissibleNormal,"normal has no finite positive-square-root point");
    if(wSquared<=R(minimumWSquared))
        return fail(Status::BelowNumericalThreshold,"positive wSquared is within the historical 1e-14 cutoff");
    const R radicand=1/wSquared;
    std::array<E,3> exactPoint;
    Eigen::Vector3d represented;
    for(int i=0;i<3;++i) {
        exactPoint[i]=E(center[i],inverseNormal[i],radicand);
        represented(i)=representedAlgebraic(exactPoint[i]);
    }
    if(!represented.allFinite())return fail(Status::NonFinitePoint,"algebraic finite point overflows binary64");
    E polynomial(coefficient(0));
    for(int i=0;i<3;++i) {
        polynomial+=E(b[i])*exactPoint[i];
        E gradient(b[i]);
        for(int j=0;j<3;++j) {
            polynomial+=E(a[i][j])*exactPoint[i]*exactPoint[j];
            gradient+=E(2*a[i][j])*exactPoint[j];
        }
        if(gradient!=E(R(0),2*n[i],radicand))
            return fail(Status::NonFinitePoint,"exact algebraic normal residual is nonzero");
    }
    if(polynomial!=E(0))return fail(Status::NonFinitePoint,"exact algebraic surface residual is nonzero");
    result.status=Status::Valid;result.point=represented;
    result.diagnostic="exact full-rank positive-root point rounded as algebraic coordinates";
    if(exactOutput)*exactOutput=exactPoint;
    return result;
}

std::pair<CGAL::Gmpq,CGAL::Gmpq> encloseRationalSourceCoordinate(
    const CGAL::Sqrt_extension<CGAL::Gmpq,CGAL::Gmpq>& value)
{
    using Q=CGAL::Gmpq;
    using Z=CGAL::Gmpz;
    if(!value.is_extended() || value.a1()==0)
        return {value.a0(),value.a0()};
    if(value.root()<0)throw std::runtime_error("negative source coordinate radicand");
    Z toleranceDenominator(1);
    mpz_mul_2exp(toleranceDenominator.mpz(),toleranceDenominator.mpz(),384);
    const Q tolerance(Z(1),toleranceDenominator);
    // Bound the algebraic value, including cancellation. The working sqrt
    // precision is not itself asserted as the final-coordinate precision.
    for(unsigned bits=384;bits<=6144;bits*=2) {
        Z denominator(1),scaled,squareRoot;
        mpz_mul_2exp(denominator.mpz(),denominator.mpz(),bits);
        mpz_mul_2exp(scaled.mpz(),value.root().numerator().mpz(),2*bits);
        mpz_fdiv_q(scaled.mpz(),scaled.mpz(),value.root().denominator().mpz());
        mpz_sqrt(squareRoot.mpz(),scaled.mpz());
        const Q rootLow(squareRoot,denominator);
        const Q rootHigh=rootLow*rootLow==value.root()
            ?rootLow:Q(squareRoot+1,denominator);
        Q lo=value.a0()+value.a1()*rootLow;
        Q hi=value.a0()+value.a1()*rootHigh;
        if(hi<lo)std::swap(lo,hi);
        if(hi-lo<=tolerance)return {lo,hi};
    }
    throw std::runtime_error("source coordinate rational enclosure resource limit");
}

std::string serializeSourceRational(const CGAL::Gmpq& value)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());stream<<value;
    return stream.str();
}

} // namespace


namespace spin
{


const char* quadricInverseGaussPointStatusToken(
    QuadricInverseGaussPointStatus status
) noexcept
{
    switch (status)
    {
    case QuadricInverseGaussPointStatus::Valid:
        return "VALID";
    case QuadricInverseGaussPointStatus::InvalidNormal:
        return "INVALID_NORMAL";
    case QuadricInverseGaussPointStatus::SingularQuadraticMatrix:
        return "SINGULAR_QUADRATIC_MATRIX";
    case QuadricInverseGaussPointStatus::DegenerateQuadricScale:
        return "DEGENERATE_QUADRIC_SCALE";
    case QuadricInverseGaussPointStatus::IdealOrInadmissibleNormal:
        return "IDEAL_OR_INADMISSIBLE_NORMAL";
    case QuadricInverseGaussPointStatus::BelowNumericalThreshold:
        return "BELOW_NUMERICAL_THRESHOLD";
    case QuadricInverseGaussPointStatus::NonFinitePoint:
        return "NONFINITE_POINT";
    }
    return "UNKNOWN";
}


bool QuadricInverseGaussPointResult::valid() const noexcept
{
    return status == QuadricInverseGaussPointStatus::Valid &&
        point.has_value();
}


const char* quadricInverseGaussDifferentialStatusToken(
    QuadricInverseGaussDifferentialStatus status
) noexcept
{
    switch (status)
    {
    case QuadricInverseGaussDifferentialStatus::Valid:
        return "VALID";
    case QuadricInverseGaussDifferentialStatus::PointUndefined:
        return "POINT_UNDEFINED";
    case QuadricInverseGaussDifferentialStatus::BorderedSystemRankDeficient:
        return "BORDERED_SYSTEM_RANK_DEFICIENT";
    case QuadricInverseGaussDifferentialStatus::BorderedSystemIllConditioned:
        return "BORDERED_SYSTEM_ILL_CONDITIONED";
    case QuadricInverseGaussDifferentialStatus::UncertifiedNumericalSolution:
        return "UNCERTIFIED_NUMERICAL_SOLUTION";
    }
    return "UNKNOWN";
}


bool QuadricInverseGaussDifferentialResult::valid() const noexcept
{
    return point.valid() &&
        status == QuadricInverseGaussDifferentialStatus::Valid &&
        coefficientJacobian.has_value();
}


QuadricInverseGaussMap::QuadricInverseGaussMap(
    const QuadricCoefficients10& quadric
)
    : quadric_(quadric)
{
    const Eigen::Matrix3d a = quadric_.quadraticMatrix();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(a);
    const double rankTolerance = 1e-12 * std::max(1.0, a.norm());
    if (
        eigen.info() != Eigen::Success ||
        (eigen.eigenvalues().array().abs() <= rankTolerance).any()
    )
    {
        preparationStatus_ =
            QuadricInverseGaussPointStatus::SingularQuadraticMatrix;
        preparationDiagnostic_ =
            "inverse Gauss map requires the triangulator's full-rank "
            "quadratic matrix";
        return;
    }

    inverseA_ = a.inverse();
    if (!inverseA_.allFinite())
    {
        preparationStatus_ =
            QuadricInverseGaussPointStatus::DegenerateQuadricScale;
        preparationDiagnostic_ =
            "inverse Gauss map quadratic inverse is nonfinite";
        return;
    }
    if (!exactInverseGaussDenominator(quadric_, denominator_))
    {
        preparationStatus_ =
            QuadricInverseGaussPointStatus::SingularQuadraticMatrix;
        preparationDiagnostic_ =
            "inverse Gauss map quadratic matrix is exactly singular";
        return;
    }
    if (!std::isfinite(denominator_) || denominator_ == 0.0)
    {
        preparationStatus_ =
            QuadricInverseGaussPointStatus::DegenerateQuadricScale;
        preparationDiagnostic_ = denominator_ == 0.0
            ? "inverse Gauss map centered level is exactly zero or "
              "underflows binary64"
            : "inverse Gauss map denominator is nonfinite";
    }
}


QuadricInverseGaussPointResult QuadricInverseGaussMap::point(
    const Eigen::Vector3d& normal
) const
{
    QuadricInverseGaussPointResult result;
    const double normalNorm = normal.norm();
    if (
        !normal.allFinite() || !std::isfinite(normalNorm) ||
        normalNorm <= minimumNormalNorm
    )
    {
        result.status = QuadricInverseGaussPointStatus::InvalidNormal;
        result.diagnostic =
            "inverse Gauss map requires a finite normal with norm above 1e-15";
        return result;
    }
    if (preparationStatus_ != QuadricInverseGaussPointStatus::Valid)
    {
        result.status = preparationStatus_;
        result.diagnostic = preparationDiagnostic_;
        return result;
    }

    const double wSquared =
        4.0 * normal.dot(inverseA_ * normal) / denominator_;
    if (!std::isfinite(wSquared))
    {
        result.status = QuadricInverseGaussPointStatus::NonFinitePoint;
        result.diagnostic = "inverse Gauss map produced nonfinite wSquared";
        return result;
    }
    if (wSquared <= 0.0)
    {
        result.status =
            QuadricInverseGaussPointStatus::IdealOrInadmissibleNormal;
        result.diagnostic =
            "normal has no finite point on the positive-square-root branch";
        return result;
    }
    if (wSquared <= minimumWSquared)
    {
        result.status = QuadricInverseGaussPointStatus::BelowNumericalThreshold;
        result.diagnostic =
            "positive wSquared is within the historical 1e-14 cutoff; "
            "normal magnitude is part of this numerical contract";
        return result;
    }

    const double w = std::sqrt(wSquared);
    const Eigen::Vector3d mapped = inverseA_ *
        (normal / w - 0.5 * quadric_.linearVector());
    if (!mapped.allFinite())
    {
        result.status = QuadricInverseGaussPointStatus::NonFinitePoint;
        result.diagnostic = "inverse Gauss point is nonfinite";
        return result;
    }
    result.status = QuadricInverseGaussPointStatus::Valid;
    result.point = mapped;
    return result;
}


QuadricInverseGaussDifferentialResult
QuadricInverseGaussMap::coefficientDifferential(
    const Eigen::Vector3d& fixedNormal
) const
{
    return coefficientDifferentialAtPoint(fixedNormal, point(fixedNormal));
}

QuadricInverseGaussDifferentialResult
QuadricInverseGaussMap::algebraicCoefficientDifferential(
    const Eigen::Vector3d& fixedNormal
) const
{
    return coefficientDifferentialAtPoint(fixedNormal,algebraicPoint(fixedNormal));
}

QuadricInverseGaussPointResult QuadricInverseGaussMap::algebraicPoint(
    const Eigen::Vector3d& normal
) const
{
    return algebraicRegularPoint(quadric_,normal);
}

QuadricRationalPointResult QuadricInverseGaussMap::rationalPoint(
    const Eigen::Vector3d& normal, bool compactRationalCoordinates
) const
{
    return rationalPointImpl(normal,compactRationalCoordinates,nullptr);
}

QuadricRationalPointResult QuadricInverseGaussMap::rationalPointFromAuthority(
    const QuadricCoefficientAuthority10& authority,
    const Eigen::Vector3d& normal, bool compactRationalCoordinates)
{
    QuadricCoefficients10::Vector display;
    const auto values=authority.geometryDisplay();
    for(int i=0;i<10;++i)display(i)=values[i];
    const QuadricInverseGaussMap map(QuadricCoefficients10::fromNormalized(display));
    return map.rationalPointImpl(normal,compactRationalCoordinates,&authority);
}

QuadricRationalPointResult QuadricInverseGaussMap::rationalPointImpl(
    const Eigen::Vector3d& normal, bool compactRationalCoordinates,
    const QuadricCoefficientAuthority10* authority
) const
{
    using Q=CGAL::Gmpq;
    using E=CGAL::Sqrt_extension<Q,Q>;
    QuadricRationalPointResult result;
    const auto a=quadric_.quadraticMatrix();
    const Q determinant=Q(a(0,0))*(Q(a(1,1))*Q(a(2,2))-Q(a(1,2))*Q(a(2,1)))
        -Q(a(0,1))*(Q(a(1,0))*Q(a(2,2))-Q(a(1,2))*Q(a(2,0)))
        +Q(a(0,2))*(Q(a(1,0))*Q(a(2,1))-Q(a(1,1))*Q(a(2,0)));
    if(!authority && determinant==0) {
        result.diagnostic="rational source point requires exact full rank";
        return result;
    }
    std::array<E,3> exactPoint;
    const auto admission=algebraicRegularPoint(quadric_,normal,&exactPoint,authority);
    if(!admission.valid()) {
        result.diagnostic=admission.diagnostic;return result;
    }
    try {
        for(int axis=0;axis<3;++axis) {
            const auto interval=encloseRationalSourceCoordinate(exactPoint[axis]);
            result.lower[axis]=serializeSourceRational(interval.first);
            result.upper[axis]=serializeSourceRational(interval.second);
            result.coordinates[axis]=serializeSourceRational(rationalSourceRepresentative(
                interval.first,interval.second,compactRationalCoordinates));
        }
    } catch(const std::runtime_error& error) {
        result.diagnostic=error.what();return result;
    }
    result.defined=true;
    result.diagnostic="explicit rational full-rank finite source; exact enclosure width <= 2^-384; not an exact analytic point";
    return result;
}

QuadricRationalFiniteDifferentialResult QuadricInverseGaussMap::rationalDifferentialFromAuthority(
    const QuadricCoefficientAuthority10& authority,const Eigen::Vector3d& normal,
    bool compactRationalCoordinates)
{
    using Q=CGAL::Gmpq;
    using E=CGAL::Sqrt_extension<Q,Q>;
    using Vector=std::array<Q,3>;
    using Matrix=std::array<Vector,3>;
    QuadricRationalFiniteDifferentialResult result;
    const auto display=authority.geometryDisplay();
    QuadricCoefficients10::Vector represented;
    for(int i=0;i<10;++i)represented(i)=display[i];
    const auto view=QuadricCoefficients10::fromNormalized(represented);
    std::array<E,3> point;
    const auto admission=algebraicRegularPoint(view,normal,&point,&authority);
    if(!admission.valid()) { result.diagnostic=admission.diagnostic;return result; }
    const std::array<std::array<int,3>,3> index{{{{4,9,8}},{{9,5,7}},{{8,7,6}}}};
    Matrix a{},inverse{};Vector b{},n{},u{},inverseB{};
    std::array<Q,10> q;
    for(int i=0;i<10;++i)q[i]=Q(authority.geometry()[i]);
    for(int i=0;i<3;++i) {
        b[i]=q[i+1];n[i]=Q(normal(i));
        for(int j=0;j<3;++j)a[i][j]=q[index[i][j]];
    }
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        inverse[i][j]=a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]-
                      a[(j+1)%3][(i+2)%3]*a[(j+2)%3][(i+1)%3];
    Q determinant=0;for(int i=0;i<3;++i)determinant+=a[0][i]*inverse[i][0];
    // Admission has already certified full rank, finite normal and positive
    // radicand. No rank-two or represented-coefficient fallback is allowed.
    Q denominator=-4*q[0],normalForm=0;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
        inverse[i][j]/=determinant;
        u[i]+=inverse[i][j]*n[j];inverseB[i]+=inverse[i][j]*b[j];
        denominator+=b[i]*inverse[i][j]*b[j];
    }
    for(int i=0;i<3;++i)normalForm+=n[i]*u[i];
    const Q radicand=denominator/(4*normalForm);
    std::array<std::array<E,10>,3> jacobian;
    try {
        for(int column=0;column<10;++column) {
            Matrix da{};Vector db{},dInverseB{},du{};
            const Q dc=column==0 ? Q(1) : Q(0);
            for(int i=0;i<3;++i) {
                db[i]=column==i+1 ? Q(1) : Q(0);
                for(int j=0;j<3;++j)da[i][j]=column==index[i][j] ? Q(1) : Q(0);
            }
            Q dDenominator=-4*dc,dForm=0;
            for(int i=0;i<3;++i) {
                dDenominator+=2*db[i]*inverseB[i];
                for(int j=0;j<3;++j) {
                    dDenominator-=inverseB[i]*da[i][j]*inverseB[j];
                    dForm-=u[i]*da[i][j]*u[j];
                    for(int k=0;k<3;++k) {
                        dInverseB[i]-=inverse[i][j]*da[j][k]*inverseB[k];
                        du[i]-=inverse[i][j]*da[j][k]*u[k];
                    }
                    dInverseB[i]+=inverse[i][j]*db[j];
                }
            }
            const Q logarithmicRadicand=dDenominator/denominator-dForm/normalForm;
            for(int axis=0;axis<3;++axis)
                jacobian[axis][column]=E(-dInverseB[axis]/2,
                    du[axis]+u[axis]*logarithmicRadicand/2,radicand);
            // Independently certify the implicit surface and normal equations
            // at this exact algebraic point, including cross coefficient 2.
            E surface(dc);std::array<E,3> gradientMotion;
            for(int i=0;i<3;++i) {
                surface+=E(db[i])*point[i]+E(b[i])*jacobian[i][column];
                gradientMotion[i]=E(db[i]);
                for(int j=0;j<3;++j) {
                    surface+=E(da[i][j])*point[i]*point[j]+
                             E(2*a[i][j])*point[j]*jacobian[i][column];
                    gradientMotion[i]+=E(2*da[i][j])*point[j]+
                                       E(2*a[i][j])*jacobian[j][column];
                }
            }
            if(surface!=E(0))throw std::runtime_error("authority finite differential has nonzero exact surface residual");
            for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)
                if(E(n[i])*gradientMotion[j]-E(n[j])*gradientMotion[i]!=E(0))
                    throw std::runtime_error("authority finite differential has nonzero exact normal residual");
        }
        for(int axis=0;axis<3;++axis) {
            E radial(0);for(int column=0;column<10;++column)radial+=E(q[column])*jacobian[axis][column];
            if(radial!=E(0))throw std::runtime_error("authority finite differential failed radial invariance");
            // Radial invariance proves that coefficient renormalization leaves
            // this local point differential unchanged, without a float tangent.
            const auto pInterval=encloseRationalSourceCoordinate(point[axis]);
            result.point.lower[axis]=serializeSourceRational(pInterval.first);
            result.point.upper[axis]=serializeSourceRational(pInterval.second);
            result.point.coordinates[axis]=serializeSourceRational(rationalSourceRepresentative(
                pInterval.first,pInterval.second,compactRationalCoordinates));
            for(int column=0;column<10;++column) {
                const auto interval=encloseRationalSourceCoordinate(jacobian[axis][column]);
                result.lower[axis][column]=serializeSourceRational(interval.first);
                result.upper[axis][column]=serializeSourceRational(interval.second);
                result.coefficientJacobian[axis][column]=serializeSourceRational(rationalSourceRepresentative(
                    interval.first,interval.second,compactRationalCoordinates));
            }
        }
    } catch(const std::runtime_error& error) {
        result.diagnostic=error.what();return result;
    }
    result.point.defined=true;
    result.point.diagnostic="same authority finite-point algebraic enclosure";
    result.defined=true;
    result.diagnostic="authority full-rank finite differential; exact algebraic residuals and <=2^-384 enclosures; not rounding or Boolean derivative";
    return result;
}

QuadricInverseGaussDifferentialResult
QuadricInverseGaussMap::coefficientDifferentialAtPoint(
    const Eigen::Vector3d& fixedNormal,
    const QuadricInverseGaussPointResult& suppliedPoint
) const
{
    QuadricInverseGaussDifferentialResult result;
    result.point = suppliedPoint;
    if (!result.point.valid())
    {
        result.status = QuadricInverseGaussDifferentialStatus::PointUndefined;
        result.diagnostic = "point is undefined: " + result.point.diagnostic;
        return result;
    }

    const Eigen::Vector3d unitNormal = fixedNormal / fixedNormal.norm();
    Eigen::Matrix4d system = Eigen::Matrix4d::Zero();
    system.topLeftCorner<3, 3>() = quadric_.quadraticMatrix();
    system.topRightCorner<3, 1>() = -unitNormal;
    system.bottomLeftCorner<1, 3>() =
        quadric_.gradient(*result.point.point).transpose();

    Eigen::FullPivLU<Eigen::Matrix4d> factorization(system);
    result.rawReciprocalCondition = factorization.rcond();

    // The bordered equations mix the unit normal, the quadric's quadratic
    // coefficients, and its gradient.  Their physical scales can differ by
    // many orders of magnitude near a cone even when the exact binary64
    // system is nonsingular.  Reversibly equilibrate rows and columns for the
    // numerical rank/condition decision, while retaining the original system
    // for the final residual and backward-error certificates below.
    Eigen::Vector4d rowScale;
    for (int row = 0; row < 4; ++row)
    {
        const double rowMaximum = system.row(row).cwiseAbs().maxCoeff();
        if (!std::isfinite(rowMaximum) || rowMaximum <= 0.0)
        {
            result.status =
                QuadricInverseGaussDifferentialStatus::BorderedSystemRankDeficient;
            result.diagnostic =
                "fixed-normal bordered system has a zero or nonfinite row";
            return result;
        }
        rowScale(row) = 1.0 / rowMaximum;
    }
    const Eigen::Matrix4d rowScaledSystem = rowScale.asDiagonal() * system;
    Eigen::Vector4d columnScale;
    for (int column = 0; column < 4; ++column)
    {
        const double columnMaximum =
            rowScaledSystem.col(column).cwiseAbs().maxCoeff();
        if (!std::isfinite(columnMaximum) || columnMaximum <= 0.0)
        {
            result.status =
                QuadricInverseGaussDifferentialStatus::BorderedSystemRankDeficient;
            result.diagnostic =
                "fixed-normal bordered system has a zero or nonfinite column";
            return result;
        }
        columnScale(column) = 1.0 / columnMaximum;
    }
    const Eigen::Matrix4d equilibratedSystem =
        rowScaledSystem * columnScale.asDiagonal();
    Eigen::FullPivLU<Eigen::Matrix4d> equilibratedFactorization(
        equilibratedSystem
    );
    if (equilibratedFactorization.rank() < 4)
    {
        result.status =
            QuadricInverseGaussDifferentialStatus::BorderedSystemRankDeficient;
        std::ostringstream diagnostic;
        diagnostic.imbue(std::locale::classic());
        diagnostic << std::setprecision(std::numeric_limits<double>::max_digits10)
                   << "reversibly equilibrated fixed-normal bordered system is "
                   << "rank deficient; raw_rank=" << factorization.rank()
                   << "; equilibrated_rank=" << equilibratedFactorization.rank()
                   << "; raw_rcond=" << result.rawReciprocalCondition;
        result.diagnostic = diagnostic.str();
        return result;
    }
    result.equilibratedReciprocalCondition =
        equilibratedFactorization.rcond();
    const double reciprocalCondition =
        result.equilibratedReciprocalCondition;
    if (
        !std::isfinite(reciprocalCondition) ||
        reciprocalCondition <= minimumBorderedReciprocalCondition
    )
    {
        result.status =
            QuadricInverseGaussDifferentialStatus::BorderedSystemIllConditioned;
        std::ostringstream diagnostic;
        diagnostic.imbue(std::locale::classic());
        diagnostic << std::setprecision(std::numeric_limits<double>::max_digits10)
                   << "fixed-normal bordered differential system has reciprocal "
                   << "condition at or below "
                   << minimumBorderedReciprocalCondition
                   << "; equilibrated_rcond=" << reciprocalCondition
                   << "; raw_rcond=" << result.rawReciprocalCondition;
        result.diagnostic = diagnostic.str();
        return result;
    }

    Eigen::Matrix<double, 4, 10> rightHandSide =
        Eigen::Matrix<double, 4, 10>::Zero();
    const Eigen::Vector3d& mapped = *result.point.point;
    rightHandSide(3, 0) = -1.0;
    for (int coordinate = 0; coordinate < 3; ++coordinate)
    {
        rightHandSide(coordinate, 1 + coordinate) = -0.5;
        rightHandSide(3, 1 + coordinate) = -mapped(coordinate);
        rightHandSide(coordinate, 4 + coordinate) = -mapped(coordinate);
        rightHandSide(3, 4 + coordinate) =
            -mapped(coordinate) * mapped(coordinate);
    }
    const std::array<std::array<int, 2>, 3> crossCoordinates = {{
        {{1, 2}}, {{0, 2}}, {{0, 1}}
    }};
    for (int cross = 0; cross < 3; ++cross)
    {
        const int first = crossCoordinates[cross][0];
        const int second = crossCoordinates[cross][1];
        rightHandSide(first, 7 + cross) = -mapped(second);
        rightHandSide(second, 7 + cross) = -mapped(first);
        rightHandSide(3, 7 + cross) =
            -2.0 * mapped(first) * mapped(second);
    }

    // The bordered system can be moderately ill-conditioned while remaining
    // backward stable.  In that regime a double-precision factorization may
    // produce a normalized residual just above the certification threshold
    // solely because a large solution is rounded.  Solve and certify the
    // small system in the platform's extended precision, then expose the
    // representable double Jacobian used by the rest of the pipeline.
    using ExtendedSystem = Eigen::Matrix<long double, 4, 4>;
    using ExtendedRightHandSide = Eigen::Matrix<long double, 4, 10>;
    const ExtendedSystem extendedSystem = system.cast<long double>();
    const ExtendedRightHandSide extendedRightHandSide =
        rightHandSide.cast<long double>();
    const ExtendedSystem extendedEquilibratedSystem =
        equilibratedSystem.cast<long double>();
    const ExtendedRightHandSide extendedEquilibratedRightHandSide =
        (rowScale.asDiagonal() * rightHandSide).cast<long double>();
    Eigen::FullPivLU<ExtendedSystem> extendedFactorization(
        extendedEquilibratedSystem
    );
    if (extendedFactorization.rank() < 4)
    {
        result.status =
            QuadricInverseGaussDifferentialStatus::BorderedSystemRankDeficient;
        result.diagnostic =
            "extended-precision fixed-normal bordered system is rank deficient";
        return result;
    }
    const ExtendedRightHandSide equilibratedSolution =
        extendedFactorization.solve(extendedEquilibratedRightHandSide);
    const ExtendedRightHandSide extendedSolution =
        columnScale.cast<long double>().asDiagonal() * equilibratedSolution;
    long double extendedResidualNorm =
        (extendedSystem * extendedSolution - extendedRightHandSide).norm();
    long double extendedRelativeResidual = extendedResidualNorm /
        std::max(1.0L, extendedRightHandSide.norm());
    Eigen::Matrix<double, 4, 10> solution =
        extendedSolution.cast<double>();
    bool exactRationalSolveUsed = false;
    if (extendedRelativeResidual > 1e-10)
    {
        // For near-cone strata, the exact binary64 system can be full rank and
        // backward stable while its very large solution makes long-double
        // residual evaluation round above the unchanged 1e-10 gate.  Solve
        // this 4x4 system exactly as rationals only on that fallback path.
        if (!solveExactBorderedSystem(system, rightHandSide, solution))
        {
            result.status =
                QuadricInverseGaussDifferentialStatus::UncertifiedNumericalSolution;
            result.diagnostic =
                "exact-rational fixed-normal bordered solve failed";
            return result;
        }
        exactRationalSolveUsed = true;
        // solveExactBorderedSystem verifies zero exact residual for the input
        // binary64 system before returning the rounded Jacobian. Keep the
        // extended-precision solution only for the diagnostic scale; the
        // separate rounded-backward-error check below audits the returned
        // binary64 matrix.
        extendedRelativeResidual = 0.0L;
        extendedResidualNorm = 0.0L;
    }
    const long double extendedBackwardScale = extendedSystem.norm() *
        extendedSolution.norm() + extendedRightHandSide.norm();
    const long double extendedBackwardError =
        std::isfinite(extendedBackwardScale) && extendedBackwardScale > 0.0L
            ? extendedResidualNorm / extendedBackwardScale
            : std::numeric_limits<long double>::infinity();
    const ExtendedRightHandSide roundedSolution = solution.cast<long double>();
    const long double roundedResidualNorm =
        (extendedSystem * roundedSolution - extendedRightHandSide).norm();
    const long double roundedBackwardScale = extendedSystem.norm() *
        roundedSolution.norm() + extendedRightHandSide.norm();
    const long double roundedBackwardError =
        std::isfinite(roundedBackwardScale) && roundedBackwardScale > 0.0L
            ? roundedResidualNorm / roundedBackwardScale
            : std::numeric_limits<long double>::infinity();
    const double relativeResidual =
        static_cast<double>(extendedRelativeResidual);
    const double backwardError =
        static_cast<double>(extendedBackwardError);
    const double roundedRelativeResidual = static_cast<double>(
        roundedResidualNorm / std::max(1.0L, extendedRightHandSide.norm()));
    const double roundedBackwardErrorDouble =
        static_cast<double>(roundedBackwardError);
    const double roundedBackwardErrorLimit =
        64.0 * std::numeric_limits<double>::epsilon();
    if (!solution.allFinite() || !std::isfinite(relativeResidual) ||
        !std::isfinite(backwardError) ||
        !std::isfinite(roundedRelativeResidual) ||
        !std::isfinite(roundedBackwardErrorDouble) ||
        relativeResidual > 1e-10 ||
        roundedBackwardErrorDouble > roundedBackwardErrorLimit)
    {
        result.status =
            QuadricInverseGaussDifferentialStatus::UncertifiedNumericalSolution;
        std::ostringstream diagnostic;
        diagnostic.imbue(std::locale::classic());
        diagnostic << std::setprecision(17)
                   << "fixed-normal bordered differential solve did not "
                      "certify a finite 1e-10 relative residual"
                   << "; relative_residual=" << relativeResidual
                   << "; reciprocal_condition=" << reciprocalCondition
                   << "; system_norm=" << system.norm()
                   << "; rhs_norm=" << rightHandSide.norm()
                   << "; solution_norm=" << solution.norm()
                   << "; backward_error=" << backwardError
                   << "; rounded_relative_residual="
                   << roundedRelativeResidual
                   << "; rounded_backward_error="
                   << roundedBackwardErrorDouble
                   << "; rounded_backward_error_limit="
                   << roundedBackwardErrorLimit
                   << "; exact_rational_solve="
                   << (exactRationalSolveUsed ? "true" : "false")
                   << "; fixed_normal=" << unitNormal.x() << ','
                   << unitNormal.y() << ',' << unitNormal.z();
        result.diagnostic = diagnostic.str();
        return result;
    }

    const QuadricCoefficients10::Vector& coefficients = quadric_.values();
    const Eigen::Matrix<double, 10, 10> normalizationTangent =
        Eigen::Matrix<double, 10, 10>::Identity() -
        coefficients * coefficients.transpose();
    const CoefficientJacobian normalizedDerivative =
        solution.topRows<3>() * normalizationTangent;
    if (!normalizedDerivative.allFinite())
    {
        result.status =
            QuadricInverseGaussDifferentialStatus::UncertifiedNumericalSolution;
        result.diagnostic = "normalized coefficient differential is nonfinite";
        return result;
    }

    result.status = QuadricInverseGaussDifferentialStatus::Valid;
    result.coefficientJacobian = normalizedDerivative;
    return result;
}


QuadricInverseGaussDifferentialResult
QuadricInverseGaussMap::paraboloidCoefficientDifferential(const Eigen::Vector3d& normal) const
{
    return coefficientDifferentialAtPoint(normal, regularParaboloidPoint(quadric_,normal));
}

const Eigen::Matrix3d& QuadricInverseGaussMap::inverseQuadraticMatrix() const
    noexcept
{
    return inverseA_;
}


QuadricIdealEdgeResult QuadricInverseGaussMap::idealEdgeDifferential(
    const Eigen::Vector3d& first, const Eigen::Vector3d& second) const
{
    const auto fail=[](const std::string& why) {
        QuadricIdealEdgeResult out;out.diagnostic=why;return out;
    };
    if(preparationStatus_!=QuadricInverseGaussPointStatus::Valid)
        return fail(preparationDiagnostic_);
    if(!first.allFinite() || !second.allFinite() ||
       first.norm()<=minimumNormalNorm || second.norm()<=minimumNormalNorm ||
       first.cross(second).norm()<=minimumNormalNorm)
        return fail("ideal edge requires independent finite nonzero normals");
    using Rational=CGAL::Gmpq;
    using Algebraic=CGAL::Sqrt_extension<Rational,Rational>;
    const auto& p=quadric_.values();
    const Rational a[3][3]={{Rational(p(4)),Rational(p(9)),Rational(p(8))},
        {Rational(p(9)),Rational(p(5)),Rational(p(7))},
        {Rational(p(8)),Rational(p(7)),Rational(p(6))}};
    Rational adj[3][3];
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        adj[i][j]=a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]-
                  a[(j+1)%3][(i+2)%3]*a[(j+2)%3][(i+1)%3];
    Rational det(0);for(int j=0;j<3;++j)det+=a[0][j]*adj[j][0];
    if(det==0)return fail("exact singular quadratic matrix at ideal edge");
    std::array<Rational,3> n,d;
    for(int i=0;i<3;++i){n[i]=Rational(first(i));d[i]=Rational(second(i))-n[i];}
    Rational quadratic(0),linear(0),constant(0);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
        quadratic+=d[i]*adj[i][j]*d[j];
        linear+=2*n[i]*adj[i][j]*d[j];
        constant+=n[i]*adj[i][j]*n[j];
    }
    if(constant==0 || quadratic+linear+constant==0)
        return fail("ideal normal lies at an edge endpoint");
    std::vector<Algebraic> roots;
    if(quadratic==0) {
        if(linear!=0)roots.emplace_back(-constant/linear);
    } else {
        const Rational discriminant=linear*linear-4*quadratic*constant;
        if(discriminant==0) {
            const Rational t=-linear/(2*quadratic);
            if(t>0 && t<1)return fail("tangent ideal-normal event");
        } else if(discriminant>0) {
            roots.emplace_back(-linear/(2*quadratic),-1/(2*quadratic),discriminant);
            roots.emplace_back(-linear/(2*quadratic),1/(2*quadratic),discriminant);
        }
    }
    std::sort(roots.begin(),roots.end());
    QuadricIdealEdgeResult out;
    const Eigen::Vector3d delta=second-first;
    const Eigen::Matrix<double,10,10> tangent=
        Eigen::Matrix<double,10,10>::Identity()-p*p.transpose();
    for(const auto& root:roots) {
        if(root<=Algebraic(0) || root>=Algebraic(1))continue;
        QuadricIdealEdgePoint point;point.parameter=representedAlgebraic(root);
        if(!std::isfinite(point.parameter) || point.parameter<=0 || point.parameter>=1)
            return fail("ideal root has no separated binary64 interior representation");
        for(int i=0;i<3;++i) {
            Algebraic coordinate(0);
            for(int j=0;j<3;++j)
                coordinate+=Algebraic(adj[i][j]/det)*(Algebraic(n[j])+Algebraic(d[j])*root);
            point.direction(i)=representedAlgebraic(coordinate);
        }
        const Eigen::Vector3d normal=first+point.parameter*delta;
        const double slope=2*delta.dot(point.direction);
        if(!point.direction.allFinite() || !std::isfinite(slope) ||
           point.direction.norm()<=minimumNormalNorm ||
           std::abs(normal.dot(point.direction))>1e-12*normal.norm()*point.direction.norm() ||
           std::abs(slope)<=1e-10*delta.norm()*point.direction.norm())
            return fail("ideal root representation or differential is ill-conditioned");
        // n^T B n=0 gives dt=(D^T dA D)/(2 delta^T D), and
        // dD=-B dA D+B delta dt. b,c do not enter before normalization.
        Eigen::Matrix<double,3,10> derivative=Eigen::Matrix<double,3,10>::Zero();
        for(int j=4;j<10;++j) {
            Eigen::Matrix3d da=Eigen::Matrix3d::Zero();
            if(j<7)da(j-4,j-4)=1;
            else {
                const int row=j==7?1:0,col=j==9?1:2;
                da(row,col)=da(col,row)=1;
            }
            const double dt=point.direction.dot(da*point.direction)/slope;
            derivative.col(j)=inverseA_*(-da*point.direction+delta*dt);
        }
        point.coefficientJacobian=derivative*tangent;
        if(!point.coefficientJacobian.allFinite())return fail("nonfinite ideal-point differential");
        out.points.push_back(point);
    }
    out.defined=true;
    out.diagnostic="paper ideal-edge points; face stitching and truncation are separate";
    return out;
}

QuadricProjectiveRayResult truncateQuadricProjectiveRay(
    const Eigen::Vector3d& origin, const Eigen::Vector3d& direction, double radius)
{
    const auto fail=[](const std::string& why) {
        QuadricProjectiveRayResult out;out.diagnostic=why;return out;
    };
    if(!origin.allFinite() || !direction.allFinite() || !std::isfinite(radius) || radius<=0)
        return fail("projective ray requires finite coordinates and positive radius");
    // Even extended precision can reverse the sign near the sphere. The
    // strict-domain decision concerns the exact represented input coordinates.
    CGAL::Gmpq exactOffset=-CGAL::Gmpq(radius)*CGAL::Gmpq(radius);
    for(int i=0;i<3;++i)exactOffset+=CGAL::Gmpq(origin(i))*CGAL::Gmpq(origin(i));
    if(exactOffset>=0)
        return fail("projective ray origin is not exactly strictly interior");
    using Vector=Eigen::Matrix<long double,3,1>;
    const Vector p=origin.cast<long double>(),d=direction.cast<long double>();
    const long double r=radius,a=d.squaredNorm(),b=p.dot(d),c=p.squaredNorm()-r*r;
    if(!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || a<=0 || c>=0)
        return fail("projective ray requires nonzero direction and strictly interior finite origin");
    const long double discriminant=b*b-a*c;
    if(!std::isfinite(discriminant) || discriminant<=0)
        return fail("invalid projective ray/sphere discriminant");
    const long double root=std::sqrt(discriminant);
    const long double t=b>=0 ? -c/(b+root) : (-b+root)/a;
    const Vector x=p+t*d;
    const long double denominator=x.dot(d);
    if(!std::isfinite(t) || t<=0 || !x.allFinite() ||
       denominator<=1e-12L*r*std::sqrt(a))
        return fail("ill-conditioned projective ray truncation");
    const Eigen::Matrix<long double,3,3> projection=
        Eigen::Matrix<long double,3,3>::Identity()-d*x.transpose()/denominator;
    QuadricProjectiveRayResult out;
    out.point=x.cast<double>();
    out.jacobian.leftCols<3>()=projection.cast<double>();
    out.jacobian.middleCols<3>(3)=(t*projection).cast<double>();
    out.jacobian.col(6)=(r*d/denominator).cast<double>();
    if(!out.point.allFinite() || !out.jacobian.allFinite() ||
       std::abs(out.point.cast<long double>().norm()-r)>1e-12L*r)
        return fail("uncertified represented projective ray output");
    out.defined=true;
    out.diagnostic="projective straight-edge truncation; not an analytic quadric intersection";
    return out;
}

} // namespace spin
