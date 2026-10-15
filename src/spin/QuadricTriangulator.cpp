#include "spin/QuadricTriangulator.h"

#include <igl/bfs_orient.h>
#include <igl/piecewise_constant_winding_number.h>

#include "manufacturing/SlicerValidation.h"
#include "spin/QuadricInverseGaussMap.h"
#include "spin/TruncationRootCertificate.h"
#include "spin/detail/ProjectedDomainClearance.h"
#include "detail/FixedNormalSphere.h"
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactBoundaryAudit.h"
#include "verification/ExactRationalCoordinates.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <igl/copyleft/cgal/assign_scalar.h>
#endif

#include <CGAL/Gmpq.h>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>


namespace
{


using Face = std::array<int, 3>;
using Edge = std::pair<int, int>;

#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
using RationalPoint=std::array<CGAL::Gmpq,3>;
RationalPoint rationalPoint(const spin::QuadricRationalPointResult& point) {
    RationalPoint out;
    for(int a=0;a<3;++a)out[a]=verification::parseExactRationalCoordinate(point.coordinates[a]);
    return out;
}
Eigen::Vector3d rationalDisplay(const spin::QuadricRationalPointResult& point) {
    Eigen::Vector3d out;const auto p=rationalPoint(point);
    for(int a=0;a<3;++a)igl::copyleft::cgal::assign_scalar(CGAL::Epeck::FT(p[a]),out(a));
    if(!out.allFinite())throw std::runtime_error("rational source display is nonfinite");
    return out;
}
void orientRationalInterfaceFace(const spin::QuadricCoefficients10& quadric,
    const std::vector<spin::QuadricRationalPointResult>& points,Face* face,
    const spin::QuadricCoefficientAuthority10* authority=nullptr) {
    const auto a=rationalPoint(points.at((*face)[0])),b=rationalPoint(points.at((*face)[1])),
        c=rationalPoint(points.at((*face)[2]));
    RationalPoint u,v,center;
    for(int axis=0;axis<3;++axis){u[axis]=b[axis]-a[axis];v[axis]=c[axis]-a[axis];center[axis]=(a[axis]+b[axis]+c[axis])/3;}
    const RationalPoint cross{{u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]}};
    CGAL::Gmpq dot(0);
    const std::array<std::array<int,3>,3> indices{{{{4,9,8}},{{9,5,7}},{{8,7,6}}}};
    for(int axis=0;axis<3;++axis) {
        CGAL::Gmpq gradient=authority ? CGAL::Gmpq(authority->geometry()[1+axis]) : CGAL::Gmpq(quadric.linearVector()(axis));
        for(int k=0;k<3;++k)gradient+=2*(authority ? CGAL::Gmpq(authority->geometry()[indices[axis][k]]) : CGAL::Gmpq(quadric.quadraticMatrix()(axis,k)))*center[k];
        dot-=cross[axis]*gradient;
    }
    if(dot==0)throw std::runtime_error("rational source interface orientation is ambiguous");
    if(dot<0)std::swap((*face)[1],(*face)[2]);
}
void assembleRationalAuthority(const std::vector<spin::QuadricRationalPointResult>& points,
                              spin::QuadricTriangulation* result) {
    for(const auto& point:points) {
        result->exactVertices.push_back(point.coordinates);
        result->coordinateLower.push_back(point.lower);result->coordinateUpper.push_back(point.upper);
    }
    Eigen::MatrixXi oriented;Eigen::VectorXi components;
    igl::bfs_orient(result->closedBoundary.faces,oriented,components);
    if(oriented.rows()!=result->closedBoundary.faces.rows() || components.size()!=oriented.rows() || components.size()==0)
        throw std::runtime_error("rational boundary could not assign all face components");
    result->closedBoundary.faces=std::move(oriented);
    std::vector<CGAL::Gmpq> volumes(components.maxCoeff()+1,CGAL::Gmpq(0));
    for(int f=0;f<result->closedBoundary.faces.rows();++f) {
        const auto a=rationalPoint(points.at(result->closedBoundary.faces(f,0))),
            b=rationalPoint(points.at(result->closedBoundary.faces(f,1))),
            c=rationalPoint(points.at(result->closedBoundary.faces(f,2)));
        volumes[components(f)]+=a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]);
    }
    for(int f=0;f<result->closedBoundary.faces.rows();++f) {
        if(volumes[components(f)]==0)throw std::runtime_error("rational source component has zero exact volume");
        if(volumes[components(f)]<0)std::swap(result->closedBoundary.faces(f,1),result->closedBoundary.faces(f,2));
    }
    const auto audit=verification::auditExactTriangleBoundaryRational(result->exactVertices,
        result->closedBoundary.faces,verification::ExactBoundaryOrientation::OutwardPositive);
    if(!audit.facts().boundaryValid)throw std::runtime_error("rational source failed exact authority boundary audit");
}
#endif

class BoundaryLoopTopologyError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

Eigen::Matrix3d canonicalEigenvectorSigns(const Eigen::Matrix3d& input)
{
    Eigen::Matrix3d result=input;
    for(int column=0;column<3;++column)
    {
        Eigen::Index pivot=0;
        result.col(column).cwiseAbs().maxCoeff(&pivot);
        if(result(pivot,column)<0.0)result.col(column)=-result.col(column);
    }
    return result;
}


using SphereMesh = spin::detail::FixedNormalSphere;

int shiftedLevelSign(const spin::QuadricCoefficients10& quadric)
{
    const Eigen::Matrix3d a = quadric.quadraticMatrix();
    const Eigen::Vector3d b = quadric.linearVector();
    const double shifted = quadric.constant() -
        0.25 * b.dot(a.inverse() * b);
    const double scale = std::max(1.0, a.norm());
    if (std::abs(shifted) > 1e-12 * scale)
    {
        return shifted > 0.0 ? 1 : -1;
    }

    // The double-precision tolerance is useful for rank decisions, but it
    // must not turn a regular, nonzero level into an exactly degenerate cone.
    // Near that boundary, determine the sign of 4*c*det(A)-b^T adj(A)b from
    // the exact binary64 input values using rational arithmetic.
    using Rational = CGAL::Gmpq;
    const auto& coefficients = quadric.values();
    const Rational a00(coefficients(4));
    const Rational a11(coefficients(5));
    const Rational a22(coefficients(6));
    const Rational a01(coefficients(9));
    const Rational a02(coefficients(8));
    const Rational a12(coefficients(7));
    const Rational determinant =
        a00 * (a11 * a22 - a12 * a12) -
        a01 * (a01 * a22 - a12 * a02) +
        a02 * (a01 * a12 - a11 * a02);
    if (determinant == 0)
    {
        throw std::invalid_argument(
            "shifted-level sign requires a full-rank quadratic matrix"
        );
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
    const Rational numerator =
        4 * Rational(coefficients(0)) * determinant - quadraticForm;
    if (numerator == 0)
    {
        return 0;
    }
    return (numerator > 0) == (determinant > 0) ? 1 : -1;
}


Edge undirectedEdge(int first, int second)
{
    return std::minmax(first, second);
}

Eigen::Vector3d refineTruncationCrossing(
    const spin::QuadricCoefficients10& quadric,
    const Eigen::Vector3d& firstNormal,
    const Eigen::Vector3d& secondNormal,
    double radius,
    const Eigen::Vector3d& initialPoint,
    bool allowSphereReseed = true
)
{
    const Eigen::Vector3d edgeNormal =
        firstNormal.cross(secondNormal).normalized();
    const Eigen::Matrix3d quadratic = quadric.quadraticMatrix();
    Eigen::Vector3d point = initialPoint;
    const double radiusScale = std::max(1.0, radius);

    const auto residual = [&](const Eigen::Vector3d& x) {
        Eigen::Vector3d value;
        const double squaredScale = std::max(1.0, x.squaredNorm());
        const Eigen::Vector3d gradient = quadric.gradient(x);
        const double gradientScale = std::max(1.0, gradient.norm());
        value <<
            quadric.evaluate(x) / squaredScale,
            (x.squaredNorm() - radius * radius) / (radiusScale * radiusScale),
            edgeNormal.dot(gradient) / gradientScale;
        return value;
    };

    for (int iteration = 0; iteration < 32; ++iteration)
    {
        const Eigen::Vector3d current = residual(point);
        if (std::abs(point.norm() - radius) <= 1e-10 * radiusScale &&
            std::abs(current[0]) <= 1e-12 && std::abs(current[2]) <= 1e-12)
        {
            // These geometric residual guards admit coordinate noise of
            // order 1e-9 at distant caps. Differentiating the implicit root
            // while returning that early seed makes the actual discrete map
            // fail finite differences. Polish the same root in extended
            // arithmetic, then reapply the unchanged geometric guards.
            using Real = long double;
            using Vector = Eigen::Matrix<Real, 3, 1>;
            using Matrix = Eigen::Matrix<Real, 3, 3>;
            const Matrix A = quadratic.cast<Real>();
            const Vector b = quadric.linearVector().cast<Real>();
            const Vector edge = firstNormal.cast<Real>().cross(
                secondNormal.cast<Real>()).normalized();
            Vector polished = point.cast<Real>();
            for (int polish = 0; polish < 8; ++polish)
            {
                const Vector gradient = b + 2 * A * polished;
                Vector value;
                value << Real(quadric.constant()) + b.dot(polished) +
                             polished.dot(A * polished),
                    polished.squaredNorm() - Real(radius) * Real(radius),
                    edge.dot(gradient);
                Matrix system;
                system.row(0) = gradient.transpose();
                system.row(1) = 2 * polished.transpose();
                system.row(2) = 2 * edge.transpose() * A;
                const auto factor = system.fullPivLu();
                if (!factor.isInvertible())
                    throw std::runtime_error("singular truncation root precision refinement");
                const Vector correction = factor.solve(-value);
                if (!correction.allFinite())
                    throw std::runtime_error("nonfinite truncation root precision correction");
                polished += correction;
                if (!polished.allFinite())
                    throw std::runtime_error("nonfinite truncation root precision result");
                if (correction.norm() <= 8 * std::numeric_limits<Real>::epsilon() *
                        std::max(Real(1), polished.norm()))
                {
                    break;
                }
            }
            const Eigen::Vector3d rounded = polished.cast<double>();
            const Eigen::Vector3d verified = residual(rounded);
            // A correction at the arithmetic floor can oscillate rather than
            // become stationary. The bounded polishing loop is not a new
            // nonlinear convergence gate: the original residual checks below
            // still decide geometric admission, and derivative consumers keep
            // their separate conditioning and consistency checks.
            if (!rounded.allFinite() || !verified.allFinite() ||
                (rounded - point).norm() > 1e-6 * radiusScale ||
                std::abs(rounded.norm() - radius) > 1e-10 * radiusScale ||
                std::abs(verified[0]) > 1e-12 || std::abs(verified[2]) > 1e-12)
                throw std::runtime_error(
                    "hyperboloid truncation root precision refinement failed");
            return rounded;
        }

        const Eigen::Vector3d gradient = quadric.gradient(point);
        Eigen::Matrix3d jacobian;
        jacobian.row(0) = gradient.transpose() /
            std::max(1.0, point.squaredNorm());
        jacobian.row(1) = 2.0 * point.transpose() /
            (radiusScale * radiusScale);
        jacobian.row(2) = 2.0 * edgeNormal.transpose() * quadratic /
            std::max(1.0, gradient.norm());
        if (!jacobian.allFinite() || !current.allFinite())
        {
            break;
        }
        const Eigen::Vector3d step = jacobian.fullPivLu().solve(-current);
        if (!step.allFinite())
        {
            break;
        }

        const double oldNorm = current.lpNorm<Eigen::Infinity>();
        bool accepted = false;
        double scale = 1.0;
        for (int lineSearch = 0; lineSearch < 24; ++lineSearch)
        {
            const Eigen::Vector3d candidate = point + scale * step;
            const Eigen::Vector3d candidateResidual = residual(candidate);
            if (candidate.allFinite() && candidateResidual.allFinite() &&
                candidateResidual.lpNorm<Eigen::Infinity>() < oldNorm)
            {
                point = candidate;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if (!accepted)
        {
            break;
        }
    }

    // Near an ideal normal, the normal-edge search can supply a tiny interior
    // point, for which the sphere equation defeats the original local Newton
    // initialization. Keep that solve first; only after its failure attempt
    // one sphere-scale seed. Both solves retain the same residual guards, and
    // source-motion consumers still require their exact root certificate.
    const double failedNorm = point.norm();
    if (allowSphereReseed && std::isfinite(failedNorm) && failedNorm > 0.0 &&
        failedNorm < 0.5 * radius)
    {
        return refineTruncationCrossing(
            quadric, firstNormal, secondNormal, radius,
            point * (radius / failedNorm), false
        );
    }

    throw std::runtime_error(
        "hyperboloid truncation crossing could not satisfy the certified "
        "quadric/sphere/normal-plane constraints"
    );
}


SphereMesh icosphere(int subdivisions)
{
    return spin::detail::fixedNormalSphere(subdivisions);
}


double automaticRadius(const Eigen::MatrixXd& vertices)
{
    if (vertices.rows() == 0 || vertices.cols() != 3 || !vertices.allFinite())
    {
        throw std::invalid_argument(
            "quadric triangulation requires finite design-domain vertices"
        );
    }
    const Eigen::Vector3d minimum = vertices.colwise().minCoeff();
    const Eigen::Vector3d maximum = vertices.colwise().maxCoeff();
    const double diameter = (maximum - minimum).norm();
    double maximumNorm = 0.0;
    for (Eigen::Index row = 0; row < vertices.rows(); ++row)
    {
        maximumNorm = std::max(maximumNorm, vertices.row(row).norm());
    }
    return std::max({1.0, 2.0 * diameter, 2.0 * maximumNorm});
}


spin::QuadricSurfaceType classifyAuthority(
    const spin::QuadricCoefficientAuthority10& authority,
    const spin::QuadricCoefficients10& display,
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>* eigen)
{
    using Q=CGAL::Gmpq;
    std::array<Q,10> q;
    for(int i=0;i<10;++i)q[i]=Q(authority.geometry()[i]);
    Q a[3][3]={{q[4],q[9],q[8]},{q[9],q[5],q[7]},{q[8],q[7],q[6]}};
    Q work[3][3];for(int i=0;i<3;++i)for(int j=0;j<3;++j)work[i][j]=a[i][j];
    int positive=0,negative=0;
    // Exact symmetric congruence: 1x1 pivots, or a zero-diagonal 2x2
    // block with one positive and one negative eigenvalue. No eigenvector
    // rounding determines the authoritative signature.
    for(int k=0;k<3;) {
        int pivot=k;while(pivot<3 && work[pivot][pivot]==0)++pivot;
        if(pivot<3) {
            if(pivot!=k) {
                for(int j=0;j<3;++j)std::swap(work[k][j],work[pivot][j]);
                for(int i=0;i<3;++i)std::swap(work[i][k],work[i][pivot]);
            }
            const Q d=work[k][k];if(d>0)++positive;else ++negative;
            for(int i=k+1;i<3;++i)for(int j=k+1;j<3;++j)
                work[i][j]-=work[i][k]*work[k][j]/d;
            ++k;
        } else {
            int i=k,j=k+1;
            while(i<3) {
                j=i+1;while(j<3 && work[i][j]==0)++j;
                if(j<3)break;
                ++i;
            }
            if(i==3)throw std::invalid_argument("coefficient authority requires exact full rank");
            if(i!=k) {
                for(int c=0;c<3;++c)std::swap(work[k][c],work[i][c]);
                for(int r=0;r<3;++r)std::swap(work[r][k],work[r][i]);
            }
            if(j!=k+1) {
                for(int c=0;c<3;++c)std::swap(work[k+1][c],work[j][c]);
                for(int r=0;r<3;++r)std::swap(work[r][k+1],work[r][j]);
            }
            const Q d=work[k][k+1];++positive;++negative;
            for(int r=k+2;r<3;++r)for(int c=k+2;c<3;++c)
                work[r][c]-=(work[r][k]*work[k+1][c]+work[r][k+1]*work[k][c])/d;
            k+=2;
        }
    }
    const Q det=a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-
        a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
    if(det==0)throw std::invalid_argument("coefficient authority requires exact full rank");
    Q shifted=q[0];
    for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
        const Q adj=a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]-
                    a[(j+1)%3][(i+2)%3]*a[(j+2)%3][(i+1)%3];
        shifted-=q[i+1]*adj*q[j+1]/(4*det);
    }
    eigen->compute(display.quadraticMatrix());
    const double tolerance=1e-12*std::max(1.,display.quadraticMatrix().norm());
    if(eigen->info()!=Eigen::Success || (eigen->eigenvalues().array().abs()<=tolerance).any() ||
       (eigen->eigenvalues().array()>0).count()!=positive)
        throw std::invalid_argument("authority signature cannot be admitted by the existing numerical rank gate");
    if(positive==3 || negative==3) {
        if((positive==3 && shifted>=0) || (negative==3 && shifted<=0))
            throw std::invalid_argument("authoritative ellipsoid has no real surface");
        return spin::QuadricSurfaceType::Ellipsoid;
    }
    if(shifted==0)throw std::invalid_argument("authoritative cone is outside the regular full-rank chart");
    return spin::QuadricSurfaceType::Hyperboloid;
}

spin::QuadricSurfaceType classify(
    const spin::QuadricCoefficients10& quadric,
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>* eigensolver
)
{
    const Eigen::Matrix3d a = quadric.quadraticMatrix();
    eigensolver->compute(a);
    if (eigensolver->info() != Eigen::Success)
    {
        throw std::invalid_argument("quadric eigen-decomposition failed");
    }
    const double scale = std::max(1.0, a.norm());
    const double tolerance = 1e-12 * scale;
    int rank = 0;
    for (int i = 0; i < 3; ++i)
    {
        rank += std::abs(eigensolver->eigenvalues()(i)) > tolerance ? 1 : 0;
    }
    if (rank == 0 && quadric.linearVector().norm() > tolerance)
    {
        return spin::QuadricSurfaceType::Plane;
    }
    if (rank == 1)
    {
        int curvedIndex = 0;
        eigensolver->eigenvalues().cwiseAbs().maxCoeff(&curvedIndex);
        const Eigen::Vector3d curvedAxis =
            eigensolver->eigenvectors().col(curvedIndex);
        const Eigen::Vector3d linear = quadric.linearVector();
        const Eigen::Vector3d transverseLinear =
            linear - curvedAxis * curvedAxis.dot(linear);
        if (transverseLinear.norm() > tolerance)
        {
            return spin::QuadricSurfaceType::ParabolicCylinder;
        }
        throw std::invalid_argument(
            "unsupported rank-one quadric with no transverse linear term "
            "(parallel, repeated, or empty planes)"
        );
    }
    if (rank == 2)
    {
        int nullIndex = 0;
        eigensolver->eigenvalues().cwiseAbs().minCoeff(&nullIndex);
        const Eigen::Vector3d nullVector =
            eigensolver->eigenvectors().col(nullIndex);
        if (std::abs(nullVector.dot(quadric.linearVector())) > tolerance)
        {
            return spin::QuadricSurfaceType::Paraboloid;
        }
        std::array<int, 2> nonzero{};
        int cursor = 0;
        for (int index = 0; index < 3; ++index)
        {
            if (index != nullIndex)
            {
                nonzero[cursor++] = index;
            }
        }
        const double first = eigensolver->eigenvalues()(nonzero[0]);
        const double second = eigensolver->eigenvalues()(nonzero[1]);
        if (first * second > 0.0)
        {
            return spin::QuadricSurfaceType::EllipticCylinder;
        }
        const double shiftedConstant = quadric.constant() -
            quadric.linearVector().dot(
                (eigensolver->eigenvectors().col(nonzero[0]) *
                     (eigensolver->eigenvectors().col(nonzero[0]).dot(
                          quadric.linearVector()) / first) +
                 eigensolver->eigenvectors().col(nonzero[1]) *
                     (eigensolver->eigenvectors().col(nonzero[1]).dot(
                          quadric.linearVector()) / second)) / 4.0);
        if (std::abs(shiftedConstant) <= tolerance)
        {
            return spin::QuadricSurfaceType::IntersectingPlanes;
        }
        return spin::QuadricSurfaceType::HyperbolicCylinder;
    }
    if (rank != 3)
    {
        throw std::invalid_argument(
            "unsupported constant/zero quadric with no nonzero quadratic or linear term"
        );
    }

    const Eigen::Vector3d eigenvalues = eigensolver->eigenvalues();
    const bool positive = (eigenvalues.array() > tolerance).all();
    const bool negative = (eigenvalues.array() < -tolerance).all();
    if (positive || negative)
    {
        const int centeredLevelSign = shiftedLevelSign(quadric);
        const bool nonempty =
            (positive && centeredLevelSign < 0) ||
            (negative && centeredLevelSign > 0);
        if (!nonempty)
        {
            throw std::invalid_argument("ellipsoid has no real surface");
        }
        return spin::QuadricSurfaceType::Ellipsoid;
    }

    // A full-rank indefinite quadratic whose translated level is zero is a
    // cone, not a regular one-/two-sheet hyperboloid. The inverse Gauss map
    // used below is singular at its apex, so reject this stratum here rather
    // than letting it masquerade as a hyperboloid and fail later while
    // searching for truncation loops.
    if (shiftedLevelSign(quadric) == 0)
    {
        return spin::QuadricSurfaceType::Cone;
    }
    return spin::QuadricSurfaceType::Hyperboloid;
}


void orientInterfaceFace(
    const spin::QuadricCoefficients10& quadric,
    const std::vector<Eigen::Vector3d>& vertices,
    Face* face
)
{
    const Eigen::Vector3d& a = vertices[(*face)[0]];
    const Eigen::Vector3d& b = vertices[(*face)[1]];
    const Eigen::Vector3d& c = vertices[(*face)[2]];
    const Eigen::Vector3d normal = (b - a).cross(c - a);
    const Eigen::Vector3d desired =
        -quadric.gradient((a + b + c) / 3.0);
    if (normal.dot(desired) < 0.0)
    {
        std::swap((*face)[1], (*face)[2]);
    }
}

std::vector<std::vector<int>> boundaryLoops(
    const std::vector<Face>& faces
)
{
    std::map<Edge, int> counts;
    std::map<Edge, std::pair<int, int>> directions;
    for (const Face& face : faces)
    {
        for (int corner = 0; corner < 3; ++corner)
        {
            const int first = face[corner];
            const int second = face[(corner + 1) % 3];
            const Edge key = undirectedEdge(first, second);
            ++counts[key];
            directions[key] = {first, second};
        }
    }
    std::map<int, std::vector<int>> outgoing;
    std::map<int, int> incoming;
    for (const auto& entry : counts)
    {
        if (entry.second == 1)
        {
            const std::pair<int, int> directed = directions.at(entry.first);
            ++incoming[directed.second];
            outgoing[directed.first].push_back(directed.second);
        }
    }
    std::map<int, int> next;
    for (const auto& [vertex, destinations] : outgoing)
    {
        if (destinations.size() != 1)
        {
            throw BoundaryLoopTopologyError(
                "quadric truncation boundary has " +
                std::to_string(destinations.size()) + " outgoing edges at vertex " +
                std::to_string(vertex) + " and " +
                std::to_string(incoming[vertex]) + " incoming edges; targets=" +
                std::to_string(destinations.front()) +
                (destinations.size() > 1 ? "," + std::to_string(destinations[1]) : ""));
        }
        next.emplace(vertex, destinations.front());
    }
    std::vector<std::vector<int>> loops;
    while (!next.empty())
    {
        std::vector<int> loop;
        const int start = next.begin()->first;
        int current = start;
        do
        {
            loop.push_back(current);
            const auto found = next.find(current);
            if (found == next.end())
            {
                throw BoundaryLoopTopologyError("quadric truncation loop is open");
            }
            current = found->second;
            next.erase(found);
        } while (current != start);
        if (loop.size() < 3)
        {
            throw BoundaryLoopTopologyError("quadric truncation loop is degenerate");
        }
        loops.push_back(std::move(loop));
    }
    return loops;
}


int idealNormalEdgeIntersections(
    const SphereMesh& sphere,
    const Eigen::Matrix3d& inverseA
)
{
    std::set<Edge> edges;
    for (const Face& face : sphere.faces)
    {
        edges.insert(undirectedEdge(face[0], face[1]));
        edges.insert(undirectedEdge(face[1], face[2]));
        edges.insert(undirectedEdge(face[2], face[0]));
    }
    int count = 0;
    for (const Edge& edge : edges)
    {
        const Eigen::Vector3d& first = sphere.vertices[edge.first];
        const Eigen::Vector3d difference =
            sphere.vertices[edge.second] - first;
        const double quadratic = difference.dot(inverseA * difference);
        const double linear = 2.0 * first.dot(inverseA * difference);
        const double constant = first.dot(inverseA * first);
        if (std::abs(quadratic) <= 1e-15)
        {
            if (std::abs(linear) > 1e-15)
            {
                const double root = -constant / linear;
                count += root > 0.0 && root < 1.0 ? 1 : 0;
            }
            continue;
        }
        const double discriminant =
            linear * linear - 4.0 * quadratic * constant;
        if (discriminant < 0.0)
        {
            continue;
        }
        const double squareRoot = std::sqrt(std::max(0.0, discriminant));
        const double firstRoot = (-linear - squareRoot) / (2.0 * quadratic);
        const double secondRoot = (-linear + squareRoot) / (2.0 * quadratic);
        count += firstRoot > 0.0 && firstRoot < 1.0 ? 1 : 0;
        if (squareRoot > 1e-15)
        {
            count += secondRoot > 0.0 && secondRoot < 1.0 ? 1 : 0;
        }
    }
    return count;
}


spin::TriangleMesh matrixMesh(
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<Face>& faces
)
{
    spin::TriangleMesh result;
    result.vertices.resize(static_cast<Eigen::Index>(vertices.size()), 3);
    for (std::size_t row = 0; row < vertices.size(); ++row)
    {
        result.vertices.row(static_cast<Eigen::Index>(row)) = vertices[row];
    }
    result.faces.resize(static_cast<Eigen::Index>(faces.size()), 3);
    for (std::size_t row = 0; row < faces.size(); ++row)
    {
        result.faces.row(static_cast<Eigen::Index>(row)) <<
            faces[row][0], faces[row][1], faces[row][2];
    }
    return result;
}


void orientComponentsOutward(spin::TriangleMesh* mesh)
{
    Eigen::MatrixXi orientedFaces;
    Eigen::VectorXi faceComponents;
    igl::bfs_orient(mesh->faces, orientedFaces, faceComponents);
    if (orientedFaces.rows() != mesh->faces.rows() ||
        faceComponents.size() != mesh->faces.rows() ||
        faceComponents.size() == 0)
    {
        throw std::runtime_error(
            "quadric boundary orientation could not assign every face component");
    }
    mesh->faces = std::move(orientedFaces);
    const int componentCount = faceComponents.maxCoeff() + 1;
    std::vector<double> componentSixVolumes(
        static_cast<std::size_t>(componentCount), 0.0);
    for (Eigen::Index face = 0; face < mesh->faces.rows(); ++face)
    {
        const Eigen::Vector3d a = mesh->vertices.row(mesh->faces(face, 0));
        const Eigen::Vector3d b = mesh->vertices.row(mesh->faces(face, 1));
        const Eigen::Vector3d c = mesh->vertices.row(mesh->faces(face, 2));
        componentSixVolumes[static_cast<std::size_t>(faceComponents(face))] +=
            a.dot(b.cross(c));
    }
    for (int component = 0; component < componentCount; ++component)
    {
        const double sixVolume =
            componentSixVolumes[static_cast<std::size_t>(component)];
        if (!std::isfinite(sixVolume) || sixVolume == 0.0)
        {
            throw std::runtime_error(
                "quadric boundary component has zero or nonfinite signed volume");
        }
        if (sixVolume < 0.0)
        {
            for (Eigen::Index face = 0; face < mesh->faces.rows(); ++face)
            {
                if (faceComponents(face) == component)
                    std::swap(mesh->faces(face, 1), mesh->faces(face, 2));
            }
        }
    }
    if (!igl::piecewise_constant_winding_number(mesh->faces))
    {
        throw std::runtime_error(
            "quadric boundary orientation normalization did not produce PWN faces");
    }
}


spin::QuadricTriangulation closeExtrudedInterface(
    spin::QuadricSurfaceType type,
    const std::vector<Eigen::Vector3d>& interfaceVertices,
    const std::vector<Face>& interfaceFaces,
    const Eigen::Vector3d& extrusion,
    double radius
)
{
    const std::vector<std::vector<int>> loops =
        boundaryLoops(interfaceFaces);
    if (loops.size() != 1)
    {
        throw std::runtime_error(
            "plane/paraboloid interface must have one boundary loop"
        );
    }
    std::vector<Eigen::Vector3d> vertices = interfaceVertices;
    const int offset = static_cast<int>(vertices.size());
    for (const Eigen::Vector3d& vertex : interfaceVertices)
    {
        vertices.push_back(vertex + extrusion);
    }
    std::vector<Face> faces = interfaceFaces;
    for (const Face& face : interfaceFaces)
    {
        faces.push_back({
            face[2] + offset,
            face[1] + offset,
            face[0] + offset
        });
    }
    const std::vector<int>& loop = loops.front();
    for (std::size_t edge = 0; edge < loop.size(); ++edge)
    {
        const int first = loop[edge];
        const int second = loop[(edge + 1) % loop.size()];
        faces.push_back({second, first, first + offset});
        faces.push_back({second, first + offset, second + offset});
    }

    spin::QuadricTriangulation result;
    result.type = type;
    result.closedBoundary = matrixMesh(vertices, faces);
    orientComponentsOutward(&result.closedBoundary);
    result.interfaceVertexCount = offset;
    result.interfaceFaceCount = static_cast<int>(interfaceFaces.size());
    result.capFaceCount =
        static_cast<int>(faces.size() - interfaceFaces.size());
    result.boundaryLoopCount = 1;
    result.truncationRadius = radius;
    result.faceIsCap = Eigen::VectorXi::Ones(
        static_cast<Eigen::Index>(faces.size())
    );
    result.faceIsCap.head(result.interfaceFaceCount).setZero();
    return result;
}


std::string triangulationHash(const spin::QuadricTriangulation& result,
                              bool algebraicFiniteMap=false)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << "QuadricTriangulation-v1;"
           << spin::quadricSurfaceTypeToken(result.type) << ';'
           << result.interfaceVertexCount << ';'
           << result.interfaceFaceCount << ';'
           << result.capFaceCount << ';'
           << result.boundaryLoopCount << ';'
           << result.idealNormalEdgeIntersectionCount << ';'
           << std::hexfloat << result.truncationRadius << ';';
    if(result.paperProjectiveHyperboloid)stream<<"paper-projective-triangles-v1;";
    if(std::any_of(result.sourceProjectiveRayByVertex.begin(),result.sourceProjectiveRayByVertex.end(),
        [](const auto& ray){return ray[0]==-2;}))stream<<"remote-planar-cap-extrusion-v1;";
    if(std::any_of(result.sourceProjectiveRayByVertex.begin(),result.sourceProjectiveRayByVertex.end(),
        [](const auto& ray){return ray[0]==-4;}))stream<<"remote-radial-cap-v1;";
    if(algebraicFiniteMap)stream<<"finite-map=algebraic-v1;";
    if(result.fixedNormalParaboloid)stream<<"finite-normal-paraboloid-v1;cosine-cutoff=0.01;";
    if (!result.exactVertices.empty())
    {
        stream << "source-coordinates=rational-enclosed-v1;absolute-width=2^-384;";
        if(result.compactRationalSourceCoordinates)
            stream << "representative=enclosed-dyadic-1024-v1;";
        for (std::size_t vertex = 0; vertex < result.exactVertices.size(); ++vertex)
            for (int axis = 0; axis < 3; ++axis)
                stream << result.exactVertices[vertex][axis] << ';'
                       << result.coordinateLower[vertex][axis] << ';'
                       << result.coordinateUpper[vertex][axis] << ';';
    }
    if (result.initialIcosphereSubdivisions >= 0 &&
        result.effectiveIcosphereSubdivisions >
            result.initialIcosphereSubdivisions)
    {
        stream << "adaptive-normal-resolution-v1;"
               << result.initialIcosphereSubdivisions << ';'
               << result.effectiveIcosphereSubdivisions << ';';
    }
    for (Eigen::Index row = 0; row < result.closedBoundary.vertices.rows(); ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            stream << result.closedBoundary.vertices(row, column) << ';';
        }
    }
    for (Eigen::Index row = 0; row < result.closedBoundary.faces.rows(); ++row)
    {
        stream << result.closedBoundary.faces(row, 0) << ';'
               << result.closedBoundary.faces(row, 1) << ';'
               << result.closedBoundary.faces(row, 2) << ';'
               << result.faceIsCap(row) << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}

std::string rationalConstructionIdentity(const spin::QuadricCoefficients10& q,
                                        const spin::QuadricTriangulation& t,
                                        const spin::QuadricCoefficientAuthority10* authority=nullptr) {
    const auto count=static_cast<std::size_t>(t.closedBoundary.vertices.rows());
    if(count<4 || t.closedBoundary.vertices.cols()!=3 ||
       t.closedBoundary.faces.rows()<4 || t.closedBoundary.faces.cols()!=3 ||
       t.faceIsCap.size()!=t.closedBoundary.faces.rows() ||
       t.exactVertices.size()!=count || t.coordinateLower.size()!=count ||
       t.coordinateUpper.size()!=count || t.fixedNormals.cols()!=3 ||
       !t.closedBoundary.vertices.allFinite() || !t.fixedNormals.allFinite())return {};
    const auto geometry=triangulationHash(t);
    if(t.hash!=geometry)return {};
    std::ostringstream stream;stream.imbue(std::locale::classic());
    stream<<"NativeRationalSourceConstruction-v1;"<<q.hash()<<';'<<geometry<<';'
          <<t.initialIcosphereSubdivisions<<';'<<t.effectiveIcosphereSubdivisions<<';'
          <<std::hexfloat;
    if(authority) {
        stream<<authority->policy<<';'<<authority->precisionBits<<';';
        for(const auto& coefficient:authority->geometry())stream<<coefficient<<';';
    }
    for(int i=0;i<t.fixedNormals.rows();++i)for(int a=0;a<3;++a)
        stream<<t.fixedNormals(i,a)<<';';
    for(const auto& edge:t.sourceNormalEdgeByVertex)stream<<edge[0]<<','<<edge[1]<<';';
    stream<<"normal-vertices;";
    for(int vertex:t.sourceNormalVertexByVertex)stream<<vertex<<';';
    stream<<"projective-rays;";
    for(const auto& ray:t.sourceProjectiveRayByVertex)for(int index:ray)stream<<index<<';';
    return manufacturing::computeContentSha256(stream.str());
}

// The sign of this homogeneous quadratic is the sign of w^2 in the
// supplemental inverse Gauss map. Exact extrema on the planar normal triangle
// detect components invisible at its three vertices and boundary crossings.
struct ExactNormalConic
{
    using R=CGAL::Gmpq;
    using V=std::array<R,3>;
    R h[3][3];
    int curvatureSign=0;
    explicit ExactNormalConic(const spin::QuadricCoefficients10& q) {
        const auto& p=q.values();
        const R a[3][3]={{R(p(4)),R(p(9)),R(p(8))},
                        {R(p(9)),R(p(5)),R(p(7))},
                        {R(p(8)),R(p(7)),R(p(6))}};
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            h[i][j]=a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]-
                    a[(j+1)%3][(i+2)%3]*a[(j+2)%3][(i+1)%3];
        R det(0);for(int j=0;j<3;++j)det+=a[0][j]*h[j][0];
        R denominator=-4*R(p(0))*det;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            denominator+=R(p(i+1))*h[i][j]*R(p(j+1));
        if(det==0 || denominator==0)throw std::runtime_error("singular paper normal conic");
        curvatureSign=denominator>0?1:-1;
        if(denominator<0)for(auto& row:h)for(auto& value:row)value=-value;
    }
    static V exact(const Eigen::Vector3d& p) {return {R(p.x()),R(p.y()),R(p.z())};}
    R bilinear(const V& a,const V& b) const {
        R value(0);for(int i=0;i<3;++i)for(int j=0;j<3;++j)value+=a[i]*h[i][j]*b[j];
        return value;
    }
    int sign(const Eigen::Vector3d& p) const {
        const auto v=exact(p);const R value=bilinear(v,v);return value>0?1:value<0?-1:0;
    }
    bool twoPositiveCones() const {
        const R determinant=
            h[0][0]*(h[1][1]*h[2][2]-h[1][2]*h[2][1])-
            h[0][1]*(h[1][0]*h[2][2]-h[1][2]*h[2][0])+
            h[0][2]*(h[1][0]*h[2][1]-h[1][1]*h[2][0]);
        // For a nonsingular symmetric 3x3 matrix det>0 permits only (+++)
        // or (+--). Sylvester's criterion separates the former exactly.
        return determinant>0 && !(h[0][0]>0 &&
            h[0][0]*h[1][1]-h[0][1]*h[1][0]>0);
    }
    int coneSide(const Eigen::Vector3d& p,const Eigen::Vector3d& positiveReference) const {
        const R value=bilinear(exact(p),exact(positiveReference));
        return value>0?1:value<0?-1:0;
    }
    bool uniform(const std::array<Eigen::Vector3d,3>& points,bool inside) const {
        std::array<V,3> n{exact(points[0]),exact(points[1]),exact(points[2])};
        std::vector<R> values;
        for(const auto& v:n)values.push_back(bilinear(v,v));
        for(int edge=0;edge<3;++edge) {
            V d;for(int i=0;i<3;++i)d[i]=n[(edge+1)%3][i]-n[edge][i];
            const R a=bilinear(d,d),b=bilinear(n[edge],d);
            if(a!=0) {
                const R t=-b/a;
                if(t>0 && t<1)values.push_back(values[edge]+2*t*b+t*t*a);
            }
        }
        V u,v;for(int i=0;i<3;++i){u[i]=n[1][i]-n[0][i];v[i]=n[2][i]-n[0][i];}
        const R a=bilinear(u,u),b=bilinear(u,v),c=bilinear(v,v);
        const R d=bilinear(n[0],u),e=bilinear(n[0],v),det=a*c-b*b;
        if(det!=0) {
            const R x=(b*e-c*d)/det,y=(b*d-a*e)/det;
            if(x>0 && y>0 && x+y<1)
                values.push_back(values[0]+2*x*d+2*y*e+x*x*a+2*x*y*b+y*y*c);
        }
        // A singular Hessian's interior stationary set reaches the boundary;
        // its value is therefore already represented by a boundary extremum.
        return std::all_of(values.begin(),values.end(),[&](const R& value) {
            return inside?value>0:value<0;
        });
    }
};

std::vector<Face> triangulateProjectedCap(const std::vector<int>& loop,
    const std::vector<Eigen::Vector3d>& vertices,const Eigen::Vector3d& axis)
{
    using Q=CGAL::Gmpq;
    Eigen::Index omitted=0;axis.cwiseAbs().maxCoeff(&omitted);
    std::array<int,2> axes{};int cursor=0;
    for(int i=0;i<3;++i)if(i!=omitted)axes[cursor++]=i;
    const auto turn=[&](int a,int b,int c) {
        const Q ax(vertices[a](axes[0])),ay(vertices[a](axes[1]));
        const Q bx(vertices[b](axes[0])),by(vertices[b](axes[1]));
        const Q cx(vertices[c](axes[0])),cy(vertices[c](axes[1]));
        return (bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
    };
    Q area(0);
    for(std::size_t k=0;k<loop.size();++k) {
        const auto& a=vertices[loop[k]];const auto& b=vertices[loop[(k+1)%loop.size()]];
        area+=Q(a(axes[0]))*Q(b(axes[1]))-Q(a(axes[1]))*Q(b(axes[0]));
    }
    if(area==0)throw std::runtime_error("paper cap projection has zero exact area");
    const int sign=area>0?1:-1;
    std::vector<int> remaining=loop;std::vector<Face> triangles;
    while(remaining.size()>3) {
        bool found=false;
        for(std::size_t k=0;k<remaining.size();++k) {
            const int a=remaining[(k+remaining.size()-1)%remaining.size()],b=remaining[k],
                c=remaining[(k+1)%remaining.size()];
            if(sign*turn(a,b,c)<=0)continue;
            bool occupied=false;
            for(int point:remaining)if(point!=a && point!=b && point!=c)
                occupied|=sign*turn(a,b,point)>=0 && sign*turn(b,c,point)>=0 && sign*turn(c,a,point)>=0;
            if(occupied)continue;
            triangles.push_back({c,b,a});remaining.erase(remaining.begin()+k);found=true;break;
        }
        if(!found)throw std::runtime_error("paper cap projection has no certified ear");
    }
    triangles.push_back({remaining[2],remaining[1],remaining[0]});return triangles;
}

spin::QuadricTriangulation triangulateProjectiveHyperboloid(
    const spin::QuadricCoefficients10& q,const SphereMesh& sphere,
    const Eigen::MatrixXd& design,const spin::QuadricTriangulationSettings& settings,
    double radius)
{
    if(settings.maximumIcosphereSubdivisions!=-1 &&
       settings.maximumIcosphereSubdivisions!=settings.icosphereSubdivisions)
        throw std::invalid_argument("paper projective construction requires fixed normal discretization");
    const Eigen::Vector3d lo=design.colwise().minCoeff(),hi=design.colwise().maxCoeff();
    if(radius<=(hi-lo).norm() || radius<=design.rowwise().norm().maxCoeff())
        throw std::invalid_argument("paper initial radius must exceed design diameter and origin distance");
    const spin::QuadricInverseGaussMap map(q);
    const ExactNormalConic conic(q);
    struct Node {int normal=-1;Edge edge{-1,-1};int ordinal=-1;Eigen::Vector3d point;};
    std::vector<Node> nodes;
    std::vector<bool> inside;
    for(int i=0;i<int(sphere.vertices.size());++i) {
        const int sign=conic.sign(sphere.vertices[i]);
        if(sign==0)throw std::runtime_error("paper normal vertex is ideal: derivative branch event");
        inside.push_back(sign>0);
        const auto finite=settings.algebraicFiniteMap
            ? map.algebraicPoint(sphere.vertices[i]) : map.point(sphere.vertices[i]);
        if(inside.back()!=finite.valid())
            throw std::runtime_error("paper finite-normal status differs from exact conic membership");
        nodes.push_back({i,{-1,-1},-1,finite.valid()?*finite.point:Eigen::Vector3d::Zero()});
    }
    std::map<Edge,std::vector<int>> edgeNodes;
    for(const Face& face:sphere.faces)for(int k=0;k<3;++k) {
        const Edge edge=undirectedEdge(face[k],face[(k+1)%3]);
        if(edgeNodes.count(edge))continue;
        const auto ideal=map.idealEdgeDifferential(sphere.vertices[edge.first],sphere.vertices[edge.second]);
        if(!ideal.defined)throw std::runtime_error("paper ideal edge "+std::to_string(edge.first)+","+
            std::to_string(edge.second)+": "+ideal.diagnostic);
        auto& ids=edgeNodes[edge];
        for(int j=0;j<int(ideal.points.size());++j) {
            ids.push_back(int(nodes.size()));nodes.push_back({-1,edge,j,ideal.points[j].direction});
        }
    }
    std::vector<Face> projectiveTriangles;
    std::map<Edge,int> intervalAnchors;
    std::map<int,Eigen::Vector3d> intervalNormals;
    bool hasBeltMultipleCrossings=false;
    for(int f=0;f<int(sphere.faces.size());++f) {
        const auto& face=sphere.faces[f];
        int roots=0;bool singleCrossings=true;
        for(int k=0;k<3;++k) {
            const auto count=edgeNodes.at(undirectedEdge(face[k],face[(k+1)%3])).size();
            roots+=int(count);singleCrossings&=count<=1;
        }
        if(roots==0) {
            if(!conic.uniform({sphere.vertices[face[0]],sphere.vertices[face[1]],sphere.vertices[face[2]]},inside[face[0]]))
                throw std::runtime_error("UNSUPPORTED_PAPER_IDEAL_FACE: hidden conic at face "+std::to_string(f));
            if(inside[face[0]])projectiveTriangles.push_back(face);
            continue;
        }
        if(roots!=2 || !singleCrossings) {
            if(!conic.twoPositiveCones()) {
                hasBeltMultipleCrossings=true;
                // The finite normal domain can instead be the connected
                // positive belt, complementary to two convex negative cones.
                // Retain its ordered projective boundary, including both roots
                // of an edge whose endpoints are finite. Never discard the
                // ideal interval just because no normal corner lies within it.
                const auto anchor=std::find(inside.begin(),inside.end(),false);
                if(anchor==inside.end())
                    throw std::runtime_error("UNSUPPORTED_PAPER_IDEAL_FACE: no negative cone reference");
                const auto& reference=sphere.vertices[std::distance(inside.begin(),anchor)];
                std::map<int,int> rootSide;
                std::map<int,int> rootsBySide;
                for(int k=0;k<3;++k) {
                    const Edge edge=undirectedEdge(face[k],face[(k+1)%3]);
                    const auto& ids=edgeNodes.at(edge);
                    if(ids.empty())continue;
                    Eigen::Vector3d negative;
                    if(ids.size()==1 && inside[edge.first]!=inside[edge.second])
                        negative=sphere.vertices[inside[edge.first]?edge.second:edge.first];
                    else if(ids.size()==2 && inside[edge.first] && inside[edge.second]) {
                        const auto ideal=map.idealEdgeDifferential(sphere.vertices[edge.first],sphere.vertices[edge.second]);
                        const double t=(ideal.points[0].parameter+ideal.points[1].parameter)/2;
                        negative=(1-t)*sphere.vertices[edge.first]+t*sphere.vertices[edge.second];
                        if(conic.sign(negative)>=0)
                            throw std::runtime_error("paper belt ideal interval has no negative interior");
                    } else throw std::runtime_error("UNSUPPORTED_PAPER_IDEAL_FACE: unanchored belt crossing");
                    const int side=conic.coneSide(negative,reference);
                    if(side==0)throw std::runtime_error("paper negative cone separator is singular");
                    for(int id:ids){rootSide[id]=side;++rootsBySide[side];}
                }
                std::map<int,int> next;
                std::set<int> incoming;
                for(int k=0;k<3;++k) {
                    const int a=face[k],b=face[(k+1)%3];
                    auto ids=edgeNodes.at(undirectedEdge(a,b));
                    if(a>b)std::reverse(ids.begin(),ids.end());
                    ids.insert(ids.begin(),a);ids.push_back(b);
                    bool positive=inside[a];
                    for(std::size_t j=0;j+1<ids.size();++j) {
                        if(positive){next.emplace(ids[j],ids[j+1]);incoming.insert(ids[j+1]);}
                        if(j+1<ids.size()-1)positive=!positive;
                    }
                    if(positive!=inside[b])throw std::runtime_error("paper belt interval membership mismatch");
                }
                std::vector<std::vector<int>> paths;
                for(const auto& entry:rootSide)if(next.count(entry.first) && !incoming.count(entry.first)) {
                    std::vector<int> path{entry.first};
                    while(next.count(path.back())) {
                        if(path.size()>nodes.size())throw std::runtime_error("paper belt boundary path cycles");
                        path.push_back(next.at(path.back()));
                    }
                    if(nodes[path.back()].normal>=0)
                        throw std::runtime_error("paper belt path does not end at an ideal root");
                    paths.push_back(std::move(path));
                }
                if(paths.empty())throw std::runtime_error("UNSUPPORTED_PAPER_IDEAL_FACE: no finite belt paths");
                std::vector<std::vector<int>> polygons;
                if(rootsBySide.size()==1) {
                    // Removing one convex negative patch leaves one positive
                    // component per positive boundary run; four/six crossings
                    // can therefore retain two/three distinct finite corners.
                    polygons=paths;
                } else {
                    for(const auto& entry:rootsBySide)if(entry.second!=2)
                        throw std::runtime_error("UNSUPPORTED_PAPER_IDEAL_FACE: ambiguous belt cone arcs");
                    std::set<int> used;
                    for(int start=0;start<int(paths.size());++start)if(!used.count(start)) {
                        std::vector<int> polygon;int current=start;
                        do {
                            if(!used.insert(current).second)throw std::runtime_error("paper belt arc connects different cycles");
                            polygon.insert(polygon.end(),paths[current].begin(),paths[current].end());
                            const int side=rootSide.at(paths[current].back());
                            int following=-1;
                            for(int j=0;j<int(paths.size());++j)if(rootSide.at(paths[j].front())==side) {
                                if(following!=-1)throw std::runtime_error("paper belt arc has multiple starts");
                                following=j;
                            }
                            if(following<0)throw std::runtime_error("paper belt arc has no matching start");
                            current=following;
                        } while(current!=start);
                        polygons.push_back(std::move(polygon));
                    }
                }
                for(auto& polygon:polygons) {
                    const auto finite=std::find_if(polygon.begin(),polygon.end(),
                        [&](int i){return nodes[i].normal>=0;});
                    if(polygon.size()<3 || finite==polygon.end())
                        throw std::runtime_error("UNSUPPORTED_PAPER_IDEAL_FACE: unanchored belt polygon");
                    std::rotate(polygon.begin(),finite,polygon.end());
                    for(std::size_t k=1;k+1<polygon.size();++k)
                        projectiveTriangles.push_back({polygon[0],polygon[k],polygon[k+1]});
                }
                continue;
            }
            // Each positive Lorentz cone intersects this planar normal
            // triangle in a convex patch. Its boundary roots must be paired
            // by cone, not by assuming only one root on each triangle edge.
            const auto anchor=std::find(inside.begin(),inside.end(),true);
            if(anchor==inside.end())
                throw std::runtime_error("UNSUPPORTED_PAPER_IDEAL_FACE: no positive cone reference");
            const auto& reference=sphere.vertices[std::distance(inside.begin(),anchor)];
            std::map<int,int> sides;
            for(int corner:face)if(inside[corner]) {
                const int side=conic.coneSide(sphere.vertices[corner],reference);
                if(side==0)throw std::runtime_error("paper positive cone separator is singular");
                sides[corner]=side;
            }
            for(int k=0;k<3;++k) {
                const Edge edge=undirectedEdge(face[k],face[(k+1)%3]);
                const auto& ids=edgeNodes.at(edge);
                if(ids.empty())continue;
                if(ids.size()==1 && inside[edge.first]!=inside[edge.second])
                    sides[ids[0]]=sides.at(inside[edge.first]?edge.first:edge.second);
                else if(ids.size()==2 && inside[edge.first] && inside[edge.second] &&
                        sides.at(edge.first)!=sides.at(edge.second)) {
                    // Stored roots are ordered from edge.first to edge.second.
                    sides[ids[0]]=sides.at(edge.first);
                    sides[ids[1]]=sides.at(edge.second);
                } else if(ids.size()==2 && !inside[edge.first] && !inside[edge.second]) {
                    const auto roots=map.idealEdgeDifferential(sphere.vertices[edge.first],sphere.vertices[edge.second]);
                    const double mid=(roots.points[0].parameter+roots.points[1].parameter)*0.5;
                    const Eigen::Vector3d normal=(1-mid)*sphere.vertices[edge.first]+mid*sphere.vertices[edge.second];
                    if(conic.sign(normal)<=0)throw std::runtime_error("paper double-root edge has no positive interior");
                    const int side=conic.coneSide(normal,reference);
                    if(side==0)throw std::runtime_error("paper interval cone separator is singular");
                    sides[ids[0]]=sides[ids[1]]=side;
                } else throw std::runtime_error(
                    "UNSUPPORTED_PAPER_IDEAL_FACE: unanchored cone crossings at face "+std::to_string(f));
            }
            for(int side:{-1,1}) {
                std::vector<int> polygon;
                for(int k=0;k<3;++k) {
                    const int a=face[k],b=face[(k+1)%3];
                    if(inside[a] && sides.at(a)==side)polygon.push_back(a);
                    auto ids=edgeNodes.at(undirectedEdge(a,b));
                    if(a>b)std::reverse(ids.begin(),ids.end());
                    for(int id:ids)if(sides.at(id)==side)polygon.push_back(id);
                }
                if(polygon.empty())continue;
                std::set<Edge> insertedIntervals;
                for(int position=0;position<int(polygon.size());++position) {
                    const auto a=nodes[polygon[position]],b=nodes[polygon[(position+1)%polygon.size()]];
                    if(a.normal>=0 || b.normal>=0 || a.edge!=b.edge || a.edge.first<0 ||
                       a.ordinal==b.ordinal || inside[a.edge.first] || inside[a.edge.second] ||
                       insertedIntervals.count(a.edge))continue;
                    insertedIntervals.insert(a.edge);
                    if(settings.rationalSourceCoordinates)
                        throw std::runtime_error("UNSUPPORTED_PAPER_INTERVAL_ANCHOR: rational enclosure not implemented");
                    // Insert a finite normal on a retained original edge interval.
                    // Its coefficient-dependent midpoint is a construction point,
                    // not a different normal grid or an answer-dependent sample.
                    const Edge selected=a.edge;
                    auto found=intervalAnchors.find(selected);
                    if(found==intervalAnchors.end()) {
                        const auto roots=map.idealEdgeDifferential(sphere.vertices[selected.first],sphere.vertices[selected.second]);
                        if(!roots.defined || roots.points.size()!=2)throw std::runtime_error("paper interval roots unavailable");
                        const double mid=(roots.points[0].parameter+roots.points[1].parameter)*0.5;
                        const Eigen::Vector3d normal=(1-mid)*sphere.vertices[selected.first]+mid*sphere.vertices[selected.second];
                        const auto point=settings.algebraicFiniteMap?map.algebraicPoint(normal):map.point(normal);
                        if(!point.valid())throw std::runtime_error("paper finite interval anchor: "+point.diagnostic);
                        const int id=int(nodes.size());nodes.push_back({id,selected,-2,*point.point});
                        intervalNormals.emplace(id,normal);found=intervalAnchors.emplace(selected,id).first;
                    }
                    polygon.insert(polygon.begin()+position+1,found->second);
                    ++position;
                    hasBeltMultipleCrossings=true; // same certified remote cap path
                }
                const auto finite=std::find_if(polygon.begin(),polygon.end(),
                    [&](int i){return nodes[i].normal>=0;});
                if(polygon.size()<3 || finite==polygon.end())
                    throw std::runtime_error("paper positive cone patch has no finite interval anchor");
                std::rotate(polygon.begin(),finite,polygon.end());
                for(std::size_t k=1;k+1<polygon.size();++k)
                    projectiveTriangles.push_back({polygon[0],polygon[k],polygon[k+1]});
            }
            continue;
        }
        std::vector<int> polygon;
        for(int k=0;k<3;++k) {
            const int a=face[k],b=face[(k+1)%3];
            const auto& ids=edgeNodes.at(undirectedEdge(a,b));
            if(inside[a])polygon.push_back(a);
            if(!ids.empty()) {
                if(inside[a]==inside[b])throw std::runtime_error("paper crossing membership mismatch");
                polygon.push_back(ids.front());
            }
        }
        auto finite=std::find_if(polygon.begin(),polygon.end(),[&](int i){return nodes[i].normal>=0;});
        if(finite==polygon.end())throw std::runtime_error("UNSUPPORTED_PAPER_IDEAL_FACE: no finite anchor");
        std::rotate(polygon.begin(),finite,polygon.end());
        for(std::size_t k=1;k+1<polygon.size();++k)
            projectiveTriangles.push_back({polygon[0],polygon[k],polygon[k+1]});
    }
    std::string reason="no certified paper closure";
    for(int doubling=0;doubling<=settings.maximumRadiusDoublings;++doubling,radius*=2) {
        if(!std::isfinite(radius))throw std::runtime_error("paper truncation radius overflow");
        bool allFiniteInside=true;
        for(int i=0;i<int(nodes.size());++i)if(i<int(inside.size())?inside[i]:nodes[i].normal>=0) {
            CGAL::Gmpq offset=-CGAL::Gmpq(radius)*CGAL::Gmpq(radius);
            for(int a=0;a<3;++a)offset+=CGAL::Gmpq(nodes[i].point(a))*CGAL::Gmpq(nodes[i].point(a));
            allFiniteInside&=offset<0;
        }
        if(!allFiniteInside){reason="finite mapped vertices not strictly within truncation sphere";continue;}
        std::vector<Eigen::Vector3d> vertices;
        std::vector<Face> faces;
        std::vector<int> normalIds;
        std::vector<std::array<int,4>> rays;
        std::map<int,int> finiteIds;
        std::map<std::pair<int,int>,int> rayIds;
        const auto addFinite=[&](int id) {
            auto found=finiteIds.find(id);if(found!=finiteIds.end())return found->second;
            const int index=int(vertices.size());vertices.push_back(nodes[id].point);
            normalIds.push_back(id);
            rays.push_back(intervalNormals.count(id)?std::array<int,4>{-3,nodes[id].edge.first,nodes[id].edge.second,0}:
                std::array<int,4>{-1,-1,-1,-1});finiteIds[id]=index;return index;
        };
        const auto addRay=[&](int finite,int ideal) {
            const auto key=std::make_pair(finite,ideal);
            auto found=rayIds.find(key);if(found!=rayIds.end())return found->second;
            const auto cut=spin::truncateQuadricProjectiveRay(nodes[finite].point,nodes[ideal].point,radius);
            if(!cut.defined)throw std::runtime_error("paper projective ray: "+cut.diagnostic);
            const int index=int(vertices.size());vertices.push_back(cut.point);normalIds.push_back(-1);
            rays.push_back({finite,nodes[ideal].edge.first,nodes[ideal].edge.second,nodes[ideal].ordinal});
            rayIds[key]=index;return index;
        };
        for(const auto& triangle:projectiveTriangles) {
            std::vector<int> polygon;
            for(int k=0;k<3;++k) {
                const int a=triangle[k],b=triangle[(k+1)%3];
                const bool fa=nodes[a].normal>=0,fb=nodes[b].normal>=0;
                if(fa)polygon.push_back(addFinite(a));
                if(fa!=fb)polygon.push_back(addRay(fa?a:b,fa?b:a));
            }
            for(std::size_t k=1;k+1<polygon.size();++k) {
                Face face{polygon[0],polygon[k],polygon[k+1]};
                // Inverse Gauss orientation is the sign of Gaussian curvature
                // (-4*c*det(A)+b^T*adj(A)*b). Use the fixed normal chart's
                // orientation globally. A coarse saddle chord can disagree
                // with the analytic gradient at its centroid, so per-chord
                // flips can reverse just one side of a shared edge.
                if(conic.curvatureSign>0)std::swap(face[1],face[2]);
                faces.push_back(face);
            }
        }
        std::map<Edge,std::vector<Edge>> incidence;
        for(const auto& face:faces)for(int k=0;k<3;++k)
            incidence[undirectedEdge(face[k],face[(k+1)%3])].push_back({face[k],face[(k+1)%3]});
        for(const auto& entry:incidence) {
            const auto& uses=entry.second;
            if(uses.size()>2 || (uses.size()==2 && uses[0]==uses[1]))
                throw std::runtime_error("paper projective source has inconsistent edge incidence");
        }
        const auto loops=boundaryLoops(faces);
        if(loops.size()!=2){reason="paper truncation requires exactly two boundary loops";continue;}
        if(hasBeltMultipleCrossings) {
            // A non-star-shaped truncation ring cannot use its mean as a cap
            // fan: that fan can cross the interface. Extend the same ring to
            // a remote plane, triangulate its projected polygon exactly, and
            // certify the complete represented source. No interface point,
            // fixed normal or design-domain operand changes here.
            const int interfaceVertices=int(vertices.size()),interfaceFaces=int(faces.size());
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> capEigen(q.quadraticMatrix());
            if(capEigen.info()!=Eigen::Success)
                throw std::runtime_error("paper radial cap eigensystem failed");
            const int capAxisIndex=capEigen.eigenvalues()(1)<0?2:0;
            const Eigen::Vector3d capCenter=q.quadraticMatrix().fullPivLu().solve(-0.5*q.linearVector());
            if(!capCenter.allFinite())throw std::runtime_error("paper radial cap center is nonfinite");
            bool exterior=true;
            for(int l=0;l<int(loops.size());++l) {
                const auto& loop=loops[l];
                Eigen::Vector3d area=Eigen::Vector3d::Zero(),mean=Eigen::Vector3d::Zero();
                for(std::size_t k=0;k<loop.size();++k) {
                    area+=vertices[loop[k]].cross(vertices[loop[(k+1)%loop.size()]]);
                    mean+=vertices[loop[k]];
                }
                if(!area.allFinite() || area.norm()==0)throw std::runtime_error("paper cap axis is singular");
                Eigen::Vector3d axis=capEigen.eigenvectors().col(capAxisIndex);
                if(axis.dot(mean-double(loop.size())*capCenter)<0)axis=-axis;
                using Q=CGAL::Gmpq;
                Q support(0);
                for(int a=0;a<3;++a)support+=Q(axis(a))*Q(axis(a)>=0?hi(a):lo(a));
                const auto beyond=[&](const Eigen::Vector3d& point) {
                    Q dot(0);for(int a=0;a<3;++a)dot+=Q(axis(a))*Q(point(a));
                    return dot>support;
                };
                std::vector<int> upper;
                for(int k=0;k<int(loop.size());++k) {
                    const int original=loop[k];
                    const Eigen::Vector3d delta=vertices[original]-capCenter;
                    const double denominator=axis.dot(delta);
                    if(!(denominator>0) || !std::isfinite(denominator))
                        throw std::runtime_error("paper radial cap projection is not separated");
                    const Eigen::Vector3d point=capCenter+delta*(2*radius/denominator);
                    exterior&=beyond(vertices[original]) && beyond(point);
                    upper.push_back(int(vertices.size()));vertices.push_back(point);
                    normalIds.push_back(-1);rays.push_back({-4,original,l,k});
                }
                for(std::size_t k=0;k<loop.size();++k) {
                    const int next=(k+1)%loop.size();
                    faces.push_back({loop[next],loop[k],upper[k]});
                    faces.push_back({loop[next],upper[k],upper[next]});
                }
                const auto cap=triangulateProjectedCap(upper,vertices,axis);
                faces.insert(faces.end(),cap.begin(),cap.end());
            }
            if(!exterior){reason="paper extruded caps not separated from design AABB";continue;}
            spin::QuadricTriangulation result;
            result.type=spin::QuadricSurfaceType::Hyperboloid;result.paperProjectiveHyperboloid=true;
            result.interfaceVertexCount=interfaceVertices;result.interfaceFaceCount=interfaceFaces;
            result.closedBoundary=matrixMesh(vertices,faces);orientComponentsOutward(&result.closedBoundary);
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
            const auto audit=verification::auditExactTriangleBoundary(result.closedBoundary.vertices,
                result.closedBoundary.faces,verification::ExactBoundaryOrientation::OutwardPositive);
            if(!audit.facts().boundaryValid){
                const auto& facts=audit.facts();
                std::ostringstream diagnostic;
                diagnostic << "paper extruded source failed exact boundary audit"
                    << ";closed=" << facts.closed
                    << ";edge_manifold=" << facts.edgeManifold
                    << ";vertex_manifold=" << facts.vertexManifold
                    << ";oriented=" << facts.consistentlyOriented
                    << ";component_orientations=" << facts.componentOrientationsValid
                    << ";degenerate_faces=" << facts.degenerateFaces.size()
                    << ";unreferenced_vertices=" << facts.unreferencedVertices.size()
                    << ";intersection_completed=" << facts.intersectionCompleted
                    << ";intersecting_pairs=" << facts.intersectingPairs.size()
                    << ";intersection_error=" << facts.intersectionError;
                if(!facts.intersectingPairs.empty())
                {
                    diagnostic << ";first_intersecting_pair="
                        << facts.intersectingPairs.front()[0] << ','
                        << facts.intersectingPairs.front()[1];
                    diagnostic << ";interface_face_count=" << interfaceFaces
                               << ";truncation_radius=" << radius;
                    int interfacePairs=0, mixedPairs=0, capPairs=0;
                    for(const auto& pair:facts.intersectingPairs) {
                        const int count=int(pair[0]<interfaceFaces)+int(pair[1]<interfaceFaces);
                        if(count==2)++interfacePairs;
                        else if(count==1)++mixedPairs;
                        else ++capPairs;
                    }
                    diagnostic << ";interface_pairs=" << interfacePairs
                               << ";mixed_pairs=" << mixedPairs << ";cap_pairs=" << capPairs;
                    diagnostic.precision(17);
                    for(const int face:facts.intersectingPairs.front()) {
                        diagnostic << ";face_" << face << '=';
                        for(int corner=0;corner<3;++corner) {
                            const int vertex=result.closedBoundary.faces(face,corner);
                            diagnostic << vertex << ':';
                            for(int axis=0;axis<3;++axis)
                                diagnostic << result.closedBoundary.vertices(vertex,axis) << ',';
                        }
                    }
                }
                if(!facts.degenerateFaces.empty())
                    diagnostic << ";first_degenerate_face=" << facts.degenerateFaces.front();
                // Preserve the actual rejected construction for deterministic
                // source-only diagnosis.  These fields never authorize a mesh.
                diagnostic.precision(17);
                diagnostic << ";witness_vertices=";
                for(int i=0;i<result.closedBoundary.vertices.rows();++i) {
                    for(int a=0;a<3;++a)diagnostic << result.closedBoundary.vertices(i,a) << ',';
                    diagnostic << '|';
                }
                diagnostic << ";witness_faces=";
                for(int i=0;i<result.closedBoundary.faces.rows();++i) {
                    for(int a=0;a<3;++a)diagnostic << result.closedBoundary.faces(i,a) << ',';
                    diagnostic << '|';
                }
                diagnostic << ";witness_loops=";
                for(const auto& loop:loops) {
                    for(int vertex:loop)diagnostic << vertex << ',';
                    diagnostic << '|';
                }
                diagnostic << ";witness_pairs=";
                for(const auto& pair:facts.intersectingPairs)
                    diagnostic << pair[0] << ',' << pair[1] << '|';
                reason=diagnostic.str();continue;
            }
#else
            throw std::runtime_error("paper extruded cap certification requires exact geometry checks");
#endif
            result.capFaceCount=int(faces.size())-interfaceFaces;result.boundaryLoopCount=2;
            result.faceIsCap=Eigen::VectorXi::Zero(faces.size());result.faceIsCap.tail(result.capFaceCount).setOnes();
            result.idealNormalEdgeIntersectionCount=int(nodes.size()-sphere.vertices.size()-intervalNormals.size());result.truncationRadius=radius;
            result.initialIcosphereSubdivisions=result.effectiveIcosphereSubdivisions=settings.icosphereSubdivisions;
            result.sourceNormalVertexByVertex=std::move(normalIds);
            result.sourceNormalEdgeByVertex.assign(vertices.size(),{-1,-1});
            result.sourceProjectiveRayByVertex=std::move(rays);
            result.fixedNormals=Eigen::MatrixXd::Zero(intervalNormals.empty()?sphere.vertices.size():nodes.size(),3);
            for(int i=0;i<int(sphere.vertices.size());++i)result.fixedNormals.row(i)=sphere.vertices[i];
            for(const auto& entry:intervalNormals)result.fixedNormals.row(entry.first)=entry.second;
            result.hash=triangulationHash(result,settings.algebraicFiniteMap);return result;
        }
        std::vector<Eigen::Vector3d> centers;
        bool separated=true;
        const double margin=256*std::numeric_limits<double>::epsilon()*std::max(1.0,radius);
        for(const auto& loop:loops) {
            Eigen::Vector3d center=Eigen::Vector3d::Zero();
            for(int i:loop)center+=vertices[i];center/=loop.size();centers.push_back(center);
            for(std::size_t k=0;k<loop.size();++k) {
                const std::array<Eigen::Vector3d,3> points{center,vertices[loop[k]],vertices[loop[(k+1)%loop.size()]]};
                bool outside=false;
                for(int a=0;a<3;++a)outside|=
                    std::all_of(points.begin(),points.end(),[&](const auto& v){return v(a)<lo(a)-margin;}) ||
                    std::all_of(points.begin(),points.end(),[&](const auto& v){return v(a)>hi(a)+margin;});
                separated&=outside;
            }
        }
        if(!separated){reason="paper caps not separated from design AABB";continue;}
        spin::QuadricTriangulation result;
        result.type=spin::QuadricSurfaceType::Hyperboloid;
        result.paperProjectiveHyperboloid=true;
        result.interfaceVertexCount=int(vertices.size());result.interfaceFaceCount=int(faces.size());
        for(std::size_t l=0;l<loops.size();++l) {
            const int center=int(vertices.size());vertices.push_back(centers[l]);
            normalIds.push_back(-1);rays.push_back({-1,-1,-1,-1});
            const auto& loop=loops[l];
            for(std::size_t k=0;k<loop.size();++k)faces.push_back({loop[(k+1)%loop.size()],loop[k],center});
        }
        result.closedBoundary=matrixMesh(vertices,faces);orientComponentsOutward(&result.closedBoundary);
        result.capFaceCount=int(faces.size())-result.interfaceFaceCount;result.boundaryLoopCount=2;
        result.faceIsCap=Eigen::VectorXi::Zero(faces.size());result.faceIsCap.tail(result.capFaceCount).setOnes();
        result.idealNormalEdgeIntersectionCount=int(nodes.size()-sphere.vertices.size());
        result.truncationRadius=radius;
        result.initialIcosphereSubdivisions=result.effectiveIcosphereSubdivisions=settings.icosphereSubdivisions;
        result.sourceNormalVertexByVertex=std::move(normalIds);
        result.sourceNormalEdgeByVertex.assign(vertices.size(),{-1,-1});
        result.sourceProjectiveRayByVertex=std::move(rays);
        result.fixedNormals.resize(sphere.vertices.size(),3);
        for(int i=0;i<int(sphere.vertices.size());++i)result.fixedNormals.row(i)=sphere.vertices[i];
        result.hash=triangulationHash(result,settings.algebraicFiniteMap);return result;
    }
    throw std::runtime_error("paper projective truncation budget exhausted: "+reason);
}


spin::QuadricTriangulation triangulatePlane(
    const spin::QuadricCoefficients10& quadric,
    double radius
)
{
    const Eigen::Vector3d normal = quadric.linearVector();
    const Eigen::Vector3d point =
        -quadric.constant() * normal / normal.squaredNorm();
    const auto frame=spin::makeQuadricPlaneFrame(normal.normalized());
    const Eigen::Vector3d first=frame.first;
    const Eigen::Vector3d second=frame.second;
    std::vector<Eigen::Vector3d> vertices = {
        point - radius * first - radius * second,
        point + radius * first - radius * second,
        point + radius * first + radius * second,
        point - radius * first + radius * second
    };
    std::vector<Face> faces = {{0, 1, 2}, {0, 2, 3}};
    for (Face& face : faces)
    {
        orientInterfaceFace(quadric, vertices, &face);
    }
    spin::QuadricTriangulation result = closeExtrudedInterface(
        spin::QuadricSurfaceType::Plane,
        vertices,
        faces,
        4.0 * radius * normal.normalized(),
        radius
    );
    result.hash = triangulationHash(result);
    return result;
}


spin::QuadricTriangulation triangulateParaboloid(
    const spin::QuadricCoefficients10& quadric,
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& eigen,
    double radius,
    int resolution
)
{
    const Eigen::Matrix3d eigenvectors=
        canonicalEigenvectorSigns(eigen.eigenvectors());
    int nullIndex = 0;
    eigen.eigenvalues().cwiseAbs().minCoeff(&nullIndex);
    std::array<int, 2> nonzero{};
    int cursor = 0;
    for (int index = 0; index < 3; ++index)
    {
        if (index != nullIndex)
        {
            nonzero[cursor++] = index;
        }
    }
    Eigen::Vector3d n = eigenvectors.col(nullIndex);
    Eigen::Vector3d u = eigenvectors.col(nonzero[0]);
    Eigen::Vector3d v = eigenvectors.col(nonzero[1]);
    double bn = quadric.linearVector().dot(n);
    if (bn < 0.0)
    {
        n = -n;
        bn = -bn;
    }
    const double bu = quadric.linearVector().dot(u);
    const double bv = quadric.linearVector().dot(v);
    const double au = eigen.eigenvalues()(nonzero[0]);
    const double av = eigen.eigenvalues()(nonzero[1]);

    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(static_cast<std::size_t>((resolution + 1) * (resolution + 1)));
    for (int row = 0; row <= resolution; ++row)
    {
        const double t = -radius + 2.0 * radius * row / resolution;
        for (int column = 0; column <= resolution; ++column)
        {
            const double s = -radius + 2.0 * radius * column / resolution;
            const double w = -(
                quadric.constant() + bu * s + bv * t +
                au * s * s + av * t * t
            ) / bn;
            vertices.push_back(u * s + v * t + n * w);
        }
    }
    std::vector<Face> faces;
    faces.reserve(static_cast<std::size_t>(2 * resolution * resolution));
    for (int row = 0; row < resolution; ++row)
    {
        for (int column = 0; column < resolution; ++column)
        {
            const int lower = row * (resolution + 1) + column;
            Face first = {lower, lower + 1, lower + resolution + 2};
            Face second = {lower, lower + resolution + 2, lower + resolution + 1};
            orientInterfaceFace(quadric, vertices, &first);
            orientInterfaceFace(quadric, vertices, &second);
            faces.push_back(first);
            faces.push_back(second);
        }
    }
    spin::QuadricTriangulation result = closeExtrudedInterface(
        spin::QuadricSurfaceType::Paraboloid,
        vertices,
        faces,
        4.0 * radius * n,
        radius
    );
    result.hash = triangulationHash(result);
    return result;
}


spin::QuadricTriangulation triangulateFixedNormalParaboloid(
    const spin::QuadricCoefficients10& q,
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& eigen,
    const Eigen::MatrixXd& domain,
    double radius,
    const spin::QuadricTriangulationSettings& settings)
{
    const auto vectors=canonicalEigenvectorSigns(eigen.eigenvectors());
    int nullIndex=0;eigen.eigenvalues().cwiseAbs().minCoeff(&nullIndex);
    Eigen::Vector3d axis=vectors.col(nullIndex);
    if(axis.dot(q.linearVector())<0)axis=-axis;
    const double bn=axis.dot(q.linearVector());
    const auto sphere=icosphere(settings.icosphereSubdivisions);
    std::vector<int> index(sphere.vertices.size(),-1),normalIds;
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Face> faces;
    // This explicit finite chart deliberately excludes near-ideal normals.
    // Missing remote faces are admissible only under the coverage gate below.
    for(const auto& face:sphere.faces) {
        bool finite=true;
        for(int id:face)finite&=sphere.vertices[id].dot(axis)>0.01;
        if(!finite)continue;
        Face mapped;
        for(int k=0;k<3;++k) {
            const int id=face[k];
            if(index[id]<0) {
                const auto& normal=sphere.vertices[id];
                const double multiplier=bn/normal.dot(axis);
                Eigen::Vector3d point=Eigen::Vector3d::Zero();
                for(int j=0;j<3;++j)if(j!=nullIndex) {
                    const Eigen::Vector3d direction=vectors.col(j);
                    point+=direction*(multiplier*normal.dot(direction)-
                        q.linearVector().dot(direction))/(2*eigen.eigenvalues()(j));
                }
                point-=axis*q.evaluate(point)/bn;
                if(!point.allFinite())throw std::runtime_error("nonfinite fixed-normal paraboloid point");
                index[id]=static_cast<int>(vertices.size());
                vertices.push_back(point);normalIds.push_back(id);
            }
            mapped[k]=index[id];
        }
        // This finite chart is a height graph along +axis, with q increasing
        // in that direction. Its positive-side extrusion therefore needs a
        // consistently -axis-facing lower boundary. A coarse saddle chord's
        // normal can oppose the analytic gradient at its centroid; orienting
        // each chord by that gradient can break otherwise shared edge winding.
        // Certify the projected sign of the represented vertices, not a
        // cancellation-prone floating triple product or an analytic tangent.
        const auto& a=vertices[mapped[0]];
        const auto& b=vertices[mapped[1]];
        const auto& c=vertices[mapped[2]];
        CGAL::Gmpq projectedArea(0);
        for(int i=0;i<3;++i) {
            const int j=(i+1)%3,k=(i+2)%3;
            projectedArea+=CGAL::Gmpq(axis(i))*(
                (CGAL::Gmpq(b(j))-CGAL::Gmpq(a(j)))*(CGAL::Gmpq(c(k))-CGAL::Gmpq(a(k)))-
                (CGAL::Gmpq(b(k))-CGAL::Gmpq(a(k)))*(CGAL::Gmpq(c(j))-CGAL::Gmpq(a(j))));
        }
        if(projectedArea==0)
            throw std::runtime_error("degenerate projected fixed-normal paraboloid face");
        if(projectedArea>0)std::swap(mapped[1],mapped[2]);
        faces.push_back(mapped);
    }
    const auto loops=boundaryLoops(faces);
    if(loops.size()!=1 || faces.empty())
        throw std::runtime_error("fixed-normal paraboloid finite patch is not one disk");
    // The finite rim may be nonconvex. Its kernel is unnecessarily restrictive:
    // require the whole projected design hull to lie inside the actual simple
    // polygon, including all hull edges and the original clearance margin.
    const auto& loop=loops.front();
    const double margin=1e-10*std::max(1.0,domain.cwiseAbs().maxCoeff());
    const spin::detail::ProjectedDomainClearance clearance(domain,axis,margin);
    std::vector<Eigen::Vector3d> rim;for(int id:loop)rim.push_back(vertices[id]);
    if(!clearance.inside(rim))
        throw std::runtime_error("fixed-normal paraboloid finite patch does not cover design hull");
    const Eigen::Vector3d lo=domain.colwise().minCoeff(),hi=domain.colwise().maxCoeff();
    for(int attempt=0;attempt<=settings.maximumRadiusDoublings;++attempt) {
        auto result=closeExtrudedInterface(spin::QuadricSurfaceType::Paraboloid,
            vertices,faces,4*radius*axis,radius);
        bool outside=true;
        // Disjointness alone is not sufficient: an offset interface and its
        // translated copy can BOTH precede the domain along +axis. Since q
        // increases along +axis, the translated graph must lie beyond the
        // whole domain to represent the selected positive side there.
        double domainAxisMaximum=-std::numeric_limits<double>::infinity();
        for(int v=0;v<domain.rows();++v)
            domainAxisMaximum=std::max(domainAxisMaximum,axis.dot(domain.row(v)));
        for(int v=result.interfaceVertexCount;v<result.closedBoundary.vertices.rows();++v)
            outside&=axis.dot(result.closedBoundary.vertices.row(v))>domainAxisMaximum+margin;
        if(!outside) {
            // Remote saddle vertices may precede the domain while the entire
            // far graph above its projection is safely beyond it. Certify
            // that local condition on actual rounded cap vertices, retaining
            // all subsequent cap-separation checks and the radius budget.
            std::vector<Eigen::Vector3d> far;
            for(int v=result.interfaceVertexCount;v<result.closedBoundary.vertices.rows();++v)
                far.push_back(result.closedBoundary.vertices.row(v));
            outside=clearance.graphAbove(far,faces,loop);
        }
        for(int f=result.interfaceFaceCount;f<result.closedBoundary.faces.rows();++f) {
            bool separated=false;
            for(int a=0;a<3;++a) {
                bool below=true,above=true;
                for(int k=0;k<3;++k) {
                    double x=result.closedBoundary.vertices(result.closedBoundary.faces(f,k),a);
                    below&=x<lo(a)-margin;above&=x>hi(a)+margin;
                }
                separated|=below||above;
            }
            if(!separated) {
                // A slanted lateral face can have an overlapping AABB while
                // lying strictly outside the domain hull. Certify separation
                // of the ACTUAL represented vertices by a linear functional.
                // Its floating normal is only a candidate: all interval bounds
                // and the positive gap are recomputed with exact rationals.
                const auto& v=result.closedBoundary.vertices;
                const auto& faces=result.closedBoundary.faces;
                const Eigen::Vector3d a=v.row(faces(f,0)),b=v.row(faces(f,1)),c=v.row(faces(f,2));
                const Eigen::Vector3d normal=(b-a).cross(c-a);
                if(normal.allFinite() && normal.norm()>0) {
                    const auto project=[&](const Eigen::Vector3d& point) {
                        CGAL::Gmpq sum(0);
                        for(int k=0;k<3;++k)sum+=CGAL::Gmpq(normal(k))*CGAL::Gmpq(point(k));
                        return sum;
                    };
                    auto tlo=project(a),thi=tlo,dlo=project(domain.row(0)),dhi=dlo;
                    for(int k=1;k<3;++k) {
                        const auto x=project(v.row(faces(f,k)));tlo=std::min(tlo,x);thi=std::max(thi,x);
                    }
                    for(int k=1;k<domain.rows();++k) {
                        const auto x=project(domain.row(k));dlo=std::min(dlo,x);dhi=std::max(dhi,x);
                    }
                    const CGAL::Gmpq gap=CGAL::Gmpq(margin)*CGAL::Gmpq(normal.norm());
                    separated=thi+gap<dlo || dhi+gap<tlo;
                }
            }
            // A nonconvex rim can be disjoint from the hull even when the
            // infinite cap plane cuts it. Use the ACTUAL cap vertices, since
            // rounded extrusion need not project to an exactly straight edge.
            if(!separated) {
                const auto& v=result.closedBoundary.vertices;
                const auto& face=result.closedBoundary.faces;
                separated=clearance.separated({v.row(face(f,0)),v.row(face(f,1)),v.row(face(f,2))});
            }
            outside&=separated;
        }
        if(outside) {
            result.fixedNormalParaboloid=true;
            result.initialIcosphereSubdivisions=result.effectiveIcosphereSubdivisions=settings.icosphereSubdivisions;
            result.fixedNormals.resize(sphere.vertices.size(),3);
            for(int i=0;i<result.fixedNormals.rows();++i)result.fixedNormals.row(i)=sphere.vertices[i];
            result.sourceNormalVertexByVertex=normalIds;
            result.sourceNormalVertexByVertex.insert(result.sourceNormalVertexByVertex.end(),normalIds.begin(),normalIds.end());
            result.hash=triangulationHash(result);return result;
        }
        radius*=2;
        if(!std::isfinite(radius))break;
    }
    throw std::runtime_error("fixed-normal paraboloid caps not certified outside design AABB");
}

spin::QuadricTriangulation triangulateEllipticCylinder(
    const spin::QuadricCoefficients10& quadric,
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& eigen,
    double radius,
    int resolution
)
{
    const Eigen::Matrix3d eigenvectors=
        canonicalEigenvectorSigns(eigen.eigenvectors());
    int nullIndex = 0;
    eigen.eigenvalues().cwiseAbs().minCoeff(&nullIndex);
    std::array<int, 2> nonzero{};
    int cursor = 0;
    for (int index = 0; index < 3; ++index)
    {
        if (index != nullIndex)
        {
            nonzero[cursor++] = index;
        }
    }
    // Eigen's self-adjoint solver sorts by signed eigenvalue. Reversing an
    // implicit quadric's sign therefore swaps the two radial directions for
    // an elliptic cylinder. The resulting point set is the same, but the
    // side quads choose different diagonals, so Boolean results for q and -q
    // no longer share an identical polyhedral interface. Order by magnitude
    // instead; this is invariant under the q -> -q representation change.
    if (std::abs(eigen.eigenvalues()(nonzero[1])) <
        std::abs(eigen.eigenvalues()(nonzero[0])))
    {
        std::swap(nonzero[0], nonzero[1]);
    }
    const Eigen::Vector3d n = eigenvectors.col(nullIndex);
    const Eigen::Vector3d u = eigenvectors.col(nonzero[0]);
    const Eigen::Vector3d v = eigenvectors.col(nonzero[1]);
    const double au = eigen.eigenvalues()(nonzero[0]);
    const double av = eigen.eigenvalues()(nonzero[1]);
    const double bu = quadric.linearVector().dot(u);
    const double bv = quadric.linearVector().dot(v);
    const double centerU = -bu / (2.0 * au);
    const double centerV = -bv / (2.0 * av);
    const double shifted = quadric.constant() -
        bu * bu / (4.0 * au) - bv * bv / (4.0 * av);
    const double radiusSquaredU = -shifted / au;
    const double radiusSquaredV = -shifted / av;
    if (radiusSquaredU <= 0.0 || radiusSquaredV <= 0.0)
    {
        throw std::invalid_argument("elliptic cylinder has no real surface");
    }
    // The interface winding is always selected against -gradient. For the
    // negative interior this points inward; negating the implicit function
    // makes it point outward while leaving the geometric cylinder unchanged.
    // Cap winding must follow that same side choice or the complement mesh is
    // not consistently oriented (and CGAL's PWN propagation rejects it).
    const double interiorSign = shifted > 0.0 ? 1.0 : -1.0;
    const double radiusU = std::sqrt(radiusSquaredU);
    const double radiusV = std::sqrt(radiusSquaredV);

    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(static_cast<std::size_t>(2 * resolution + 2));
    for (int layer = 0; layer < 2; ++layer)
    {
        const double height = layer == 0 ? -2.0 * radius : 2.0 * radius;
        for (int segment = 0; segment < resolution; ++segment)
        {
            const double angle =
                2.0 * std::acos(-1.0) * (segment + 0.5) / resolution;
            vertices.push_back(
                u * (centerU + radiusU * std::cos(angle)) +
                v * (centerV + radiusV * std::sin(angle)) +
                n * height
            );
        }
    }
    std::vector<Face> faces;
    faces.reserve(static_cast<std::size_t>(4 * resolution));
    for (int segment = 0; segment < resolution; ++segment)
    {
        const int next = (segment + 1) % resolution;
        Face first = {segment, next, next + resolution};
        Face second = {segment, next + resolution, segment + resolution};
        orientInterfaceFace(quadric, vertices, &first);
        orientInterfaceFace(quadric, vertices, &second);
        faces.push_back(first);
        faces.push_back(second);
    }
    const int interfaceVertexCount = static_cast<int>(vertices.size());
    const int interfaceFaceCount = static_cast<int>(faces.size());
    const int lowerCenter = static_cast<int>(vertices.size());
    vertices.push_back(u * centerU + v * centerV - 2.0 * radius * n);
    const int upperCenter = static_cast<int>(vertices.size());
    vertices.push_back(u * centerU + v * centerV + 2.0 * radius * n);
    for (int segment = 0; segment < resolution; ++segment)
    {
        const int next = (segment + 1) % resolution;
        Face lowerFace={next,segment,lowerCenter};
        const Eigen::Vector3d lowerNormal=
            (vertices[lowerFace[1]]-vertices[lowerFace[0]]).cross(
             vertices[lowerFace[2]]-vertices[lowerFace[0]]);
        const Eigen::Vector3d desiredLowerNormal = -interiorSign * n;
        if (lowerNormal.dot(desiredLowerNormal) < 0.0)
            std::swap(lowerFace[1], lowerFace[2]);
        faces.push_back(lowerFace);
        Face upperFace={
            segment + resolution,
            next + resolution,
            upperCenter
        };
        const Eigen::Vector3d upperNormal=
            (vertices[upperFace[1]]-vertices[upperFace[0]]).cross(
             vertices[upperFace[2]]-vertices[upperFace[0]]);
        const Eigen::Vector3d desiredUpperNormal = interiorSign * n;
        if (upperNormal.dot(desiredUpperNormal) < 0.0)
            std::swap(upperFace[1], upperFace[2]);
        faces.push_back(upperFace);
    }

    spin::QuadricTriangulation result;
    result.type = spin::QuadricSurfaceType::EllipticCylinder;
    result.closedBoundary = matrixMesh(vertices, faces);
    orientComponentsOutward(&result.closedBoundary);
    result.interfaceVertexCount = interfaceVertexCount;
    result.interfaceFaceCount = interfaceFaceCount;
    result.capFaceCount = static_cast<int>(faces.size()) - interfaceFaceCount;
    result.boundaryLoopCount = 2;
    result.truncationRadius = radius;
    result.faceIsCap = Eigen::VectorXi::Ones(
        static_cast<Eigen::Index>(faces.size())
    );
    result.faceIsCap.head(interfaceFaceCount).setZero();
    result.hash = triangulationHash(result);
    return result;
}


spin::QuadricTriangulation triangulateHyperbolicCylinder(
    const spin::QuadricCoefficients10& quadric,
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& eigen,
    double radius,
    int resolution)
{
    const Eigen::Vector3d values=eigen.eigenvalues();
    const Eigen::Matrix3d vectors=canonicalEigenvectorSigns(eigen.eigenvectors());
    int nullIndex=0;
    values.cwiseAbs().minCoeff(&nullIndex);
    std::array<int,2> nonzero{};
    int cursor=0;
    for(int i=0;i<3;++i)if(i!=nullIndex)nonzero[cursor++]=i;
    const double tolerance=1e-12*std::max(1.0,values.cwiseAbs().maxCoeff());
    if(std::abs(values(nonzero[0]))<=tolerance ||
       std::abs(values(nonzero[1]))<=tolerance ||
       values(nonzero[0])*values(nonzero[1])>=0.0)
        throw std::invalid_argument("hyperbolic-cylinder triangulation requires a simple rank-two indefinite quadratic");
    const Eigen::Vector3d b=quadric.linearVector();
    const Eigen::Vector3d n=vectors.col(nullIndex);
    const double bn=b.dot(n);
    if(std::abs(bn)>tolerance*std::max(1.0,b.norm()))
        throw std::invalid_argument("rank-two quadric is a paraboloid, not a hyperbolic cylinder");
    const double b0=b.dot(vectors.col(nonzero[0]));
    const double b1=b.dot(vectors.col(nonzero[1]));
    const double c0=-b0/(2.0*values(nonzero[0]));
    const double c1=-b1/(2.0*values(nonzero[1]));
    const double shifted=quadric.constant()-b0*b0/(4.0*values(nonzero[0]))-
        b1*b1/(4.0*values(nonzero[1]));
    if(!std::isfinite(shifted) || std::abs(shifted)<=tolerance)
        throw std::invalid_argument("hyperbolic-cylinder level is degenerate or numerically zero");

    int longitudinalPosition=-1;
    for(int position=0;position<2;++position)
        if(-shifted/values(nonzero[position])>0.0)longitudinalPosition=position;
    if(longitudinalPosition<0)
        throw std::invalid_argument("hyperbolic-cylinder has no real branch scale");
    const int transversePosition=1-longitudinalPosition;
    const int longitudinalIndex=nonzero[longitudinalPosition];
    const int transverseIndex=nonzero[transversePosition];
    const Eigen::Vector3d longitudinal=vectors.col(longitudinalIndex);
    const Eigen::Vector3d transverse=vectors.col(transverseIndex);
    const double longitudinalCenter=longitudinalPosition==0?c0:c1;
    const double transverseCenter=transversePosition==0?c0:c1;
    const double longitudinalEigenvalue=values(longitudinalIndex);
    const double transverseEigenvalue=values(transverseIndex);
    const double longitudinalScaleSquared=-shifted/longitudinalEigenvalue;
    const double transverseScaleSquared=shifted/transverseEigenvalue;
    if(!(longitudinalScaleSquared>0.0) || !(transverseScaleSquared>0.0) ||
       !std::isfinite(longitudinalScaleSquared) || !std::isfinite(transverseScaleSquared))
        throw std::invalid_argument("hyperbolic-cylinder real branch scales are invalid");
    const double longitudinalScale=std::sqrt(longitudinalScaleSquared);
    const double transverseScale=std::sqrt(transverseScaleSquared);
    const int side=resolution+1;
    const std::size_t perBranch=static_cast<std::size_t>(side)*static_cast<std::size_t>(side);
    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(2*perBranch+2);
    std::vector<Face> faces;
    faces.reserve(static_cast<std::size_t>(4*resolution*resolution)+8U*static_cast<std::size_t>(resolution));
    for(int branch=0;branch<2;++branch) {
        const double branchSign=branch==0?-1.0:1.0;
        const int offset=static_cast<int>(vertices.size());
        for(int row=0;row<=resolution;++row) {
            const double y=-2.0*radius+4.0*radius*row/resolution;
            const double longitudinalOffset=branchSign*longitudinalScale*
                std::sqrt(1.0+(y/transverseScale)*(y/transverseScale));
            for(int column=0;column<=resolution;++column) {
                const double z=-2.0*radius+4.0*radius*column/resolution;
                vertices.push_back(longitudinal*(longitudinalCenter+longitudinalOffset)+
                    transverse*(transverseCenter+y)+n*z);
            }
        }
        for(int row=0;row<resolution;++row)for(int column=0;column<resolution;++column) {
            const int lower=offset+row*side+column;
            Face first{lower,lower+1,lower+side+1};
            Face second{lower,lower+side+1,lower+side};
            orientInterfaceFace(quadric,vertices,&first);
            orientInterfaceFace(quadric,vertices,&second);
            // Both compact patches use the same material-side convention
            // (`-grad(q)`) at the quadric boundary. The closed components are
            // subsequently oriented independently; retaining this convention
            // keeps the interface orientation tied to the same q-side on both
            // branches.
            faces.push_back(first);
            faces.push_back(second);
        }
    }
    const int interfaceVertexCount=static_cast<int>(vertices.size());
    const int interfaceFaceCount=static_cast<int>(faces.size());
    const auto loops=boundaryLoops(faces);
    if(loops.size()!=2)
        throw std::runtime_error("hyperbolic-cylinder grid did not produce two branch boundary loops");
    for(const auto& loop:loops) {
        Eigen::Vector3d center=Eigen::Vector3d::Zero();
        for(const int vertex:loop)center+=vertices[static_cast<std::size_t>(vertex)];
        center/=static_cast<double>(loop.size());
        const double coordinate=(center-longitudinal*longitudinalCenter-
            transverse*transverseCenter).dot(longitudinal);
        const double branchSign=coordinate<0.0?-1.0:1.0;
        const int capCenter=static_cast<int>(vertices.size());
        vertices.push_back(center+branchSign*4.0*radius*longitudinal);
        for(std::size_t i=0;i<loop.size();++i) {
            const int current=loop[i];
            const int next=loop[(i+1)%loop.size()];
            faces.push_back({next,current,capCenter});
        }
    }

    spin::QuadricTriangulation result;
    result.type=spin::QuadricSurfaceType::HyperbolicCylinder;
    result.closedBoundary=matrixMesh(vertices,faces);
    orientComponentsOutward(&result.closedBoundary);
    result.interfaceVertexCount=interfaceVertexCount;
    result.interfaceFaceCount=interfaceFaceCount;
    result.capFaceCount=static_cast<int>(faces.size())-interfaceFaceCount;
    result.boundaryLoopCount=static_cast<int>(loops.size());
    result.truncationRadius=radius;
    result.faceIsCap=Eigen::VectorXi::Ones(static_cast<Eigen::Index>(faces.size()));
    result.faceIsCap.head(interfaceFaceCount).setZero();
    result.hash=triangulationHash(result);
    return result;
}

spin::QuadricTriangulation triangulateParabolicCylinder(
    const spin::QuadricCoefficients10& quadric,
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& eigen,
    double radius,
    int resolution)
{
    const Eigen::Vector3d values=eigen.eigenvalues();
    const Eigen::Matrix3d vectors=canonicalEigenvectorSigns(eigen.eigenvectors());
    int curvedIndex=0;
    values.cwiseAbs().maxCoeff(&curvedIndex);
    const double lambda=values(curvedIndex);
    const double scale=std::max(1.0,values.cwiseAbs().maxCoeff());
    if(std::abs(lambda)<=1e-12*scale)
        throw std::invalid_argument("parabolic-cylinder curvature is singular");
    for(int i=0;i<3;++i)if(i!=curvedIndex &&
       std::abs(values(i))>1e-12*scale)
        throw std::invalid_argument("parabolic-cylinder quadratic block is not rank one");
    const Eigen::Vector3d u=vectors.col(curvedIndex);
    const Eigen::Vector3d b=quadric.linearVector();
    const double bu=b.dot(u);
    const Eigen::Vector3d transverseLinear=b-bu*u;
    const double transverseNorm=transverseLinear.norm();
    if(!std::isfinite(transverseNorm) || transverseNorm<=1e-12*std::max(1.0,b.norm()))
        throw std::invalid_argument("rank-one quadric is not a real parabolic cylinder");
    const Eigen::Vector3d v=transverseLinear/transverseNorm;
    const Eigen::Vector3d w=u.cross(v).normalized();
    const int side=resolution+1;
    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(static_cast<std::size_t>(side)*static_cast<std::size_t>(side));
    for(int row=0;row<=resolution;++row) {
        const double r=-radius+2.0*radius*row/resolution;
        // p=a+b.x+x^T A x, hence p(u*r+v*s+w*z)=a+bu*r+beta*s+lambda*r*r.
        const double s=-(quadric.constant()+bu*r+lambda*r*r)/transverseNorm;
        for(int column=0;column<=resolution;++column) {
            const double z=-radius+2.0*radius*column/resolution;
            vertices.push_back(u*r+v*s+w*z);
        }
    }
    std::vector<Face> faces;
    faces.reserve(static_cast<std::size_t>(2*resolution*resolution));
    for(int row=0;row<resolution;++row)for(int column=0;column<resolution;++column) {
        const int lower=row*side+column;
        Face first{lower,lower+1,lower+side+1};
        Face second{lower,lower+side+1,lower+side};
        orientInterfaceFace(quadric,vertices,&first);
        orientInterfaceFace(quadric,vertices,&second);
        faces.push_back(first);
        faces.push_back(second);
    }
    auto result=closeExtrudedInterface(
        spin::QuadricSurfaceType::ParabolicCylinder,
        vertices,faces,4.0*radius*v,radius);
    result.hash=triangulationHash(result);
    return result;
}


} // namespace


namespace spin
{

QuadricPlaneFrame makeQuadricPlaneFrame(const Eigen::Vector3d& normal)
{
    QuadricPlaneFrame frame;
    if(!normal.allFinite() || std::abs(normal.norm()-1.0)>1e-10) {
        frame.diagnostic="plane frame requires a finite unit normal";
        return frame;
    }
    // Use the x-axis reference except in a cap around its singular direction,
    // where the y-axis reference is regular. This keeps both coordinate-axis
    // normals and the symmetric (1,1,1) direction away from chart seams.
    constexpr double switchThreshold=0.9;
    constexpr double chartSeparation=1e-10;
    const double xMagnitude=std::abs(normal.x());
    frame.chart=xMagnitude<switchThreshold?0:1;
    frame.chartStable=std::abs(xMagnitude-switchThreshold)>chartSeparation;
    const int referenceAxis=frame.chart==0?0:1;
    Eigen::Vector3d reference=Eigen::Vector3d::Zero();
    reference(referenceAxis)=1.0;
    const Eigen::Vector3d tangent=normal.cross(reference);
    const double length=tangent.norm();
    if(!std::isfinite(length) || length<=1e-12) {
        frame.chartStable=false;
        frame.diagnostic="plane frame chart has a singular reference axis";
        return frame;
    }
    frame.first=tangent/length;
    frame.second=normal.cross(frame.first);
    frame.diagnostic=frame.chartStable
        ? "stable plane tangent-frame chart"
        : "plane frame is on a chart boundary";
    return frame;
}


const char* quadricSurfaceTypeToken(QuadricSurfaceType type) noexcept
{
    switch (type)
    {
    case QuadricSurfaceType::Ellipsoid:
        return "ELLIPSOID";
    case QuadricSurfaceType::Hyperboloid:
        return "HYPERBOLOID";
    case QuadricSurfaceType::Paraboloid:
        return "PARABOLOID";
    case QuadricSurfaceType::EllipticCylinder:
        return "ELLIPTIC_CYLINDER";
    case QuadricSurfaceType::Plane:
        return "PLANE";
    case QuadricSurfaceType::Cone:
        return "CONE";
    case QuadricSurfaceType::HyperbolicCylinder:
        return "HYPERBOLIC_CYLINDER";
    case QuadricSurfaceType::ParabolicCylinder:
        return "PARABOLIC_CYLINDER";
    case QuadricSurfaceType::IntersectingPlanes:
        return "INTERSECTING_PLANES";
    }
    return "UNKNOWN";
}


const char* quadricSurfaceSubtypeToken(QuadricSurfaceSubtype type) noexcept
{
    switch (type)
    {
    case QuadricSurfaceSubtype::Plane:
        return "PLANE";
    case QuadricSurfaceSubtype::EllipticCylinder:
        return "ELLIPTIC_CYLINDER";
    case QuadricSurfaceSubtype::Ellipsoid:
        return "ELLIPSOID";
    case QuadricSurfaceSubtype::OneSheetHyperboloid:
        return "ONE_SHEET_HYPERBOLOID";
    case QuadricSurfaceSubtype::TwoSheetHyperboloid:
        return "TWO_SHEET_HYPERBOLOID";
    case QuadricSurfaceSubtype::Paraboloid:
        return "PARABOLOID";
    case QuadricSurfaceSubtype::Cone:
        return "CONE";
    case QuadricSurfaceSubtype::HyperbolicCylinder:
        return "HYPERBOLIC_CYLINDER";
    case QuadricSurfaceSubtype::ParabolicCylinder:
        return "PARABOLIC_CYLINDER";
    case QuadricSurfaceSubtype::IntersectingPlanes:
        return "INTERSECTING_PLANES";
    }
    return "UNKNOWN";
}


QuadricSurfaceSubtype classifyQuadricSurfaceSubtype(
    const QuadricCoefficients10& quadric
)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen;
    const QuadricSurfaceType coarse = classify(quadric, &eigen);
    switch (coarse)
    {
    case QuadricSurfaceType::Plane:
        return QuadricSurfaceSubtype::Plane;
    case QuadricSurfaceType::Cone:
        return QuadricSurfaceSubtype::Cone;
    case QuadricSurfaceType::HyperbolicCylinder:
        return QuadricSurfaceSubtype::HyperbolicCylinder;
    case QuadricSurfaceType::ParabolicCylinder:
        return QuadricSurfaceSubtype::ParabolicCylinder;
    case QuadricSurfaceType::IntersectingPlanes:
        return QuadricSurfaceSubtype::IntersectingPlanes;
    case QuadricSurfaceType::EllipticCylinder:
        return QuadricSurfaceSubtype::EllipticCylinder;
    case QuadricSurfaceType::Ellipsoid:
        return QuadricSurfaceSubtype::Ellipsoid;
    case QuadricSurfaceType::Paraboloid:
        return QuadricSurfaceSubtype::Paraboloid;
    case QuadricSurfaceType::Hyperboloid:
        break;
    }

    const int rightHandSideSign = -shiftedLevelSign(quadric);
    int matchingSigns = 0;
    for (int index = 0; index < 3; ++index)
    {
        if ((eigen.eigenvalues()(index) > 0.0 && rightHandSideSign > 0) ||
            (eigen.eigenvalues()(index) < 0.0 && rightHandSideSign < 0))
        {
            ++matchingSigns;
        }
    }
    if (matchingSigns == 2)
    {
        return QuadricSurfaceSubtype::OneSheetHyperboloid;
    }
    if (matchingSigns == 1)
    {
        return QuadricSurfaceSubtype::TwoSheetHyperboloid;
    }
    throw std::invalid_argument("quadric has no real hyperboloid subtype");
}


static QuadricTriangulation triangulateQuadricAtResolution(
    const QuadricCoefficients10& quadric,
    const Eigen::MatrixXd& designDomainVertices,
    const QuadricTriangulationSettings& settings,
    int initialIcosphereSubdivisions,
    const QuadricCoefficientAuthority10* authority=nullptr
)
{
    if (
        settings.icosphereSubdivisions < 0 ||
        settings.icosphereSubdivisions > 7 ||
        (settings.maximumIcosphereSubdivisions != -1 &&
         (settings.maximumIcosphereSubdivisions <
              settings.icosphereSubdivisions ||
          settings.maximumIcosphereSubdivisions > 7)) ||
        settings.paraboloidGridResolution < 1 ||
        settings.parabolicCylinderGridResolution < 1 ||
        settings.cylinderAngularResolution < 8 ||
        settings.hyperbolicCylinderGridResolution < 2 ||
        settings.maximumRadiusDoublings < 0 ||
        !std::isfinite(settings.truncationRadius) ||
        settings.truncationRadius < 0.0
    )
    {
        throw std::invalid_argument("invalid quadric triangulation settings");
    }
    double radius = settings.truncationRadius;
    if (radius == 0.0)
    {
        radius = automaticRadius(designDomainVertices);
    }
    else
    {
        (void)automaticRadius(designDomainVertices);
    }
    const Eigen::Vector3d designDomainMinimum =
        designDomainVertices.colwise().minCoeff().transpose();
    const Eigen::Vector3d designDomainMaximum =
        designDomainVertices.colwise().maxCoeff().transpose();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen;
    const QuadricSurfaceType type = authority ? classifyAuthority(*authority,quadric,&eigen) : classify(quadric, &eigen);
    if(settings.fixedNormalParaboloid &&
       (type!=QuadricSurfaceType::Paraboloid || settings.algebraicFiniteMap ||
        settings.paperProjectiveHyperboloid || settings.maximumIcosphereSubdivisions!=-1))
        throw std::invalid_argument("fixed-normal paraboloid requires its explicit fixed-resolution profile only");
    if(settings.algebraicFiniteMap && type!=QuadricSurfaceType::Ellipsoid &&
       type!=QuadricSurfaceType::Hyperboloid)
        throw std::invalid_argument("algebraic finite-map source requires full rank; rank-deficient source integration is not implemented");
    if (type == QuadricSurfaceType::Plane)
    {
        return triangulatePlane(quadric, radius);
    }
    if (type == QuadricSurfaceType::Paraboloid)
    {
        if(settings.fixedNormalParaboloid)
            return triangulateFixedNormalParaboloid(quadric,eigen,designDomainVertices,radius,settings);
        return triangulateParaboloid(
            quadric,
            eigen,
            radius,
            settings.paraboloidGridResolution
        );
    }
    if (type == QuadricSurfaceType::ParabolicCylinder)
    {
        return triangulateParabolicCylinder(
            quadric,eigen,radius,settings.parabolicCylinderGridResolution);
    }
    if (type == QuadricSurfaceType::EllipticCylinder)
    {
        // The finite cylinder boundary must be bitwise sign-invariant, not
        // merely geometrically close: p and -p are evaluated separately and
        // their Boolean outputs are later required to form an exact
        // polyhedral partition. Eigen can return the two radial eigenspaces in
        // a different order (and with slightly different binary64 bases) for
        // Q and -Q. Canonicalize the polynomial before solving that basis;
        // the Boolean side is determined independently from the original
        // polynomial by quadricSolidInteriorIsPositiveRegion().
        const QuadricCoefficients10 cylinderGeometry =
            quadric.quadraticMatrix().trace() < 0.0
                ? quadric.inverted()
                : quadric;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> canonicalEigen;
        const QuadricSurfaceType canonicalType =
            classify(cylinderGeometry, &canonicalEigen);
        if (canonicalType != QuadricSurfaceType::EllipticCylinder)
            throw std::logic_error(
                "sign-canonical elliptic-cylinder polynomial changed type");
        return triangulateEllipticCylinder(
            cylinderGeometry,
            canonicalEigen,
            radius,
            settings.cylinderAngularResolution
        );
    }
    if (type == QuadricSurfaceType::Cone)
    {
        throw std::invalid_argument(
            "unsupported degenerate cone: inverse-Gauss hyperboloid "
            "triangulation requires a nonzero centered level"
        );
    }
    if (type == QuadricSurfaceType::HyperbolicCylinder)
    {
        return triangulateHyperbolicCylinder(
            quadric,eigen,radius,settings.hyperbolicCylinderGridResolution);
    }
    if (type == QuadricSurfaceType::IntersectingPlanes)
    {
        throw std::invalid_argument(
            "unsupported intersecting planes: singular union is outside the "
            "regular quadric source triangulation"
        );
    }

    const QuadricInverseGaussMap inverseGaussMap(quadric);
    const Eigen::Matrix3d& inverseA =
        inverseGaussMap.inverseQuadraticMatrix();
    const SphereMesh sphere = icosphere(settings.icosphereSubdivisions);
    if (settings.rationalSourceCoordinates && type != QuadricSurfaceType::Ellipsoid && type != QuadricSurfaceType::Hyperboloid)
        throw std::invalid_argument("rational source coordinates require a supported full-rank source");
    if (type == QuadricSurfaceType::Ellipsoid)
    {
        std::vector<Eigen::Vector3d> vertices;
        std::vector<QuadricRationalPointResult> rationalPoints;
        vertices.reserve(sphere.vertices.size());
        for (const Eigen::Vector3d& normal : sphere.vertices)
        {
            if (settings.rationalSourceCoordinates)
            {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
                auto mapped = authority ? QuadricInverseGaussMap::rationalPointFromAuthority(*authority,normal,settings.compactRationalSourceCoordinates)
                                        : inverseGaussMap.rationalPoint(normal,settings.compactRationalSourceCoordinates);
                if (!mapped.defined)
                    throw std::runtime_error("rational ellipsoid finite source undefined: " + mapped.diagnostic);
                Eigen::Vector3d display;
                for (int axis = 0; axis < 3; ++axis)
                    igl::copyleft::cgal::assign_scalar(
                        CGAL::Epeck::FT(verification::parseExactRationalCoordinate(mapped.coordinates[axis])), display(axis));
                if (!display.allFinite()) throw std::runtime_error("rational source display is nonfinite");
                vertices.push_back(display);
                rationalPoints.push_back(std::move(mapped));
                continue;
#else
                throw std::invalid_argument("rational source coordinates require exact geometry checks");
#endif
            }
            const QuadricInverseGaussPointResult mapped =
                settings.algebraicFiniteMap ? inverseGaussMap.algebraicPoint(normal)
                                            : inverseGaussMap.point(normal);
            if (!mapped.valid())
            {
                throw std::runtime_error(
                    "ellipsoid inverse Gauss map unexpectedly failed"
                );
            }
            vertices.push_back(*mapped.point);
        }
        std::vector<Face> faces = sphere.faces;
        for (Face& face : faces)
        {
            // Icosphere connectivity is consistently oriented. Authority
            // orientation is decided below from its exact signed volume,
            // never from a possibly collapsed display triangle.
            if(!settings.rationalSourceCoordinates)orientInterfaceFace(quadric, vertices, &face);
        }
        QuadricTriangulation result;
        result.type = type;
        result.initialIcosphereSubdivisions = initialIcosphereSubdivisions;
        result.effectiveIcosphereSubdivisions =
            settings.icosphereSubdivisions;
        result.closedBoundary = matrixMesh(vertices, faces);
        if (settings.rationalSourceCoordinates)
        {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
            for (auto& point : rationalPoints)
            {
                result.exactVertices.push_back(std::move(point.coordinates));
                result.coordinateLower.push_back(std::move(point.lower));
                result.coordinateUpper.push_back(std::move(point.upper));
            }
            CGAL::Gmpq sixVolume(0);
            for (int f = 0; f < result.closedBoundary.faces.rows(); ++f)
            {
                std::array<std::array<CGAL::Gmpq,3>,3> p;
                for (int k = 0; k < 3; ++k) for (int a = 0; a < 3; ++a)
                    p[k][a] = verification::parseExactRationalCoordinate(result.exactVertices[result.closedBoundary.faces(f,k)][a]);
                sixVolume += p[0][0]*(p[1][1]*p[2][2]-p[1][2]*p[2][1])
                    -p[0][1]*(p[1][0]*p[2][2]-p[1][2]*p[2][0])
                    +p[0][2]*(p[1][0]*p[2][1]-p[1][1]*p[2][0]);
            }
            if (sixVolume < 0) result.closedBoundary.faces.col(1).swap(result.closedBoundary.faces.col(2));
            const auto proof = verification::auditExactTriangleBoundaryRational(
                result.exactVertices, result.closedBoundary.faces,
                verification::ExactBoundaryOrientation::OutwardPositive);
            if (!proof.facts().boundaryValid)
                throw std::runtime_error("rational ellipsoid source failed exact authority boundary audit");
#endif
        }
        else orientComponentsOutward(&result.closedBoundary);
        result.faceIsCap = Eigen::VectorXi::Zero(
            static_cast<Eigen::Index>(faces.size())
        );
        result.interfaceVertexCount = static_cast<int>(vertices.size());
        result.interfaceFaceCount = static_cast<int>(faces.size());
        result.truncationRadius = radius;
        result.compactRationalSourceCoordinates = settings.compactRationalSourceCoordinates;
        result.hash = triangulationHash(result,settings.algebraicFiniteMap);
        if (settings.retainFixedNormalProvenance || settings.rationalSourceCoordinates)
        {
            result.fixedNormals.resize(static_cast<Eigen::Index>(sphere.vertices.size()), 3);
            for (int i = 0; i < result.fixedNormals.rows(); ++i)
                result.fixedNormals.row(i) = sphere.vertices[static_cast<std::size_t>(i)];
        }
        if(settings.rationalSourceCoordinates)
            result.rationalConstructionIdentity=rationalConstructionIdentity(quadric,result,authority);
        return result;
    }

    if(settings.paperProjectiveHyperboloid)
        return triangulateProjectiveHyperboloid(quadric,sphere,designDomainVertices,settings,radius);

    std::string lastBoundaryLoopError;
    for (int doubling = 0; doubling <= settings.maximumRadiusDoublings;
         ++doubling)
    {
        std::vector<bool> inside(sphere.vertices.size(), false);
        std::vector<Eigen::Vector3d> mapped(sphere.vertices.size());
        std::vector<QuadricRationalPointResult> rationalMapped(sphere.vertices.size());
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        const TruncationRootCertifier rationalRoots=authority ? TruncationRootCertifier(*authority,radius) : TruncationRootCertifier(quadric,radius);
#endif
        for (std::size_t index = 0; index < sphere.vertices.size(); ++index)
        {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
            if(settings.rationalSourceCoordinates) {
                rationalMapped[index]=authority ? QuadricInverseGaussMap::rationalPointFromAuthority(*authority,sphere.vertices[index],settings.compactRationalSourceCoordinates)
                    : inverseGaussMap.rationalPoint(sphere.vertices[index],settings.compactRationalSourceCoordinates);
                if(rationalMapped[index].defined) {
                    mapped[index]=rationalDisplay(rationalMapped[index]);
                    const auto point=rationalPoint(rationalMapped[index]);
                    CGAL::Gmpq squared(0);for(const auto& x:point)squared+=x*x;
                    inside[index]=squared<=CGAL::Gmpq(radius)*CGAL::Gmpq(radius);
                }
                continue;
            }
#endif
            const QuadricInverseGaussPointResult point =
                settings.algebraicFiniteMap ? inverseGaussMap.algebraicPoint(sphere.vertices[index])
                                            : inverseGaussMap.point(sphere.vertices[index]);
            inside[index] = point.valid();
            if (inside[index])
            {
                mapped[index] = *point.point;
                inside[index] = mapped[index].norm() <= radius;
            }
        }

        std::vector<Eigen::Vector3d> vertices;
        std::vector<QuadricRationalPointResult> rationalVertices;
        std::vector<Face> faces;
        std::vector<std::array<int, 2>> sourceNormalEdgeByVertex;
        std::vector<int> sourceNormalVertexByVertex;
        std::map<int, int> originalVertex;
        std::map<Edge, int> crossingVertex;
        const auto addOriginal = [&](int index) {
            const auto found = originalVertex.find(index);
            if (found != originalVertex.end())
            {
                return found->second;
            }
            const int output = static_cast<int>(vertices.size());
            vertices.push_back(mapped[index]);
            if(settings.rationalSourceCoordinates)rationalVertices.push_back(rationalMapped[index]);
            sourceNormalEdgeByVertex.push_back({{-1, -1}});
            sourceNormalVertexByVertex.push_back(index);
            originalVertex.emplace(index, output);
            return output;
        };
        const auto addCrossing = [&](int first, int second) {
            const Edge edge = undirectedEdge(first, second);
            const auto found = crossingVertex.find(edge);
            if (found != crossingVertex.end())
            {
                return found->second;
            }
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
            if(settings.rationalSourceCoordinates) {
                auto point=rationalRoots.rationalPoint(sphere.vertices[edge.first],sphere.vertices[edge.second],settings.compactRationalSourceCoordinates);
                if(!point.defined)throw std::runtime_error("rational hyperboloid truncation undefined: "+point.diagnostic);
                const int output=static_cast<int>(vertices.size());
                vertices.push_back(rationalDisplay(point));rationalVertices.push_back(std::move(point));
                sourceNormalEdgeByVertex.push_back({{edge.first,edge.second}});
                sourceNormalVertexByVertex.push_back(-1);crossingVertex.emplace(edge,output);
                return output;
            }
#endif
            int insideIndex = inside[first] ? first : second;
            int outsideIndex = inside[first] ? second : first;
            Eigen::Vector3d insideNormal = sphere.vertices[insideIndex];
            Eigen::Vector3d outsideNormal = sphere.vertices[outsideIndex];
            Eigen::Vector3d insidePoint = mapped[insideIndex];
            for (int iteration = 0; iteration < 80; ++iteration)
            {
                Eigen::Vector3d middleNormal =
                    0.5 * (insideNormal + outsideNormal);
                middleNormal.normalize();
                const QuadricInverseGaussPointResult middle =
                    settings.algebraicFiniteMap ? inverseGaussMap.algebraicPoint(middleNormal)
                                                : inverseGaussMap.point(middleNormal);
                const bool middleInside = middle.valid() &&
                    middle.point->norm() <= radius;
                if (middleInside)
                {
                    insideNormal = middleNormal;
                    insidePoint = *middle.point;
                }
                else
                {
                    outsideNormal = middleNormal;
                }
            }
            insidePoint = refineTruncationCrossing(
                quadric,
                sphere.vertices[std::min(first, second)],
                sphere.vertices[std::max(first, second)],
                radius,
                insidePoint
            );
            const int output = static_cast<int>(vertices.size());
            vertices.push_back(insidePoint);
            sourceNormalEdgeByVertex.push_back({{edge.first, edge.second}});
            sourceNormalVertexByVertex.push_back(-1);
            crossingVertex.emplace(edge, output);
            return output;
        };

        for (const Face& input : sphere.faces)
        {
            std::vector<int> polygon;
            for (int corner = 0; corner < 3; ++corner)
            {
                const int current = input[corner];
                const int next = input[(corner + 1) % 3];
                if (inside[current])
                {
                    polygon.push_back(addOriginal(current));
                }
                if (inside[current] != inside[next])
                {
                    polygon.push_back(addCrossing(current, next));
                }
            }
            for (std::size_t corner = 1; corner + 1 < polygon.size(); ++corner)
            {
                Face face = {polygon[0], polygon[corner], polygon[corner + 1]};
                if(settings.rationalSourceCoordinates) {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
                    orientRationalInterfaceFace(quadric,rationalVertices,&face,authority);
#endif
                } else orientInterfaceFace(quadric, vertices, &face);
                faces.push_back(face);
            }
        }
        std::vector<std::vector<int>> loops;
        try
        {
            loops = boundaryLoops(faces);
        }
        catch (const BoundaryLoopTopologyError& error)
        {
            // A sampled sphere may pass through a mapped quadric vertex, making
            // the cap boundary combinatorially ambiguous at this radius. The
            // existing bounded radius-refinement protocol can move that cap
            // farther away; retain the last concrete topology witness if no
            // refinement produces a proper closed loop.
            lastBoundaryLoopError = error.what();
            radius *= 2.0;
            continue;
        }
        // The truncation sphere is centered at the model origin, not at the
        // canonical center of the quadric. An offset regular hyperboloid can
        // therefore meet this sphere in one loop (or more than two) even
        // though the centered canonical form has the usual two asymptotic
        // ends. Every recovered directed loop can be capped, provided each
        // fan triangle is certified outside the design-domain AABB below.
        if (loops.empty())
        {
            lastBoundaryLoopError = "hyperboloid truncation produced no "
                "boundary loops at radius " +
                std::to_string(radius);
            radius *= 2.0;
            continue;
        }

        // The cap is an implementation closure and must never alter the
        // Boolean material inside the requested design domain. AABB
        // separation is a conservative certificate: because each cap face
        // lies in the convex hull of its three vertices, one coordinate
        // half-space containing all three vertices proves the entire face is
        // outside the box. If any face lacks such a witness, increase the
        // truncation radius instead of accepting a cap that touches or enters
        // the design domain.
        const double capSeparationTolerance = 256.0 *
            std::numeric_limits<double>::epsilon() *
            std::max({1.0, radius, designDomainMinimum.cwiseAbs().maxCoeff(),
                      designDomainMaximum.cwiseAbs().maxCoeff()});
        std::vector<Eigen::Vector3d> capCenters;
        std::vector<QuadricRationalPointResult> rationalCapCenters;
        capCenters.reserve(loops.size());
        bool capsCertifiedOutsideDesignDomain = true;
        for (const std::vector<int>& loop : loops)
        {
            Eigen::Vector3d center = Eigen::Vector3d::Zero();
            for (int vertex : loop) center += vertices[vertex];
            center /= static_cast<double>(loop.size());
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
            if(settings.rationalSourceCoordinates) {
                QuadricRationalPointResult cap;cap.defined=true;
                for(int axis=0;axis<3;++axis) {
                    CGAL::Gmpq x(0),lo(0),hi(0);
                    for(int vertex:loop) {
                        x+=verification::parseExactRationalCoordinate(rationalVertices.at(vertex).coordinates[axis]);
                        lo+=verification::parseExactRationalCoordinate(rationalVertices.at(vertex).lower[axis]);
                        hi+=verification::parseExactRationalCoordinate(rationalVertices.at(vertex).upper[axis]);
                    }
                    const auto text=[&](CGAL::Gmpq value){std::ostringstream out;out.imbue(std::locale::classic());out<<value/CGAL::Gmpq(int(loop.size()));return out.str();};
                    cap.coordinates[axis]=text(x);cap.lower[axis]=text(lo);cap.upper[axis]=text(hi);
                }
                center=rationalDisplay(cap);rationalCapCenters.push_back(std::move(cap));
            }
#endif
            capCenters.push_back(center);
            for (std::size_t edge = 0; edge < loop.size(); ++edge)
            {
                const std::array<Eigen::Vector3d, 3> triangle = {
                    vertices[loop[(edge + 1) % loop.size()]],
                    vertices[loop[edge]],
                    center};
                bool separated = false;
                for (int axis = 0; axis < 3 && !separated; ++axis)
                {
                    const bool below = std::all_of(
                        triangle.begin(), triangle.end(), [&](const auto& point) {
                            return point(axis) <
                                designDomainMinimum(axis) - capSeparationTolerance;
                        });
                    const bool above = std::all_of(
                        triangle.begin(), triangle.end(), [&](const auto& point) {
                            return point(axis) >
                                designDomainMaximum(axis) + capSeparationTolerance;
                        });
                    separated = below || above;
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
                    if(authority) {
                        const std::array<RationalPoint,3> exactTriangle{{
                            rationalPoint(rationalVertices.at(loop[(edge+1)%loop.size()])),
                            rationalPoint(rationalVertices.at(loop[edge])),
                            rationalPoint(rationalCapCenters.back())}};
                        const CGAL::Gmpq lower=CGAL::Gmpq(designDomainMinimum(axis))-CGAL::Gmpq(capSeparationTolerance);
                        const CGAL::Gmpq upper=CGAL::Gmpq(designDomainMaximum(axis))+CGAL::Gmpq(capSeparationTolerance);
                        const bool exactBelow=std::all_of(exactTriangle.begin(),exactTriangle.end(),
                            [&](const auto& p){return p[axis]<lower;});
                        const bool exactAbove=std::all_of(exactTriangle.begin(),exactTriangle.end(),
                            [&](const auto& p){return p[axis]>upper;});
                        separated=exactBelow || exactAbove;
                    }
#endif
                }
                if (!separated)
                {
                    capsCertifiedOutsideDesignDomain = false;
                    break;
                }
            }
            if (!capsCertifiedOutsideDesignDomain) break;
        }
        if (!capsCertifiedOutsideDesignDomain)
        {
            lastBoundaryLoopError =
                "hyperboloid cap triangle could not be certified outside the "
                "design-domain AABB at radius " + std::to_string(radius);
            radius *= 2.0;
            continue;
        }

        const int interfaceVertexCount = static_cast<int>(vertices.size());
        const int interfaceFaceCount = static_cast<int>(faces.size());
        for (std::size_t loopIndex = 0; loopIndex < loops.size(); ++loopIndex)
        {
            const auto& loop = loops[loopIndex];
            const int centerIndex = static_cast<int>(vertices.size());
            vertices.push_back(capCenters[loopIndex]);
            if(settings.rationalSourceCoordinates)rationalVertices.push_back(rationalCapCenters[loopIndex]);
            sourceNormalEdgeByVertex.push_back({{-1, -1}});
            sourceNormalVertexByVertex.push_back(-1);
            for (std::size_t edge = 0; edge < loop.size(); ++edge)
            {
                faces.push_back({
                    loop[(edge + 1) % loop.size()],
                    loop[edge],
                    centerIndex
                });
            }
        }

        TriangleMesh boundary = matrixMesh(vertices, faces);
        if(!settings.rationalSourceCoordinates)orientComponentsOutward(&boundary);

        QuadricTriangulation result;
        result.type = type;
        result.closedBoundary = std::move(boundary);
        result.interfaceVertexCount = interfaceVertexCount;
        result.interfaceFaceCount = interfaceFaceCount;
        result.capFaceCount =
            static_cast<int>(faces.size()) - interfaceFaceCount;
        result.boundaryLoopCount = static_cast<int>(loops.size());
        result.initialIcosphereSubdivisions = initialIcosphereSubdivisions;
        result.effectiveIcosphereSubdivisions =
            settings.icosphereSubdivisions;
        result.idealNormalEdgeIntersectionCount =
            idealNormalEdgeIntersections(sphere, inverseA);
        result.truncationRadius = radius;
        if (settings.retainFixedNormalProvenance || settings.rationalSourceCoordinates)
        {
            result.sourceNormalEdgeByVertex = std::move(sourceNormalEdgeByVertex);
            result.sourceNormalVertexByVertex =
                std::move(sourceNormalVertexByVertex);
        }
        result.faceIsCap = Eigen::VectorXi::Zero(
            static_cast<Eigen::Index>(faces.size())
        );
        result.faceIsCap.tail(result.capFaceCount).setOnes();
        result.compactRationalSourceCoordinates = settings.compactRationalSourceCoordinates;
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        if(settings.rationalSourceCoordinates) {
            assembleRationalAuthority(rationalVertices,&result);
            result.fixedNormals.resize(sphere.vertices.size(),3);
            for(std::size_t i=0;i<sphere.vertices.size();++i)result.fixedNormals.row(i)=sphere.vertices[i];
        }
#endif
        result.hash = triangulationHash(result,settings.algebraicFiniteMap);
        if(settings.rationalSourceCoordinates)
            result.rationalConstructionIdentity=rationalConstructionIdentity(quadric,result,authority);
        return result;
    }
    if (settings.maximumIcosphereSubdivisions >
        settings.icosphereSubdivisions)
    {
        auto refinedSettings = settings;
        ++refinedSettings.icosphereSubdivisions;
        try
        {
            return triangulateQuadricAtResolution(
                quadric, designDomainVertices, refinedSettings,
                initialIcosphereSubdivisions,authority);
        }
        catch (const std::exception& error)
        {
            throw std::runtime_error(
                "bounded normal-sphere refinement failed after level " +
                std::to_string(settings.icosphereSubdivisions) +
                ": " + error.what());
        }
    }
    throw std::runtime_error(lastBoundaryLoopError.empty()
        ? "hyperboloid truncation did not produce any closed loops"
        : "hyperboloid truncation failed after bounded radius refinement: " +
              lastBoundaryLoopError);
}


QuadricTriangulation triangulateQuadric(
    const QuadricCoefficients10& quadric,
    const Eigen::MatrixXd& designDomainVertices,
    const QuadricTriangulationSettings& settings
)
{
    if(settings.rankTwoArrangementSource)
        throw std::invalid_argument(
            "rank-two arrangement requires PreparedRankTwoSource; legacy entry cannot construct it");
    if(settings.compactRationalSourceCoordinates && !settings.rationalSourceCoordinates)
        throw std::invalid_argument("compact rational source requires rational source coordinates");
    if (settings.rationalSourceCoordinates)
    {
#if !MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        throw std::invalid_argument("rational source coordinates require exact geometry checks");
#endif
        const auto subtype=classifyQuadricSurfaceSubtype(quadric);
        if (settings.paperProjectiveHyperboloid || settings.fixedNormalParaboloid ||
            (subtype != QuadricSurfaceSubtype::Ellipsoid && subtype != QuadricSurfaceSubtype::OneSheetHyperboloid &&
             subtype != QuadricSurfaceSubtype::TwoSheetHyperboloid))
            throw std::invalid_argument("rational source coordinates require a supported full-rank source");
    }
    return triangulateQuadricAtResolution(
        quadric, designDomainVertices, settings,
        settings.icosphereSubdivisions);
}

static bool constructionIdentityMatches(
    const QuadricCoefficients10& q,const QuadricTriangulation& t,
    const Eigen::MatrixXd& design,const QuadricTriangulationSettings& settings,
    const QuadricCoefficientAuthority10* authority) {
    if(t.rationalConstructionIdentity.empty())return false;
    if(t.compactRationalSourceCoordinates!=settings.compactRationalSourceCoordinates)return false;
    if(!std::isfinite(settings.truncationRadius) || settings.truncationRadius<0 ||
       settings.maximumRadiusDoublings<0)return false;
    double radius;
    try {
        const auto automatic=automaticRadius(design);
        radius=settings.truncationRadius==0?automatic:settings.truncationRadius;
    } catch(const std::invalid_argument&) {return false;}
    // Reuse only a radius obtainable under this exact request's doubling
    // budget. A larger/smaller or unrelated recipe is not the same mapping.
    for(int doubling=0;doubling<settings.maximumRadiusDoublings &&
        radius<t.truncationRadius && std::isfinite(radius);++doubling)radius*=2;
    if(radius!=t.truncationRadius)return false;
    const auto recomputed=rationalConstructionIdentity(q,t,authority);
    return !recomputed.empty() && recomputed==t.rationalConstructionIdentity;
}

bool rationalQuadricConstructionIdentityMatches(
    const QuadricCoefficients10& q,const QuadricTriangulation& t,
    const Eigen::MatrixXd& design,const QuadricTriangulationSettings& settings) {
    return constructionIdentityMatches(q,t,design,settings,nullptr);
}

static QuadricCoefficients10 coefficientDisplay(const QuadricCoefficientAuthority10& authority) {
    const auto values=authority.geometryDisplay();
    QuadricCoefficients10::Vector q;
    for(int i=0;i<10;++i)q(i)=values[i];
    return QuadricCoefficients10::fromNormalized(q);
}

QuadricTriangulation triangulateQuadric(
    const QuadricCoefficientAuthority10& authority,const Eigen::MatrixXd& design,
    const QuadricTriangulationSettings& settings) {
#if !MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    throw std::invalid_argument("coefficient authority requires exact geometry checks");
#endif
    if(!settings.rationalSourceCoordinates || settings.rankTwoArrangementSource ||
       settings.paperProjectiveHyperboloid || settings.fixedNormalParaboloid ||
       settings.algebraicFiniteMap ||
       (settings.maximumIcosphereSubdivisions!=-1 &&
        settings.maximumIcosphereSubdivisions!=settings.icosphereSubdivisions))
        throw std::invalid_argument("coefficient authority requires the explicit fixed full-rank rational source profile");
    return triangulateQuadricAtResolution(coefficientDisplay(authority),design,settings,
                                         settings.icosphereSubdivisions,&authority);
}

bool rationalQuadricConstructionIdentityMatches(
    const QuadricCoefficientAuthority10& authority,const QuadricTriangulation& t,
    const Eigen::MatrixXd& design,const QuadricTriangulationSettings& settings) {
    if(!settings.rationalSourceCoordinates || settings.rankTwoArrangementSource ||
       settings.paperProjectiveHyperboloid || settings.fixedNormalParaboloid ||
       settings.algebraicFiniteMap ||
       (settings.maximumIcosphereSubdivisions!=-1 &&
        settings.maximumIcosphereSubdivisions!=settings.icosphereSubdivisions) ||
       settings.icosphereSubdivisions!=t.initialIcosphereSubdivisions ||
       settings.icosphereSubdivisions!=t.effectiveIcosphereSubdivisions)
        return false;
    return constructionIdentityMatches(coefficientDisplay(authority),t,design,settings,&authority);
}


} // namespace spin
