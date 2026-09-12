#include "spin2014/SurfaceRawMomentDifferential.h"
#include "physics/MassProperties.h"
#include "SurfaceTriangleMoments.h"

#include <Eigen/Geometry>
#include <array>
#include <cmath>
#include <stdexcept>

namespace
{
Eigen::Matrix3d skew(const Eigen::Vector3d& v)
{
    Eigen::Matrix3d result;
    result << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
    return result;
}

Eigen::Array3d cycle(const Eigen::Array3d& v)
{
    return {v.y(), v.z(), v.x()};
}
} // namespace

namespace spin2014
{
detail::TriangleMomentValues detail::surfaceTriangleMoments(
    const Eigen::Vector3d& av, const Eigen::Vector3d& bv,
    const Eigen::Vector3d& cv, TriangleMomentJacobian* jacobian,
    NormalOverflowCategory normalOverflow)
{
    const Eigen::Array3d a=av.array(), b=bv.array(), c=cv.array();
    if (!a.isFinite().all() || !b.isFinite().all() || !c.isFinite().all())
        throw std::overflow_error("surface triangle has nonfinite coordinates");
    constexpr std::array<int,10> nativeRow={0,1,2,3,9,7,8,4,5,6};
    constexpr std::array<double,10> factor={1./6,1./24,1./24,1./24,
        1./120,1./120,1./120,1./60,1./60,1./60};
    const Eigen::Vector3d u=bv-av, v=cv-av;
    const Eigen::Array3d n=u.cross(v).array();
    if (!n.isFinite().all()) {
        if (normalOverflow==NormalOverflowCategory::LegacyInvalidArgument)
            throw std::invalid_argument("surface triangle is degenerate or nonfinite");
        throw std::overflow_error("surface triangle normal overflow");
    }
    if (n.abs().maxCoeff()==0) throw std::invalid_argument("surface triangle is degenerate");
    const Eigen::Array3d h1=a+b+c;
    const Eigen::Array3d h2=a*a+b*(a+b);
    const Eigen::Array3d h3=h2+c*h1;
    const Eigen::Array3d h4=a*a*a+b*h2+c*h3;
    const Eigen::Array3d h5=h3+a*(h1+a), h6=h3+b*(h1+b), h7=h3+c*(h1+c);
    const Eigen::Array3d h8=cycle(a)*h5+cycle(b)*h6+cycle(c)*h7;
    TriangleMomentValues local, values;
    local << n.x()*h1.x(), (n*h3).matrix(), (n*h8).matrix(), (n*h4).matrix();
    for (int row=0; row<10; ++row) values(nativeRow[row])=factor[row]*local(row);
    if (!values.allFinite()) throw std::overflow_error("surface triangle moment overflow");
    if (!jacobian) return values;
    const std::array<Eigen::Matrix3d,3> dn={skew(v)-skew(u),-skew(v),skew(u)};
    const std::array<Eigen::Array3d,3> d3={2*a+b+c,a+2*b+c,a+b+2*c};
    const Eigen::Array3d d5a=6*a+2*b+2*c, d6b=2*a+6*b+2*c, d7c=2*a+2*b+6*c;
    const Eigen::Array3d d5b=2*a+2*b+c, d5c=2*a+b+2*c, d6c=a+2*b+2*c;
    const std::array<Eigen::Array3d,3> d8={
        cycle(a)*d5a+cycle(b)*d5b+cycle(c)*d5c,
        cycle(a)*d5b+cycle(b)*d6b+cycle(c)*d6c,
        cycle(a)*d5c+cycle(b)*d6c+cycle(c)*d7c};
    const std::array<Eigen::Array3d,3> d4={h5,h6,h7};
    for (int vertex=0; vertex<3; ++vertex) {
        Eigen::Matrix<double,10,3> block;
        block.row(0)=h1.x()*dn[vertex].row(0);
        block(0,0)+=n.x();
        for (int axis=0; axis<3; ++axis) {
            block.row(1+axis)=h3(axis)*dn[vertex].row(axis);
            block(1+axis,axis)+=n(axis)*d3[vertex](axis);
            block.row(4+axis)=h8(axis)*dn[vertex].row(axis);
            block(4+axis,axis)+=n(axis)*d8[vertex](axis);
            block(4+axis,(axis+1)%3)+=n(axis)*d4[vertex](axis);
            block.row(7+axis)=h4(axis)*dn[vertex].row(axis);
            block(7+axis,axis)+=n(axis)*d4[vertex](axis);
        }
        for (int row=0; row<10; ++row)
            jacobian->block<1,3>(nativeRow[row],3*vertex)=factor[row]*block.row(row);
    }
    if (!jacobian->allFinite()) throw std::overflow_error("surface triangle derivative overflow");
    return values;
}
namespace {
physics::RawMassMoments10 integrate(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density,
    Eigen::MatrixXd* jacobian
)
{
    if (vertices.cols() != 3 || vertices.rows() == 0 || !vertices.allFinite() ||
        faces.cols() != 3 || faces.rows() == 0 ||
        faces.minCoeff() < 0 || faces.maxCoeff() >= vertices.rows() ||
        !std::isfinite(density) || density <= 0.0 ||
        !isClosedOrientedTriangleChain(faces))
    {
        throw std::invalid_argument(
            "raw-moment differential requires a finite closed triangle surface"
        );
    }

    Eigen::Matrix<double, 10, 1> values = Eigen::Matrix<double, 10, 1>::Zero();
    if (jacobian) *jacobian = Eigen::MatrixXd::Zero(10, 3 * vertices.rows());
    for (Eigen::Index face = 0; face < faces.rows(); ++face)
    {
        detail::TriangleMomentJacobian localJacobian;
        values += density*detail::surfaceTriangleMoments(
            vertices.row(faces(face,0)).transpose(),vertices.row(faces(face,1)).transpose(),
            vertices.row(faces(face,2)).transpose(),jacobian ? &localJacobian : nullptr,
            detail::NormalOverflowCategory::LegacyInvalidArgument);
        if (!jacobian) continue;
        for (int vertex = 0; vertex < 3; ++vertex)
            jacobian->block<10,3>(0,3*faces(face,vertex)) += density*localJacobian.block<10,3>(0,3*vertex);
    }
    MassMoments mass;
    mass.mass = values(0);
    mass.firstMoment = values.segment<3>(1);
    mass.secondMoment << values(4), values(9), values(8),
                         values(9), values(5), values(7),
                         values(8), values(7), values(6);
    const auto result = physics::RawMassMoments10::fromMassMoments(mass);
    if (!result.isFinite() || (jacobian && !jacobian->allFinite()))
    {
        throw std::overflow_error("raw-moment differential overflowed");
    }
    return result;
}
} // namespace
physics::RawMassMoments10 integrateSurfaceRawMoments(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, double density)
{
    return integrate(vertices, faces, density, nullptr);
}
SurfaceRawMomentDifferential differentiateSurfaceRawMoments(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, double density)
{
    SurfaceRawMomentDifferential result;
    result.moments = integrate(vertices, faces, density, &result.jacobian);
    return result;
}
} // namespace spin2014
