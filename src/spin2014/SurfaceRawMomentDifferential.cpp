#include "spin2014/SurfaceRawMomentDifferential.h"
#include "physics/MassProperties.h"

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
SurfaceRawMomentDifferential differentiateSurfaceRawMoments(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double density
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

    // Spin-It (2014) supplement, Algorithms 1 and 2. This surface-normal
    // polynomial route is deliberately separate from signed tetrahedra.
    // Printed order: [1,x,y,z,xy,yz,xz,x^2,y^2,z^2].
    // Native order:  [M,LX,LY,LZ,QXX,QYY,QZZ,QYZ,QXZ,QXY].
    constexpr std::array<int, 10> nativeRow = {0, 1, 2, 3, 9, 7, 8, 4, 5, 6};
    constexpr std::array<double, 10> factor = {
        1./6, 1./24, 1./24, 1./24, 1./120, 1./120, 1./120,
        1./60, 1./60, 1./60
    };
    Eigen::Matrix<double, 10, 1> values = Eigen::Matrix<double, 10, 1>::Zero();
    SurfaceRawMomentDifferential result;
    result.jacobian = Eigen::MatrixXd::Zero(10, 3 * vertices.rows());
    for (Eigen::Index face = 0; face < faces.rows(); ++face)
    {
        const Eigen::Array3d a = vertices.row(faces(face, 0)).transpose();
        const Eigen::Array3d b = vertices.row(faces(face, 1)).transpose();
        const Eigen::Array3d c = vertices.row(faces(face, 2)).transpose();
        const Eigen::Vector3d u = (b-a).matrix(), v = (c-a).matrix();
        const Eigen::Array3d n = u.cross(v).array();
        if (!n.isFinite().all() || n.abs().maxCoeff() == 0.0)
        {
            throw std::invalid_argument("surface triangle is degenerate or nonfinite");
        }
        const Eigen::Array3d h1 = a+b+c;
        const Eigen::Array3d h2 = a*a + b*(a+b);
        const Eigen::Array3d h3 = h2 + c*h1;
        const Eigen::Array3d h4 = a*a*a + b*h2 + c*h3;
        const Eigen::Array3d h5 = h3 + a*(h1+a);
        const Eigen::Array3d h6 = h3 + b*(h1+b);
        const Eigen::Array3d h7 = h3 + c*(h1+c);
        const Eigen::Array3d h8 = cycle(a)*h5 + cycle(b)*h6 + cycle(c)*h7;

        Eigen::Matrix<double, 10, 1> local;
        local << n.x()*h1.x(), (n*h3).matrix(), (n*h8).matrix(), (n*h4).matrix();
        const std::array<Eigen::Matrix3d, 3> dn = {skew(v)-skew(u), -skew(v), skew(u)};
        const std::array<Eigen::Array3d, 3> d3 = {2*a+b+c, a+2*b+c, a+b+2*c};
        const Eigen::Array3d d5a = 6*a+2*b+2*c, d6b = 2*a+6*b+2*c;
        const Eigen::Array3d d7c = 2*a+2*b+6*c;
        const Eigen::Array3d d5b = 2*a+2*b+c, d5c = 2*a+b+2*c, d6c = a+2*b+2*c;
        const std::array<Eigen::Array3d, 3> d8 = {
            cycle(a)*d5a + cycle(b)*d5b + cycle(c)*d5c,
            cycle(a)*d5b + cycle(b)*d6b + cycle(c)*d6c,
            cycle(a)*d5c + cycle(b)*d6c + cycle(c)*d7c
        };
        const std::array<Eigen::Array3d, 3> d4 = {h5, h6, h7};
        for (int printed = 0; printed < 10; ++printed)
        {
            values(nativeRow[printed]) += density*factor[printed]*local(printed);
        }
        for (int vertex = 0; vertex < 3; ++vertex)
        {
            Eigen::Matrix<double, 10, 3> block;
            // Row-expanded Algorithm 2: d(n_i h_i) = h_i dn_i + n_i dh_i.
            // The printed H_i replication multiplies dn componentwise; it
            // must not be read as an ordinary dense matrix product.
            block.row(0) = h1.x()*dn[vertex].row(0);
            block(0, 0) += n.x();
            for (int axis = 0; axis < 3; ++axis)
            {
                block.row(1+axis) = h3(axis)*dn[vertex].row(axis);
                block(1+axis, axis) += n(axis)*d3[vertex](axis);
                block.row(4+axis) = h8(axis)*dn[vertex].row(axis);
                block(4+axis, axis) += n(axis)*d8[vertex](axis);
                // Published [n*h5/6/7]_* accounts for d(cycle(vertex)).
                block(4+axis, (axis+1)%3) += n(axis)*d4[vertex](axis);
                block.row(7+axis) = h4(axis)*dn[vertex].row(axis);
                block(7+axis, axis) += n(axis)*d4[vertex](axis);
            }
            for (int printed = 0; printed < 10; ++printed)
            {
                result.jacobian.block<1, 3>(nativeRow[printed], 3*faces(face, vertex)) +=
                    density*factor[printed]*block.row(printed);
            }
        }
    }
    MassMoments mass;
    mass.mass = values(0);
    mass.firstMoment = values.segment<3>(1);
    mass.secondMoment << values(4), values(9), values(8),
                         values(9), values(5), values(7),
                         values(8), values(7), values(6);
    result.moments = physics::RawMassMoments10::fromMassMoments(mass);
    if (!result.moments.isFinite() || !result.jacobian.allFinite())
    {
        throw std::overflow_error("raw-moment differential overflowed");
    }
    return result;
}
} // namespace spin2014
