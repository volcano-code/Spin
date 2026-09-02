#include "spin2014/SurfaceRawMomentDifferential.h"

#include "physics/MassProperties.h"

#include <Eigen/Core>

#include <array>
#include <cmath>
#include <stdexcept>


namespace
{


using LocalGradient = Eigen::Matrix<double, 9, 1>;


struct Dual
{
    double value = 0.0;
    LocalGradient gradient = LocalGradient::Zero();
};


Dual variable(double value, int coordinate)
{
    Dual result;
    result.value = value;
    result.gradient(coordinate) = 1.0;
    return result;
}


Dual operator+(const Dual& left, const Dual& right)
{
    return {left.value + right.value, left.gradient + right.gradient};
}


Dual operator-(const Dual& left, const Dual& right)
{
    return {left.value - right.value, left.gradient - right.gradient};
}


Dual operator*(const Dual& left, const Dual& right)
{
    return {
        left.value * right.value,
        left.gradient * right.value + right.gradient * left.value
    };
}


Dual operator*(double scale, const Dual& value)
{
    return {scale * value.value, scale * value.gradient};
}


using DualVector3 = std::array<Dual, 3>;


DualVector3 cross(const DualVector3& left, const DualVector3& right)
{
    return {
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0]
    };
}


Dual dot(const DualVector3& left, const DualVector3& right)
{
    return left[0] * right[0] +
        left[1] * right[1] +
        left[2] * right[2];
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
    if (
        vertices.cols() != 3 || !vertices.allFinite() ||
        faces.cols() != 3 ||
        (faces.size() > 0 &&
         (faces.minCoeff() < 0 || faces.maxCoeff() >= vertices.rows())) ||
        !std::isfinite(density) || density <= 0.0 ||
        !isClosedOrientedTriangleChain(faces)
    )
    {
        throw std::invalid_argument(
            "raw-moment differential requires a finite closed triangle surface"
        );
    }

    Eigen::Matrix<double, 10, 1> values =
        Eigen::Matrix<double, 10, 1>::Zero();
    SurfaceRawMomentDifferential result;
    result.jacobian = Eigen::MatrixXd::Zero(10, 3 * vertices.rows());
    for (Eigen::Index face = 0; face < faces.rows(); ++face)
    {
        DualVector3 point[3];
        for (int local = 0; local < 3; ++local)
        {
            const int vertex = faces(face, local);
            for (int axis = 0; axis < 3; ++axis)
            {
                point[local][axis] = variable(
                    vertices(vertex, axis),
                    3 * local + axis
                );
            }
        }
        const Dual volume = (density / 6.0) * dot(
            point[0], cross(point[1], point[2])
        );
        std::array<Dual, 10> moments;
        moments[0] = volume;
        DualVector3 sum;
        for (int axis = 0; axis < 3; ++axis)
        {
            sum[axis] = point[0][axis] + point[1][axis] + point[2][axis];
            moments[1 + axis] = 0.25 * (volume * sum[axis]);
        }
        const int diagonalRaw[3] = {4, 5, 6};
        for (int axis = 0; axis < 3; ++axis)
        {
            Dual vertexProducts =
                point[0][axis] * point[0][axis] +
                point[1][axis] * point[1][axis] +
                point[2][axis] * point[2][axis];
            moments[diagonalRaw[axis]] = 0.05 * (
                volume * (sum[axis] * sum[axis] + vertexProducts)
            );
        }
        const int firstAxis[3] = {1, 0, 0};
        const int secondAxis[3] = {2, 2, 1};
        const int crossRaw[3] = {7, 8, 9};
        for (int term = 0; term < 3; ++term)
        {
            const int first = firstAxis[term];
            const int second = secondAxis[term];
            Dual vertexProducts =
                point[0][first] * point[0][second] +
                point[1][first] * point[1][second] +
                point[2][first] * point[2][second];
            moments[crossRaw[term]] = 0.05 * (
                volume * (sum[first] * sum[second] + vertexProducts)
            );
        }

        for (int moment = 0; moment < 10; ++moment)
        {
            values(moment) += moments[static_cast<std::size_t>(moment)].value;
            for (int local = 0; local < 3; ++local)
            {
                const int vertex = faces(face, local);
                result.jacobian.block<1, 3>(moment, 3 * vertex) +=
                    moments[static_cast<std::size_t>(moment)]
                        .gradient.segment<3>(3 * local).transpose();
            }
        }
    }
    MassMoments mass;
    mass.mass = values(0);
    mass.firstMoment = values.segment<3>(1);
    mass.secondMoment <<
        values(4), values(9), values(8),
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
