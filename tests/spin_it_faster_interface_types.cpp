#include "spin/QuadricTriangulator.h"

#include <cstdlib>
#include <iostream>
#include <string>


namespace
{


using Vector10 = spin::QuadricCoefficients10::Vector;


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-it-faster-interface-types "
              << message << '\n';
    return EXIT_FAILURE;
}


spin::QuadricCoefficients10 quadric(
    double constant,
    const Eigen::Vector3d& linear,
    const Eigen::Vector3d& diagonal
)
{
    Vector10 values = Vector10::Zero();
    values(0) = constant;
    values.segment<3>(1) = linear;
    values.segment<3>(4) = diagonal;
    return spin::QuadricCoefficients10::fromUnnormalized(values);
}


} // namespace


int main()
{
    struct Case
    {
        spin::QuadricCoefficients10 value;
        const char* expected;
    };
    const Case cases[] = {
        {quadric(0.1, Eigen::Vector3d(1.0, 0.0, 0.0),
                 Eigen::Vector3d::Zero()), "PLANE"},
        {quadric(-1.0, Eigen::Vector3d::Zero(),
                 Eigen::Vector3d(1.0, 1.0, 0.0)), "ELLIPTIC_CYLINDER"},
        {quadric(-1.0, Eigen::Vector3d::Zero(),
                 Eigen::Vector3d::Ones()), "ELLIPSOID"},
        {quadric(-1.0, Eigen::Vector3d::Zero(),
                 Eigen::Vector3d(1.0, 1.0, -1.0)),
                 "ONE_SHEET_HYPERBOLOID"},
        {quadric(1.0, Eigen::Vector3d::Zero(),
                 Eigen::Vector3d(1.0, 1.0, -1.0)),
                 "TWO_SHEET_HYPERBOLOID"},
        {quadric(0.0, Eigen::Vector3d(0.0, 0.0, -1.0),
                 Eigen::Vector3d(1.0, 1.0, 0.0)), "PARABOLOID"},
    };
    for (const Case& item : cases)
    {
        const std::string actual = spin::quadricSurfaceSubtypeToken(
            spin::classifyQuadricSurfaceSubtype(item.value)
        );
        if (actual != item.expected)
        {
            return fail(
                "expected " + std::string(item.expected) + ", got " + actual
            );
        }
    }
    std::cout << "[PASS] scenario=spin-it-faster-interface-types\n";
    return EXIT_SUCCESS;
}
