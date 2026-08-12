#include "manufacturing/ManufacturingGeometry.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>


namespace
{


using manufacturing::BoundaryAudit;
using manufacturing::BoundaryOrientation;


struct TriangleMesh
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
};


int fail(const std::string& scenario, const std::string& reason)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << " reason=" << reason << '\n';
    return EXIT_FAILURE;
}


int pass(const std::string& scenario)
{
    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


TriangleMesh makeBox(
    const Eigen::Vector3d& minimum,
    const Eigen::Vector3d& maximum,
    bool outward = true
)
{
    TriangleMesh result;
    result.vertices.resize(8, 3);
    result.vertices <<
        minimum(0), minimum(1), minimum(2),
        maximum(0), minimum(1), minimum(2),
        maximum(0), maximum(1), minimum(2),
        minimum(0), maximum(1), minimum(2),
        minimum(0), minimum(1), maximum(2),
        maximum(0), minimum(1), maximum(2),
        maximum(0), maximum(1), maximum(2),
        minimum(0), maximum(1), maximum(2);

    result.faces.resize(12, 3);
    result.faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;

    if (!outward)
    {
        result.faces.col(1).swap(result.faces.col(2));
    }
    return result;
}


TriangleMesh join(const TriangleMesh& first, const TriangleMesh& second)
{
    TriangleMesh result;
    result.vertices.resize(
        first.vertices.rows() + second.vertices.rows(),
        3
    );
    result.vertices.topRows(first.vertices.rows()) = first.vertices;
    result.vertices.bottomRows(second.vertices.rows()) = second.vertices;

    result.faces.resize(first.faces.rows() + second.faces.rows(), 3);
    result.faces.topRows(first.faces.rows()) = first.faces;
    result.faces.bottomRows(second.faces.rows()) =
        (second.faces.array() + first.vertices.rows()).matrix();
    return result;
}


BoundaryAudit auditOutward(const TriangleMesh& mesh)
{
    return manufacturing::auditTriangleBoundary(
        mesh.vertices,
        mesh.faces,
        BoundaryOrientation::OutwardPositive
    );
}


int runCube()
{
    const std::string scenario = "cube";
    const TriangleMesh outward = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    const BoundaryAudit outerAudit = auditOutward(outward);

    const TriangleMesh inward = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),
        false
    );
    const BoundaryAudit innerAudit = manufacturing::auditTriangleBoundary(
        inward.vertices,
        inward.faces,
        BoundaryOrientation::InwardNegative
    );

    if (
        !outerAudit.strictBoundaryValid() ||
        outerAudit.faceComponentCount != 1 ||
        std::abs(outerAudit.signedVolume - 1.0) > 1.0e-12 ||
        !innerAudit.strictBoundaryValid() ||
        innerAudit.faceComponentCount != 1 ||
        std::abs(innerAudit.signedVolume + 1.0) > 1.0e-12
    )
    {
        return fail(scenario, "oriented analytic cube contract failed");
    }
    return pass(scenario);
}


int runMissing()
{
    const std::string scenario = "missing";
    TriangleMesh mesh = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    mesh.faces.conservativeResize(mesh.faces.rows() - 1, 3);
    const BoundaryAudit audit = auditOutward(mesh);
    if (audit.strictBoundaryValid() || audit.closed || audit.boundaryEdgeCount == 0)
    {
        return fail(scenario, "open boundary was not rejected");
    }
    return pass(scenario);
}


int runReversed()
{
    const std::string scenario = "reversed";
    const TriangleMesh mesh = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),
        false
    );
    const BoundaryAudit audit = auditOutward(mesh);
    if (
        audit.strictBoundaryValid() ||
        !audit.consistentlyOriented ||
        audit.componentOrientationsValid ||
        !(audit.signedVolume < 0.0)
    )
    {
        return fail(scenario, "globally reversed shell was not diagnosed");
    }
    return pass(scenario);
}


int runDuplicate()
{
    const std::string scenario = "duplicate";
    TriangleMesh mesh = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    const Eigen::Index oldRows = mesh.faces.rows();
    mesh.faces.conservativeResize(oldRows + 1, 3);
    mesh.faces.row(oldRows) = mesh.faces.row(0);
    const BoundaryAudit audit = auditOutward(mesh);
    if (
        audit.strictBoundaryValid() ||
        audit.nonManifoldEdgeCount == 0 ||
        audit.edgeManifold
    )
    {
        return fail(scenario, "duplicate facet was not rejected");
    }
    return pass(scenario);
}


int runDegenerate()
{
    const std::string scenario = "degenerate";
    TriangleMesh mesh = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    mesh.faces(0, 2) = mesh.faces(0, 1);
    const BoundaryAudit audit = auditOutward(mesh);
    if (
        audit.strictBoundaryValid() ||
        audit.nondegenerateFaces ||
        audit.degenerateFaceCount == 0
    )
    {
        return fail(scenario, "zero-area facet was not rejected");
    }
    return pass(scenario);
}


int runDoubleOuter()
{
    const std::string scenario = "double-outer";
    const TriangleMesh mesh = join(
        makeBox(Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones()),
        makeBox(
            Eigen::Vector3d(2.0, 0.0, 0.0),
            Eigen::Vector3d(3.0, 1.0, 1.0)
        )
    );
    const BoundaryAudit audit = auditOutward(mesh);
    const bool acceptedAsSingleOuter =
        audit.strictBoundaryValid() && audit.faceComponentCount == 1;
    if (
        acceptedAsSingleOuter ||
        audit.faceComponentCount != 2 ||
        !audit.selfIntersections.intersectionFree
    )
    {
        return fail(scenario, "multiple outer components were not diagnosed");
    }
    return pass(scenario);
}


int runSelfIntersection()
{
    const std::string scenario = "self-intersection";
    const TriangleMesh mesh = join(
        makeBox(Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones()),
        makeBox(
            Eigen::Vector3d(0.5, 0.5, 0.5),
            Eigen::Vector3d(1.5, 1.5, 1.5)
        )
    );
    const BoundaryAudit audit = auditOutward(mesh);
    if (
        audit.strictBoundaryValid() ||
        !audit.selfIntersections.completed ||
        audit.selfIntersections.intersectionFree ||
        audit.selfIntersections.intersectingPairCount == 0
    )
    {
        return fail(scenario, "crossing closed shells escaped intersection audit");
    }
    return pass(scenario);
}


} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr
            << "usage: manufacturing-geometry-regression "
            << "<cube|missing|reversed|duplicate|degenerate|double-outer|"
            << "self-intersection>\n";
        return EXIT_FAILURE;
    }

    const std::string scenario = argv[1];
    if (scenario == "cube") return runCube();
    if (scenario == "missing") return runMissing();
    if (scenario == "reversed") return runReversed();
    if (scenario == "duplicate") return runDuplicate();
    if (scenario == "degenerate") return runDegenerate();
    if (scenario == "double-outer") return runDoubleOuter();
    if (scenario == "self-intersection") return runSelfIntersection();

    return fail(scenario, "unknown scenario");
}
