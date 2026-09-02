#include "viewer/ResultView.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

int fail(const std::string& message)
{
    std::cerr << "[FAIL] result-view: " << message << '\n';
    return 1;
}

bool exactlyEqual(const Eigen::MatrixXd& lhs, const Eigen::MatrixXd& rhs)
{
    return lhs.rows() == rhs.rows() && lhs.cols() == rhs.cols() &&
        (lhs.array() == rhs.array()).all();
}

} // namespace

int main()
{
    const char* spinArguments[] = {
        "make-it-stand", "--spin-result", "outer.obj", "design.obj",
        "material.obj", "cavity.obj", "diagnostics.json", "Fidget result"
    };
    const result_view::ResultViewArguments parsedSpin =
        result_view::parseResultViewArguments(8, spinArguments);
    if (!parsedSpin.enabled || !parsedSpin.spinEnabled ||
        parsedSpin.designPath != "design.obj" ||
        parsedSpin.materialPath != "material.obj" ||
        parsedSpin.innerPath != "cavity.obj" ||
        parsedSpin.diagnosticsPath != "diagnostics.json" ||
        parsedSpin.title != "Fidget result")
    {
        return fail("--spin-result arguments were not preserved");
    }

    const char* arguments[] = {
        "make-it-stand",
        "--result",
        "sealed/outer-flattened.obj",
        "sealed/inner-cavity.obj",
        "Head-down bunny"
    };
    const result_view::ResultViewArguments parsed =
        result_view::parseResultViewArguments(5, arguments);
    if (!parsed.enabled ||
        parsed.outerPath != "sealed/outer-flattened.obj" ||
        parsed.innerPath != "sealed/inner-cavity.obj" ||
        parsed.title != "Head-down bunny")
    {
        return fail("--result arguments were not preserved");
    }

    bool rejectedIncomplete = false;
    try
    {
        const char* incomplete[] = {"make-it-stand", "--result", "outer.obj"};
        (void)result_view::parseResultViewArguments(3, incomplete);
    }
    catch (const std::invalid_argument&)
    {
        rejectedIncomplete = true;
    }
    if (!rejectedIncomplete)
    {
        return fail("incomplete --result arguments were accepted");
    }

    Eigen::MatrixXd outer(4, 3);
    outer <<
        -4.0, -3.0,  2.0,
        -1.0, -2.0,  2.5,
        -3.0,  5.0,  8.0,
         2.0,  1.0, -6.0;
    Eigen::MatrixXi outerFaces(2, 3);
    outerFaces << 0, 1, 2, 0, 3, 1;

    Eigen::MatrixXd inner(3, 3);
    inner <<
        -2.5, -1.5,  1.0,
        -2.0,  0.0,  2.0,
        -1.5,  1.5, -1.0;
    Eigen::MatrixXi innerFaces(1, 3);
    innerFaces << 0, 1, 2;

    const result_view::ResultViewGeometry geometry =
        result_view::prepareResultViewGeometry(
            outer,
            outerFaces,
            inner,
            innerFaces,
            0.10
        );

    // This catches the original viewer's alignToAxis path: a sealed result is
    // already in its experiment posture and both surfaces must remain in the
    // exact same coordinate frame.
    if (!exactlyEqual(geometry.outerVertices, outer) ||
        !exactlyEqual(geometry.innerVertices, inner))
    {
        return fail("sealed experiment posture or cavity registration changed");
    }
    if (geometry.groundVertices.rows() != 4 ||
        geometry.groundFaces.rows() != 2 ||
        std::abs(geometry.groundHeight + 3.0) > 1.0e-15 ||
        (geometry.groundVertices.col(1).array() != -3.0).any())
    {
        return fail("ground plane does not coincide with the sealed support height");
    }
    if (!(geometry.groundVertices.col(0).minCoeff() < outer.col(0).minCoeff()) ||
        !(geometry.groundVertices.col(0).maxCoeff() > outer.col(0).maxCoeff()) ||
        !(geometry.groundVertices.col(2).minCoeff() < outer.col(2).minCoeff()) ||
        !(geometry.groundVertices.col(2).maxCoeff() > outer.col(2).maxCoeff()))
    {
        return fail("ground plane does not extend beyond the outer silhouette");
    }
    const Eigen::Vector3d edgeA =
        (geometry.groundVertices.row(geometry.groundFaces(0, 1)) -
            geometry.groundVertices.row(geometry.groundFaces(0, 0))).transpose();
    const Eigen::Vector3d edgeB =
        (geometry.groundVertices.row(geometry.groundFaces(0, 2)) -
            geometry.groundVertices.row(geometry.groundFaces(0, 0))).transpose();
    const Eigen::Vector3d groundNormal = edgeA.cross(edgeB);
    if (!(groundNormal(1) > 0.0))
    {
        return fail("ground plane faces away from the visible upper side");
    }

    result_view::SpinResultViewMetadata metadata;
    metadata.mode = "STATIC_BALANCE";
    metadata.centerOfMass = Eigen::Vector3d(0.1, 0.2, 0.3);
    metadata.targetOrigin = Eigen::Vector3d(0.0, 0.0, -1.0);
    metadata.targetAxes = Eigen::Matrix3d::Identity();
    metadata.principalAxes = Eigen::Matrix3d::Identity();
    const result_view::SpinResultViewGeometry spinGeometry =
        result_view::prepareSpinResultViewGeometry(
            outer, outerFaces, outer, outerFaces,
            inner, innerFaces, inner, innerFaces, metadata
        );
    if (!exactlyEqual(spinGeometry.outerVertices, outer) ||
        !exactlyEqual(spinGeometry.designVertices, outer) ||
        !exactlyEqual(spinGeometry.materialVertices, inner) ||
        !exactlyEqual(spinGeometry.cavityVertices, inner) ||
        spinGeometry.metadata.mode != "STATIC_BALANCE")
    {
        return fail("Spin-It layers changed coordinate frames");
    }
    const Eigen::Vector3d targetEdgeA =
        (spinGeometry.targetPlaneVertices.row(1) -
         spinGeometry.targetPlaneVertices.row(0)).transpose();
    const Eigen::Vector3d targetEdgeB =
        (spinGeometry.targetPlaneVertices.row(2) -
         spinGeometry.targetPlaneVertices.row(0)).transpose();
    if (targetEdgeA.cross(targetEdgeB).dot(metadata.targetAxes.col(2)) <= 0.0)
    {
        return fail("Spin-It target plane does not face along the target axis");
    }

    std::cout << "[PASS] result-view preserves all sealed layer registrations\n";
    return 0;
}
