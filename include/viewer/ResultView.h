#pragma once

#include <Eigen/Core>

#include <string>

namespace result_view
{

struct ResultViewArguments
{
    bool enabled = false;
    bool spinEnabled = false;
    std::string outerPath;
    std::string innerPath;
    std::string designPath;
    std::string materialPath;
    std::string diagnosticsPath;
    std::string title = "Make It Stand — sealed result";
};


struct SpinResultViewMetadata
{
    std::string artifactFamily = "Spin-It Faster";
    std::string cavityLayerName = "Cavity / signed quadric";
    std::string planeLayerName = "Target/support plane";
    std::string mode;
    bool buoyancyAvailable = false;
    double displacedLiquidVolume = 0.0;
    double buoyantMetacenterHeight = 0.0;
    double buoyantStabilityMargin = 0.0;
    Eigen::Vector3d centerOfMass = Eigen::Vector3d::Zero();
    Eigen::Vector3d targetOrigin = Eigen::Vector3d::Zero();
    Eigen::Matrix3d targetAxes = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d principalAxes = Eigen::Matrix3d::Identity();
};


struct SpinResultViewGeometry
{
    Eigen::MatrixXd outerVertices;
    Eigen::MatrixXi outerFaces;
    Eigen::MatrixXd designVertices;
    Eigen::MatrixXi designFaces;
    Eigen::MatrixXd materialVertices;
    Eigen::MatrixXi materialFaces;
    Eigen::MatrixXd cavityVertices;
    Eigen::MatrixXi cavityFaces;
    Eigen::MatrixXd targetPlaneVertices;
    Eigen::MatrixXi targetPlaneFaces;
    SpinResultViewMetadata metadata;
    double axisDisplayLength = 1.0;
};


SpinResultViewMetadata loadSpinResultViewMetadata(
    const std::string& diagnosticsPath
);


SpinResultViewGeometry prepareSpinResultViewGeometry(
    const Eigen::MatrixXd& outerVertices,
    const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& designVertices,
    const Eigen::MatrixXi& designFaces,
    const Eigen::MatrixXd& materialVertices,
    const Eigen::MatrixXi& materialFaces,
    const Eigen::MatrixXd& cavityVertices,
    const Eigen::MatrixXi& cavityFaces,
    const SpinResultViewMetadata& metadata,
    double displayScaleFraction = 0.35
);

ResultViewArguments parseResultViewArguments(
    int argc,
    const char* const argv[]
);

struct ResultViewGeometry
{
    Eigen::MatrixXd outerVertices;
    Eigen::MatrixXi outerFaces;
    Eigen::MatrixXd innerVertices;
    Eigen::MatrixXi innerFaces;
    Eigen::MatrixXd groundVertices;
    Eigen::MatrixXi groundFaces;
    double groundHeight = 0.0;
};

ResultViewGeometry prepareResultViewGeometry(
    const Eigen::MatrixXd& outerVertices,
    const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& innerVertices,
    const Eigen::MatrixXi& innerFaces,
    double groundMarginFraction = 0.10
);

} // namespace result_view
