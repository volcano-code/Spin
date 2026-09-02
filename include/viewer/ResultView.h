#pragma once

#include <Eigen/Core>

#include <string>

namespace result_view
{

struct ResultViewArguments
{
    bool enabled = false;
    std::string outerPath;
    std::string innerPath;
    std::string title = "Make It Stand — sealed result";
};

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
