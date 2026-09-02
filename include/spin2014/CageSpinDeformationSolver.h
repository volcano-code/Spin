#pragma once

#include "spin2014/CageSpinDeformation.h"

#include <Eigen/Core>

#include <string>
#include <vector>


namespace spin2014
{


enum class CageSpinStatus
{
    Converged,
    MaximumIterations,
    LineSearchExhausted,
    InvalidGeometry
};


const char* cageSpinStatusToken(CageSpinStatus status) noexcept;


struct CageSpinDeformationSolverSettings
{
    int maximumIterations = 120;
    int maximumLineSearchSteps = 24;
    int lbfgsMemory = 8;
    double equalityTolerance = 1e-8;
    double projectedGradientTolerance = 1e-7;
    double meritPenalty = 100.0;
    double lineSearchReduction = 0.5;
    double maximumTranslation = 0.25;
    double minimumScale = 0.8;
    double maximumScale = 1.2;
};


struct CageSpinIteration
{
    int iteration = 0;
    double objective = 0.0;
    double equalityInfinityNorm = 0.0;
    double projectedGradientInfinityNorm = 0.0;
    double acceptedStep = 0.0;
    int activeVariables = 0;
    int lbfgsHistorySize = 0;
    bool geometryAccepted = false;
};


struct CageSpinDeformationResult
{
    CageSpinStatus status = CageSpinStatus::InvalidGeometry;
    std::string diagnostic;
    Eigen::VectorXd beta;
    Eigen::VectorXd handleParameters;
    CageSpinDeformationEvaluation evaluation;
    bool geometryValid = false;
    bool silhouetteChanged = false;
    std::vector<CageSpinIteration> trace;
    std::string traceHash;
    std::string outputHash;
};


CageSpinDeformationResult solveCageSpinDeformation(
    const CageSpinDeformationFunctional& functional,
    const Eigen::VectorXd& initialBeta,
    const CageSpinDeformationSolverSettings& settings = {}
);


} // namespace spin2014
