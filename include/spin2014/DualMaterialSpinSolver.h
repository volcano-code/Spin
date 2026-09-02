#pragma once

#include "physics/RawMassMoments10.h"
#include "spin2014/AdaptiveVoxelSpinSolver.h"

#include <string>
#include <vector>


namespace spin2014
{


struct DualMaterialSpinRequest
{
    AdaptiveVoxelSpinRequest hollowing;
    double boundaryDensity = 1.0;
    double interiorDensity = 8.0;
};


struct DualMaterialSpinResult
{
    AdaptiveVoxelSpinResult hollowing;
    physics::RawMassMoments10 boundaryContribution;
    physics::RawMassMoments10 interiorContribution;
    /** 0=void/exterior, 1=boundary material, 2=dense retained interior. */
    std::vector<unsigned char> materialLabels;
    std::string outputHash;
};


DualMaterialSpinResult solveDualMaterialSpin(
    const DualMaterialSpinRequest& request
);


} // namespace spin2014
