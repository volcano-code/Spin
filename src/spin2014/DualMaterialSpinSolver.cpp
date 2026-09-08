#include "spin2014/DualMaterialSpinSolver.h"

#include "manufacturing/SlicerValidation.h"
#include "spin2014/AdaptiveVoxelDomain.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>


namespace spin2014
{


DualMaterialSpinResult solveDualMaterialSpin(
    const DualMaterialSpinRequest& request
)
{
    if (request.hollowing.formulation == Spin2014Formulation::PaperGivens)
    {
        throw std::invalid_argument(
            "paper Givens profile is not validated for the dual-material product/export wrapper"
        );
    }
    if (
        !std::isfinite(request.boundaryDensity) ||
        request.boundaryDensity <= 0.0 ||
        !std::isfinite(request.interiorDensity) ||
        request.interiorDensity <= 0.0
    )
    {
        throw std::invalid_argument(
            "dual-material densities must be finite and positive"
        );
    }
    AdaptiveVoxelSpinRequest weighted = request.hollowing;
    weighted.boundaryDensity = request.boundaryDensity;
    weighted.interiorDensity = request.interiorDensity;

    DualMaterialSpinResult result;
    result.hollowing = solveAdaptiveVoxelSpin(weighted);
    if (result.hollowing.status != AdaptiveSpinStatus::Converged)
    {
        result.outputHash = manufacturing::computeContentSha256(
            std::string("spin-it-2014-dual-material-failure-v1") +
            result.hollowing.outputHash
        );
        return result;
    }

    AdaptiveVoxelDomain original(
        weighted.domain,
        weighted.maximumLevel,
        weighted.maximumLevel
    );
    const physics::RawMassMoments10 boundaryUnit =
        original.solidMoments() - original.variableMoments();
    result.boundaryContribution =
        boundaryUnit.scaledDensity(request.boundaryDensity);
    result.interiorContribution =
        result.hollowing.finalMoments - result.boundaryContribution;

    result.materialLabels.resize(weighted.domain->cellCount(), 0U);
    for (std::size_t index = 0; index < result.materialLabels.size(); ++index)
    {
        if (!weighted.domain->material(index))
        {
            continue;
        }
        if (!weighted.domain->carveable(index))
        {
            result.materialLabels[index] = 1U;
        }
        else if (result.hollowing.binaryVoidMask[index] == 0U)
        {
            result.materialLabels[index] = 2U;
        }
    }
    std::string bytes("spin-it-2014-dual-material-result-v1");
    bytes += result.hollowing.outputHash;
    bytes.append(
        reinterpret_cast<const char*>(result.materialLabels.data()),
        result.materialLabels.size()
    );
    result.outputHash = manufacturing::computeContentSha256(bytes);
    return result;
}


} // namespace spin2014
