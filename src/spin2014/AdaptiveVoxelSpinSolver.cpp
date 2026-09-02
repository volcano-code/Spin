#include "spin2014/AdaptiveVoxelSpinSolver.h"

#include "manufacturing/SlicerValidation.h"
#include "spin2014/AdaptiveFillFunctional.h"
#include "spin2014/AdaptiveVoxelDomain.h"
#include "spin2014/SpinIt2014MomentProblem.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace
{


void appendInt(std::string& bytes, int value)
{
    bytes.append(reinterpret_cast<const char*>(&value), sizeof(value));
}


void appendDouble(std::string& bytes, double value)
{
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    for (int shift = 0; shift < 64; shift += 8)
    {
        bytes.push_back(static_cast<char>((bits >> shift) & 0xffU));
    }
}


spin2014::SpinIt2014MomentProblem momentProblem(
    const spin2014::AdaptiveVoxelSpinRequest& request
)
{
    switch (request.mode)
    {
    case spin::SpinMode::StaticBalance:
        return spin2014::SpinIt2014MomentProblem::staticBalance(
            request.centerOfMassWeight
        );
    case spin::SpinMode::YoYo:
        return spin2014::SpinIt2014MomentProblem::yoYo();
    case spin::SpinMode::SpinningTop:
        return spin2014::SpinIt2014MomentProblem::spinningTop(
            request.centerOfMassWeight
        );
    case spin::SpinMode::FidgetSpinner:
        break;
    }
    throw std::invalid_argument(
        "the 2014 paper does not define a fidget-spinner objective"
    );
}


std::string hashTrace(const std::vector<spin2014::AdaptiveSpinPass>& trace)
{
    std::string bytes("spin-it-2014-adaptive-trace-v1");
    for (const spin2014::AdaptiveSpinPass& pass : trace)
    {
        appendInt(bytes, pass.pass);
        appendInt(bytes, pass.leafCountBefore);
        appendInt(bytes, pass.fractionalLeafCount);
        appendInt(bytes, pass.splitLeafCount);
        appendInt(bytes, pass.mergedParentCount);
        bytes += pass.partitionHash;
        bytes += pass.slqpTraceHash;
    }
    return manufacturing::computeContentSha256(bytes);
}


} // namespace


namespace spin2014
{


const char* adaptiveSpinStatusToken(AdaptiveSpinStatus status) noexcept
{
    switch (status)
    {
    case AdaptiveSpinStatus::Converged: return "CONVERGED";
    case AdaptiveSpinStatus::RelaxationDidNotConverge:
        return "RELAXATION_DID_NOT_CONVERGE";
    case AdaptiveSpinStatus::AdaptivePassLimit: return "ADAPTIVE_PASS_LIMIT";
    case AdaptiveSpinStatus::BinaryValidationFailed:
        return "BINARY_VALIDATION_FAILED";
    case AdaptiveSpinStatus::InvalidRequest: return "INVALID_REQUEST";
    }
    return "INVALID_REQUEST";
}


void validateAdaptiveVoxelSpinRequest(const AdaptiveVoxelSpinRequest& request)
{
    if (!request.domain)
    {
        throw std::invalid_argument("Spin-It 2014 domain must not be null");
    }
    if (request.mode == spin::SpinMode::FidgetSpinner)
    {
        throw std::invalid_argument(
            "Spin-It 2014 supports static, yo-yo and top modes"
        );
    }
    if (
        request.initialLevel < 0 ||
        request.maximumLevel < request.initialLevel ||
        request.maximumLevel >= 30 ||
        request.maximumAdaptivePasses <= 0 ||
        !std::isfinite(request.initialFill) || request.initialFill < 0.0 ||
        request.initialFill > 1.0 ||
        !std::isfinite(request.fractionalEpsilon) ||
        request.fractionalEpsilon <= 0.0 ||
        request.fractionalEpsilon >= 0.5 ||
        !std::isfinite(request.laplacianWeight) ||
        request.laplacianWeight < 0.0 ||
        !std::isfinite(request.centerOfMassWeight) ||
        request.centerOfMassWeight < 0.0 ||
        !std::isfinite(request.boundaryDensity) ||
        request.boundaryDensity <= 0.0 ||
        !std::isfinite(request.interiorDensity) ||
        request.interiorDensity <= 0.0 ||
        !std::isfinite(request.binaryConstraintTolerance) ||
        request.binaryConstraintTolerance <= 0.0
    )
    {
        throw std::invalid_argument("invalid Spin-It 2014 request settings");
    }
    const int rootSide = 1 << request.maximumLevel;
    if ((request.domain->grid().dimensions.array() > rootSide).any())
    {
        throw std::invalid_argument(
            "maximum octree level cannot contain the source grid"
        );
    }
    (void)momentProblem(request);
}


std::string computeAdaptiveVoxelSpinRequestHash(
    const AdaptiveVoxelSpinRequest& request
)
{
    validateAdaptiveVoxelSpinRequest(request);
    std::string bytes("spin-it-2014-adaptive-request-v1");
    bytes += request.domain->artifactId();
    appendInt(bytes, static_cast<int>(request.mode));
    appendInt(bytes, request.initialLevel);
    appendInt(bytes, request.maximumLevel);
    appendInt(bytes, request.maximumAdaptivePasses);
    appendDouble(bytes, request.initialFill);
    appendDouble(bytes, request.fractionalEpsilon);
    appendDouble(bytes, request.laplacianWeight);
    appendDouble(bytes, request.centerOfMassWeight);
    appendDouble(bytes, request.boundaryDensity);
    appendDouble(bytes, request.interiorDensity);
    appendDouble(bytes, request.binaryConstraintTolerance);
    appendInt(bytes, request.slqp.maximumIterations);
    appendInt(bytes, request.slqp.maximumLineSearchSteps);
    appendInt(bytes, request.slqp.lbfgsMemory);
    appendDouble(bytes, request.slqp.equalityTolerance);
    appendDouble(bytes, request.slqp.projectedGradientTolerance);
    appendDouble(bytes, request.slqp.activeSetTolerance);
    appendDouble(bytes, request.slqp.meritPenalty);
    appendDouble(bytes, request.slqp.lineSearchReduction);
    appendDouble(bytes, request.slqp.armijo);
    return manufacturing::computeContentSha256(bytes);
}


AdaptiveVoxelSpinResult solveAdaptiveVoxelSpin(
    const AdaptiveVoxelSpinRequest& request
)
{
    validateAdaptiveVoxelSpinRequest(request);
    AdaptiveVoxelSpinResult result;
    result.inputHash = computeAdaptiveVoxelSpinRequestHash(request);
    AdaptiveVoxelDomain domain(
        request.domain,
        request.initialLevel,
        request.maximumLevel
    );
    Eigen::VectorXd beta = Eigen::VectorXd::Constant(
        static_cast<Eigen::Index>(domain.leaves().size()),
        request.initialFill
    );
    bool adaptiveConverged = false;

    for (int passIndex = 0;
         passIndex < request.maximumAdaptivePasses;
         ++passIndex)
    {
        const AdaptiveFillFunctional functional(
            domain,
            momentProblem(request),
            request.laplacianWeight,
            request.boundaryDensity,
            request.interiorDensity
        );
        const SlqpLbfgsResult relaxed = solveSlqpLbfgs(
            functional,
            beta,
            request.slqp
        );
        beta = relaxed.beta;

        AdaptiveSpinPass pass;
        pass.pass = passIndex;
        pass.leafCountBefore = static_cast<int>(domain.leaves().size());
        pass.partitionHash = domain.partitionHash();
        pass.slqpTraceHash = relaxed.traceHash;
        for (std::size_t index = 0; index < domain.leaves().size(); ++index)
        {
            const double value = beta(static_cast<Eigen::Index>(index));
            if (
                value >= request.fractionalEpsilon &&
                value <= 1.0 - request.fractionalEpsilon
            )
            {
                ++pass.fractionalLeafCount;
            }
        }
        if (relaxed.status != SlqpStatus::Converged)
        {
            result.adaptiveTrace.push_back(pass);
            result.relaxedFill = beta;
            result.status = AdaptiveSpinStatus::RelaxationDidNotConverge;
            result.diagnostic = std::string("adaptive pass failed: ") +
                relaxed.diagnostic;
            result.adaptiveTraceHash = hashTrace(result.adaptiveTrace);
            result.outputHash = manufacturing::computeContentSha256(
                result.inputHash + result.adaptiveTraceHash + result.diagnostic
            );
            return result;
        }

        pass.splitLeafCount = static_cast<int>(
            domain.splitFractionalLeaves(beta, request.fractionalEpsilon)
        );
        pass.mergedParentCount = static_cast<int>(
            domain.mergeBinarySiblings(beta, request.fractionalEpsilon)
        );
        result.adaptiveTrace.push_back(pass);
        if (pass.splitLeafCount == 0)
        {
            adaptiveConverged = true;
            break;
        }
    }
    result.relaxedFill = beta;
    result.adaptiveTraceHash = hashTrace(result.adaptiveTrace);
    if (!adaptiveConverged)
    {
        result.status = AdaptiveSpinStatus::AdaptivePassLimit;
        result.diagnostic = "adaptive refinement reached its pass limit";
        result.outputHash = manufacturing::computeContentSha256(
            result.inputHash + result.adaptiveTraceHash + result.diagnostic
        );
        return result;
    }

    Eigen::VectorXd binary(beta.size());
    result.binaryFill.resize(static_cast<std::size_t>(beta.size()));
    for (Eigen::Index index = 0; index < beta.size(); ++index)
    {
        binary(index) = beta(index) >= 0.5 ? 1.0 : 0.0;
        result.binaryFill[static_cast<std::size_t>(index)] = binary(index);
    }
    const AdaptiveFillFunctional finalFunctional(
        domain,
        momentProblem(request),
        request.laplacianWeight,
        request.boundaryDensity,
        request.interiorDensity
    );
    const AdaptiveFillEvaluation finalEvaluation =
        finalFunctional.evaluate(binary);
    result.finalMoments = finalEvaluation.finalMoments;
    result.binaryConstraints = finalEvaluation.constraints;
    result.binaryConstraintInfinityNorm = result.binaryConstraints.size() > 0
        ? result.binaryConstraints.lpNorm<Eigen::Infinity>()
        : 0.0;
    result.binaryVoidMask = domain.binaryVoidMask(binary);

    std::vector<unsigned char> finalMaterial = request.domain->materialMask();
    for (std::size_t index = 0; index < finalMaterial.size(); ++index)
    {
        if (result.binaryVoidMask[index] != 0U)
        {
            finalMaterial[index] = 0U;
        }
    }
    const std::shared_ptr<const voxel::VoxelDomain> finalDomain =
        voxel::VoxelDomain::create(
            request.domain->grid(),
            finalMaterial,
            finalMaterial,
            finalMaterial
        );
    result.materialBoundary = voxel::extractMaterialBoundaries(finalDomain);
    AdaptiveVoxelDomain finalMomentDomain(
        finalDomain,
        request.maximumLevel,
        request.maximumLevel
    );
    const physics::RawMassMoments10 boundaryMoments =
        domain.solidMoments() - domain.variableMoments();
    const physics::RawMassMoments10 retainedInteriorMoments =
        finalMomentDomain.solidMoments() - boundaryMoments;
    const physics::RawMassMoments10 independentlyWeighted =
        boundaryMoments.scaledDensity(request.boundaryDensity) +
        retainedInteriorMoments.scaledDensity(request.interiorDensity);
    if (!independentlyWeighted.values().isApprox(
            result.finalMoments.values(), 1e-10
        ))
    {
        result.status = AdaptiveSpinStatus::BinaryValidationFailed;
        result.diagnostic = "binary geometry and analytic moments disagree";
    }
    else if (
        result.binaryConstraintInfinityNorm >
        request.binaryConstraintTolerance
    )
    {
        result.status = AdaptiveSpinStatus::BinaryValidationFailed;
        result.diagnostic = "binary rounding violates the equality constraints";
    }
    else
    {
        result.status = AdaptiveSpinStatus::Converged;
        result.diagnostic =
            "adaptive relaxation, binary rounding and geometry audit passed";
    }

    result.geometryHash = finalDomain->artifactId();
    std::string outputBytes = result.inputHash + result.adaptiveTraceHash +
        result.geometryHash + result.diagnostic;
    for (double value : result.binaryFill)
    {
        appendDouble(outputBytes, value);
    }
    for (int coordinate = 0;
         coordinate < physics::RawMassMoments10::Count;
         ++coordinate)
    {
        appendDouble(outputBytes, result.finalMoments.values()(coordinate));
    }
    result.outputHash = manufacturing::computeContentSha256(outputBytes);
    return result;
}


} // namespace spin2014
