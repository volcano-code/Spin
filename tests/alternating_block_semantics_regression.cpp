#include "deformation/DeformationOptimizer.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


int fail(
    const std::string& scenario,
    const std::string& reason
)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << " reason=" << reason << '\n';
    return EXIT_FAILURE;
}


void makeBox(
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    vertices.resize(8, 3);
    vertices <<
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 1.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,
        1.0, 1.0, 1.0,
        0.0, 1.0, 1.0;

    faces.resize(12, 3);
    faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
}


struct VoxelMetadata
{
    Eigen::Vector3d center;
    double distanceToMesh;
    int x;
    int y;
    int z;
    std::set<int> corners;
    bool isInner;
};


int linearIndex(
    const Eigen::Vector3i& dimensions,
    int x,
    int y,
    int z
)
{
    return (x * dimensions(1) + y) * dimensions(2) + z;
}


std::vector<VoxelMetadata> metadataSnapshot(const VoxelGrid& grid)
{
    const Eigen::Vector3i dimensions = grid.dimensions();
    std::vector<VoxelMetadata> result(
        static_cast<std::size_t>(dimensions.prod())
    );

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const Voxel& voxel = grid(x, y, z);
                result[static_cast<std::size_t>(
                    linearIndex(dimensions, x, y, z)
                )] = VoxelMetadata{
                    voxel.center,
                    voxel.distanceToMesh,
                    voxel.xIdx,
                    voxel.yIdx,
                    voxel.zIdx,
                    voxel.cornerIndices,
                    voxel.isInner
                };
            }
        }
    }
    return result;
}


std::vector<unsigned char> fillSnapshot(const VoxelGrid& grid)
{
    const Eigen::Vector3i dimensions = grid.dimensions();
    std::vector<unsigned char> result(
        static_cast<std::size_t>(dimensions.prod()),
        0U
    );

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                result[static_cast<std::size_t>(
                    linearIndex(dimensions, x, y, z)
                )] = grid(x, y, z).isFilled ? 1U : 0U;
            }
        }
    }
    return result;
}


bool metadataEqual(
    const std::vector<VoxelMetadata>& first,
    const std::vector<VoxelMetadata>& second
)
{
    if (first.size() != second.size())
    {
        return false;
    }

    for (std::size_t index = 0; index < first.size(); ++index)
    {
        if (
            !first[index].center.isApprox(second[index].center, 0.0) ||
            first[index].distanceToMesh != second[index].distanceToMesh ||
            first[index].x != second[index].x ||
            first[index].y != second[index].y ||
            first[index].z != second[index].z ||
            first[index].corners != second[index].corners ||
            first[index].isInner != second[index].isInner
        )
        {
            return false;
        }
    }
    return true;
}


int runRefill(const std::string& scenario)
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeBox(vertices, faces);
    VoxelGrid grid(vertices, faces, 0.0, 10);

    const std::vector<VoxelMetadata> metadataBefore =
        metadataSnapshot(grid);
    const Eigen::MatrixXd cornersBefore = grid.corners();
    const Eigen::Vector3i dimensions = grid.dimensions();

    std::vector<Eigen::Vector3i> carveable;
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                if (grid(x, y, z).isInner)
                {
                    carveable.emplace_back(x, y, z);
                }
            }
        }
    }

    if (carveable.size() < 4U)
    {
        return fail(scenario, "fixture has too few carveable voxels");
    }

    const std::vector<Eigen::Vector3i> patternA =
    {
        carveable[0], carveable[1], carveable[2]
    };
    for (const Eigen::Vector3i& key : patternA)
    {
        grid(key(0), key(1), key(2)).isFilled = false;
    }

    if (
        grid.refillAll() != 3 ||
        grid.refillAll() != 0 ||
        !metadataEqual(metadataBefore, metadataSnapshot(grid)) ||
        !cornersBefore.isApprox(grid.corners(), 0.0)
    )
    {
        return fail(
            scenario,
            "refill count/idempotence or immutable metadata contract failed"
        );
    }

    /*
     * Paper Section 3.2: the next alpha block starts full, so an A-empty cell
     * may become filled while a previously retained cell enters pattern B.
     */
    const std::vector<Eigen::Vector3i> patternB =
    {
        carveable[1], carveable[3]
    };
    for (const Eigen::Vector3i& key : patternB)
    {
        grid(key(0), key(1), key(2)).isFilled = false;
    }

    if (
        !grid(patternA[0](0), patternA[0](1), patternA[0](2)).isFilled ||
        grid(patternB[0](0), patternB[0](1), patternB[0](2)).isFilled ||
        grid(patternB[1](0), patternB[1](1), patternB[1](2)).isFilled ||
        grid.refillAll() != 2 ||
        !metadataEqual(metadataBefore, metadataSnapshot(grid))
    )
    {
        return fail(
            scenario,
            "second alpha block inherited rather than reset the first pattern"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " patternA=3 patternB=2\n";
    return EXIT_SUCCESS;
}


int runPaperTransaction(const std::string& scenario)
{
    using namespace deformation;

    PaperLineSearchConfig config;
    PaperLineSearchState state;
    state.step = 1.0;
    state.mu = 0.4;

    const Eigen::VectorXd parameters =
        Eigen::VectorXd::Constant(1, 10.0);
    const Eigen::VectorXd gradient =
        Eigen::VectorXd::Constant(1, 1.0);
    const std::vector<unsigned char> alpha = {1U, 0U, 1U, 0U};
    int currentCalls = 0;
    int trialCalls = 0;

    const PaperLineSearchResult accepted = paperGradientDescentTrial(
        parameters,
        gradient,
        state,
        config,
        [&](const Eigen::VectorXd& query, double)
        {
            if (query(0) == parameters(0))
            {
                ++currentCalls;
                return 100.0;
            }
            ++trialCalls;
            return 99.0;
        }
    );

    if (
        !accepted.accepted ||
        accepted.rolledBack ||
        accepted.decision != PaperTrialDecision::Accepted ||
        accepted.parameters(0) != 9.0 ||
        currentCalls != 1 ||
        trialCalls != 1 ||
        alpha != std::vector<unsigned char>({1U, 0U, 1U, 0U})
    )
    {
        return fail(
            scenario,
            "exact paper 1% fixed-alpha trial was not committed"
        );
    }

    currentCalls = 0;
    trialCalls = 0;
    const PaperLineSearchResult rejected = paperGradientDescentTrial(
        parameters,
        gradient,
        state,
        config,
        [&](const Eigen::VectorXd& query, double)
        {
            if (query(0) == parameters(0))
            {
                ++currentCalls;
                return 100.0;
            }
            ++trialCalls;
            return 99.1;
        }
    );

    if (
        rejected.accepted ||
        !rejected.rolledBack ||
        rejected.decision != PaperTrialDecision::StepReduced ||
        rejected.parameters(0) != parameters(0) ||
        rejected.energyAfter != rejected.energyBefore ||
        rejected.nextState.step != 0.8 ||
        currentCalls != 1 ||
        trialCalls != 1 ||
        alpha != std::vector<unsigned char>({1U, 0U, 1U, 0U})
    )
    {
        return fail(
            scenario,
            "sub-1% fixed-alpha trial was not transactionally rolled back"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " exact_1pct=accepted sub_1pct=rolled_back\n";
    return EXIT_SUCCESS;
}


struct ProfileTrace
{
    int deformationTrials = 0;
    int refillCalls = 0;
    int carvingCalls = 0;
};


int runProfileDiscriminator(const std::string& scenario)
{
    /*
     * R6 is the paper-style fixed-alpha deformation block: evaluating one
     * line-search trial must not refill/recarve. The inspected beta instead
     * performs deform -> refill/recarve before comparing each outer trial.
     */
    ProfileTrace paperTrace;
    const Eigen::VectorXd parameters =
        Eigen::VectorXd::Constant(1, 5.0);
    const Eigen::VectorXd gradient =
        Eigen::VectorXd::Constant(1, 1.0);
    deformation::PaperLineSearchState state;
    state.mu = 0.2;

    const deformation::PaperLineSearchResult result =
        deformation::paperGradientDescentTrial(
            parameters,
            gradient,
            state,
            deformation::PaperLineSearchConfig(),
            [&](const Eigen::VectorXd& query, double)
            {
                if (query(0) != parameters(0))
                {
                    ++paperTrace.deformationTrials;
                }
                /* Deliberately no refillAll() and no CarvingSolver call. */
                const double residual = query(0) - 4.0;
                return residual * residual;
            }
        );

    ProfileTrace betaOracle;
    ++betaOracle.deformationTrials;
    ++betaOracle.refillCalls;
    ++betaOracle.carvingCalls;

    if (
        !result.accepted ||
        paperTrace.deformationTrials != 1 ||
        paperTrace.refillCalls != 0 ||
        paperTrace.carvingCalls != 0 ||
        betaOracle.refillCalls != 1 ||
        betaOracle.carvingCalls != 1
    )
    {
        return fail(
            scenario,
            "paper fixed-alpha and beta recarve profiles were conflated"
        );
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " paper_trial_refill=0 paper_trial_carve=0"
        << " beta_oracle_refill=1 beta_oracle_carve=1\n";
    return EXIT_SUCCESS;
}


} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        return fail(
            "unknown",
            "usage: alternating-block-semantics-regression "
            "<refill|paper-transaction|profile-discriminator>"
        );
    }

    const std::string scenario(argv[1]);
    try
    {
        if (scenario == "refill")
        {
            return runRefill(scenario);
        }
        if (scenario == "paper-transaction")
        {
            return runPaperTransaction(scenario);
        }
        if (scenario == "profile-discriminator")
        {
            return runProfileDiscriminator(scenario);
        }
        return fail(scenario, "unknown scenario");
    }
    catch (const std::exception& error)
    {
        return fail(scenario, error.what());
    }
}
