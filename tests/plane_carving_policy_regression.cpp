#include "carving/PlaneCarvingPolicy.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>


namespace
{

bool close(double actual, double expected)
{
    const double scale = std::max(
        std::abs(actual),
        std::abs(expected)
    );

    return
        std::isfinite(actual) &&
        std::isfinite(expected) &&
        std::abs(actual - expected) <=
            1e-12 + 1e-10 * scale;
}


int fail(
    const std::string& scenario,
    const std::string& message
)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << ' ' << message << '\n';
    return EXIT_FAILURE;
}


PlaneCarvingCandidate candidate(
    int x,
    const Eigen::Vector3d& center,
    double mass,
    double score
)
{
    PlaneCarvingCandidate result;
    result.key = {x, 0, 0};
    result.center = center;
    result.massProperties = MassProperties(
        mass,
        mass * center
    );
    result.score = score;
    return result;
}


PlaneCarvingCandidate keyedCandidate(
    int x,
    int y,
    int z,
    double score
)
{
    PlaneCarvingCandidate result = candidate(
        x,
        Eigen::Vector3d::UnitX(),
        1.0,
        score
    );
    result.key = {x, y, z};
    return result;
}


int runPlanes()
{
    const std::string scenario = "planes";
    BalanceObjectives objectives;
    objectives.emplace_back(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(0.0, 0.0, -7.0)
    );
    objectives.emplace_back(
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, -1.0)
    );

    const CuttingPlanes planes = makeCuttingPlanes(
        Eigen::Vector3d(1.0, 0.0, 4.0),
        objectives
    );

    if (
        planes.size() != 2 ||
        !planes[0].active ||
        !close(planes[0].normal(0), 1.0) ||
        !close(planes[0].normal(1), 0.0) ||
        !close(planes[0].normal.dot(planes[0].gravity), 0.0) ||
        planes[1].active ||
        !close(
            aggregatePositivePlaneScore(
                Eigen::Vector3d(3.0, 10.0, -9.0),
                planes
            ),
            3.0
        ) ||
        !close(
            aggregatePositivePlaneScore(
                Eigen::Vector3d(-3.0, 10.0, -9.0),
                planes
            ),
            0.0
        )
    )
    {
        return fail(scenario, "analytic Equation (8) plane contract failed");
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


int runMultiScore()
{
    const std::string scenario = "multi-score";
    BalanceObjectives objectives;
    objectives.emplace_back(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(0.0, 1.0, 1.0)
    );
    objectives.emplace_back(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(1.0, 0.0, 1.0)
    );
    objectives.emplace_back(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(1.0, 1.0, 0.0)
    );

    const CuttingPlanes planes = makeCuttingPlanes(
        Eigen::Vector3d(1.0, 1.0, 1.0),
        objectives
    );

    const double a = aggregatePositivePlaneScore(
        Eigen::Vector3d(3.0, -2.0, -2.0),
        planes
    );
    const double b = aggregatePositivePlaneScore(
        Eigen::Vector3d(-1.0, 4.0, -3.0),
        planes
    );
    const double c = aggregatePositivePlaneScore(
        Eigen::Vector3d(2.0, 1.0, 1.0),
        planes
    );
    const double d = aggregatePositivePlaneScore(
        Eigen::Vector3d(-1.0, -1.0, -1.0),
        planes
    );

    if (
        !close(a, 3.0) ||
        !close(b, 4.0) ||
        !close(c, 4.0) ||
        !close(d, 0.0)
    )
    {
        return fail(
            scenario,
            "score is not the sum of independently positive distances"
        );
    }

    /*
     * Printed Equation (8) uses the raw projected residual as its normal.
     * Normalizing each plane would reverse this two-candidate order:
     * raw(A,B)=(30,100), normalized(A,B)=(21,10).
     */
    BalanceObjectives rawObjectives;
    rawObjectives.emplace_back(
        Eigen::Vector3d(-10.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, -1.0)
    );
    rawObjectives.emplace_back(
        Eigen::Vector3d(0.0, -1.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, -1.0)
    );

    const CuttingPlanes rawPlanes = makeCuttingPlanes(
        Eigen::Vector3d::Zero(),
        rawObjectives
    );
    const double rawA = aggregatePositivePlaneScore(
        Eigen::Vector3d(-9.0, 19.0, 0.0),
        rawPlanes
    );
    const double rawB = aggregatePositivePlaneScore(
        Eigen::Vector3d(0.0, -1.0, 0.0),
        rawPlanes
    );

    if (!close(rawA, 30.0) || !close(rawB, 100.0) || rawB <= rawA)
    {
        return fail(
            scenario,
            "multi-target score normalized the raw Equation (8) residual"
        );
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


int runOrdering()
{
    const std::string scenario = "ordering";
    PlaneCarvingCandidates candidates;
    candidates.push_back(keyedCandidate(9, 0, 0, 2.0));
    candidates.push_back(keyedCandidate(3, 5, 9, 2.0));
    candidates.push_back(keyedCandidate(8, 0, 0, 4.0));
    candidates.push_back(keyedCandidate(3, 4, 9, 2.0));
    candidates.push_back(keyedCandidate(3, 4, 2, 2.0));

    sortPlaneCarvingCandidates(candidates);

    if (
        candidates[0].key.x != 8 ||
        candidates[1].key.x != 3 ||
        candidates[1].key.y != 4 ||
        candidates[1].key.z != 2 ||
        candidates[2].key.x != 3 ||
        candidates[2].key.y != 4 ||
        candidates[2].key.z != 9 ||
        candidates[3].key.x != 3 ||
        candidates[3].key.y != 5 ||
        candidates[4].key.x != 9
    )
    {
        return fail(scenario, "score/key deterministic order failed");
    }

    bool rejected = false;
    candidates[0].score =
        std::numeric_limits<double>::quiet_NaN();

    try
    {
        sortPlaneCarvingCandidates(candidates);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }

    if (!rejected)
    {
        return fail(scenario, "non-finite score was accepted");
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


int runTransient()
{
    const std::string scenario = "transient-prefix";
    BalanceObjectives objectives;
    objectives.emplace_back(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(0.0, 0.0, -1.0)
    );

    PlaneCarvingCandidates candidates;
    candidates.push_back(candidate(
        0,
        Eigen::Vector3d(2.0, -10.0, 0.0),
        1.0,
        2.0
    ));
    candidates.push_back(candidate(
        1,
        Eigen::Vector3d(1.0, 10.0, 0.0),
        1.0,
        1.0
    ));

    const PlaneCarvingPrefixResult result =
        selectBestCarvingPrefix(
            MassProperties(
                3.0,
                Eigen::Vector3d(3.0, 0.0, 0.0)
            ),
            objectives,
            candidates
        );

    if (
        result.evaluatedPrefixCount != 3 ||
        result.bestPrefixLength != 2 ||
        !close(result.initialEnergy, 0.5) ||
        !close(result.bestEnergy, 0.0) ||
        !close(result.bestMassProperties.mass, 1.0) ||
        result.bestCenterOfMass.norm() > 1e-12
    )
    {
        return fail(
            scenario,
            "full scan stopped at a transient energy increase"
        );
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


int runMiddle()
{
    const std::string scenario = "middle-prefix";
    BalanceObjectives objectives;
    objectives.emplace_back(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(0.0, -1.0, 0.0)
    );

    PlaneCarvingCandidates candidates;
    candidates.push_back(candidate(0, Eigen::Vector3d(5.0, 0.0, 0.0), 1.0, 5.0));
    candidates.push_back(candidate(1, Eigen::Vector3d(4.0, 0.0, 0.0), 1.0, 4.0));
    candidates.push_back(candidate(2, Eigen::Vector3d(0.1, 0.0, 0.0), 1.0, 0.1));

    const PlaneCarvingPrefixResult result =
        selectBestCarvingPrefix(
            MassProperties(
                10.0,
                Eigen::Vector3d(10.0, 0.0, 0.0)
            ),
            objectives,
            candidates
        );

    if (
        result.evaluatedPrefixCount != 4 ||
        result.bestPrefixLength != 2 ||
        !close(result.bestEnergy, 0.0078125) ||
        !close(result.bestMassProperties.mass, 8.0) ||
        !close(result.bestMassProperties.firstMoment(0), 1.0) ||
        !close(result.bestCenterOfMass(0), 0.125)
    )
    {
        return fail(scenario, "global optimum was not the middle prefix");
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


int runNoImprovement()
{
    const std::string scenario = "no-improvement";
    BalanceObjectives objectives;
    objectives.emplace_back(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(0.0, -1.0, 0.0)
    );

    PlaneCarvingCandidates candidates;
    candidates.push_back(candidate(
        0,
        Eigen::Vector3d(0.5, 0.0, 0.0),
        1.0,
        0.5
    ));

    const MassProperties initial(
        3.0,
        Eigen::Vector3d(3.0, 0.0, 0.0)
    );
    const PlaneCarvingPrefixResult result =
        selectBestCarvingPrefix(
            initial,
            objectives,
            candidates
        );

    if (
        result.bestPrefixLength != 0 ||
        !close(result.bestEnergy, 0.5) ||
        !close(result.bestMassProperties.mass, initial.mass) ||
        (result.bestMassProperties.firstMoment - initial.firstMoment).norm() > 0.0
    )
    {
        return fail(scenario, "prefix zero did not remain eligible");
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


int runTiePrefix()
{
    const std::string scenario = "tie-prefix";
    BalanceObjectives objectives;
    objectives.emplace_back(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(0.0, -1.0, 0.0)
    );

    PlaneCarvingCandidates candidates;
    candidates.push_back(candidate(
        0,
        Eigen::Vector3d(5.0, 0.0, 0.0),
        1.0,
        1.0
    ));

    const PlaneCarvingPrefixResult result =
        selectBestCarvingPrefix(
            MassProperties(
                3.0,
                Eigen::Vector3d(3.0, 0.0, 0.0)
            ),
            objectives,
            candidates
        );

    if (
        result.evaluatedPrefixCount != 2 ||
        result.bestPrefixLength != 0 ||
        !close(result.bestEnergy, 0.5)
    )
    {
        return fail(scenario, "equal energy did not retain the shorter prefix");
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


int runScaleInvariantPrefix()
{
    const std::string scenario = "scale-prefix";
    const double scales[3] = {1e-8, 1.0, 1e8};

    for (double scale : scales)
    {
        BalanceObjectives objectives;
        objectives.emplace_back(
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d(0.0, -1.0, 0.0)
        );

        PlaneCarvingCandidates candidates;
        candidates.push_back(candidate(
            0,
            scale * Eigen::Vector3d(2.0, 0.0, 0.0),
            1.0,
            1.0
        ));

        const PlaneCarvingPrefixResult result =
            selectBestCarvingPrefix(
                MassProperties(
                    3.0,
                    scale * Eigen::Vector3d(3.0, 0.0, 0.0)
                ),
                objectives,
                candidates
            );

        if (
            result.bestPrefixLength != 1 ||
            !close(result.bestEnergy, 0.125 * scale * scale)
        )
        {
            return fail(
                scenario,
                "best-prefix decision changed under uniform coordinate scaling"
            );
        }
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


int runOverThousand()
{
    const std::string scenario = "over-thousand";
    BalanceObjectives objectives;
    objectives.emplace_back(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(0.0, -1.0, 0.0)
    );

    PlaneCarvingCandidates candidates;
    candidates.reserve(1001);

    for (int index = 0; index < 1001; ++index)
    {
        candidates.push_back(candidate(
            index,
            Eigen::Vector3d(2.0, 0.0, 0.0),
            1.0,
            2.0
        ));
    }

    const PlaneCarvingPrefixResult result =
        selectBestCarvingPrefix(
            MassProperties(
                2000.0,
                Eigen::Vector3d(2002.0, 0.0, 0.0)
            ),
            objectives,
            candidates
        );

    if (
        result.evaluatedPrefixCount != 1002 ||
        result.bestPrefixLength != 1001 ||
        !close(result.bestEnergy, 0.0) ||
        !close(result.bestMassProperties.mass, 999.0)
    )
    {
        return fail(scenario, "legacy 1000-candidate cap remains active");
    }

    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}

} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        return fail("unknown", "usage: plane-carving-policy-regression <scenario>");
    }

    const std::string scenario = argv[1];

    if (scenario == "planes")
    {
        return runPlanes();
    }
    if (scenario == "multi-score")
    {
        return runMultiScore();
    }
    if (scenario == "ordering")
    {
        return runOrdering();
    }
    if (scenario == "transient-prefix")
    {
        return runTransient();
    }
    if (scenario == "middle-prefix")
    {
        return runMiddle();
    }
    if (scenario == "no-improvement")
    {
        return runNoImprovement();
    }
    if (scenario == "tie-prefix")
    {
        return runTiePrefix();
    }
    if (scenario == "scale-prefix")
    {
        return runScaleInvariantPrefix();
    }
    if (scenario == "over-thousand")
    {
        return runOverThousand();
    }

    return fail(scenario, "unknown scenario");
}
