#include "spin2014/CageSpinDeformationSolver.h"

#include "manufacturing/ManufacturingGeometry.h"
#include "manufacturing/SlicerValidation.h"

#include <Eigen/Core>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


struct Evaluation
{
    spin2014::CageSpinDeformationEvaluation physical;
    Eigen::VectorXd gradient;
    Eigen::MatrixXd jacobian;
};


Evaluation evaluate(
    const spin2014::CageSpinDeformationFunctional& functional,
    const Eigen::VectorXd& state,
    Eigen::Index fillCount
)
{
    Evaluation result;
    result.physical = functional.evaluate(
        state.head(fillCount), state.tail(state.size() - fillCount)
    );
    result.gradient.resize(state.size());
    result.gradient <<
        result.physical.fillGradient,
        result.physical.handleGradient;
    result.jacobian.resize(result.physical.constraints.size(), state.size());
    result.jacobian <<
        result.physical.fillConstraintJacobian,
        result.physical.handleConstraintJacobian;
    return result;
}


Eigen::VectorXd solve(const Eigen::MatrixXd& matrix, const Eigen::VectorXd& rhs)
{
    if (matrix.rows() == 0) return Eigen::VectorXd();
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition(matrix);
    decomposition.setThreshold(1e-12);
    return decomposition.solve(rhs);
}


struct History
{
    int capacity = 8;
    std::deque<Eigen::VectorXd> s;
    std::deque<Eigen::VectorXd> y;
    std::deque<double> rho;

    Eigen::VectorXd apply(const Eigen::VectorXd& vector) const
    {
        if (s.empty()) return vector;
        Eigen::VectorXd q = vector;
        std::vector<double> alpha(s.size());
        for (std::size_t offset = 0; offset < s.size(); ++offset)
        {
            const std::size_t index = s.size() - 1 - offset;
            alpha[index] = rho[index] * s[index].dot(q);
            q -= alpha[index] * y[index];
        }
        q *= std::clamp(
            s.back().dot(y.back()) / y.back().squaredNorm(), 1e-6, 1e6
        );
        for (std::size_t index = 0; index < s.size(); ++index)
        {
            const double coefficient = rho[index] * y[index].dot(q);
            q += s[index] * (alpha[index] - coefficient);
        }
        return q;
    }

    void update(const Eigen::VectorXd& step, const Eigen::VectorXd& change)
    {
        const double curvature = step.dot(change);
        if (curvature <= 1e-12 * step.norm() * change.norm()) return;
        if (static_cast<int>(s.size()) == capacity)
        {
            s.pop_front(); y.pop_front(); rho.pop_front();
        }
        s.push_back(step); y.push_back(change); rho.push_back(1.0 / curvature);
    }
};


Eigen::VectorXd multipliers(const Evaluation& value)
{
    return solve(
        value.jacobian * value.jacobian.transpose(),
        -value.jacobian * value.gradient
    );
}


Eigen::VectorXd lagrangian(
    const Evaluation& value,
    const Eigen::VectorXd& lambda
)
{
    return value.gradient + value.jacobian.transpose() * lambda;
}


bool validGeometry(
    const spin2014::CageSpinDeformationFunctional& functional,
    const Evaluation& value,
    std::string* diagnostic = nullptr
)
{
    const manufacturing::BoundaryAudit audit =
        manufacturing::auditTriangleBoundary(
            value.physical.deformedOuterVertices,
            functional.outerFaces(),
            manufacturing::BoundaryOrientation::OutwardPositive
        );
    if (!audit.strictBoundaryValid())
    {
        if (diagnostic) *diagnostic = "outer boundary audit failed";
        return false;
    }
    const std::vector<spin::TriangleMesh>& cells = functional.voidCells();
    if (cells.size() != value.physical.deformedVoidCellVertices.size())
    {
        if (diagnostic) *diagnostic = "void-cell count mismatch";
        return false;
    }
    for (std::size_t index = 0; index < cells.size(); ++index)
    {
        Eigen::MatrixXi inward = cells[index].faces;
        inward.col(1).swap(inward.col(2));
        const manufacturing::BoundaryAudit cellAudit =
            manufacturing::auditTriangleBoundary(
                value.physical.deformedVoidCellVertices[index],
                cells[index].faces,
                manufacturing::BoundaryOrientation::OutwardPositive
            );
        if (!cellAudit.strictBoundaryValid())
        {
            if (diagnostic) *diagnostic = "void-cell boundary audit failed";
            return false;
        }
        const manufacturing::MaterialShellAudit shell =
            manufacturing::auditMaterialShells(
                value.physical.deformedOuterVertices,
                functional.outerFaces(),
                value.physical.deformedVoidCellVertices[index],
                inward,
                functional.minimumShellThickness()
            );
        if (
            !shell.wallClearance.completed ||
            !shell.wallClearance.intersectionFree ||
            !shell.wallClearance.meetsRequestedThickness ||
            !shell.containment.completed ||
            !shell.containment.allCavitiesInsideOuter
        )
        {
            if (diagnostic)
            {
                if (!shell.wallClearance.completed)
                    *diagnostic = "shell clearance audit did not complete";
                else if (!shell.wallClearance.intersectionFree)
                    *diagnostic = "void cell intersects the outer boundary"
                        " (minimum_clearance=" +
                        std::to_string(
                            shell.wallClearance.measuredMinimumClearance
                        ) + ", outer_face=" +
                        std::to_string(shell.wallClearance.outerFace) +
                        ", inner_face=" +
                        std::to_string(shell.wallClearance.innerFace) + ')';
                else if (!shell.wallClearance.meetsRequestedThickness)
                    *diagnostic = "void cell violates shell thickness";
                else if (!shell.containment.completed)
                    *diagnostic = "shell containment audit did not complete";
                else
                    *diagnostic = "void cell lies outside the outer boundary";
            }
            return false;
        }
    }
    return true;
}


double merit(const Evaluation& value, double penalty)
{
    return value.physical.objective +
        penalty * value.physical.constraints.lpNorm<1>();
}


void appendDouble(std::string& bytes, double value)
{
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    bytes.append(reinterpret_cast<const char*>(&bits), sizeof(bits));
}


std::string traceHash(const std::vector<spin2014::CageSpinIteration>& trace)
{
    std::string bytes("spin-it-2014-cage-trace-v1");
    for (const auto& entry : trace)
    {
        bytes.append(reinterpret_cast<const char*>(&entry.iteration), sizeof(int));
        appendDouble(bytes, entry.objective);
        appendDouble(bytes, entry.equalityInfinityNorm);
        appendDouble(bytes, entry.projectedGradientInfinityNorm);
        appendDouble(bytes, entry.acceptedStep);
        bytes.append(
            reinterpret_cast<const char*>(&entry.activeVariables), sizeof(int)
        );
        bytes.append(
            reinterpret_cast<const char*>(&entry.lbfgsHistorySize), sizeof(int)
        );
        bytes.push_back(entry.geometryAccepted ? '\1' : '\0');
    }
    return manufacturing::computeContentSha256(bytes);
}


} // namespace


namespace spin2014
{


const char* cageSpinStatusToken(CageSpinStatus status) noexcept
{
    switch (status)
    {
    case CageSpinStatus::Converged: return "CONVERGED";
    case CageSpinStatus::MaximumIterations: return "MAXIMUM_ITERATIONS";
    case CageSpinStatus::LineSearchExhausted: return "LINE_SEARCH_EXHAUSTED";
    case CageSpinStatus::InvalidGeometry: return "INVALID_GEOMETRY";
    }
    return "INVALID_GEOMETRY";
}


CageSpinDeformationResult solveCageSpinDeformation(
    const CageSpinDeformationFunctional& functional,
    const Eigen::VectorXd& initialBeta,
    const CageSpinDeformationSolverSettings& settings
)
{
    if (
        initialBeta.size() == 0 || !initialBeta.allFinite() ||
        (initialBeta.array() < 0.0).any() ||
        (initialBeta.array() > 1.0).any() ||
        settings.maximumIterations <= 0 || settings.maximumLineSearchSteps <= 0 ||
        settings.lbfgsMemory <= 0 || settings.equalityTolerance <= 0.0 ||
        settings.projectedGradientTolerance <= 0.0 ||
        settings.meritPenalty <= 0.0 || settings.lineSearchReduction <= 0.0 ||
        settings.lineSearchReduction >= 1.0 || settings.maximumTranslation <= 0.0 ||
        settings.minimumScale <= 0.0 ||
        settings.maximumScale <= settings.minimumScale
    )
    {
        throw std::invalid_argument("invalid cage solver request");
    }
    const Eigen::VectorXd initialHandles = functional.initialHandleParameters();
    if (initialHandles.size() % 4 != 0)
    {
        throw std::invalid_argument("cage handles must use four parameters each");
    }
    const Eigen::Index fillCount = initialBeta.size();
    Eigen::VectorXd state(initialBeta.size() + initialHandles.size());
    state << initialBeta, initialHandles;
    Eigen::VectorXd lower = Eigen::VectorXd::Zero(state.size());
    Eigen::VectorXd upper = Eigen::VectorXd::Ones(state.size());
    for (Eigen::Index index = 0; index < initialHandles.size(); ++index)
    {
        const Eigen::Index global = fillCount + index;
        if (index % 4 == 3)
        {
            lower(global) = settings.minimumScale;
            upper(global) = settings.maximumScale;
        }
        else
        {
            lower(global) = initialHandles(index) - settings.maximumTranslation;
            upper(global) = initialHandles(index) + settings.maximumTranslation;
        }
    }

    CageSpinDeformationResult result;
    Evaluation current = evaluate(functional, state, fillCount);
    result.evaluation = current.physical;
    std::string initialGeometryDiagnostic;
    result.geometryValid = validGeometry(
        functional, current, &initialGeometryDiagnostic
    );
    if (!result.geometryValid)
    {
        result.status = CageSpinStatus::InvalidGeometry;
        result.diagnostic = "initial cage geometry failed: " +
            initialGeometryDiagnostic;
        return result;
    }
    History history;
    history.capacity = settings.lbfgsMemory;

    for (int iteration = 0; iteration < settings.maximumIterations; ++iteration)
    {
        const Eigen::VectorXd lambda = multipliers(current);
        Eigen::VectorXd stationarity = lagrangian(current, lambda);
        std::vector<unsigned char> active(static_cast<std::size_t>(state.size()), 0U);
        Eigen::MatrixXd freeA = current.jacobian;
        CageSpinIteration record;
        record.iteration = iteration;
        record.objective = current.physical.objective;
        record.equalityInfinityNorm =
            current.physical.constraints.lpNorm<Eigen::Infinity>();
        record.lbfgsHistorySize = static_cast<int>(history.s.size());
        record.geometryAccepted = true;
        for (Eigen::Index index = 0; index < state.size(); ++index)
        {
            const bool fixed =
                (state(index) <= lower(index) + 1e-10 && stationarity(index) > 0.0) ||
                (state(index) >= upper(index) - 1e-10 && stationarity(index) < 0.0);
            if (fixed)
            {
                active[static_cast<std::size_t>(index)] = 1U;
                freeA.col(index).setZero();
                stationarity(index) = 0.0;
                ++record.activeVariables;
            }
        }
        const Eigen::VectorXd freeLambda = solve(
            freeA * freeA.transpose(), -freeA * current.gradient
        );
        Eigen::VectorXd projected =
            current.gradient + current.jacobian.transpose() * freeLambda;
        for (Eigen::Index index = 0; index < state.size(); ++index)
            if (active[static_cast<std::size_t>(index)]) projected(index) = 0.0;
        record.projectedGradientInfinityNorm =
            projected.lpNorm<Eigen::Infinity>();
        if (
            record.equalityInfinityNorm <= settings.equalityTolerance &&
            record.projectedGradientInfinityNorm <=
                settings.projectedGradientTolerance
        )
        {
            result.trace.push_back(record);
            result.status = CageSpinStatus::Converged;
            result.diagnostic = "joint cage SLQP tolerances met";
            break;
        }

        const auto applyFree = [&history, &active](const Eigen::VectorXd& input)
        {
            Eigen::VectorXd value = input;
            for (Eigen::Index index = 0; index < value.size(); ++index)
                if (active[static_cast<std::size_t>(index)]) value(index) = 0.0;
            value = history.apply(value);
            for (Eigen::Index index = 0; index < value.size(); ++index)
                if (active[static_cast<std::size_t>(index)]) value(index) = 0.0;
            return value;
        };
        const Eigen::VectorXd hg = applyFree(current.gradient);
        Eigen::MatrixXd hat(state.size(), current.jacobian.rows());
        for (Eigen::Index row = 0; row < current.jacobian.rows(); ++row)
            hat.col(row) = applyFree(current.jacobian.row(row).transpose());
        const Eigen::VectorXd qpLambda = solve(
            current.jacobian * hat,
            current.physical.constraints - current.jacobian * hg
        );
        Eigen::VectorXd direction = -hg - hat * qpLambda;
        const auto removeOutwardComponents = [&state, &lower, &upper](
            Eigen::VectorXd* candidateDirection
        )
        {
            for (Eigen::Index index = 0; index < state.size(); ++index)
            {
                if (
                    (state(index) <= lower(index) + 1e-12 &&
                     (*candidateDirection)(index) < 0.0) ||
                    (state(index) >= upper(index) - 1e-12 &&
                     (*candidateDirection)(index) > 0.0)
                )
                {
                    (*candidateDirection)(index) = 0.0;
                }
            }
        };
        removeOutwardComponents(&direction);
        const auto boundStep = [&state, &lower, &upper](
            const Eigen::VectorXd& candidateDirection
        )
        {
            double maximumStep = 1.0;
            for (Eigen::Index index = 0; index < state.size(); ++index)
            {
                if (candidateDirection(index) > 0.0)
                    maximumStep = std::min(maximumStep,
                        (upper(index) - state(index)) /
                        candidateDirection(index));
                else if (candidateDirection(index) < 0.0)
                    maximumStep = std::min(maximumStep,
                        (lower(index) - state(index)) /
                        candidateDirection(index));
            }
            return std::clamp(maximumStep, 0.0, 1.0);
        };
        double step = boundStep(direction);
        const double oldMerit = merit(current, settings.meritPenalty);
        bool accepted = false;
        Evaluation candidate;
        Eigen::VectorXd candidateState;
        for (int search = 0; search < settings.maximumLineSearchSteps; ++search)
        {
            candidateState = (state + step * direction).cwiseMax(lower).cwiseMin(upper);
            candidate = evaluate(functional, candidateState, fillCount);
            if (
                validGeometry(functional, candidate) &&
                merit(candidate, settings.meritPenalty) < oldMerit - 1e-12
            )
            {
                accepted = true;
                break;
            }
            step *= settings.lineSearchReduction;
        }
        if (!accepted)
        {
            Eigen::VectorXd freeGradient = current.gradient;
            for (Eigen::Index index = 0; index < state.size(); ++index)
                if (active[static_cast<std::size_t>(index)])
                    freeGradient(index) = 0.0;
            const Eigen::VectorXd identityLambda = solve(
                freeA * freeA.transpose(),
                current.physical.constraints - freeA * freeGradient
            );
            direction = -freeGradient - freeA.transpose() * identityLambda;
            removeOutwardComponents(&direction);
            step = boundStep(direction);
            for (int search = 0;
                 search < settings.maximumLineSearchSteps;
                 ++search)
            {
                candidateState = (state + step * direction)
                    .cwiseMax(lower).cwiseMin(upper);
                candidate = evaluate(functional, candidateState, fillCount);
                if (
                    validGeometry(functional, candidate) &&
                    merit(candidate, settings.meritPenalty) < oldMerit - 1e-12
                )
                {
                    accepted = true;
                    break;
                }
                step *= settings.lineSearchReduction;
            }
        }
        if (!accepted)
        {
            Eigen::VectorXd constraintSign =
                current.physical.constraints.unaryExpr([](double value)
                {
                    return value > 0.0 ? 1.0 : (value < 0.0 ? -1.0 : 0.0);
                });
            direction = -(
                current.gradient + settings.meritPenalty *
                current.jacobian.transpose() * constraintSign
            );
            for (Eigen::Index index = 0; index < direction.size(); ++index)
                if (active[static_cast<std::size_t>(index)]) direction(index) = 0.0;
            removeOutwardComponents(&direction);
            step = boundStep(direction);
            for (int search = 0;
                 search < settings.maximumLineSearchSteps;
                 ++search)
            {
                candidateState = (state + step * direction)
                    .cwiseMax(lower).cwiseMin(upper);
                candidate = evaluate(functional, candidateState, fillCount);
                if (
                    validGeometry(functional, candidate) &&
                    merit(candidate, settings.meritPenalty) < oldMerit - 1e-12
                )
                {
                    accepted = true;
                    break;
                }
                step *= settings.lineSearchReduction;
            }
        }
        record.acceptedStep = accepted ? step : 0.0;
        record.geometryAccepted = accepted;
        result.trace.push_back(record);
        if (!accepted)
        {
            result.status = CageSpinStatus::LineSearchExhausted;
            result.diagnostic = "joint cage line search exhausted";
            break;
        }
        const Eigen::VectorXd oldLagrangian = lagrangian(current, lambda);
        const Eigen::VectorXd nextLambda = multipliers(candidate);
        history.update(
            candidateState - state,
            lagrangian(candidate, nextLambda) - oldLagrangian
        );
        state = candidateState;
        current = candidate;
    }
    if (result.trace.size() == static_cast<std::size_t>(settings.maximumIterations))
    {
        result.status = CageSpinStatus::MaximumIterations;
        result.diagnostic = "joint cage solver reached its iteration limit";
    }
    result.beta = state.head(fillCount);
    result.handleParameters = state.tail(state.size() - fillCount);
    result.evaluation = current.physical;
    result.geometryValid = validGeometry(functional, current);
    result.silhouetteChanged =
        (result.handleParameters - initialHandles).lpNorm<Eigen::Infinity>() > 1e-12;
    result.traceHash = traceHash(result.trace);
    std::string bytes = result.traceHash + result.diagnostic;
    for (Eigen::Index index = 0; index < state.size(); ++index)
        appendDouble(bytes, state(index));
    result.outputHash = manufacturing::computeContentSha256(bytes);
    return result;
}


} // namespace spin2014
