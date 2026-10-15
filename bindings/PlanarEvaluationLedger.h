#pragma once

#include "spin/PlanarFidgetCircleProblem.h"
#include <pybind11/pybind11.h>
#include <cmath>
#include <exception>
#include <vector>

namespace spin_bindings
{

/** Binding-side assessments are outside the dogleg's scientific call budget.
 * The observer stores native data while the GIL is released; conversion and
 * exception attachment occur only after the binding has reacquired the GIL. */
struct PlanarEvaluationLedger
{
    std::vector<spin::PlanarDoglegEvaluationEvent> solver;
    std::vector<spin::PlanarDoglegEvaluationEvent> auxiliary;
    const char* phase = "request_geometry";

    spin::PlanarFidgetCircleEvaluation assess(
        const spin::PlanarFidgetCircleProblem& problem,
        const Eigen::Vector2d& center, const char* purpose)
    {
        phase = purpose;
        spin::PlanarDoglegEvaluationEvent event;
        event.evaluationIndex = auxiliary.size()+1;
        event.purpose = purpose;
        event.center = center;
        auxiliary.push_back(std::move(event));
        auto& recorded = auxiliary.back();
        try
        {
            auto value = problem.evaluate(center);
            recorded.normalizedResidual = value.normalizedResidual;
            recorded.hasResidual = true;
            recorded.outcome = "evaluated";
            return value;
        }
        catch (const std::exception& error)
        {
            recorded.outcome = "failed";
            recorded.diagnostic = error.what();
            throw;
        }
    }

    static pybind11::object scalar(double value)
    {
        if (std::isfinite(value)) return pybind11::float_(value);
        return pybind11::str(std::isnan(value) ? "nan" : value > 0 ? "+inf" : "-inf");
    }

    static pybind11::list records(
        const std::vector<spin::PlanarDoglegEvaluationEvent>& events)
    {
        namespace py = pybind11;
        py::list rows;
        for (const auto& event : events)
        {
            py::dict row;
            row["evaluation_index"] = event.evaluationIndex;
            row["iteration"] = event.iteration;
            row["purpose"] = event.purpose;
            row["difference_axis"] = event.differenceAxis;
            py::list center;
            center.append(scalar(event.center.x()));
            center.append(scalar(event.center.y()));
            row["center"] = std::move(center);
            row["has_residual"] = event.hasResidual;
            row["normalized_residual"] = event.hasResidual
                ? scalar(event.normalizedResidual.norm()) : py::none();
            row["normalized_residual_vector"] = py::none();
            if (event.hasResidual)
            {
                py::list residual;
                residual.append(scalar(event.normalizedResidual.x()));
                residual.append(scalar(event.normalizedResidual.y()));
                row["normalized_residual_vector"] = std::move(residual);
            }
            row["outcome"] = event.outcome;
            row["accepted_trial"] = event.acceptedTrial;
            row["diagnostic"] = event.diagnostic;
            rows.append(std::move(row));
        }
        return rows;
    }

    pybind11::dict toPython() const
    {
        pybind11::dict result;
        result["schema"] = "planar-dogleg-evaluation-ledger-v1";
        result["scope"] = "FUNCTION_EVALUATION_LEDGER_NOT_MATERIAL_OR_PAPER_APPROVAL";
        result["solver_evaluations"] = records(solver);
        result["auxiliary_evaluations"] = records(auxiliary);
        result["solver_function_evaluations"] = solver.size();
        result["auxiliary_geometry_evaluations"] = auxiliary.size();
        result["total_residual_evaluations"] = solver.size()+auxiliary.size();
        result["terminal_phase"] = phase;
        return result;
    }

    [[noreturn]] void rethrowWithLedger() const
    {
        // Use this pinned pybind11 version's ordinary native/Python exception
        // translation, retaining ValueError/OverflowError/etc. and messages.
        pybind11::detail::translate_exception(std::current_exception());
        pybind11::error_already_set error;
        error.value().attr("evaluation_ledger") = toPython();
        throw error;
    }
};

} // namespace spin_bindings
