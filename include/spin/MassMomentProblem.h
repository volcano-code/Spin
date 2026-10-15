#pragma once

#include "spin/SpinMomentProblem.h"

#include <functional>
#include <stdexcept>
#include <utility>

namespace spin
{

/** Owned objective/constraint contract shared by moment-based solvers.
 * Moments already include material density and are in the target frame.
 * No geometry, initialization, acceptance or material validation is hidden
 * here. The historical problem is copied, not captured by dangling reference.
 */
class MassMomentProblem
{
public:
    using Vector10 = SpinMomentProblem::Vector10;
    using Evaluator = std::function<SpinMomentEvaluation(const Vector10&)>;

    MassMomentProblem(SpinMomentProblem problem)
        : evaluator_([problem=std::move(problem)](const Vector10& moments) {
              return problem.evaluate(moments);
          })
    {}

    explicit MassMomentProblem(Evaluator evaluator)
        : evaluator_(std::move(evaluator))
    {
        if (!evaluator_)
            throw std::invalid_argument("mass moment problem evaluator is empty");
    }

    SpinMomentEvaluation evaluate(const Vector10& moments) const
    {
        return evaluator_(moments);
    }

private:
    Evaluator evaluator_;
};

} // namespace spin
