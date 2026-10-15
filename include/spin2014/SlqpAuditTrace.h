#pragma once
#include "spin2014/SlqpLbfgsSolver.h"
#include <functional>
#include <ostream>
namespace spin2014 {
using SlqpIterationTraversal = std::function<void(
    const std::function<void(const SlqpIteration&)>&)>;
std::string slqpMembershipHash(const std::vector<unsigned char>&);
std::string hashAuditedSlqp(const SlqpLbfgsResult&,const SlqpLbfgsSettings&);
std::string hashAuditedSlqp(const SlqpLbfgsResult&,const SlqpLbfgsSettings&,
    const SlqpIterationTraversal&);
void captureSlqpState(SlqpIteration&,const AdaptiveFillEvaluation&,const Eigen::VectorXd&);
void writeSlqpAuditJson(std::ostream&,const SlqpIteration&);
void writeSlqpAuditControlsJson(std::ostream&,const SlqpLbfgsSettings&);
}
