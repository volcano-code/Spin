#pragma once
#include "spin2014/SlqpLbfgsSolver.h"
#include <ostream>
namespace spin2014 {
std::string slqpMembershipHash(const std::vector<unsigned char>&);
std::string hashAuditedSlqp(const SlqpLbfgsResult&,const SlqpLbfgsSettings&);
void captureSlqpState(SlqpIteration&,const AdaptiveFillEvaluation&,const Eigen::VectorXd&);
void writeSlqpAuditJson(std::ostream&,const SlqpIteration&);
void writeSlqpAuditControlsJson(std::ostream&,const SlqpLbfgsSettings&);
}
