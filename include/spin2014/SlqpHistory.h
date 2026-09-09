#pragma once
#include <Eigen/Core>
#include <deque>
#include <vector>

namespace spin2014 {
// Shared L-BFGS metric; the audited policy resets this actual history before
// applying a metric to a different active membership.
struct LbfgsHistory {
    int capacity=0;
    std::deque<Eigen::VectorXd> s,y;
    std::deque<double> inverseCurvature;
    bool update(const Eigen::VectorXd&,const Eigen::VectorXd&);
    Eigen::VectorXd applyInverse(const Eigen::VectorXd&) const;
    bool resetForMembership(const std::vector<unsigned char>& previous,
                            const std::vector<unsigned char>& current);
};
}
