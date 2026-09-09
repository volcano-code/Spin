#include "spin2014/SlqpHistory.h"
#include <algorithm>
#include <cmath>
namespace spin2014 {
bool LbfgsHistory::update(const Eigen::VectorXd& step,const Eigen::VectorXd& gradientChange) {
    if(capacity<=0 || step.size()!=gradientChange.size() || !step.allFinite() || !gradientChange.allFinite())return false;
    const double curvature=step.dot(gradientChange);
    if(!std::isfinite(curvature)||curvature<=1e-12*step.norm()*gradientChange.norm())return false;
    const double reciprocal=1./curvature;
    if(!std::isfinite(reciprocal))return false;
    if(static_cast<int>(s.size())==capacity){s.pop_front();y.pop_front();inverseCurvature.pop_front();}
    s.push_back(step);y.push_back(gradientChange);inverseCurvature.push_back(reciprocal);return true;
}
Eigen::VectorXd LbfgsHistory::applyInverse(const Eigen::VectorXd& vector) const {
    if(s.empty())return vector;
    Eigen::VectorXd result=vector;std::vector<double> alpha(s.size(),0.);
    for(std::size_t offset=0;offset<s.size();++offset) {
        const std::size_t index=s.size()-1U-offset;
        alpha[index]=inverseCurvature[index]*s[index].dot(result);result-=alpha[index]*y[index];
    }
    const double yy=y.back().squaredNorm();const double scale=yy>0.?s.back().dot(y.back())/yy:1.;
    result*=std::clamp(scale,1e-6,1e6);
    for(std::size_t index=0;index<s.size();++index) {
        const double beta=inverseCurvature[index]*y[index].dot(result);result+=s[index]*(alpha[index]-beta);
    }return result;
}
bool LbfgsHistory::resetForMembership(const std::vector<unsigned char>& previous,const std::vector<unsigned char>& current) {
    if(previous==current)return false;
    s.clear();y.clear();inverseCurvature.clear();return true;
}
}
