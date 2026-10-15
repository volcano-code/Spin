#pragma once
#include <Eigen/Core>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <limits>
namespace spin2014 { namespace detail {
// Per-QP cache, never shared across different inverse metrics. A rejected fast
// solve is NOT a rank decision: the caller must use its original COD path.
// Bounds below are conservative floating-point guards, not interval certificates.
class SchurFactorCache {
    Eigen::MatrixXd matrix_,lower_,inverse_;
    std::vector<int> order_;
    double factorError_=0.,inverseError_=0.;
public:
    void clear() {matrix_.resize(0,0);lower_.resize(0,0);inverse_.resize(0,0);order_.clear();factorError_=inverseError_=0.;}
    bool trySolve(const Eigen::MatrixXd& s,const Eigen::VectorXd& rhs,
                  const std::vector<int>& keys,Eigen::VectorXd& solution) {
        const Eigen::Index n=s.rows();
        auto reject=[&](){clear();solution.resize(0);return false;};
        if(n<=0 || s.cols()!=n || rhs.size()!=n || keys.size()!=std::size_t(n) ||
           !s.allFinite() || !rhs.allFinite())return reject();
        // Symmetry is a matrix property, independent of the current RHS.
        // Permit only roundoff-sized differences in independently formed H*C'.
        if((s-s.transpose()).cwiseAbs().maxCoeff() >
           128*std::numeric_limits<double>::epsilon()*std::max(1.,s.cwiseAbs().maxCoeff()))return reject();
        std::map<int,Eigen::Index> positions;
        for(Eigen::Index i=0;i<n;++i)if(!positions.emplace(keys[i],i).second)return reject();
        bool reuse=order_.size()<=keys.size();
        std::vector<Eigen::Index> indices;indices.reserve(n);
        for(int k:order_)if(!positions.count(k))reuse=false;
        if(reuse)for(int k:order_)indices.push_back(positions.at(k));
        if(reuse)for(Eigen::Index j=0;j<matrix_.cols();++j)
            for(Eigen::Index i=0;i<matrix_.rows();++i)
                if(matrix_(i,j)!=s(indices[i],indices[j]))reuse=false;
        if(!reuse)clear(); // Release or changed metric: rebuild, never unstable downdate.
        const Eigen::Index old=order_.size();
        std::map<int,bool> present;for(int k:order_)present[k]=true;
        for(int k:keys)if(!present.count(k))order_.push_back(k);
        indices.clear();for(int k:order_)indices.push_back(positions.at(k));
        matrix_.resize(n,n);
        for(Eigen::Index j=0;j<n;++j)for(Eigen::Index i=0;i<n;++i)
            matrix_(i,j)=s(indices[i],indices[j]);
        lower_.conservativeResize(n,n);inverse_.conservativeResize(n,n);
        if(old<n){lower_.bottomRows(n-old).setZero();lower_.rightCols(n-old).setZero();
            inverse_.bottomRows(n-old).setZero();inverse_.rightCols(n-old).setZero();}
        const double eps=std::numeric_limits<double>::epsilon();
        for(Eigen::Index i=old;i<n;++i) {
            for(Eigen::Index j=0;j<i;++j) {
                const double dot=lower_.row(i).head(j).dot(lower_.row(j).head(j));
                const double absdot=lower_.row(i).head(j).cwiseAbs().dot(lower_.row(j).head(j).cwiseAbs());
                const double v=(matrix_(i,j)-dot)/lower_(j,j);
                lower_(i,j)=v;
                const double round=64*eps*(j+2)*(std::abs(matrix_(i,j))+absdot+std::abs(v*lower_(j,j)));
                const double err=std::abs(matrix_(i,j)-(dot+v*lower_(j,j)))+round;
                factorError_+=2*err+std::abs(matrix_(i,j)-matrix_(j,i));
            }
            const double sum=lower_.row(i).head(i).squaredNorm();
            const double diag=matrix_(i,i)-sum;
            if(!std::isfinite(diag) || diag<=0)return reject();
            lower_(i,i)=std::sqrt(diag);
            factorError_+=std::abs(matrix_(i,i)-(sum+lower_(i,i)*lower_(i,i)))+
                64*eps*(i+2)*(std::abs(matrix_(i,i))+sum+diag);
            if(i) {
                const Eigen::RowVectorXd product=lower_.row(i).head(i)*inverse_.topLeftCorner(i,i);
                const Eigen::RowVectorXd absolute=lower_.row(i).head(i).cwiseAbs()*inverse_.topLeftCorner(i,i).cwiseAbs();
                inverse_.row(i).head(i)=-product/lower_(i,i);
                for(Eigen::Index j=0;j<i;++j)
                    inverseError_+=std::abs(product[j]+lower_(i,i)*inverse_(i,j))+
                        64*eps*(i+2)*(absolute[j]+std::abs(lower_(i,i)*inverse_(i,j)));
            }
            inverse_(i,i)=1./lower_(i,i);
            inverseError_+=std::abs(lower_(i,i)*inverse_(i,i)-1.)+128*eps;
        }
        // If ||I-LB||<1, ||L^-1||_2 <= ||B||_F/(1-error). Subtract
        // factor reconstruction/asymmetry error, and require a wide separation
        // from the existing COD 1e-12 rank threshold (not a new science tolerance).
        const double invSquared=inverse_.squaredNorm(),norm=matrix_.norm();
        const double lowerBound=(1.-inverseError_)*(1.-inverseError_)/invSquared-factorError_;
        if(!lower_.allFinite() || !inverse_.allFinite() || !std::isfinite(lowerBound) ||
           !std::isfinite(norm) || inverseError_>=1e-6 || lowerBound<=1e-8*std::max(1.,norm))return reject();
        Eigen::VectorXd ordered(n);for(Eigen::Index i=0;i<n;++i)ordered[i]=rhs[indices[i]];
        Eigen::VectorXd y=lower_.triangularView<Eigen::Lower>().solve(ordered);
        const Eigen::VectorXd answer=lower_.transpose().triangularView<Eigen::Upper>().solve(y);
        solution.resize(n);for(Eigen::Index i=0;i<n;++i)solution[indices[i]]=answer[i];
        if(!solution.allFinite() || (s*solution-rhs).cwiseAbs().maxCoeff()>
            1e-12*std::max(1.,rhs.cwiseAbs().maxCoeff()))return reject();
        return true;
    }
};
}}
