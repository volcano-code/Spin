#pragma once
#include <Eigen/Core>
#include <cmath>
#include <stdexcept>
#include <string>

namespace surface_test {
struct Mesh { Eigen::MatrixXd v; Eigen::MatrixXi f; };
inline void require(bool b, const char* message) { if (!b) throw std::runtime_error(message); }
template<class F> void rejects(F f, const char* message) {
    bool rejected=false; try { f(); } catch(const std::exception&) { rejected=true; }
    require(rejected,message);
}
inline Mesh box(Eigen::Vector3d lo, Eigen::Vector3d hi) {
    Mesh m; m.v.resize(8,3); m.f.resize(12,3);
    for(int i=0;i<8;++i) for(int a=0;a<3;++a) m.v(i,a)=(i&(1<<a))?hi(a):lo(a);
    m.f << 0,2,1,1,2,3,4,5,6,5,7,6,0,1,4,1,5,4,
           2,6,3,3,6,7,0,4,2,2,4,6,1,3,5,3,7,5;
    return m;
}
inline Mesh cube(double lo,double hi) { return box(Eigen::Vector3d::Constant(lo),Eigen::Vector3d::Constant(hi)); }
inline Mesh join(const Mesh& a, const Mesh& b, bool reverse=false) {
    Mesh m; m.v.resize(a.v.rows()+b.v.rows(),3); m.v<<a.v,b.v;
    m.f.resize(a.f.rows()+b.f.rows(),3); m.f.topRows(a.f.rows())=a.f;
    Eigen::MatrixXi f=b.f; if(reverse) f=b.f.rowwise().reverse();
    m.f.bottomRows(b.f.rows())=(f.array()+a.v.rows()).matrix(); return m;
}
using Ten=Eigen::Matrix<double,10,1>;
inline Ten boxMoments(Eigen::Vector3d lo, Eigen::Vector3d hi) {
    const double volume=(hi-lo).prod(); const Eigen::Vector3d mean=(lo+hi)/2;
    Ten r; r<<volume,volume*mean,volume*(lo.array().square()+lo.array()*hi.array()+hi.array().square()).matrix()/3,
        volume*mean.y()*mean.z(),volume*mean.x()*mean.z(),volume*mean.x()*mean.y(); return r;
}
inline void close(const Ten& got,const Ten& expected,const char* label) {
    require((got-expected).cwiseAbs().maxCoeff()/std::max(1e-30,expected.cwiseAbs().maxCoeff())<=1e-10,label);
}
// Uniform simplex barycentric expectations: E(lambda_i)=1/4,
// E(lambda_i lambda_j)=(1+delta_ij)/20. Independent of surface integration.
inline Ten tetraMoments(const Eigen::MatrixXd& v) {
    Eigen::Matrix3d edges; for(int i=0;i<3;++i) edges.col(i)=(v.row(i+1)-v.row(0)).transpose();
    const double volume=std::abs(edges.determinant())/6;
    const Eigen::Vector3d sum=v.colwise().sum().transpose();
    const Eigen::Matrix3d q=volume*(sum*sum.transpose()+v.transpose()*v)/20;
    Ten r; r<<volume,volume*sum/4,q(0,0),q(1,1),q(2,2),q(1,2),q(0,2),q(0,1); return r;
}
}
