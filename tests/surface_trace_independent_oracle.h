#pragma once
#include "spin2014/SurfaceRawMomentDifferential.h"
#include "spin2014/SurfaceTriangleMoments.h"
#include <array>
#include <complex>
#include <cmath>
#include <stdexcept>

namespace surface_trace_oracle {
using C=std::complex<double>;
using Point=std::array<C,3>;
using Triangle=std::array<Point,3>;
// Four-point degree-three rule on the reference triangle, area 1/2.
constexpr double bary[4][3]={{1./3,1./3,1./3},{.6,.2,.2},{.2,.6,.2},{.2,.2,.6}};
constexpr double weight[4]={-27./96,25./96,25./96,25./96};
inline void verifyCubature() {
    const double factorial[6]={1,1,2,6,24,120};
    for(int i=0;i<=3;++i) for(int j=0;j<=3-i;++j) {
        double value=0;
        for(int q=0;q<4;++q) value+=weight[q]*std::pow(bary[q][1],i)*std::pow(bary[q][2],j);
        const double exact=factorial[i]*factorial[j]/factorial[i+j+2];
        if(std::abs(value-exact)>1e-15) throw std::runtime_error("independent cubature calibration failed");
    }
}
// Divergence-field integration, no h1-h8 or native derivative formulas.
inline std::array<C,10> integrate(const Triangle& t) {
    Point u,v,n;
    for(int k=0;k<3;++k) { u[k]=t[1][k]-t[0][k]; v[k]=t[2][k]-t[0][k]; }
    for(int k=0;k<3;++k) n[k]=u[(k+1)%3]*v[(k+2)%3]-u[(k+2)%3]*v[(k+1)%3];
    std::array<C,10> out{};
    for(int q=0;q<4;++q) {
        Point p{};
        for(int k=0;k<3;++k) for(int a=0;a<3;++a) p[k]+=bary[q][a]*t[a][k];
        const C x=p[0],y=p[1],z=p[2];
        const std::array<C,10> field={n[0]*x,n[0]*x*x/2.,n[1]*y*y/2.,n[2]*z*z/2.,
            n[0]*x*x*x/3.,n[1]*y*y*y/3.,n[2]*z*z*z/3.,
            n[1]*y*y*z/2.,n[2]*z*z*x/2.,n[0]*x*x*y/2.};
        for(int row=0;row<10;++row) out[row]+=weight[q]*field[row];
    }
    return out;
}
struct PrefixErrors { double moment=0, derivative=0; int faces=0; };
inline PrefixErrors verifyPrefixes(const Eigen::MatrixXd& vertices,const Eigen::MatrixXi& faces,double density,double length=1.) {
    if(!std::isfinite(length) || length<=0 || !std::isfinite(density) || density<=0)
        throw std::runtime_error("invalid fixed normalization");
    PrefixErrors errors;
    Eigen::Matrix<double,10,1> expected=Eigen::Matrix<double,10,1>::Zero();
    Eigen::MatrixXd derivative=Eigen::MatrixXd::Zero(10,3*vertices.rows());
    int count=0;
    spin2014::differentiateSurfaceRawMomentsObserved(vertices,faces,density,
        [&](Eigen::Index face,const spin2014::detail::TriangleMomentTrace& trace,
            const physics::RawMassMoments10::Vector& moments,const Eigen::MatrixXd& jacobian) {
            if(face!=count++ || !trace.completed) throw std::runtime_error("missing prefix");
            Triangle triangle;
            for(int a=0;a<3;++a) for(int k=0;k<3;++k) triangle[a][k]=vertices(faces(face,a),k);
            const auto value=integrate(triangle);
            for(int row=0;row<10;++row) expected[row]+=density*value[row].real();
            for(int a=0;a<3;++a) for(int k=0;k<3;++k) {
                auto perturbed=triangle;
                perturbed[a][k]+=C(0,1e-30);
                const auto shifted=integrate(perturbed);
                for(int row=0;row<10;++row)
                    derivative(row,3*faces(face,a)+k)+=density*shifted[row].imag()/1e-30;
            }
            for(int row=0;row<10;++row) {
                const int degree=row==0?0:(row<4?1:2);
                errors.moment=std::max(errors.moment,std::abs(moments[row]-expected[row])/(density*std::pow(length,3+degree)));
                errors.derivative=std::max(errors.derivative,(jacobian.row(row)-derivative.row(row)).cwiseAbs().maxCoeff()/(density*std::pow(length,2+degree)));
            }
            if(!moments.allFinite() || !jacobian.allFinite() || !expected.allFinite() || !derivative.allFinite() ||
                errors.moment>1e-10 || errors.derivative>1e-10)
                throw std::runtime_error("independent face-prefix moment/derivative mismatch");
        });
    if(count!=faces.rows()) throw std::runtime_error("truncated surface trace");
    errors.faces=count;
    return errors;
}
}
