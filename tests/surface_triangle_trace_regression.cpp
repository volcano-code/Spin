#include "spin2014/SurfaceTriangleMoments.h"
#include "spin2014/SurfaceRawMomentDifferential.h"
#include <iostream>
#include <stdexcept>
#include "surface_trace_independent_oracle.h"

int main() {
    using namespace spin2014::detail;
    const Eigen::Vector3d a(0,0,0), b(1,0,0), c(0,2,3);
    TriangleMomentJacobian plain, traced;
    const auto baseline = surfaceTriangleMoments(a,b,c,&plain);
    TriangleMomentTrace trace;
    const auto observed = traceSurfaceTriangleMoments(a,b,c,&traced,trace);
    if (!(baseline.array()==observed.array()).all() || !(plain.array()==traced.array()).all())
        throw std::runtime_error("trace changes scientific output");
    const double expected[8][3] = {{1,2,3},{1,0,0},{1,4,9},{1,8,27},
                                 {1,4,9},{3,4,9},{1,12,27},{2,36,9}};
    for (int i=0;i<8;++i) for (int j=0;j<3;++j)
        if(trace.h[i](j)!=expected[i][j]) throw std::runtime_error("supplement h mismatch");
    if (trace.normal != Eigen::Vector3d(0,-3,2)) throw std::runtime_error("normal mismatch");
    if (!trace.derivativesRecorded) throw std::runtime_error("missing derivatives");
    const double d8[3][3]={{2,12,3},{4,12,6},{4,36,6}};
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        if(trace.d8[i](j)!=d8[i][j]) throw std::runtime_error("d8 mismatch");
    TriangleMomentValues local;
    local << 0,0,-12,18,0,-108,18,0,-24,54;
    if(!(trace.unscaledPaperOrder.array()==local.array()).all())
        throw std::runtime_error("paper/native moment order mismatch");
    TriangleMomentTrace valuesOnly;
    const auto second = traceSurfaceTriangleMoments(a,b,c,nullptr,valuesOnly);
    if(valuesOnly.derivativesRecorded || !(second.array()==baseline.array()).all())
        throw std::runtime_error("value-only trace contract");
    try {
        traceSurfaceTriangleMoments(a,a,c,&traced,trace);
        throw std::runtime_error("degenerate triangle accepted");
    } catch(const std::invalid_argument&) {}
    if(trace.completed || trace.derivativesRecorded)
        throw std::runtime_error("failed call retained success flags");
    Eigen::MatrixXd vertices(4,3);
    vertices << 0,0,0, 1,0,0, 0,2,0, 0,0,3;
    Eigen::MatrixXi faces(4,3);
    faces << 0,2,1, 0,1,3, 0,3,2, 1,2,3;
    const auto reference=spin2014::differentiateSurfaceRawMoments(vertices,faces,2.5);
    int count=0;
    TriangleMomentValues accumulated=TriangleMomentValues::Zero();
    Eigen::MatrixXd lastJacobian;
    const int rows[10]={0,1,2,3,9,7,8,4,5,6};
    const double factors[10]={1./6,1./24,1./24,1./24,1./120,1./120,1./120,1./60,1./60,1./60};
    const auto measured=spin2014::differentiateSurfaceRawMomentsObserved(vertices,faces,2.5,
        [&](Eigen::Index face,const TriangleMomentTrace& local,
            const TriangleMomentValues& sum,const Eigen::MatrixXd& jacobian) {
            if(face!=count++ || !local.completed || !local.derivativesRecorded)
                throw std::runtime_error("incomplete or unordered face trace");
            for(int row=0;row<10;++row)
                accumulated(rows[row])+=2.5*(factors[row]*local.unscaledPaperOrder(row));
            if(!(sum.array()==accumulated.array()).all())
                throw std::runtime_error("face accumulation mismatch");
            lastJacobian=jacobian;
        });
    if(count!=4 || !(reference.moments.values().array()==measured.moments.values().array()).all()
        || !(reference.jacobian.array()==measured.jacobian.array()).all()
        || !(lastJacobian.array()==measured.jacobian.array()).all())
        throw std::runtime_error("whole surface trace changed result");
    TriangleMomentValues analytic;
    analytic << 2.5,.625,1.25,1.875,.25,1,2.25,.75,.375,.25;
    if((measured.moments.values()-analytic).cwiseAbs().maxCoeff()>1e-14)
        throw std::runtime_error("analytic tetrahedron density/moments mismatch");
    bool emptyRejected=false, interrupted=false;
    try { spin2014::differentiateSurfaceRawMomentsObserved(vertices,faces,2.5,{}); }
    catch(const std::invalid_argument&) { emptyRejected=true; }
    count=0;
    try {
        spin2014::differentiateSurfaceRawMomentsObserved(vertices,faces,2.5,
            [&](Eigen::Index,const TriangleMomentTrace&,const TriangleMomentValues&,const Eigen::MatrixXd&) {
                ++count; throw std::runtime_error("observer stopped");
            });
    } catch(const std::runtime_error& error) { interrupted=std::string(error.what())=="observer stopped"; }
    if(!emptyRejected || !interrupted || count!=1)
        throw std::runtime_error("observer interruption not propagated");
    surface_trace_oracle::verifyCubature();
    surface_trace_oracle::verifyPrefixes(vertices,faces,2.5);
    Eigen::Matrix3d transform;
    transform << 1,.2,.1,-.4,1,.3,.2,-.1,1;
    Eigen::MatrixXd shifted=vertices*transform.transpose();
    shifted.rowwise()+=Eigen::RowVector3d(.75,-.5,1.25);
    surface_trace_oracle::verifyPrefixes(shifted,faces,.7);
    Eigen::MatrixXi reversed=faces;
    reversed.col(1).swap(reversed.col(2));
    surface_trace_oracle::verifyPrefixes(shifted,reversed,3.1);
    std::cout << "local and complete trace checks passed; independent cubature/complex-step checked 12 face prefixes\n";
}
