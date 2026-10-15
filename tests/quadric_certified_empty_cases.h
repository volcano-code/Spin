#pragma once

#include "spin/QuadricMomentEvaluator.h"
#include "spin/QuadricHollowingSolver.h"
#include "../src/spin/QuadricMaterialExportGuard.h"
#include <Eigen/Geometry>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <CGAL/Gmpq.h>
#include <vector>

namespace certified_empty_tests {
using V = Eigen::Matrix<double,10,1>;

inline void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class F> inline void mustReject(F&& call, const char* message) {
    bool rejected=false;
    try { call(); } catch(const std::exception&) { rejected=true; }
    require(rejected,message);
}

inline spin::TriangleMesh box(double lo, double hi) {
    spin::TriangleMesh m;
    m.vertices.resize(8,3);
    m.vertices << lo,lo,lo,hi,lo,lo,hi,hi,lo,lo,hi,lo,
        lo,lo,hi,hi,lo,hi,hi,hi,hi,lo,hi,hi;
    m.faces.resize(12,3);
    m.faces << 0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,
        3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
    return m;
}

inline spin::SpinDesignRequest request() {
    spin::SpinDesignRequest r;
    r.designDomain=box(-1,1); r.outerSurface=box(-2,2);
    r.materialDensity=1.7; r.minimumWallThickness=.5;
    MassMoments fixed;
    fixed.mass=95.2; fixed.firstMoment.setZero();
    fixed.secondMoment=Eigen::Matrix3d::Identity()*(248./3.*1.7);
    r.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(fixed);
    return r;
}

inline spin::QuadricCoefficients10 farSphere() {
    V p=V::Zero(); p(0)=-99; p(1)=20; p.segment<3>(4).setConstant(-1);
    return spin::QuadricCoefficients10::fromUnnormalized(p);
}

// Construct an exact coplanar contact with one native tessellation face,
// while the whole small tetrahedron lies strictly inside the analytic sphere.
// Rational barycentric coordinates are retained only when exactly binary64
// representable, so this is not a tolerance-based near-contact fixture.
inline spin::SpinDesignRequest discreteFaceContact(
    const spin::QuadricCoefficients10& exterior,
    const spin::QuadricMomentEvaluatorSettings& settings)
{
    using Q=CGAL::Gmpq;
    auto r=request();
    const auto tri=spin::triangulateQuadric(exterior,r.designDomain.vertices,settings.triangulation);
    for(int f=0;f<tri.closedBoundary.faces.rows();++f) {
        std::vector<Eigen::Vector3d> points;
        for(int i=0;i<=16;++i) for(int j=0;j<=16;++j) {
            const int k=16-i-j; if(k<0 || k>16)continue;
            const int weights[3]={i,j,k};
            Eigen::Vector3d p; bool exact=true;
            for(int axis=0;axis<3;++axis) {
                Q value(0);
                for(int n=0;n<3;++n)
                    value+=Q(weights[n])*Q(tri.closedBoundary.vertices(tri.closedBoundary.faces(f,n),axis))/Q(16);
                p(axis)=CGAL::to_double(value);
                exact=exact && Q(p(axis))==value;
            }
            if(exact)points.push_back(p);
        }
        for(std::size_t a=0;a<points.size();++a)
            for(std::size_t b=a+1;b<points.size();++b)
                for(std::size_t c=b+1;c<points.size();++c) {
                    if((points[b]-points[a]).cross(points[c]-points[a]).norm()<1e-5)continue;
                    spin::TriangleMesh tetra;
                    tetra.vertices.resize(4,3);
                    tetra.vertices.row(0)=points[a];tetra.vertices.row(1)=points[b];tetra.vertices.row(2)=points[c];
                    tetra.vertices.row(3)=points[a]*.5;
                    Q boxUpper(-1);
                    for(int axis=0;axis<3;++axis) {
                        Q lo(tetra.vertices.col(axis).minCoeff()),hi(tetra.vertices.col(axis).maxCoeff());
                        boxUpper+=std::max(lo*lo,hi*hi);
                    }
                    if(!(boxUpper<Q(0)))continue;
                    tetra.faces.resize(4,3);tetra.faces<<0,2,1,0,1,3,1,2,3,2,0,3;
                    double volume6=0;
                    for(int face=0;face<4;++face)
                        volume6+=tetra.vertices.row(tetra.faces(face,0)).dot(
                            Eigen::Vector3d(tetra.vertices.row(tetra.faces(face,1))).cross(
                            Eigen::Vector3d(tetra.vertices.row(tetra.faces(face,2)))));
                    if(volume6<0)tetra.faces.col(1).swap(tetra.faces.col(2));
                    r.designDomain=tetra;
                    return r;
                }
    }
    throw std::runtime_error("could not construct exact discrete-contact fixture");
}

inline int run() {
    const auto r=request();
    spin::QuadricMomentEvaluatorSettings s;
    s.triangulation.icosphereSubdivisions=1;
    s.emptyRegionPolicy=spin::QuadricEmptyRegionPolicy::CertifyStrictlySeparatedV1;
    // Mutation caught: rejecting every empty Boolean result before adding the
    // fixed shell, or accepting empty but dropping fixed moments.
    const auto e=spin::evaluateQuadricMaterialMoments(r,farSphere(),s);
    require(e.materialRegion==spin::QuadricMaterialRegion::CertifiedEmpty,"missing certified empty state");
    require(e.variableMomentsTarget.values().isZero(0),"empty variable moments must be exactly zero");
    require((e.totalMomentsTarget.values()-r.fixedShellMoments.values()).norm()<1e-12,"fixed shell lost");
    require(e.derivativeDefined && e.momentJacobian.isZero(0),"strictly separated state must have zero derivative");
    require(e.materialBoundary.faces.rows()==0 && e.materialBoundaryTarget.faces.rows()==0,"empty geometry fabricated");
    require(!e.emptyRegionCertificate.empty() && e.geometryHash.size()==64 && e.derivativeScientificHash.size()==64,"missing empty evidence identity");
    std::cout<<"[PASS] certified-empty basic fixed-shell state\n";

    const auto again=spin::evaluateQuadricMaterialMoments(r,farSphere(),s);
    require(e.geometryHash==again.geometryHash && e.derivativeScientificHash==again.derivativeScientificHash,
        "certificate identity must be deterministic");
    for (int coordinate=0;coordinate<10;++coordinate) for (int sign : {-1,1}) {
        V p=farSphere().values(); p(coordinate)+=sign*1e-6;
        const auto near=spin::evaluateQuadricMaterialMoments(r,spin::QuadricCoefficients10::fromUnnormalized(p),s);
        require(near.materialRegion==spin::QuadricMaterialRegion::CertifiedEmpty &&
            near.variableMomentsTarget.values().isZero(0) && near.momentJacobian.isZero(0),
            "strict empty neighborhood not locally constant");
        require(near.geometryHash!=e.geometryHash,"empty identity not bound to candidate parameters");
    }
    std::cout<<"[PASS] certified-empty deterministic identity and 20 parameter perturbations\n";

    auto discrete=s;
    discrete.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteEllipsoid;
    const auto de=spin::evaluateQuadricMaterialMoments(r,farSphere(),discrete);
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    require(de.derivativeDefined && de.momentJacobian.isZero(0),"discrete empty zero derivative missing");
#else
    require(!de.derivativeDefined,"discrete derivative must be withheld without exact clearance backend");
#endif
    require(de.derivativeScientificHash!=e.derivativeScientificHash,"derivative method hash collision");
    V negative=V::Zero(); negative(0)=-100; negative.segment<3>(4).setOnes();
    const auto outside=spin::evaluateQuadricMaterialMoments(r,spin::QuadricCoefficients10::fromUnnormalized(negative),s);
    require(outside.materialRegion==spin::QuadricMaterialRegion::CertifiedEmpty,"enclosing sphere negative region must be empty");
    std::cout<<"[PASS] certified-empty both methods and negative-region Boolean\n";

    V unitExterior=V::Zero();unitExterior(0)=-1;unitExterior.segment<3>(4).setOnes();
    const auto exterior=spin::QuadricCoefficients10::fromUnnormalized(unitExterior);
    const auto contactRequest=discreteFaceContact(exterior,discrete);
    const auto contact=spin::evaluateQuadricMaterialMoments(contactRequest,exterior,discrete);
    require(contact.materialRegion==spin::QuadricMaterialRegion::CertifiedEmpty && !contact.derivativeDefined,
        "analytic separation must not certify a discrete coplanar-contact derivative");
    std::cout<<"[PASS] certified-empty exact discrete face-contact derivative withheld\n";

    auto transformed=r;
    const Eigen::Matrix3d R=Eigen::AngleAxisd(.4,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
    const Eigen::Vector3d t(3,-2,1);
    for(auto* mesh : {&transformed.outerSurface,&transformed.designDomain})
        for(int i=0;i<mesh->vertices.rows();++i)
            mesh->vertices.row(i)=(R*mesh->vertices.row(i).transpose()+t).transpose();
    transformed.targetFrame.axes=R; transformed.targetFrame.origin=t;
    MassMoments fm;
    fm.mass=95.2; fm.firstMoment=95.2*t;
    fm.secondMoment=Eigen::Matrix3d::Identity()*(248./3.*1.7)+95.2*t*t.transpose();
    transformed.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(fm);
    const auto frame=spin::evaluateQuadricMaterialMoments(transformed,farSphere(),s);
    require((frame.totalMomentsTarget.values()-r.fixedShellMoments.values()).norm()<1e-10,
        "empty fixed-shell target transform incorrect");
    require(frame.geometryHash!=e.geometryHash,"empty identity ignores frame/request");
    std::cout<<"[PASS] certified-empty nonidentity frame and fixed-shell transform\n";

    auto legacy=s; legacy.emptyRegionPolicy=spin::QuadricEmptyRegionPolicy::RejectEmpty;
    mustReject([&]{spin::evaluateQuadricMaterialMoments(r,farSphere(),legacy);},"legacy empty behavior changed");
    require(spin::computeQuadricEvaluatorRequestHash(r,legacy)!=spin::computeQuadricEvaluatorRequestHash(r,s),
        "empty policy not bound to request hash");
    spin::QuadricMomentEvaluatorSettings defaults;
    require(spin::computeQuadricEvaluatorRequestHash(r,defaults)!=spin::computeSpinDesignRequestHash(r),
        "current cap-radius default omitted from evaluator request identity");
    auto historicalDefaults=defaults;
    historicalDefaults.triangulation.maximumRadiusDoublings=8;
    require(spin::computeQuadricEvaluatorRequestHash(r,historicalDefaults)==spin::computeSpinDesignRequestHash(r),
        "historical evaluator request hash compatibility changed");
    V ball=V::Zero();ball(0)=.01;ball.segment<3>(4).setConstant(-1);
    const auto small=spin::QuadricCoefficients10::fromUnnormalized(ball);
    const auto n=spin::evaluateQuadricMaterialMoments(r,small,s);
    const auto old=spin::evaluateQuadricMaterialMoments(r,small,legacy);
    require(n.materialRegion==spin::QuadricMaterialRegion::Nonempty && n.variableMomentsTarget.mass()>0,
        "negative vertex samples incorrectly declared interior empty");
    require(n.geometryHash==old.geometryHash && n.derivativeScientificHash!=old.derivativeScientificHash &&
        (n.variableMomentsTarget.values()-old.variableMomentsTarget.values()).isZero(0),
        "opt-in altered nonempty evaluation");
    ball(0)=100;
    const auto full=spin::evaluateQuadricMaterialMoments(r,spin::QuadricCoefficients10::fromUnnormalized(ball),s);
    require(full.materialRegion==spin::QuadricMaterialRegion::Nonempty &&
        std::abs(full.variableMomentsTarget.mass()-13.6)<1e-10,"full domain was lost");
    std::cout<<"[PASS] certified-empty legacy compatibility, nonempty and full states\n";

    // Removing either production role's rejection must fail this test. Real
    // evaluator outputs are used; no Newton/evaluator return value is mocked.
    for (const auto role : {spin::detail::QuadricExportRole::VariableMaterial,
                            spin::detail::QuadricExportRole::Cavity}) {
        bool emptyRejected=false;
        try { spin::detail::requireNonemptyQuadricExportRegion(e,role); }
        catch (const std::runtime_error&) { emptyRejected=true; }
        require(emptyRejected,"export gate accepted real certified-empty evaluation");
        spin::detail::requireNonemptyQuadricExportRegion(n,role);
        spin::detail::requireNonemptyQuadricExportRegion(full,role);
        std::cout<<"[PASS] certified-empty direct export gate role="
                 <<(role==spin::detail::QuadricExportRole::VariableMaterial?"material":"cavity")<<'\n';
    }

    V tangent=V::Zero();tangent(0)=-3;tangent(1)=4;tangent.segment<3>(4).setConstant(-1);
    mustReject([&]{spin::evaluateQuadricMaterialMoments(r,spin::QuadricCoefficients10::fromUnnormalized(tangent),s);},
        "contact must not be certified as strictly separated");
    // A small real cap may be unresolved by tessellation. Either a nonempty
    // result or a geometric rejection is acceptable; false empty is not.
    const double center=1.999;
    tangent(0)=1-center*center;tangent(1)=2*center;
    bool capReturned=false; spin::QuadricMomentEvaluation cap;
    try { cap=spin::evaluateQuadricMaterialMoments(r,spin::QuadricCoefficients10::fromUnnormalized(tangent),s); capReturned=true; }
    catch(const std::runtime_error&) {}
    if(capReturned)require(cap.materialRegion==spin::QuadricMaterialRegion::Nonempty && cap.variableMomentsTarget.mass()>0,
        "small positive material cap falsely certified empty");
    std::cout<<"[PASS] certified-empty contact and small-cap no-false-empty checks\n";

    auto bad=s;bad.emptyRegionPolicy=static_cast<spin::QuadricEmptyRegionPolicy>(99);
    mustReject([&]{spin::validateQuadricMomentEvaluatorSettings(bad);},"unknown empty policy accepted");
    auto invalid=r;invalid.materialDensity=std::numeric_limits<double>::quiet_NaN();
    mustReject([&]{spin::evaluateQuadricMaterialMoments(invalid,farSphere(),s);},"nonfinite density accepted");
    invalid=r; invalid.designDomain.faces.conservativeResize(11,3);
    mustReject([&]{spin::evaluateQuadricMaterialMoments(invalid,farSphere(),s);},"open design accepted");
    invalid=r; invalid.fixedShellMoments=physics::RawMassMoments10();
    mustReject([&]{spin::evaluateQuadricMaterialMoments(invalid,farSphere(),s);},"zero total mass accepted");
    std::cout<<"[PASS] certified-empty invalid inputs remain fail-closed\n";

    auto solveRequest=r;
    solveRequest.settings.maximumIterations=1;
    solveRequest.settings.maximumRestarts=1;
    const auto solved=spin::solveQuadricHollowing(solveRequest,s,{farSphere()});
    std::cout<<"[INFO] empty export status="<<spin::spinSolverStatusToken(solved.status)
        <<" post_solve_calls="<<solved.postSolveMaterialEvaluationCount<<" diagnostic="<<solved.diagnostic<<'\n';
    require(solved.status!=spin::SpinSolverStatus::Converged && solved.materialBoundary.faces.rows()==0 && solved.geometryHash.empty(),
        "moment-only empty state incorrectly delivered as a completed geometry");
    std::cout<<"[PASS] certified-empty downstream remains fail-closed (no geometry delivery)\n";
    return 0;
}
} // namespace certified_empty_tests

inline int certifiedEmptyRegression() {
    try { return certified_empty_tests::run(); }
    catch(const std::exception& e) {
        std::cerr<<"[FAIL] certified-empty: "<<e.what()<<'\n'; return 1;
    }
}
