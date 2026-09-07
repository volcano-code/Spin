#include "spin/QuadricMomentEvaluator.h"
#include "spin/QuadricHollowingSolver.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>

namespace {
using Vector = Eigen::Matrix<double, 10, 1>;
using Matrix = Eigen::Matrix<double, 10, 10>;
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
spin::TriangleMesh cube(double lo, double hi) {
    spin::TriangleMesh m;
    m.vertices.resize(8,3);
    m.vertices << lo,lo,lo, hi,lo,lo, hi,hi,lo, lo,hi,lo,
                  lo,lo,hi, hi,lo,hi, hi,hi,hi, lo,hi,hi;
    m.faces.resize(12,3);
    m.faces << 0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,
               3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
    return m;
}
spin::SpinDesignRequest request(double half, bool frame) {
    spin::SpinDesignRequest r;
    r.designDomain = cube(-half,half);
    r.outerSurface = cube(-half-1,half+1);
    r.minimumWallThickness = .5;
    r.materialDensity = 1.7;
    if (frame) {
        r.targetFrame.origin << .31,-.24,.18;
        r.targetFrame.axes = Eigen::AngleAxisd(.43,
            Eigen::Vector3d(1,2,-1).normalized()).toRotationMatrix();
        for (auto* m : {&r.designDomain,&r.outerSurface})
            for (int v=0; v<m->vertices.rows(); ++v)
                m->vertices.row(v) = r.targetFrame.axes *
                    m->vertices.row(v).transpose() + r.targetFrame.origin;
    }
    return r;
}
// Independent signed-tetrahedron integration, not the production surface kernel.
Vector tetraMoments(const spin::TriangleMesh& mesh, double density) {
    Vector out = Vector::Zero();
    for (int f=0; f<mesh.faces.rows(); ++f) {
        const Eigen::Vector3d a=mesh.vertices.row(mesh.faces(f,0));
        const Eigen::Vector3d b=mesh.vertices.row(mesh.faces(f,1));
        const Eigen::Vector3d c=mesh.vertices.row(mesh.faces(f,2));
        const double v = density*a.dot(b.cross(c))/6;
        const Eigen::Vector3d s=a+b+c;
        const Eigen::Matrix3d second =
            v*(s*s.transpose()+a*a.transpose()+b*b.transpose()+c*c.transpose())/20;
        out(0)+=v;
        out.segment<3>(1)+=v*s/4;
        out.segment<3>(4)+=second.diagonal();
        out(7)+=second(1,2); out(8)+=second(0,2); out(9)+=second(0,1);
    }
    return out;
}

std::set<std::string> correspondence(const spin::QuadricMomentEvaluation& e) {
    std::set<std::string> keys(e.vertexProvenance.begin(),e.vertexProvenance.end());
    require(keys.size()==e.vertexProvenance.size() &&
            keys.size()==static_cast<std::size_t>(e.materialBoundaryTarget.vertices.rows()),
            "missing or ambiguous output vertex identities");
    return keys;
}

int verifyProvenance(const spin::QuadricMomentEvaluation& e,
                     const spin::SpinDesignRequest& r) {
    auto design=r.designDomain;
    for(int i=0;i<design.vertices.rows();++i)
        design.vertices.row(i)=r.targetFrame.axes.transpose()*
            (design.vertices.row(i).transpose()-r.targetFrame.origin);
    int diagonalCount=0;
    for(int v=0;v<e.materialBoundaryTarget.vertices.rows();++v) {
        const Eigen::Vector3d x=e.materialBoundaryTarget.vertices.row(v);
        const auto& key=e.vertexProvenance[v];
        const auto separator=key.find('|');
        require(separator!=std::string::npos,"missing source separator");
        for(const auto& token : {key.substr(0,separator),key.substr(separator+1)}) {
            if(token.empty()) continue;
            const auto& mesh=token[0]=='d'?design:e.quadricTriangulation.closedBoundary;
            const int first=std::stoi(token.substr(3));
            if(token[2]=='v') {
                require((x-mesh.vertices.row(first).transpose()).norm()<1e-10,"source vertex mismatch");
            } else if(token[2]=='e') {
                const int second=std::stoi(token.substr(token.find(',')+1));
                const Eigen::Vector3d a=mesh.vertices.row(first),b=mesh.vertices.row(second);
                const double t=(x-a).dot(b-a)/(b-a).squaredNorm();
                require(t>0 && t<1 && (x-((1-t)*a+t*b)).norm()<1e-10,"source edge mismatch");
                if(token[0]=='d') {
                    std::vector<Eigen::Vector3d> normals;
                    for(int f=0;f<mesh.faces.rows();++f) {
                        bool hasFirst=false,hasSecond=false;
                        for(int k=0;k<3;++k) {hasFirst|=mesh.faces(f,k)==first;hasSecond|=mesh.faces(f,k)==second;}
                        if(hasFirst && hasSecond) {
                            const Eigen::Vector3d aa=mesh.vertices.row(mesh.faces(f,0));
                            const Eigen::Vector3d bb=mesh.vertices.row(mesh.faces(f,1));
                            const Eigen::Vector3d cc=mesh.vertices.row(mesh.faces(f,2));
                            normals.push_back((bb-aa).cross(cc-aa).normalized());
                        }
                    }
                    if(normals.size()==2 && normals[0].cross(normals[1]).norm()<1e-12) ++diagonalCount;
                }
            } else {
                require(token[2]=='f',"unknown provenance feature");
                const Eigen::Vector3d a=mesh.vertices.row(mesh.faces(first,0));
                const Eigen::Vector3d b=mesh.vertices.row(mesh.faces(first,1));
                const Eigen::Vector3d c=mesh.vertices.row(mesh.faces(first,2));
                const double area=(b-a).cross(c-a).norm();
                const double pieces=(a-x).cross(b-x).norm()+(b-x).cross(c-x).norm()+(c-x).cross(a-x).norm();
                require(std::abs(pieces-area)<1e-10 &&
                        std::abs((x-a).dot((b-a).cross(c-a)))<1e-10,"source face mismatch");
            }
        }
    }
    return diagonalCount;
}

void check(const std::string& name, const spin::SpinDesignRequest& r,
           const spin::QuadricCoefficients10& q, int subdivisions, bool clipped) {
    const auto before=spin::computeSpinDesignRequestHash(r);
    spin::QuadricMomentEvaluatorSettings settings;
    settings.triangulation.icosphereSubdivisions=subdivisions;
    const auto legacy=spin::evaluateQuadricMaterialMoments(r,q,settings);
    settings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteEllipsoid;
    const auto base=spin::evaluateQuadricMaterialMoments(r,q,settings);
    require(base.derivativeDefined, name+": "+base.derivativeDiagnostic);
    require(base.derivativeMethod==settings.derivativeMethod, "method lost");
    require(base.geometryHash==legacy.geometryHash &&
            base.quadricTriangulation.hash==legacy.quadricTriangulation.hash,
            "opt-in changed geometry");
    require(base.derivativeScientificHash!=legacy.derivativeScientificHash,
            "scientific identity aliases legacy");
    require((base.momentJacobian*q.values()).norm()<2e-10,"radial null");
    require((tetraMoments(base.materialBoundaryTarget,r.materialDensity)-
        base.variableMomentsTarget.values()).norm()<2e-10,"tetra oracle");
    if(clipped) require(base.constructedIntersectionVertexCount>0,
                        "fixture must actually intersect");
    const auto repeated=spin::evaluateQuadricMaterialMoments(r,q,settings);
    const auto keys=correspondence(base);
    const int diagonalCount=verifyProvenance(base,r);
    if(clipped) require(diagonalCount>0,"clipped coverage lacks coplanar domain diagonal construction");
    require(base.geometryHash==repeated.geometryHash &&
            base.derivativeScientificHash==repeated.derivativeScientificHash &&
            base.provenanceTopologyHash==repeated.provenanceTopologyHash &&
            base.momentJacobian==repeated.momentJacobian,"nondeterminism");
    require(before==spin::computeSpinDesignRequestHash(r),"input mutation");
    std::cout << name << " intersections=" << base.constructedIntersectionVertexCount
              << " coplanar_diagonal_intersections=" << diagonalCount
              << " constraint_residual=" << base.maximumConstraintResidual
              << " constraint_singular_ratio=" << base.minimumConstraintSingularRatio
              << " geometry=" << base.geometryHash << " provenance="
              << base.provenanceTopologyHash << " scientific="
              << base.derivativeScientificHash << '\n';
    std::vector<Matrix> differences;
    std::vector<double> relativeErrors;
    for(double h : {3e-4,1e-4,3e-5}) {
        Matrix fd;
        for(int j=0;j<10;++j) {
            Vector p=q.values(),m=q.values(); p(j)+=h; m(j)-=h;
            const auto plus=spin::evaluateQuadricMaterialMoments(r,
                spin::QuadricCoefficients10::fromUnnormalized(p),settings);
            const auto minus=spin::evaluateQuadricMaterialMoments(r,
                spin::QuadricCoefficients10::fromUnnormalized(m),settings);
            require(plus.derivativeDefined && minus.derivativeDefined,
                    "perturbed applicability changed");
            require(keys==correspondence(plus) && keys==correspondence(minus),
                    "source feature correspondence changed");
            (void)verifyProvenance(plus,r); (void)verifyProvenance(minus,r);
            require(base.provenanceTopologyHash==plus.provenanceTopologyHash &&
                    base.provenanceTopologyHash==minus.provenanceTopologyHash,
                    name+": topology/provenance change at h="+std::to_string(h));
            fd.col(j)=(plus.variableMomentsTarget.values()-
                       minus.variableMomentsTarget.values())/(2*h);
        }
        const double rel=(fd-base.momentJacobian).norm()/fd.norm();
        differences.push_back(fd); relativeErrors.push_back(rel);
        double entryRel=0,entryAbs=0;
        for(int i=0;i<10;++i) for(int j=0;j<10;++j) {
            const double scale=std::max(std::abs(fd(i,j)),std::abs(base.momentJacobian(i,j)));
            const double error=std::abs(fd(i,j)-base.momentJacobian(i,j));
            if(scale>1e-5) entryRel=std::max(entryRel,error/scale);
            else entryAbs=std::max(entryAbs,error);
        }
        std::cout << name << " h=" << h << " relative_norm=" << rel
                  << " entry_relative=" << entryRel << " near_zero_absolute="
                  << entryAbs << std::endl;
        if(h==3e-5) {
            require(rel<=1e-6,"fine full-map 1e-6 gate");
            require(entryRel<=1e-4 && entryAbs<=1e-8,"entrywise derivative gate");
        }
    }
    // Supplemental controller-required check added after raw FD results:
    // remove O(h^2) central-difference truncation using the SAME preset steps.
    // This does not replace any raw error or the fine-step full-map gate.
    require(relativeErrors[2]<relativeErrors[1] && relativeErrors[1]<relativeErrors[0],
            "preset finite-difference errors do not decrease");
    const Matrix richardsonCoarse=(9*differences[1]-differences[0])/8;
    const Matrix richardsonFine=(100*differences[2]-9*differences[1])/91;
    double richardsonRelative=0,richardsonAbsolute=0,extrapolationAgreement=0;
    for(int i=0;i<10;++i) for(int j=0;j<10;++j) {
        const double scale=std::max(std::abs(richardsonFine(i,j)),std::abs(base.momentJacobian(i,j)));
        const double error=std::abs(richardsonFine(i,j)-base.momentJacobian(i,j));
        if(scale>1e-5) {
            richardsonRelative=std::max(richardsonRelative,error/scale);
            extrapolationAgreement=std::max(extrapolationAgreement,
                std::abs(richardsonFine(i,j)-richardsonCoarse(i,j))/scale);
        } else {
            richardsonAbsolute=std::max(richardsonAbsolute,error);
            require(std::abs(richardsonFine(i,j)-richardsonCoarse(i,j))<=1e-8,
                    "near-zero extrapolation disagreement");
        }
    }
    std::cout << name << " supplemental_richardson_entry_relative=" << richardsonRelative
              << " near_zero_absolute=" << richardsonAbsolute
              << " two_extrapolations_entry_agreement=" << extrapolationAgreement
              << " raw_error_ratio_middle_over_coarse=" << relativeErrors[1]/relativeErrors[0]
              << " raw_error_ratio_fine_over_middle=" << relativeErrors[2]/relativeErrors[1] << '\n';
    require(richardsonRelative<=1e-6 && richardsonAbsolute<=1e-8 && extrapolationAgreement<=1e-6,
            "supplemental strict1e-6 Richardson entrywise gate");
    const auto complement=spin::evaluateQuadricMaterialMoments(r,q.inverted(),settings);
    // r's domain is a rigidly transformed cube: recover it independently.
    spin::TriangleMesh target=r.designDomain;
    for(int i=0;i<target.vertices.rows();++i)
        target.vertices.row(i)=r.targetFrame.axes.transpose()*
            (target.vertices.row(i).transpose()-r.targetFrame.origin);
    require((base.variableMomentsTarget.values()+complement.variableMomentsTarget.values()-
        tetraMoments(target,r.materialDensity)).norm()<3e-10,"complement moments");
    // At -q, both coefficient negation and material complementation negate.
    require((base.momentJacobian-complement.momentJacobian).norm()<3e-9,
            "complement differential");
    auto unit=r; unit.materialDensity=1;
    const auto u=spin::evaluateQuadricMaterialMoments(unit,q,settings);
    require((base.momentJacobian-r.materialDensity*u.momentJacobian).norm()<3e-10,
            "density differential");
}
}

int main() {
    try {
        std::cout << std::setprecision(17);
        Vector p=Vector::Zero(); p(0)=.25; p.segment<3>(4).setConstant(-1);
        auto sphere=spin::QuadricCoefficients10::fromUnnormalized(p);
        auto sr=request(1,false); sr.materialDensity=1;
        check("historical-sphere",sr,sphere,3,false);
        spin::QuadricMomentEvaluatorSettings sphereSettings;
        sphereSettings.triangulation.icosphereSubdivisions=3;
        const auto sphereGeometry=spin::evaluateQuadricMaterialMoments(sr,sphere,sphereSettings);
        const double sphereTriangulationError=1-sphereGeometry.variableMomentsTarget.mass()/(std::acos(-1.)/6);
        std::cout << "sphere analytic-volume triangulation_relative_error=" << sphereTriangulationError << '\n';
        require(sphereTriangulationError>1e-5 && sphereTriangulationError<.01,
                "sphere triangulation error must remain distinct from differentiation error");
        p << -1.4,.2,-.4,.1,1.1,.8,1.3,.12,-.08,.17;
        auto generic=spin::QuadricCoefficients10::fromUnnormalized(p);
        for(bool inverted : {false,true}) {
            auto q=inverted?generic.inverted():generic;
            check(inverted?"contained-negative":"contained-positive",request(2,true),q,1,false);
            check(inverted?"clipped-negative":"clipped-positive",request(1,true),q,1,true);
        }
        spin::QuadricMomentEvaluatorSettings settings;
        settings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteEllipsoid;
        settings.triangulation.icosphereSubdivisions=1;
        p.setZero(); p(1)=1;
        auto plane=spin::QuadricCoefficients10::fromUnnormalized(p);
        auto r=request(1,false); r.settings.maximumRestarts=1;
        const auto unsupported=spin::evaluateQuadricMaterialMoments(r,plane,settings);
        require(!unsupported.derivativeDefined && !unsupported.derivativeDiagnostic.empty()
                && unsupported.materialBoundary.faces.rows()>0,"unsupported plane fallback");
        auto solved=spin::solveQuadricHollowing(r,settings,{plane});
        require(solved.status!=spin::SpinSolverStatus::Converged &&
                solved.diagnostic.find(unsupported.derivativeDiagnostic)!=std::string::npos,
                "solver lost undefined differential reason");
        require(solved.derivativeMethod==spin::quadricMomentDerivativeMethodToken(settings.derivativeMethod),
                "solver lost selected method");
        p.setZero(); p(0)=1; p.segment<3>(4).setConstant(-1);
        const auto tangent=spin::evaluateQuadricMaterialMoments(r,
            spin::QuadricCoefficients10::fromUnnormalized(p),settings);
        require(!tangent.derivativeDefined &&
                tangent.derivativeDiagnostic.find("contact")!=std::string::npos &&
                tangent.materialBoundary.faces.rows()>0,"tangent contact was certified");
        std::cout << "[PASS] discrete ellipsoid moment differential\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "[FAIL] " << e.what() << '\n'; return 1;
    }
}
