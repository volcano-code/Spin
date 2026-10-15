#include "spin/PlaneMomentEvaluator.h"
#include "manufacturing/SlicerValidation.h"
#include <CGAL/Gmpq.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace {
using Q=CGAL::Gmpq;
using Q3=std::array<Q,3>;
Q3 point(const Eigen::MatrixXd& v,int i) {
    return {{Q(v(i,0)),Q(v(i,1)),Q(v(i,2))}};
}
Q3 subtract(const Q3& a,const Q3& b) {
    return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};
}
Q3 cross(const Q3& a,const Q3& b) {
    return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
}
Q dot(const Q3& a,const Q3& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}

// The finite plane construction is a convex extrusion. Strict containment
// behind every outward cap/side face excludes truncation from the design
// domain, including nonconvex design domains. Exact represented-double signs.
std::string capCertificate(const spin::QuadricTriangulation& tri,const Eigen::MatrixXd& design) {
    if(tri.type!=spin::QuadricSurfaceType::Plane || tri.interfaceFaceCount!=2 ||
       tri.faceIsCap.size()!=tri.closedBoundary.faces.rows())
        throw std::runtime_error("unsupported native plane construction");
    const auto& mesh=tri.closedBoundary;
    std::ostringstream proof;proof.imbue(std::locale::classic());
    proof<<"PlaneCapStrictHalfspaces-v1;"<<tri.hash<<';';
    for(int f=0;f<mesh.faces.rows();++f) {
        if(!tri.faceIsCap(f))continue;
        const Q3 a=point(mesh.vertices,mesh.faces(f,0));
        const Q3 n=cross(subtract(point(mesh.vertices,mesh.faces(f,1)),a),
                         subtract(point(mesh.vertices,mesh.faces(f,2)),a));
        if(dot(n,n)==Q(0))throw std::runtime_error("degenerate plane cap");
        Q upper=dot(n,subtract(point(design,0),a));
        for(int v=1;v<design.rows();++v)upper=std::max(upper,dot(n,subtract(point(design,v),a)));
        if(!(upper<Q(0)))throw std::runtime_error("plane truncation cap intersects or fails to cover design domain");
        proof<<f<<':'<<upper<<';';
    }
    return proof.str();
}

Eigen::Matrix<double,10,1> monomials(const Eigen::Vector3d& p) {
    Eigen::Matrix<double,10,1> m;
    m<<1,p.x(),p.y(),p.z(),p.x()*p.x(),p.y()*p.y(),p.z()*p.z(),p.y()*p.z(),p.x()*p.z(),p.x()*p.y();
    return m;
}

bool nativeInterfaceSeparated(const spin::QuadricTriangulation& tri,
                              const Eigen::MatrixXd& design,bool inside) {
    const auto& mesh=tri.closedBoundary;
    for(int f=0;f<tri.interfaceFaceCount;++f) {
        const Q3 a=point(mesh.vertices,mesh.faces(f,0));
        const Q3 n=cross(subtract(point(mesh.vertices,mesh.faces(f,1)),a),
                         subtract(point(mesh.vertices,mesh.faces(f,2)),a));
        if(dot(n,n)==Q(0))return false;
        for(int v=0;v<design.rows();++v) {
            const Q d=dot(n,subtract(point(design,v),a));
            if(inside ? !(d<Q(0)) : !(d>Q(0)))return false;
        }
    }
    return true;
}

using Velocity=Eigen::Matrix<double,3,4>;
std::array<Velocity,4> sourceVelocities(const spin::QuadricCoefficients10& q,
                                      const spin::QuadricTriangulation& tri) {
    const Eigen::Vector3d b=q.linearVector(),n=b.normalized();
    const double length=b.norm(),r=tri.truncationRadius,a=q.constant();
    Velocity db=Velocity::Zero();db.rightCols<3>().setIdentity();
    Eigen::RowVector4d da;da<<1,0,0,0;
    const Velocity dn=(Eigen::Matrix3d::Identity()-n*n.transpose())*db/length;
    const Eigen::Vector3d center=-a*b/b.squaredNorm();
    const Velocity dc=-b/b.squaredNorm()*da-a/b.squaredNorm()*db+
        (2*a/(b.squaredNorm()*b.squaredNorm()))*b*b.transpose()*db;
    const auto frame=spin::makeQuadricPlaneFrame(n);
    if(!frame.chartStable)
        throw std::runtime_error("plane tangent-frame chart boundary: "+frame.diagnostic);
    const Eigen::Vector3d first=frame.first,second=frame.second;
    const int referenceAxis=frame.chart==0?0:1;
    Eigen::Vector3d reference=Eigen::Vector3d::Zero();reference(referenceAxis)=1.0;
    const double frameNorm=n.cross(reference).norm();
    Velocity df;
    for(int j=0;j<4;++j) {
        const Eigen::Vector3d raw=dn.col(j).cross(reference);
        df.col(j)=(raw-first*first.dot(raw))/frameNorm;
    }
    Velocity ds;
    for(int j=0;j<4;++j)ds.col(j)=dn.col(j).cross(first)+n.cross(df.col(j));
    const int signs[4][2]={{-1,-1},{1,-1},{1,1},{-1,1}};
    std::array<Velocity,4> velocity;
    for(int v=0;v<4;++v) {
        const Eigen::Vector3d expected=center+signs[v][0]*r*first+signs[v][1]*r*second;
        if((expected-tri.closedBoundary.vertices.row(v).transpose()).norm()>1e-12*std::max(1.,expected.norm()))
            throw std::runtime_error("plane vertex differential does not match native construction");
        velocity[v]=dc+signs[v][0]*r*df+signs[v][1]*r*ds;
        if(!velocity[v].allFinite())throw std::runtime_error("nonfinite native plane vertex derivative");
    }
    return velocity;
}

bool inFacet(const Eigen::Vector3d& x,const Eigen::Vector3d& a,
             const Eigen::Vector3d& b,const Eigen::Vector3d& c,double tolerance) {
    const Eigen::Vector3d u=b-a,v=c-a,n=u.cross(v);
    if(n.norm()==0 || std::abs(n.dot(x-a))/n.norm()>tolerance)return false;
    const double s=(x-a).cross(v).dot(n)/n.squaredNorm();
    const double t=u.cross(x-a).dot(n)/n.squaredNorm();
    return std::isfinite(s) && std::isfinite(t) && s>=-1e-12 && t>=-1e-12 && s+t<=1+1e-12;
}
}

namespace spin {
PlaneMomentEvaluation evaluatePlaneMaterialMoments(
    const SpinDesignRequest& request, const Eigen::Vector4d& coefficients,
    const PlaneMomentEvaluatorSettings& settings)
{
    validateSpinDesignRequest(request);
    const double inputNorm=coefficients.norm();
    if(!coefficients.allFinite() || !std::isfinite(inputNorm) || inputNorm<=1e-15 ||
       !std::isfinite(settings.truncationRadius) || settings.truncationRadius<0)
        throw std::invalid_argument("invalid plane coefficients or truncation radius");
    QuadricCoefficients10::Vector p=QuadricCoefficients10::Vector::Zero();
    p.head<4>()=coefficients;
    const auto q=QuadricCoefficients10::fromUnnormalized(p);
    const Eigen::Vector4d normalized=q.values().head<4>();
    const double normalNorm=q.linearVector().norm();
    if(normalNorm<=1e-12)throw std::invalid_argument("plane normal is zero or ill-conditioned");
    Eigen::MatrixXd design(request.designDomain.vertices.rows(),3);
    for(int v=0;v<design.rows();++v)
        design.row(v)=(request.targetFrame.axes.transpose()*
            (request.designDomain.vertices.row(v).transpose()-request.targetFrame.origin)).transpose();
    if(!design.allFinite())throw std::invalid_argument("nonfinite target-frame design");
    QuadricMomentEvaluatorSettings nativeSettings;
    // This existing mode returns unavailable for planes; use only its values
    // and birth-face count, never its Jacobian. Avoid the legacy10 quadrature.
    nativeSettings.derivativeMethod=QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
    nativeSettings.triangulation.truncationRadius=settings.truncationRadius;
    nativeSettings.emptyRegionPolicy=QuadricEmptyRegionPolicy::CertifyStrictlySeparatedV1;
    const auto tri=triangulateQuadric(q,design,nativeSettings.triangulation);
    PlaneMomentEvaluation out;
    out.normalizedCoefficients=normalized;
    out.capSeparationCertificate=capCertificate(tri,design);
    const double scale=std::max(1.,design.cwiseAbs().maxCoeff());
    const double tolerance=1e-10*scale;
    bool positive=true,negative=true,contact=false;
    for(int v=0;v<design.rows();++v) {
        const double d=q.evaluate(design.row(v))/normalNorm;
        if(!std::isfinite(d))throw std::invalid_argument("nonfinite plane signed distance");
        positive=positive && d>tolerance;negative=negative && d<-tolerance;
        contact=contact || std::abs(d)<=32*tolerance;
    }
    const auto geometry=evaluateQuadricMaterialMoments(request,q,nativeSettings);
    if(geometry.quadricTriangulation.hash!=tri.hash)
        throw std::runtime_error("plane triangulation changed between coverage and value evaluation");
    out.variableMomentsTarget=geometry.variableMomentsTarget;
    out.totalMomentsTarget=geometry.totalMomentsTarget;
    out.variableMaterialBoundaryTarget=geometry.materialBoundaryTarget;
    out.geometryHash=geometry.geometryHash;
    const auto finish=[&]() {
        std::ostringstream id;id.imbue(std::locale::classic());id<<std::hexfloat;
        id<<"PlaneFamilyNativeFacet4-v1;input-unnormalized-target;contact-32e-10;"
          <<computeQuadricEvaluatorRequestHash(request,nativeSettings)<<';'<<out.geometryHash<<';'
          <<out.capSeparationCertificate<<';'<<out.diagnostic<<';'<<out.interfaceTriangleCount<<';'
          <<static_cast<int>(out.materialRegion)<<';';
        id<<out.variableMaterialBoundaryTarget.vertices.rows()<<';'
          <<out.variableMaterialBoundaryTarget.faces.rows()<<';';
        for(Eigen::Index i=0;i<out.variableMaterialBoundaryTarget.vertices.rows();++i)
            for(int axis=0;axis<3;++axis)id<<out.variableMaterialBoundaryTarget.vertices(i,axis)<<';';
        for(Eigen::Index i=0;i<out.variableMaterialBoundaryTarget.faces.rows();++i)
            for(int corner=0;corner<3;++corner)id<<out.variableMaterialBoundaryTarget.faces(i,corner)<<';';
        for(int i=0;i<4;++i)id<<coefficients(i)<<';';
        id<<out.coefficientJacobian.has_value()<<';';
        if(out.coefficientJacobian)for(int i=0;i<10;++i)for(int j=0;j<4;++j)id<<(*out.coefficientJacobian)(i,j)<<';';
        out.scientificHash=manufacturing::computeContentSha256(id.str());
        return out;
    };
    if(contact) {out.diagnostic="plane/design vertex contact or near-contact; restricted derivative unavailable";return finish();}
    if(positive || negative) {
        if(negative && geometry.materialRegion!=QuadricMaterialRegion::CertifiedEmpty)
            throw std::runtime_error("empty plane state disagrees with native geometry");
        if(!nativeInterfaceSeparated(tri,design,positive) || geometry.interfaceTriangleCount!=0) {
            out.diagnostic="ideal full/empty sign disagrees with native interface separation";return finish();
        }
        out.coefficientJacobian=PlaneMomentEvaluation::Jacobian::Zero();
        out.materialRegion=positive?PlaneMaterialRegion::StrictFull:PlaneMaterialRegion::StrictEmpty;
        out.diagnostic="strict constant full/empty plane-family state";return finish();
    }
    if(geometry.materialRegion==QuadricMaterialRegion::CertifiedEmpty)
        throw std::runtime_error("mixed-sign plane unexpectedly returned empty geometry");
    if(tri.truncationRadius>1e6*scale) {
        out.diagnostic="truncation scale is too large for a conditioned plane derivative";return finish();
    }
    std::array<Velocity,4> velocity;
    try {velocity=sourceVelocities(q,tri);}
    catch(const std::runtime_error& e){out.diagnostic=e.what();return finish();}
    PlaneMomentEvaluation::Jacobian jac=PlaneMomentEvaluation::Jacobian::Zero();
    const auto& mesh=geometry.materialBoundaryTarget;
    const std::array<Eigen::Vector3d,4> bary={{Eigen::Vector3d(1./3,1./3,1./3),
        Eigen::Vector3d(.6,.2,.2),Eigen::Vector3d(.2,.6,.2),Eigen::Vector3d(.2,.2,.6)}};
    const std::array<double,4> weights={{-27./48,25./48,25./48,25./48}};
    for(int f=0;f<mesh.faces.rows();++f) {
        const Eigen::Vector3d a=mesh.vertices.row(mesh.faces(f,0)),b=mesh.vertices.row(mesh.faces(f,1)),c=mesh.vertices.row(mesh.faces(f,2));
        if(std::abs(q.evaluate(a))/normalNorm>tolerance ||
           std::abs(q.evaluate(b))/normalNorm>tolerance ||
           std::abs(q.evaluate(c))/normalNorm>tolerance)continue;
        const Eigen::Vector3d area=.5*(b-a).cross(c-a);
        const double areaNorm=area.norm();
        if(!std::isfinite(areaNorm) || areaNorm<=0 ||
           (area/areaNorm+q.linearVector()/normalNorm).norm()>1e-8) {
            out.diagnostic="plane interface orientation or coplanarity inconsistent";return finish();
        }
        std::vector<int> candidates;
        const auto& source=tri.closedBoundary;
        for(int sf=0;sf<tri.interfaceFaceCount;++sf) {
            const Eigen::Vector3d sa=source.vertices.row(source.faces(sf,0)),
                sb=source.vertices.row(source.faces(sf,1)),sc=source.vertices.row(source.faces(sf,2));
            if(inFacet(a,sa,sb,sc,tolerance) && inFacet(b,sa,sb,sc,tolerance) && inFacet(c,sa,sb,sc,tolerance))
                candidates.push_back(sf);
        }
        if(candidates.empty()) {out.diagnostic="no native interface facet contains Boolean triangle";return finish();}
        ++out.interfaceTriangleCount;
        for(int k=0;k<4;++k) {
            const auto& w=bary[k];const Eigen::Vector3d x=w(0)*a+w(1)*b+w(2)*c;
            std::optional<Eigen::RowVector4d> speed;
            for(int sf:candidates) {
                const int ia=source.faces(sf,0),ib=source.faces(sf,1),ic=source.faces(sf,2);
                const Eigen::Vector3d sa=source.vertices.row(ia),sb=source.vertices.row(ib),sc=source.vertices.row(ic);
                const Eigen::Vector3d normal=(sb-sa).cross(sc-sa);
                if(normal.dot(area)<=0){out.diagnostic="native interface orientation mismatch";return finish();}
                Eigen::RowVector4d current;
                for(int j=0;j<4;++j) {
                    const Eigen::Vector3d dnormal=(velocity[ib].col(j)-velocity[ia].col(j)).cross(sc-sa)+
                        (sb-sa).cross(velocity[ic].col(j)-velocity[ia].col(j));
                    current(j)=(normal.dot(velocity[ia].col(j))-dnormal.dot(x-sa))/normal.norm();
                }
                if(!current.allFinite() || (speed && (current-*speed).norm()>1e-8*std::max(1.,speed->norm()))) {
                    out.diagnostic="ambiguous native interface facet velocity";return finish();
                }
                speed=current;
            }
            jac+=request.materialDensity*weights[k]*areaNorm*monomials(x)**speed;
        }
    }
    if(out.interfaceTriangleCount==0 || out.interfaceTriangleCount!=geometry.interfaceTriangleCount) {
        out.diagnostic="plane interface identification disagrees with native Boolean provenance";return finish();
    }
    // Transport the native source facet, including its changing vertices and
    // normal (not the ideal plane substituted for that facet). Quadratic
    // monomials times affine normal motion are cubic. Chain normalization once.
    jac=jac*(Eigen::Matrix4d::Identity()-normalized*normalized.transpose())/inputNorm;
    if(!jac.allFinite()) {out.diagnostic="nonfinite plane-family differential";return finish();}
    out.coefficientJacobian=jac;
    out.diagnostic="restricted native plane-facet differential in real arithmetic; excludes rounding and topology/basis changes";
    return finish();
}
} // namespace spin
