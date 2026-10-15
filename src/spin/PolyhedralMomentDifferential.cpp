#include "spin/PolyhedralMomentDifferential.h"
#include "verification/ExactMaterialRegion.h"
#include "spin/QuadricMomentEvaluator.h"
#include "manufacturing/ManufacturingGeometry.h"
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactBoundaryAudit.h"
#include "verification/ExactRationalCoordinates.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <igl/copyleft/cgal/assign_scalar.h>
#endif
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <sstream>
#include <locale>
#include <stdexcept>

namespace {
using Jacobian = Eigen::MatrixXd;
using Row = Eigen::RowVectorXd;
constexpr double membershipTolerance = 2e-9;
constexpr double minimumSingularRatio = 1e-10;
constexpr double consistencyTolerance = 2e-9;

std::string classifiedFailure(const std::string& why) {
    const auto has=[&](const char* token){return why.find(token)!=std::string::npos;};
    if (has("unsupported") || has("ambiguous") || has("membership") ||
        has("contact") || has("source feature") || has("lacks both source"))
        return "TOPOLOGY_TRANSITION: "+why;
    if (has("rank deficient") || has("ill-conditioned"))
        return "RANK_DEFICIENT: "+why;
    if (has("invalid") || has("missing") || has("nonfinite"))
        return "INVALID_INPUT_OR_NONFINITE: "+why;
    return "NUMERICAL_FAILURE: "+why;
}

struct Feature {
    int dimension = 3;
    std::string key;
};

// Barycentric membership is conditioned independently of the motion solve.
bool barycentric(const spin::TriangleMesh& mesh, int f,
                 const Eigen::Vector3d& x, Eigen::Vector3d* bary) {
    const Eigen::Vector3d a=mesh.vertices.row(mesh.faces(f,0));
    const Eigen::Vector3d u=mesh.vertices.row(mesh.faces(f,1)).transpose()-a;
    const Eigen::Vector3d v=mesh.vertices.row(mesh.faces(f,2)).transpose()-a;
    const Eigen::Vector3d n=u.cross(v);
    const double scale=std::max({1.0,a.norm(),u.norm(),v.norm(),x.norm()});
    if(n.norm()<=1e-12*u.norm()*v.norm() ||
       std::abs(n.dot(x-a))>membershipTolerance*scale*n.norm()) return false;
    (*bary)(1)=(x-a).cross(v).dot(n)/n.squaredNorm();
    (*bary)(2)=u.cross(x-a).dot(n)/n.squaredNorm();
    (*bary)(0)=1-(*bary)(1)-(*bary)(2);
    return bary->allFinite() && bary->minCoeff()>=-membershipTolerance &&
           bary->maxCoeff()<=1+membershipTolerance;
}

void addFacetConstraints(const spin::TriangleMesh& mesh, int f,
    const std::vector<Jacobian>& velocity, const Eigen::Vector3d& x,
    const Eigen::Vector3d& bary, const std::string& prefix, Feature* feature,
    std::vector<Eigen::Vector3d>* lhs, std::vector<Row>* rhs) {
    std::array<int,3> ids={{mesh.faces(f,0),mesh.faces(f,1),mesh.faces(f,2)}};
    const Eigen::Vector3d a=mesh.vertices.row(ids[0]);
    const Eigen::Vector3d b=mesh.vertices.row(ids[1]);
    const Eigen::Vector3d c=mesh.vertices.row(ids[2]);
    const Eigen::Vector3d n=(b-a).cross(c-a);
    Row plane=Row::Zero(velocity[ids[0]].cols());
    for(Eigen::Index j=0;j<velocity[ids[0]].cols();++j) {
        const Eigen::Vector3d va=velocity[ids[0]].col(j);
        const Eigen::Vector3d vb=velocity[ids[1]].col(j);
        const Eigen::Vector3d vc=velocity[ids[2]].col(j);
        const Eigen::Vector3d dn=(vb-va).cross(c-a)+(b-a).cross(vc-va);
        plane(j)=(n.dot(va)-dn.dot(x-a))/n.norm();
    }
    lhs->push_back(n.normalized()); rhs->push_back(plane);
    Feature current{2,prefix+"f"+std::to_string(f)};
    for(int k=0;k<3;++k) {
        if(bary(k)>=1-membershipTolerance) {
            current={0,prefix+"v"+std::to_string(ids[k])};
            for(int axis=0;axis<3;++axis) {
                lhs->push_back(Eigen::Vector3d::Unit(axis));
                rhs->push_back(velocity[ids[k]].row(axis));
            }
            break;
        }
        if(std::abs(bary(k))<=membershipTolerance) {
            int first=ids[(k+1)%3], second=ids[(k+2)%3];
            if(first>second) std::swap(first,second);
            const Eigen::Vector3d edge=mesh.vertices.row(second)-mesh.vertices.row(first);
            const double s=edge.dot(x-mesh.vertices.row(first).transpose())/edge.squaredNorm();
            const Eigen::Matrix3d perpendicular=Eigen::Matrix3d::Identity()-
                edge*edge.transpose()/edge.squaredNorm();
            const Jacobian moving=perpendicular*((1-s)*velocity[first]+s*velocity[second]);
            for(int axis=0;axis<3;++axis) {
                lhs->push_back(perpendicular.row(axis)); rhs->push_back(moving.row(axis));
            }
            current={1,prefix+"e"+std::to_string(first)+","+std::to_string(second)};
            // Vertex detection takes precedence when two coordinates vanish.
            for(int j=0;j<3;++j) if(bary(j)>=1-membershipTolerance) {
                current={0,prefix+"v"+std::to_string(ids[j])};
                for(int axis=0;axis<3;++axis) {
                    lhs->push_back(Eigen::Vector3d::Unit(axis));
                    rhs->push_back(velocity[ids[j]].row(axis));
                }
            }
            break;
        }
    }
    if(current.dimension<feature->dimension) *feature=current;
    else if(current.dimension==feature->dimension && current.key!=feature->key)
        feature->key="ambiguous";
}

#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
using Rational=CGAL::Gmpq;
using ExactPoint=std::array<Rational,3>;
// Columns of a 3 x d velocity; no binary64 round-trip in the authority path.
using ExactJacobian=std::vector<ExactPoint>;
ExactPoint exactPoint(const Eigen::Vector3d& p) {
    return {{Rational(p(0)),Rational(p(1)),Rational(p(2))}};
}
ExactPoint subtract(const ExactPoint& a,const ExactPoint& b) {
    return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};
}
ExactPoint cross(const ExactPoint& a,const ExactPoint& b) {
    return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
}
Rational dot(const ExactPoint& a,const ExactPoint& b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
double booleanRound(const Rational& value) {
    double out;
    // The pinned Boolean uses this nearest assignment, not Gmpq::to_double's
    // possibly different approximation. The represented mesh must not change.
    igl::copyleft::cgal::assign_scalar(CGAL::Epeck::FT(value),out);
    return out;
}

bool addExactFacetConstraints(const spin::TriangleMesh& mesh,int f,
    const std::vector<Jacobian>& velocity,const ExactPoint& x,
    const std::string& prefix,Feature* feature,std::vector<Eigen::Vector3d>* lhs,
    std::vector<Row>* rhs,Jacobian* interpolated,
    const std::vector<ExactPoint>* authority=nullptr,
    const std::vector<ExactJacobian>* authorityVelocity=nullptr,
    std::vector<ExactPoint>* rationalLhs=nullptr,
    std::vector<std::vector<Rational>>* rationalRhs=nullptr,
    ExactJacobian* rationalInterpolated=nullptr,
    std::string* membershipFailure=nullptr,bool fixedFacet=false) {
    const auto rejectMembership=[&](const std::string& why) {
        if(membershipFailure)*membershipFailure=why;
        return false;
    };
    const std::array<int,3> ids={{mesh.faces(f,0),mesh.faces(f,1),mesh.faces(f,2)}};
    const Eigen::Vector3d da=mesh.vertices.row(ids[0]),db=mesh.vertices.row(ids[1]),dc=mesh.vertices.row(ids[2]);
    const auto point=[&](int index) {
        return authority?authority->at(index):exactPoint(mesh.vertices.row(index));
    };
    const auto a=point(ids[0]),b=point(ids[1]),c=point(ids[2]);
    const auto u=subtract(b,a),v=subtract(c,a),n=cross(u,v),offset=subtract(x,a);
    Eigen::Vector3d du=db-da,dv=dc-da;
    if(authority) {
        Rational scale(0);
        for(int axis=0;axis<3;++axis)scale=std::max(scale,std::max(CGAL::abs(u[axis]),CGAL::abs(v[axis])));
        if(scale==0)return rejectMembership("ZERO_EDGE_SCALE");
        for(int axis=0;axis<3;++axis) {
            du(axis)=booleanRound(u[axis]/scale);dv(axis)=booleanRound(v[axis]/scale);
        }
    }
    // Moving facets still require the established conditioning guard. Fixed
    // facets have exactly zero motion: their plane and membership are computed
    // rationally below, and no ill-conditioned barycentric velocity is needed.
    // A skinny tessellation must not invalidate an otherwise identical static
    // domain. Zero area, exact membership, and the later motion rank/consistency
    // guards remain mandatory; do not infer stationarity from small velocities.
    const bool stationaryFacet=fixedFacet && velocity[ids[0]].isZero(0) &&
        velocity[ids[1]].isZero(0) && velocity[ids[2]].isZero(0);
    if(!stationaryFacet && du.cross(dv).norm()<=1e-12*du.norm()*dv.norm()) {
        std::ostringstream detail;
        detail.precision(17);
        detail<<"SOURCE_FACET_CONDITIONING cross="<<du.cross(dv).norm()
              <<" edge_product="<<du.norm()*dv.norm();
        return rejectMembership(detail.str());
    }
    const auto squared=dot(n,n);
    if(squared==0)return rejectMembership("EXACT_ZERO_AREA");
    if(dot(n,offset)!=0)return rejectMembership("EXACT_OFF_PLANE");
    std::array<Rational,3> bary;
    bary[1]=dot(cross(offset,v),n)/squared;
    bary[2]=dot(cross(u,offset),n)/squared;
    bary[0]=Rational(1)-bary[1]-bary[2];
    for(int k=0;k<3;++k)if(bary[k]<0 || bary[k]>1) {
        std::ostringstream detail;
        detail<<"EXACT_OUTSIDE_TRIANGLE bary["<<k<<"]="<<bary[k];
        return rejectMembership(detail.str());
    }
    Eigen::Vector3d selector;
    for(int axis=0;axis<3;++axis)selector(axis)=CGAL::to_double(bary[axis]);

    const Eigen::Index dimensions=velocity[ids[0]].cols();
    Rational scale(0);
    for(const auto& value:n)scale=std::max(scale,CGAL::abs(value));
    Eigen::Vector3d normal;
    for(int axis=0;axis<3;++axis)normal(axis)=CGAL::to_double(n[axis]/scale);
    const double norm=normal.norm();
    Row plane(dimensions);
    interpolated->resize(3,dimensions);
    if(authorityVelocity)rationalInterpolated->resize(dimensions);
    std::vector<Rational> exactPlane(dimensions);
    for(Eigen::Index j=0;j<dimensions;++j) {
        ExactPoint motion;
        for(int axis=0;axis<3;++axis) {
            motion[axis]=Rational(0);
            for(int k=0;k<3;++k)motion[axis]+=bary[k]*(authorityVelocity?
                authorityVelocity->at(ids[k])[j][axis]:Rational(velocity[ids[k]](axis,j)));
            (*interpolated)(axis,j)=CGAL::to_double(motion[axis]);
        }
        if(authorityVelocity)(*rationalInterpolated)[j]=motion;
        // On the exact facet, n·v_interpolated is algebraically identical to
        // n·v_a - dn·(x-a), without injecting a rounded off-plane point.
        plane(j)=CGAL::to_double(dot(n,motion)/scale)/norm;
        exactPlane[j]=dot(n,motion)/scale;
    }
    lhs->push_back(normal/norm);rhs->push_back(plane);
    if(authorityVelocity) {
        auto row=n;for(auto& value:row)value/=scale;
        rationalLhs->push_back(row);rationalRhs->push_back(std::move(exactPlane));
    }
    const auto identityMotion=[&](int vertex) {
        if(!authorityVelocity)return;
        for(int axis=0;axis<3;++axis) {
            ExactPoint row{};row[axis]=1;rationalLhs->push_back(row);
            std::vector<Rational> values(dimensions);
            for(Eigen::Index j=0;j<dimensions;++j)values[j]=authorityVelocity->at(vertex)[j][axis];
            rationalRhs->push_back(std::move(values));
        }
    };
    Feature current{2,prefix+"f"+std::to_string(f)};
    for(int k=0;k<3;++k) {
        if(authority?bary[k]==1:selector(k)>=1-membershipTolerance) {
            current={0,prefix+"v"+std::to_string(ids[k])};
            for(int axis=0;axis<3;++axis) {
                lhs->push_back(Eigen::Vector3d::Unit(axis));rhs->push_back(velocity[ids[k]].row(axis));
            }
            identityMotion(ids[k]);
            break;
        }
        if(authority?bary[k]==0:std::abs(selector(k))<=membershipTolerance) {
            int first=ids[(k+1)%3],second=ids[(k+2)%3];
            if(first>second)std::swap(first,second);
            const auto firstPoint=point(first);
            const auto edge=subtract(point(second),firstPoint);
            const auto edgeSquared=dot(edge,edge);
            if(edgeSquared==0)return false;
            const auto t=dot(edge,subtract(x,firstPoint))/edgeSquared;
            for(int axis=0;axis<3;++axis) {
                Eigen::Vector3d perpendicular;
                Row moving(dimensions);
                ExactPoint exactPerpendicular{};
                for(int j=0;j<3;++j)perpendicular(j)=CGAL::to_double(
                    Rational(axis==j?1:0)-edge[axis]*edge[j]/edgeSquared);
                for(int j=0;j<3;++j)exactPerpendicular[j]=
                    Rational(axis==j?1:0)-edge[axis]*edge[j]/edgeSquared;
                std::vector<Rational> exactMoving(dimensions);
                for(Eigen::Index j=0;j<dimensions;++j) {
                    Rational sum(0);
                    for(int s=0;s<3;++s) {
                        const auto projection=Rational(axis==s?1:0)-edge[axis]*edge[s]/edgeSquared;
                        sum+=projection*((Rational(1)-t)*(authorityVelocity?
                            authorityVelocity->at(first)[j][s]:Rational(velocity[first](s,j)))+
                            t*(authorityVelocity?authorityVelocity->at(second)[j][s]:Rational(velocity[second](s,j))));
                    }
                    moving(j)=CGAL::to_double(sum);
                    exactMoving[j]=sum;
                }
                lhs->push_back(perpendicular);rhs->push_back(moving);
                if(authorityVelocity) {
                    rationalLhs->push_back(exactPerpendicular);
                    rationalRhs->push_back(std::move(exactMoving));
                }
            }
            current={1,prefix+"e"+std::to_string(first)+","+std::to_string(second)};
            for(int j=0;j<3;++j)if(authority?bary[j]==1:selector(j)>=1-membershipTolerance) {
                current={0,prefix+"v"+std::to_string(ids[j])};
                for(int axis=0;axis<3;++axis) {
                    lhs->push_back(Eigen::Vector3d::Unit(axis));rhs->push_back(velocity[ids[j]].row(axis));
                }
                identityMotion(ids[j]);
            }
            break;
        }
    }
    if(current.dimension<feature->dimension)*feature=current;
    else if(current.dimension==feature->dimension && current.key!=feature->key)feature->key="ambiguous";
    return interpolated->allFinite();
}
#endif

Eigen::Matrix<double,10,1> monomials(const Eigen::Vector3d& p) {
    Eigen::Matrix<double,10,1> m;
    m << 1,p.x(),p.y(),p.z(),p.x()*p.x(),p.y()*p.y(),p.z()*p.z(),
         p.y()*p.z(),p.x()*p.z(),p.x()*p.y();
    return m;
}
}

namespace spin {
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density) {
    return differentiatePolyhedralMaterialMoments(design,source,moving,material,birthFaces,density,false);
}
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection) {
    return differentiatePolyhedralMaterialMoments(design,source,moving,material,
        birthFaces,density,exactBinary64SelfIntersection,false);
}
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary) {
    return differentiatePolyhedralMaterialMoments(design,source,moving,material,
        birthFaces,density,exactBinary64SelfIntersection,exactBinary64SourceBoundary,false);
}
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary) {
    return differentiatePolyhedralMaterialMoments(design,source,moving,material,birthFaces,
        density,exactBinary64SelfIntersection,exactBinary64SourceBoundary,exactBinary64MaterialBoundary,{});
}
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary,const ExactBooleanVertexCoordinates& exactVertices) {
    return differentiatePolyhedralMaterialMoments(design,source,moving,material,birthFaces,
        density,exactBinary64SelfIntersection,exactBinary64SourceBoundary,
        exactBinary64MaterialBoundary,exactVertices,{},{});
}
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary,const ExactBooleanVertexCoordinates& exactVertices,
    const ExactBooleanVertexCoordinates& exactSource, const ExactBooleanVertexCoordinates& exactDesign,
    const std::function<void(const std::string&)>& diagnosticObserver) {
    return differentiatePolyhedralMaterialMoments(design,source,moving,material,birthFaces,density,
        exactBinary64SelfIntersection,exactBinary64SourceBoundary,exactBinary64MaterialBoundary,
        exactVertices,exactSource,exactDesign,diagnosticObserver,nullptr);
}
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary,const ExactBooleanVertexCoordinates& exactVertices,
    const ExactBooleanVertexCoordinates& exactSource, const ExactBooleanVertexCoordinates& exactDesign,
    const std::function<void(const std::string&)>& diagnosticObserver,
    const verification::ExactMaterialRegion* materialAuthority) {
    std::vector<Eigen::MatrixXd> reduced;
    reduced.reserve(moving.size());
    for(const auto& velocity:moving)reduced.emplace_back(velocity);
    const auto value=differentiateRestrictedPolyhedralMaterialMoments(
        design,source,reduced,material,birthFaces,density,exactBinary64SelfIntersection,
        exactBinary64SourceBoundary,exactBinary64MaterialBoundary,exactVertices,exactSource,exactDesign,
        diagnosticObserver,materialAuthority);
    QuadricDiscreteMomentDifferential out;
    out.defined=value.defined;out.diagnostic=value.diagnostic;
    out.vertexProvenance=value.vertexProvenance;
    out.constructedIntersectionVertexCount=value.constructedIntersectionVertexCount;
    out.maximumConstraintResidual=value.maximumConstraintResidual;
    out.minimumConstraintSingularRatio=value.minimumConstraintSingularRatio;
    if(value.jacobian.rows()==10 && value.jacobian.cols()==10)
        out.jacobian=value.jacobian;
    return out;
}

RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density) {
    return differentiateRestrictedPolyhedralMaterialMoments(design,source,moving,material,birthFaces,density,false);
}
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection) {
    return differentiateRestrictedPolyhedralMaterialMoments(design,source,moving,material,
        birthFaces,density,exactBinary64SelfIntersection,false);
}
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary) {
    return differentiateRestrictedPolyhedralMaterialMoments(design,source,moving,material,
        birthFaces,density,exactBinary64SelfIntersection,exactBinary64SourceBoundary,false);
}
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary) {
    return differentiateRestrictedPolyhedralMaterialMoments(design,source,moving,material,birthFaces,
        density,exactBinary64SelfIntersection,exactBinary64SourceBoundary,exactBinary64MaterialBoundary,{});
}
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary,const ExactBooleanVertexCoordinates& exactVertices) {
    return differentiateRestrictedPolyhedralMaterialMoments(design,source,moving,material,birthFaces,
        density,exactBinary64SelfIntersection,exactBinary64SourceBoundary,
        exactBinary64MaterialBoundary,exactVertices,{},{});
}
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary,const ExactBooleanVertexCoordinates& exactVertices,
    const ExactBooleanVertexCoordinates& exactSource, const ExactBooleanVertexCoordinates& exactDesign,
    const std::function<void(const std::string&)>& diagnosticObserver) {
    return differentiateRestrictedPolyhedralMaterialMoments(design,source,moving,material,birthFaces,
        density,exactBinary64SelfIntersection,exactBinary64SourceBoundary,exactBinary64MaterialBoundary,
        exactVertices,exactSource,exactDesign,diagnosticObserver,nullptr);
}
static RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMomentsImpl(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& moving,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary,const ExactBooleanVertexCoordinates& exactVertices,
    const ExactBooleanVertexCoordinates& exactSource, const ExactBooleanVertexCoordinates& exactDesign,
    const std::function<void(const std::string&)>& diagnosticObserver,
    const verification::ExactMaterialRegion* materialAuthority,
    const std::vector<RationalSourceVertexJacobian>* authorityMotion,
    RationalPolyhedralMomentDifferential* rationalResult) {
    RestrictedPolyhedralMomentDifferential out;
    const auto stage=[&](const std::string& name){if(diagnosticObserver)diagnosticObserver(name);};
    const auto fail=[&](const std::string& why) {
        out.defined=false; out.diagnostic=classifiedFailure(why);
        if(out.jacobian.size())out.jacobian.setZero();
        if(rationalResult) {
            rationalResult->defined=false;
            rationalResult->coefficientJacobian={};
            rationalResult->materialVertexJacobian.clear();
            rationalResult->maximumConstraintResidualExact.clear();
        }
        return out;
    };
    const auto validMesh=[exactBinary64SelfIntersection](const TriangleMesh& m, bool exactMaterial=false) {
        if(m.vertices.cols()!=3 || m.faces.cols()!=3 ||
           m.vertices.rows()<4 || m.faces.rows()<4 || !m.vertices.allFinite()) return false;
        if((m.faces.array()<0).any() || (m.faces.array()>=m.vertices.rows()).any()) return false;
        const auto audit=auditQuadricMaterialBoundary(m,exactBinary64SelfIntersection,exactMaterial);
        return audit.strictBoundaryValid() && audit.signedVolume>0;
    };
    const Eigen::Index dimensions=moving.empty()?0:moving.front().cols();
    if(!std::isfinite(density) || density<=0 || dimensions<=0 ||
       moving.size()!=static_cast<std::size_t>(source.vertices.rows()) ||
       birthFaces.size()!=material.faces.rows())
        return fail("invalid density, motion shape or provenance shape");
    for(const auto& velocity:moving)
        if(velocity.rows()!=3 || velocity.cols()!=dimensions || !velocity.allFinite())
            return fail("invalid or nonfinite source vertex motion");
#if !MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(exactBinary64MaterialBoundary)
        return fail("unsupported exact material boundary audit in this build");
#endif
    out.jacobian=Eigen::MatrixXd::Zero(10,dimensions);
    const bool rationalAuthority=!exactSource.empty() || !exactDesign.empty();
    if(authorityMotion && (!rationalAuthority || !rationalResult || dimensions!=10 ||
        authorityMotion->size()!=moving.size()))
        return fail("invalid rational motion without matching complete authority tapes");
    if(materialAuthority && !rationalAuthority)
        return fail("invalid material authority reuse without full rational tapes");
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    std::vector<ExactPoint> sourcePoints,designPoints;
    std::vector<ExactJacobian> rationalMoving,rationalFixed,rationalOutput;
    Rational maximumExactResidual(0);
    if(authorityMotion) {
        try {
            rationalMoving.resize(authorityMotion->size(),ExactJacobian(dimensions));
            for(std::size_t v=0;v<authorityMotion->size();++v)
                for(Eigen::Index j=0;j<dimensions;++j)for(int axis=0;axis<3;++axis)
                    rationalMoving[v][j][axis]=verification::parseExactRationalCoordinate(
                        authorityMotion->at(v)[axis][j]);
            rationalFixed.resize(design.vertices.rows(),ExactJacobian(dimensions));
            rationalOutput.resize(material.vertices.rows(),ExactJacobian(dimensions));
        } catch(const std::exception&) {return fail("invalid rational source velocity");}
    }
    if(rationalAuthority) {
        stage("transport-authority-audits-begin");
        const auto load=[&](const ExactBooleanVertexCoordinates& tape,const TriangleMesh& mesh,
                           std::vector<ExactPoint>* points) {
            if(tape.size()!=static_cast<std::size_t>(mesh.vertices.rows()))return false;
            points->resize(tape.size());
            for(std::size_t v=0;v<tape.size();++v)for(int axis=0;axis<3;++axis) {
                (*points)[v][axis]=verification::parseExactRationalCoordinate(tape[v][axis]);
                if(booleanRound((*points)[v][axis])!=mesh.vertices(v,axis))return false;
            }
            stage("transport-authority-region-begin:faces="+std::to_string(mesh.faces.rows()));
            const bool valid=verification::makeExactMaterialRegionFromRationals(tape,mesh.faces).facts().valid;
            stage("transport-authority-region-end:faces="+std::to_string(mesh.faces.rows()));
            return valid;
        };
        try {
            if(exactVertices.empty() || !load(exactSource,source,&sourcePoints) ||
               !load(exactDesign,design,&designPoints))
                return fail("invalid rational authority boundary or display assignment");
            if(materialAuthority) {
                if(!materialAuthority->matchesRationalBoundary(exactVertices,material.faces))
                    return fail("invalid rational material authority binding");
                stage("transport-material-authority-reused");
            } else {
                stage("transport-material-authority-audit-begin");
                const bool valid=verification::makeExactMaterialRegionFromRationals(
                    exactVertices,material.faces).facts().valid;
                stage("transport-material-authority-audit-end");
                if(!valid)return fail("invalid rational authority boundary or display assignment");
            }
        } catch(const std::exception&) {
            return fail("invalid rational authority coordinate or boundary");
        }
        stage("transport-authority-audits-end");
    }
#else
    if(rationalAuthority || authorityMotion)return fail("unsupported rational authority in this build");
#endif
    bool sourceValid=false;
    stage("transport-source-boundary-audit-begin");
    if(rationalAuthority)sourceValid=true;
    else if(exactBinary64SourceBoundary) {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        try {
            const auto proof=verification::auditExactTriangleBoundary(source.vertices,source.faces,
                verification::ExactBoundaryOrientation::OutwardPositive);
            sourceValid=proof.facts().boundaryValid;
        } catch(const std::invalid_argument&) {
            return fail("invalid exact source boundary input");
        }
#else
        return fail("unsupported exact source boundary audit in this build");
#endif
    } else sourceValid=validMesh(source);
    stage("transport-source-boundary-audit-end");
    bool designValid=true,materialValid=true;
    if(sourceValid && !rationalAuthority) {
        stage("transport-design-boundary-audit-begin");
        // Process-local immutable proof, bound to every indexed binary64
        // coordinate/face and the exact audit policy. Only successful proofs
        // are retained; moving source/material boundaries never use it.
        struct FixedDesignProof {
            std::string indexedBoundary;
        };
        static thread_local FixedDesignProof fixedProof;
        std::ostringstream key;key.imbue(std::locale::classic());
        key<<"FixedDesignPositiveBoundary-v1;"<<exactBinary64SelfIntersection<<';'
           <<exactBinary64SourceBoundary<<';'<<design.vertices.rows()<<';'
           <<design.vertices.cols()<<';'<<design.faces.rows()<<';'<<design.faces.cols()<<';'<<std::hexfloat;
        for(Eigen::Index row=0;row<design.vertices.rows();++row)
            for(Eigen::Index axis=0;axis<design.vertices.cols();++axis)
                key<<design.vertices(row,axis)<<';';
        for(Eigen::Index row=0;row<design.faces.rows();++row)
            for(Eigen::Index corner=0;corner<design.faces.cols();++corner)
                key<<design.faces(row,corner)<<';';
        const std::string indexedBoundary=key.str();
        if(fixedProof.indexedBoundary==indexedBoundary) {
            designValid=true;
            stage("transport-design-proof-cache-hit");
        } else {
            stage("transport-design-proof-cache-miss");
            designValid=validMesh(design,exactBinary64SourceBoundary);
            if(designValid)fixedProof.indexedBoundary=indexedBoundary;
        }
        stage("transport-design-boundary-audit-end");
        if(designValid) {
            stage("transport-material-boundary-audit-begin");
            materialValid=validMesh(material,exactBinary64MaterialBoundary);
            stage("transport-material-boundary-audit-end");
        }
    }
    if(!sourceValid || !designValid || !materialValid)
        return fail("input is not a finite positive closed audited boundary");
    if((birthFaces.array()<0).any() ||
       (birthFaces.array()>=design.faces.rows()+source.faces.rows()).any())
        return fail("invalid Boolean birth facet");
    if(!exactVertices.empty() && exactVertices.size()!=static_cast<std::size_t>(material.vertices.rows()))
        return fail("invalid exact Boolean coordinate shape");
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    std::vector<ExactPoint> exactPoints(exactVertices.size());
    try {
        for(std::size_t v=0;v<exactVertices.size();++v)for(int axis=0;axis<3;++axis) {
            exactPoints[v][axis]=verification::parseExactRationalCoordinate(exactVertices[v][axis]);
            if(booleanRound(exactPoints[v][axis])!=material.vertices(v,axis))
                return fail("invalid exact Boolean coordinate differs from actual material assignment");
        }
    } catch(const std::invalid_argument&) {
        return fail("invalid exact Boolean rational coordinate");
    }
#else
    if(!exactVertices.empty())return fail("unsupported exact Boolean motion coordinates in this build");
#endif
    std::vector<Jacobian> fixed(static_cast<std::size_t>(design.vertices.rows()),
        Jacobian::Zero(3,dimensions));
    std::vector<std::set<int>> incident(material.vertices.rows());
    for(int f=0;f<material.faces.rows();++f) {
        if(birthFaces(f)<0 || birthFaces(f)>=design.faces.rows()+source.faces.rows())
            return fail("invalid Boolean birth facet");
        for(int k=0;k<3;++k) incident[material.faces(f,k)].insert(birthFaces(f));
    }
    std::vector<Jacobian> outputVelocity(static_cast<std::size_t>(material.vertices.rows()),
        Jacobian::Zero(3,dimensions));
    stage("transport-feature-assembly-begin");
    for(int v=0;v<material.vertices.rows();++v) {
        const Eigen::Vector3d x=material.vertices.row(v);
        std::vector<Eigen::Vector3d> lhs;
        std::vector<Row> rhs;
        std::array<Feature,2> feature;
        Jacobian sourceFeatureMotion=Jacobian::Zero(3,dimensions);
        bool hasSourceFeatureMotion=false;
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        std::vector<ExactPoint> rationalLhs;
        std::vector<std::vector<Rational>> rationalRhs;
        ExactJacobian rationalSourceFeatureMotion;
#endif
        for(int birth:incident[v]) {
            const bool isMoving=birth>=design.faces.rows();
            const auto& mesh=isMoving?source:design;
            const int f=isMoving?birth-design.faces.rows():birth;
            Jacobian candidate;
            bool precise=false;
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
            if(!exactPoints.empty()) {
                precise=true;
                ExactJacobian rationalCandidate;
                std::string membershipFailure;
                if(!addExactFacetConstraints(mesh,f,isMoving?moving:fixed,exactPoints[v],
                    isMoving?"s:":"d:",&feature[isMoving],&lhs,&rhs,&candidate,
                    rationalAuthority?(isMoving?&sourcePoints:&designPoints):nullptr,
                    authorityMotion?(isMoving?&rationalMoving:&rationalFixed):nullptr,
                    authorityMotion?&rationalLhs:nullptr,authorityMotion?&rationalRhs:nullptr,
                    authorityMotion?&rationalCandidate:nullptr,&membershipFailure,!isMoving)) {
                    std::ostringstream detail;
                    detail<<"exact Boolean birth facet membership is ill-conditioned or inconsistent"
                          <<" [reason="<<membershipFailure<<", material_vertex="<<v
                          <<", operand="<<(isMoving?"source":"design")<<", face="<<f
                          <<", indices="<<mesh.faces(f,0)<<","<<mesh.faces(f,1)
                          <<","<<mesh.faces(f,2)<<"]";
                    return fail(detail.str());
                }
                if(authorityMotion && isMoving) {
                    if(hasSourceFeatureMotion && rationalCandidate!=rationalSourceFeatureMotion)
                        return fail("inconsistent rational source-feature interpolation at Boolean vertex");
                    rationalSourceFeatureMotion=std::move(rationalCandidate);
                }
            }
#endif
            if(!precise) {
                Eigen::Vector3d bary;
                if(!barycentric(mesh,f,x,&bary)) return fail("Boolean birth facet membership is ill-conditioned or inconsistent");
                addFacetConstraints(mesh,f,isMoving?moving:fixed,x,bary,
                                    isMoving?"s:":"d:",&feature[isMoving],&lhs,&rhs);
                const auto& ids=mesh.faces.row(f);
                candidate=bary(0)*(isMoving?moving:fixed)[static_cast<std::size_t>(ids(0))]+
                    bary(1)*(isMoving?moving:fixed)[static_cast<std::size_t>(ids(1))]+
                    bary(2)*(isMoving?moving:fixed)[static_cast<std::size_t>(ids(2))];
            }
            if(isMoving) {
                if(hasSourceFeatureMotion &&
                   (candidate-sourceFeatureMotion).norm()>consistencyTolerance*
                       std::max(1.0,candidate.norm()))
                    return fail("inconsistent source-feature interpolation at Boolean vertex");
                sourceFeatureMotion=candidate;
                hasSourceFeatureMotion=true;
            }
        }
        for(int side=0;side<2;++side) {
            if(feature[side].key=="ambiguous") return fail("ambiguous Boolean source feature");
            if(feature[side].dimension==0) {
                // Exact input contacts may have only one surviving birth side.
                const auto& other=side==0?source:design;
                for(int f=0;f<other.faces.rows();++f) {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
                    if(rationalAuthority) {
                        const auto& points=side==0?sourcePoints:designPoints;
                        const auto a=points[other.faces(f,0)],u=subtract(points[other.faces(f,1)],a),
                            w=subtract(points[other.faces(f,2)],a),n=cross(u,w),d=subtract(exactPoints[v],a);
                        const auto s=dot(n,n);
                        if(s!=0 && dot(n,d)==0) {
                            const auto b=dot(cross(d,w),n)/s,c=dot(cross(u,d),n)/s;
                            if(b>=0 && c>=0 && b+c<=1)
                                return fail("source/design input-vertex contact is unsupported");
                        }
                        continue;
                    }
#endif
                    Eigen::Vector3d bary;
                    if(barycentric(other,f,x,&bary))
                        return fail("source/design input-vertex contact is unsupported");
                }
            }
        }
        // Feature dimension 3 is the sentinel for "no source feature"; 0 is
        // a valid input vertex feature. A Boolean vertex supported by only
        // one input mesh is transported directly by that feature instead of
        // being mistaken for a newly constructed two-input intersection.
        const bool hasDesignFeature=feature[0].dimension<3;
        const bool hasMovingFeature=feature[1].dimension<3;
        const bool constructed=hasDesignFeature && hasMovingFeature;
        if(constructed) {
            if(feature[0].dimension==1 && feature[1].dimension==1)
                return fail("source/design edge-edge contact is unsupported");
            ++out.constructedIntersectionVertexCount;
        }
        if(!constructed) {
            if(hasMovingFeature) {
                if(!hasSourceFeatureMotion)
                    return fail("moving source feature has no interpolation witness");
                outputVelocity[v]=sourceFeatureMotion;
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
                if(authorityMotion)rationalOutput[v]=rationalSourceFeatureMotion;
#endif
            }
            out.vertexProvenance.push_back(feature[0].key+"|"+feature[1].key);
            continue;
        }
        Eigen::MatrixXd a(lhs.size(),3),b(rhs.size(),dimensions);
        for(int row=0;row<a.rows();++row) { a.row(row)=lhs[row]; b.row(row)=rhs[row]; }
        if(!a.allFinite() || !b.allFinite()) return fail("nonfinite Boolean motion constraints");
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);
        if(svd.singularValues().size()!=3 || svd.singularValues()(0)<=0)
            return fail("Boolean motion constraints are rank deficient");
        const double ratio=svd.singularValues()(2)/svd.singularValues()(0);
        out.minimumConstraintSingularRatio=std::min(out.minimumConstraintSingularRatio,ratio);
        if(ratio<minimumSingularRatio) return fail("Boolean motion constraints are rank deficient or ill-conditioned");
        outputVelocity[v]=svd.solve(b);
        const double residual=(a*outputVelocity[v]-b).norm()/std::max(1.0,b.norm());
        out.maximumConstraintResidual=std::max(out.maximumConstraintResidual,residual);
        if(!outputVelocity[v].allFinite() || !std::isfinite(residual) || residual>consistencyTolerance)
            return fail("Boolean motion constraints have inconsistent differential");
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        if(authorityMotion) {
            // Preserve the protected numerical conditioning gate, then solve
            // the SAME constraints without rounding the authority velocities.
            if(rationalLhs.size()!=lhs.size() || rationalRhs.size()!=rhs.size())
                return fail("invalid rational Boolean motion constraint assembly");
            std::vector<std::pair<int,ExactPoint>> basis;
            std::vector<int> selected;
            for(std::size_t row=0;row<rationalLhs.size() && basis.size()<3;++row) {
                auto reduced=rationalLhs[row];
                for(const auto& pivot:basis) {
                    const Rational factor=reduced[pivot.first];
                    for(int axis=0;axis<3;++axis)reduced[axis]-=factor*pivot.second[axis];
                }
                int column=0;while(column<3 && reduced[column]==0)++column;
                if(column==3)continue;
                const Rational scale=reduced[column];
                for(auto& value:reduced)value/=scale;
                basis.push_back({column,reduced});selected.push_back(static_cast<int>(row));
            }
            if(selected.size()!=3)return fail("rational Boolean motion constraints are rank deficient");
            const std::array<ExactPoint,3> rows{{rationalLhs[selected[0]],
                rationalLhs[selected[1]],rationalLhs[selected[2]]}};
            const auto determinant=[](const std::array<ExactPoint,3>& matrix) {
                return dot(matrix[0],cross(matrix[1],matrix[2]));
            };
            const Rational det=determinant(rows);
            if(det==0)return fail("rational Boolean motion constraints are rank deficient");
            for(Eigen::Index j=0;j<dimensions;++j) {
                for(int axis=0;axis<3;++axis) {
                    auto replaced=rows;
                    for(int row=0;row<3;++row)replaced[row][axis]=rationalRhs[selected[row]][j];
                    rationalOutput[v][j][axis]=determinant(replaced)/det;
                    outputVelocity[v](axis,j)=booleanRound(rationalOutput[v][j][axis]);
                }
                for(std::size_t row=0;row<rationalLhs.size();++row) {
                    const Rational difference=CGAL::abs(dot(rationalLhs[row],rationalOutput[v][j])-rationalRhs[row][j]);
                    const Rational normalized=difference/std::max(Rational(1),CGAL::abs(rationalRhs[row][j]));
                    maximumExactResidual=std::max(maximumExactResidual,normalized);
                    if(normalized>Rational(consistencyTolerance))
                        return fail("rational Boolean motion constraints have inconsistent differential");
                }
            }
        }
#endif
        out.vertexProvenance.push_back(feature[0].key+"|"+feature[1].key);
    }
    stage("transport-feature-assembly-end");
    stage("transport-moment-integration-begin");
    // Reynolds transport on the actual polyhedron. A degree-two monomial
    // times affine vertex velocity is degree three. This four-point rule is
    // exactly cubic: weights sum to one, area vector carries the orientation.
    const std::array<Eigen::Vector3d,4> bary={{Eigen::Vector3d(1./3,1./3,1./3),
        Eigen::Vector3d(.6,.2,.2),Eigen::Vector3d(.2,.6,.2),Eigen::Vector3d(.2,.2,.6)}};
    const std::array<double,4> weights={{-27./48,25./48,25./48,25./48}};
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(rationalAuthority) {
        const std::array<std::array<Rational,3>,4> exactBary{{
            {{Rational(1,3),Rational(1,3),Rational(1,3)}},
            {{Rational(3,5),Rational(1,5),Rational(1,5)}},
            {{Rational(1,5),Rational(3,5),Rational(1,5)}},
            {{Rational(1,5),Rational(1,5),Rational(3,5)}}}};
        const std::array<Rational,4> exactWeights{{Rational(-27,48),Rational(25,48),Rational(25,48),Rational(25,48)}};
        using MomentBlock=std::vector<std::array<Rational,10>>;
        // Exact addition is associative. Accumulate complete face blocks in a
        // binary carry tree instead of multiplying each new small fraction
        // into an ever-growing global denominator. No rounding, quadrature,
        // source velocity, facet or validity predicate is changed.
        std::vector<MomentBlock> partials;
        const Rational exactDensity(density);
        for(int f=0;f<material.faces.rows();++f) {
            MomentBlock faceSums(dimensions);
            const std::array<int,3> ids{{material.faces(f,0),material.faces(f,1),material.faces(f,2)}};
            auto area=cross(subtract(exactPoints[ids[1]],exactPoints[ids[0]]),
                            subtract(exactPoints[ids[2]],exactPoints[ids[0]]));
            for(auto& value:area)value/=2;
            for(int s=0;s<4;++s) {
                ExactPoint p{};
                for(int axis=0;axis<3;++axis)for(int k=0;k<3;++k)
                    p[axis]+=exactBary[s][k]*exactPoints[ids[k]][axis];
                const std::array<Rational,10> terms{{Rational(1),p[0],p[1],p[2],p[0]*p[0],
                    p[1]*p[1],p[2]*p[2],p[1]*p[2],p[0]*p[2],p[0]*p[1]}};
                for(Eigen::Index j=0;j<dimensions;++j) {
                    Rational flux(0);
                    for(int axis=0;axis<3;++axis)for(int k=0;k<3;++k)
                        flux+=area[axis]*exactBary[s][k]*(authorityMotion?
                            rationalOutput[ids[k]][j][axis]:Rational(outputVelocity[ids[k]](axis,j)));
                    const Rational weightedFlux=exactDensity*exactWeights[s]*flux;
                    if(weightedFlux!=0)
                        for(int m=0;m<10;++m)faceSums[j][m]+=terms[m]*weightedFlux;
                }
            }
            std::size_t level=0;
            while(level<partials.size() && !partials[level].empty()) {
                for(Eigen::Index j=0;j<dimensions;++j)for(int m=0;m<10;++m)
                    faceSums[j][m]+=partials[level][j][m];
                partials[level].clear();++level;
            }
            if(level==partials.size())partials.emplace_back();
            partials[level]=std::move(faceSums);
        }
        MomentBlock sums(dimensions);
        for(const auto& block:partials)if(!block.empty())
            for(Eigen::Index j=0;j<dimensions;++j)for(int m=0;m<10;++m)
                sums[j][m]+=block[j][m];
        for(Eigen::Index j=0;j<dimensions;++j)for(int m=0;m<10;++m) {
            out.jacobian(m,j)=booleanRound(sums[j][m]);
            if(authorityMotion) {
                std::ostringstream text;text<<sums[j][m];
                rationalResult->coefficientJacobian[m][j]=text.str();
            }
        }
        if(authorityMotion) {
            rationalResult->materialVertexJacobian.resize(rationalOutput.size());
            for(std::size_t v=0;v<rationalOutput.size();++v)
                for(Eigen::Index j=0;j<dimensions;++j)for(int axis=0;axis<3;++axis) {
                    std::ostringstream text;text<<rationalOutput[v][j][axis];
                    rationalResult->materialVertexJacobian[v][axis][j]=text.str();
                }
            std::ostringstream text;text<<maximumExactResidual;
            rationalResult->maximumConstraintResidualExact=text.str();
        }
    } else
#endif
    for(int f=0;f<material.faces.rows();++f) {
        const int ia=material.faces(f,0),ib=material.faces(f,1),ic=material.faces(f,2);
        const Eigen::Vector3d a=material.vertices.row(ia),b=material.vertices.row(ib),c=material.vertices.row(ic);
        const Eigen::Vector3d area=.5*(b-a).cross(c-a);
        for(int s=0;s<4;++s) {
            const auto& w=bary[s];
            const Jacobian velocity=w(0)*outputVelocity[ia]+w(1)*outputVelocity[ib]+w(2)*outputVelocity[ic];
            out.jacobian += density*weights[s]*monomials(w(0)*a+w(1)*b+w(2)*c)*(area.transpose()*velocity);
        }
    }
    stage("transport-moment-integration-end");
    if(!out.jacobian.allFinite()) return fail("nonfinite discrete moment differential");
    out.defined=true;
    out.diagnostic="VALID_FIXED_BOOLEAN_TOPOLOGY: supplied source motion with fixed Boolean source-feature topology";
    return out;
}

RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design,const TriangleMesh& source,const std::vector<Eigen::MatrixXd>& moving,
    const TriangleMesh& material,const Eigen::VectorXi& birthFaces,double density,
    bool exactBinary64SelfIntersection,bool exactBinary64SourceBoundary,bool exactBinary64MaterialBoundary,
    const ExactBooleanVertexCoordinates& exactVertices,const ExactBooleanVertexCoordinates& exactSource,
    const ExactBooleanVertexCoordinates& exactDesign,const std::function<void(const std::string&)>& observer,
    const verification::ExactMaterialRegion* materialAuthority) {
    return differentiateRestrictedPolyhedralMaterialMomentsImpl(design,source,moving,material,birthFaces,density,
        exactBinary64SelfIntersection,exactBinary64SourceBoundary,exactBinary64MaterialBoundary,
        exactVertices,exactSource,exactDesign,observer,materialAuthority,nullptr,nullptr);
}

RationalPolyhedralMomentDifferential differentiateRationalPolyhedralMaterialMoments(
    const TriangleMesh& design,const TriangleMesh& source,
    const std::vector<RationalSourceVertexJacobian>& sourceVelocity,const TriangleMesh& material,
    const Eigen::VectorXi& birthFaces,double density,const ExactBooleanVertexCoordinates& exactMaterial,
    const ExactBooleanVertexCoordinates& exactSource,const ExactBooleanVertexCoordinates& exactDesign,
    const std::function<void(const std::string&)>& observer,
    const verification::ExactMaterialRegion* materialAuthority) {
    RationalPolyhedralMomentDifferential out;
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    try {
        std::vector<Eigen::MatrixXd> display(sourceVelocity.size(),Eigen::MatrixXd(3,10));
        for(std::size_t v=0;v<sourceVelocity.size();++v)for(int axis=0;axis<3;++axis)for(int j=0;j<10;++j)
            display[v](axis,j)=booleanRound(verification::parseExactRationalCoordinate(sourceVelocity[v][axis][j]));
        const auto value=differentiateRestrictedPolyhedralMaterialMomentsImpl(
            design,source,display,material,birthFaces,density,true,true,true,
            exactMaterial,exactSource,exactDesign,observer,materialAuthority,&sourceVelocity,&out);
        out.defined=value.defined;out.diagnostic=value.diagnostic;
        out.summary.defined=value.defined;out.summary.diagnostic=value.diagnostic;
        out.summary.vertexProvenance=value.vertexProvenance;
        out.summary.constructedIntersectionVertexCount=value.constructedIntersectionVertexCount;
        out.summary.maximumConstraintResidual=value.maximumConstraintResidual;
        out.summary.minimumConstraintSingularRatio=value.minimumConstraintSingularRatio;
        if(value.jacobian.rows()==10 && value.jacobian.cols()==10)out.summary.jacobian=value.jacobian;
    } catch(const std::exception& error) {
        out=RationalPolyhedralMomentDifferential{};
        out.diagnostic=classifiedFailure(std::string("invalid rational motion: ")+error.what());
        out.summary.diagnostic=out.diagnostic;
    }
#else
    out.diagnostic="UNSUPPORTED: rational motion requires exact geometry build";
    out.summary.diagnostic=out.diagnostic;
#endif
    return out;
}
}
