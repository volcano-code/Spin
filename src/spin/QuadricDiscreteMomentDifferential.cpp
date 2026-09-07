#include "spin/QuadricDiscreteMomentDifferential.h"
#include "spin/QuadricInverseGaussMap.h"
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <sstream>

namespace {
using Jacobian = Eigen::Matrix<double,3,10>;
using Row = Eigen::Matrix<double,1,10>;
constexpr double membershipTolerance = 2e-9;
constexpr double minimumSingularRatio = 1e-10;
constexpr double consistencyTolerance = 2e-9;

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
    Row plane;
    for(int j=0;j<10;++j) {
        const Eigen::Vector3d dn=(velocity[ids[1]].col(j)-velocity[ids[0]].col(j)).cross(c-a)+
            (b-a).cross(velocity[ids[2]].col(j)-velocity[ids[0]].col(j));
        plane(j)=(n.dot(velocity[ids[0]].col(j))-dn.dot(x-a))/n.norm();
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

Eigen::Matrix<double,10,1> monomials(const Eigen::Vector3d& p) {
    Eigen::Matrix<double,10,1> m;
    m << 1,p.x(),p.y(),p.z(),p.x()*p.x(),p.y()*p.y(),p.z()*p.z(),
         p.y()*p.z(),p.x()*p.z(),p.x()*p.y();
    return m;
}
}

namespace spin {
QuadricDiscreteMomentDifferential differentiateDiscreteEllipsoidMoments(
    const TriangleMesh& design, const QuadricTriangulation& triangulation,
    const QuadricCoefficients10& quadric, const TriangleMesh& material,
    const Eigen::VectorXi& birthFaces, double density) {
    QuadricDiscreteMomentDifferential out;
    const auto fail=[&](const std::string& why) {
        out.defined=false; out.diagnostic=why; out.jacobian.setZero(); return out;
    };
    const auto& source=triangulation.closedBoundary;
    if(triangulation.type!=QuadricSurfaceType::Ellipsoid)
        return fail("discrete differential supports finite ellipsoids only; caps and other types unsupported");
    if(triangulation.fixedNormals.rows()!=source.vertices.rows() ||
       triangulation.fixedNormals.cols()!=3 || !triangulation.fixedNormals.allFinite())
        return fail("missing or nonfinite fixed-normal point provenance");
    if(!design.vertices.allFinite() || !source.vertices.allFinite() ||
       !material.vertices.allFinite() || !std::isfinite(density) || density<=0 ||
       birthFaces.size()!=material.faces.rows()) return fail("invalid discrete geometry/provenance data");
    std::vector<Jacobian> fixed(design.vertices.rows(),Jacobian::Zero());
    std::vector<Jacobian> moving(source.vertices.rows());
    QuadricInverseGaussMap map(quadric);
    for(int v=0;v<source.vertices.rows();++v) {
        const auto differential=map.coefficientDifferential(triangulation.fixedNormals.row(v));
        if(!differential.valid()) return fail("inverse-Gauss differential: "+differential.diagnostic);
        if((*differential.point.point-source.vertices.row(v).transpose()).norm()>
           1e-12*std::max(1.0,source.vertices.row(v).norm()))
            return fail("fixed-normal point does not match actual triangulation vertex");
        moving[v]=*differential.coefficientJacobian;
    }
    std::vector<std::set<int>> incident(material.vertices.rows());
    for(int f=0;f<material.faces.rows();++f) {
        if(birthFaces(f)<0 || birthFaces(f)>=design.faces.rows()+source.faces.rows())
            return fail("invalid Boolean birth facet");
        for(int k=0;k<3;++k) incident[material.faces(f,k)].insert(birthFaces(f));
    }
    std::vector<Jacobian> outputVelocity(material.vertices.rows());
    for(int v=0;v<material.vertices.rows();++v) {
        const Eigen::Vector3d x=material.vertices.row(v);
        std::vector<Eigen::Vector3d> lhs;
        std::vector<Row> rhs;
        std::array<Feature,2> feature;
        for(int birth:incident[v]) {
            const bool isMoving=birth>=design.faces.rows();
            const auto& mesh=isMoving?source:design;
            const int f=isMoving?birth-design.faces.rows():birth;
            Eigen::Vector3d bary;
            if(!barycentric(mesh,f,x,&bary)) return fail("Boolean birth facet membership is ill-conditioned or inconsistent");
            addFacetConstraints(mesh,f,isMoving?moving:fixed,x,bary,
                                isMoving?"q:":"d:",&feature[isMoving],&lhs,&rhs);
        }
        for(int side=0;side<2;++side) {
            if(feature[side].key=="ambiguous") return fail("ambiguous Boolean source feature");
            if(feature[side].dimension==0) {
                // Exact input contacts may have only one surviving birth side.
                const auto& other=side==0?source:design;
                for(int f=0;f<other.faces.rows();++f) {
                    Eigen::Vector3d bary;
                    if(barycentric(other,f,x,&bary))
                        return fail("quadric/design input-vertex contact is unsupported");
                }
            }
        }
        const bool constructed=feature[0].dimension!=0 && feature[1].dimension!=0;
        if(constructed) {
            if(feature[0].dimension==3 || feature[1].dimension==3)
                return fail("constructed Boolean vertex lacks both source features");
            ++out.constructedIntersectionVertexCount;
        }
        Eigen::MatrixXd a(lhs.size(),3),b(rhs.size(),10);
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
        out.vertexProvenance.push_back(feature[0].key+"|"+feature[1].key);
    }
    // Reynolds transport on the actual polyhedron. A degree-two monomial
    // times affine vertex velocity is degree three. This four-point rule is
    // exactly cubic: weights sum to one, area vector carries the orientation.
    const std::array<Eigen::Vector3d,4> bary={{Eigen::Vector3d(1./3,1./3,1./3),
        Eigen::Vector3d(.6,.2,.2),Eigen::Vector3d(.2,.6,.2),Eigen::Vector3d(.2,.2,.6)}};
    const std::array<double,4> weights={{-27./48,25./48,25./48,25./48}};
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
    if(!out.jacobian.allFinite()) return fail("nonfinite discrete moment differential");
    out.defined=true;
    out.diagnostic="supported finite ellipsoid with fixed source-feature topology";
    return out;
}
}
