#include "spin/QuadricSourceMotion.h"
#include "spin/QuadricInverseGaussMap.h"
#include "spin/QuadricStratumTangent.h"
#include "spin/TruncationRootCertificate.h"
#include "manufacturing/SlicerValidation.h"
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>

namespace spin {

namespace {
Eigen::Matrix3d quadraticDirection(
    const QuadricCoefficients10::Vector& direction) {
    Eigen::Matrix3d out=Eigen::Matrix3d::Zero();
    out.diagonal()<<direction(4),direction(5),direction(6);
    out(1,2)=out(2,1)=direction(7);
    out(0,2)=out(2,0)=direction(8);
    out(0,1)=out(1,0)=direction(9);
    return out;
}

Eigen::Matrix3d canonicalEigenvectorSigns(const Eigen::Matrix3d& input) {
    Eigen::Matrix3d result=input;
    for(int column=0;column<3;++column) {
        Eigen::Index pivot=0;
        result.col(column).cwiseAbs().maxCoeff(&pivot);
        if(result(pivot,column)<0.0)result.col(column)=-result.col(column);
    }
    return result;
}

bool eigenvectorSignPivotsAreSeparated(const Eigen::Matrix3d& vectors) {
    constexpr double minimumPivotGap=1e-10;
    for(int column=0;column<3;++column) {
        std::array<std::pair<double,int>,3> magnitudes{};
        for(int row=0;row<3;++row)
            magnitudes[static_cast<std::size_t>(row)]=
                {std::abs(vectors(row,column)),row};
        std::sort(magnitudes.begin(),magnitudes.end(),
            [](const auto& lhs,const auto& rhs) {return lhs.first>rhs.first;});
        const bool oppositeSigns=vectors(magnitudes[0].second,column)*
            vectors(magnitudes[1].second,column)<0.0;
        if(oppositeSigns && magnitudes[0].first-magnitudes[1].first<=minimumPivotGap)
            return false;
    }
    return true;
}

QuadricSourceMotion paperProjectiveMotion(const QuadricCoefficients10& q,
                                        QuadricTriangulation triangulation,
                                        bool algebraicFiniteMap) {
    QuadricSourceMotion out;out.triangulation=std::move(triangulation);
    const auto fail=[&](const std::string& message) {
        out.defined=false;out.velocity.clear();out.provenance.clear();
        out.scientificHash.clear();out.diagnostic=message;return out;
    };
    const auto& t=out.triangulation;
    const auto& v=t.closedBoundary.vertices;
    if(t.sourceProjectiveRayByVertex.size()!=std::size_t(v.rows()) ||
       t.sourceNormalVertexByVertex.size()!=std::size_t(v.rows()) || t.fixedNormals.cols()!=3)
        return fail("paper projective motion lacks complete native provenance");
    QuadricInverseGaussMap map(q);
    std::map<int,QuadricInverseGaussDifferentialResult> finite;
    std::map<std::pair<int,int>,QuadricIdealEdgeResult> ideal;
    std::map<int,std::array<int,4>> movingNormals;
    for(int i=0;i<t.interfaceVertexCount;++i)if(t.sourceProjectiveRayByVertex[i][0]==-3)
        movingNormals.emplace(t.sourceNormalVertexByVertex[i],t.sourceProjectiveRayByVertex[i]);
    const auto finitePoint=[&](int id)->const QuadricInverseGaussDifferentialResult& {
        if(id<0 || id>=t.fixedNormals.rows())throw std::runtime_error("invalid paper finite-normal identity");
        auto found=finite.find(id);
        if(found==finite.end()) {
            auto differential=algebraicFiniteMap
            ? map.algebraicCoefficientDifferential(t.fixedNormals.row(id))
            : map.coefficientDifferential(t.fixedNormals.row(id));
            if(movingNormals.count(id) && differential.valid()) {
                const auto recipe=movingNormals.at(id);const int a=recipe[1],b=recipe[2];
                if(a<0 || b<=a || b>=t.fixedNormals.rows())throw std::runtime_error("invalid moving interval normal identity");
                const Eigen::Vector3d first=t.fixedNormals.row(a),delta=t.fixedNormals.row(b)-t.fixedNormals.row(a);
                const auto roots=map.idealEdgeDifferential(first,t.fixedNormals.row(b));
                if(!roots.defined || roots.points.size()!=2)throw std::runtime_error("moving interval roots unavailable");
                const double mid=(roots.points[0].parameter+roots.points[1].parameter)*0.5;
                const Eigen::Vector3d normal=first+mid*delta;
                if((normal-t.fixedNormals.row(id).transpose()).norm()>1e-12)
                    throw std::runtime_error("stale moving interval normal");
                Eigen::Matrix<double,1,10> dt=Eigen::Matrix<double,1,10>::Zero();
                for(const auto& root:roots.points) {
                    const double slope=2*delta.dot(root.direction);
                    if(!std::isfinite(slope) || slope==0)throw std::runtime_error("singular moving interval root");
                    for(int j=4;j<10;++j) {
                        QuadricCoefficients10::Vector unit=QuadricCoefficients10::Vector::Zero();unit(j)=1;
                        dt(j)+=0.5*root.direction.dot(quadraticDirection(unit)*root.direction)/slope;
                    }
                }
                dt=dt*(Eigen::Matrix<double,10,10>::Identity()-q.values()*q.values().transpose());
                const auto& inverse=map.inverseQuadraticMatrix();const Eigen::Vector3d u=inverse*normal;
                const double s=normal.dot(u);
                if(!std::isfinite(s) || s==0 || u.squaredNorm()==0)throw std::runtime_error("singular moving finite normal");
                const Eigen::Vector3d center=-0.5*inverse*q.linearVector();
                const double factor=(*differential.point.point-center).dot(u)/u.squaredNorm();
                *differential.coefficientJacobian+=factor*(inverse-u*u.transpose()/s)*delta*dt;
            }
            found=finite.emplace(id,std::move(differential)).first;
        }
        return found->second;
    };
    out.velocity.resize(v.rows(),SourceVertexJacobian::Zero());out.provenance.resize(v.rows());
    for(int i=0;i<t.interfaceVertexCount;++i) {
        const int normal=t.sourceNormalVertexByVertex[i];
        const auto id=t.sourceProjectiveRayByVertex[i];
        Eigen::Vector3d expected;
        if(normal>=0) {
            if(id!=std::array<int,4>{-1,-1,-1,-1} && id[0]!=-3)return fail("finite paper point has a ray identity");
            const auto& f=finitePoint(normal);
            if(!f.valid())return fail(f.diagnostic);
            expected=*f.point.point;out.velocity[i]=*f.coefficientJacobian;
            out.provenance[i]=(id[0]==-3?"paper-interval-normal:":"paper-normal:")+std::to_string(normal);
        } else {
            if(id[0]<0 || id[1]<0 || id[2]<=id[1] || id[2]>=t.fixedNormals.rows() || id[3]<0)
                return fail("invalid paper finite/ideal ray identity");
            const auto& f=finitePoint(id[0]);if(!f.valid())return fail(f.diagnostic);
            const auto edge=std::make_pair(id[1],id[2]);
            auto found=ideal.find(edge);
            if(found==ideal.end())found=ideal.emplace(edge,map.idealEdgeDifferential(
                t.fixedNormals.row(id[1]),t.fixedNormals.row(id[2]))).first;
            const auto& roots=found->second;
            if(!roots.defined || id[3]>=int(roots.points.size()))return fail("paper ideal-root branch unavailable");
            const auto& point=roots.points[id[3]];
            const auto ray=truncateQuadricProjectiveRay(*f.point.point,point.direction,t.truncationRadius);
            if(!ray.defined)return fail(ray.diagnostic);
            expected=ray.point;
            out.velocity[i]=ray.jacobian.leftCols<3>()*(*f.coefficientJacobian)+
                            ray.jacobian.middleCols<3>(3)*point.coefficientJacobian;
            out.provenance[i]="paper-ray:"+std::to_string(id[0])+":"+
                std::to_string(id[1])+","+std::to_string(id[2])+":"+std::to_string(id[3]);
        }
        if(!out.velocity[i].allFinite() ||
           (expected-v.row(i).transpose()).norm()>1e-12*std::max(1.0,expected.norm()))
            return fail("paper source motion does not match actual native construction");
    }
    std::map<int,std::map<int,int>> extrusionLoops;
    for(int i=t.interfaceVertexCount;i<v.rows();++i) {
        const auto id=t.sourceProjectiveRayByVertex[i];
        if(id[0]==-2) {
            if(id[1]<0 || id[1]>=t.interfaceVertexCount || id[2]<0 || id[3]<0 ||
               !extrusionLoops[id[2]].emplace(id[3],id[1]).second)
                return fail("invalid paper cap extrusion identity");
        }
    }
    std::map<int,std::pair<Eigen::Vector3d,SourceVertexJacobian>> extrusionAxes;
    for(const auto& entry:extrusionLoops) {
        std::vector<int> loop;int ordinal=0;
        for(const auto& point:entry.second) {
            if(point.first!=ordinal++)return fail("paper extrusion loop ordinals are incomplete");
            loop.push_back(point.second);
        }
        if(loop.size()<3)return fail("paper extrusion loop has fewer than three points");
        Eigen::Vector3d area=Eigen::Vector3d::Zero(),mean=Eigen::Vector3d::Zero();
        SourceVertexJacobian areaMotion=SourceVertexJacobian::Zero();
        for(std::size_t k=0;k<loop.size();++k) {
            const int a=loop[k],b=loop[(k+1)%loop.size()];
            const Eigen::Vector3d x=v.row(a),y=v.row(b);
            area+=x.cross(y);mean+=x;
            for(int j=0;j<10;++j)areaMotion.col(j)+=out.velocity[a].col(j).cross(y)+x.cross(out.velocity[b].col(j));
        }
        const double norm=area.norm();
        if(!std::isfinite(norm) || norm<=0)return fail("paper extrusion axis has zero area");
        Eigen::Vector3d axis=area/norm;
        if(std::abs(axis.dot(mean))<=1e-12*std::max(1.0,mean.norm()))
            return fail("paper extrusion axis orientation is not separated");
        const double sign=axis.dot(mean)>0?1:-1;
        const SourceVertexJacobian motion=sign*(Eigen::Matrix3d::Identity()-axis*axis.transpose())*areaMotion/norm;
        axis*=sign;extrusionAxes.emplace(entry.first,std::make_pair(axis,motion));
    }
    for(int i=t.interfaceVertexCount;i<v.rows();++i) {
        const auto id=t.sourceProjectiveRayByVertex[i];
        if(id[0]==-2) {
            if(t.sourceNormalVertexByVertex[i]!=-1)return fail("paper extrusion point has finite normal identity");
            const auto& construction=extrusionAxes.at(id[2]);
            const auto& axis=construction.first;const auto& axisMotion=construction.second;
            const Eigen::Vector3d point=v.row(id[1]);
            const double offset=2*t.truncationRadius-axis.dot(point);
            const Eigen::Vector3d expected=point+axis*offset;
            out.velocity[i]=out.velocity[id[1]]+offset*axisMotion-
                axis*(point.transpose()*axisMotion+axis.transpose()*out.velocity[id[1]]);
            if(!out.velocity[i].allFinite() ||
               (expected-v.row(i).transpose()).norm()>1e-12*std::max(1.0,expected.norm()))
                return fail("paper extrusion motion differs from actual cap construction");
            out.provenance[i]="paper-extruded-cap-v1:"+std::to_string(id[1])+":"+
                std::to_string(id[2])+":"+std::to_string(id[3]);
            continue;
        }
        if(t.sourceNormalVertexByVertex[i]!=-1 ||
           t.sourceProjectiveRayByVertex[i]!=std::array<int,4>{-1,-1,-1,-1})
            return fail("paper cap point has interface provenance");
        std::set<int> loop;
        for(int f=t.interfaceFaceCount;f<t.closedBoundary.faces.rows();++f) {
            bool contains=false;for(int k=0;k<3;++k)contains|=t.closedBoundary.faces(f,k)==i;
            if(contains)for(int k=0;k<3;++k)if(t.closedBoundary.faces(f,k)!=i)
                loop.insert(t.closedBoundary.faces(f,k));
        }
        if(loop.size()<3)return fail("paper cap lacks a complete loop");
        Eigen::Vector3d center=Eigen::Vector3d::Zero();out.provenance[i]="paper-cap:";
        for(int vertex:loop) {
            if(vertex<0 || vertex>=t.interfaceVertexCount)return fail("invalid paper cap-loop vertex");
            center+=v.row(vertex).transpose();out.velocity[i]+=out.velocity[vertex];
            out.provenance[i]+=std::to_string(vertex)+",";
        }
        center/=loop.size();out.velocity[i]/=loop.size();
        if((center-v.row(i).transpose()).norm()>1e-12*std::max(1.0,center.norm()))
            return fail("paper cap centroid differs from native construction");
    }
    std::ostringstream identity;identity.imbue(std::locale::classic());
    identity<<"PaperProjectiveSourceMotion-v1;"<<q.hash()<<';'<<t.hash<<';'<<std::hexfloat;
    for(int i=0;i<v.rows();++i) {
        identity<<out.provenance[i]<<';';
        for(int a=0;a<3;++a)for(int b=0;b<10;++b)identity<<out.velocity[i](a,b)<<';';
    }
    out.defined=true;out.diagnostic=extrusionLoops.empty()
        ? "paper finite/ideal triangles, straight-ray truncation and cap centroids; fixed topology/radius"
        : "paper finite/ideal triangles, straight-ray truncation and remote planar cap extrusion v1; fixed topology/radius";
    out.scientificHash=manufacturing::computeContentSha256(identity.str());return out;
}
}

static QuadricSourceMotion sourceMotionImpl(const QuadricCoefficients10& q,
    const Eigen::MatrixXd& design, const QuadricTriangulationSettings& settings,
    const QuadricTriangulation* actualSource) {
    QuadricSourceMotion out;
    const auto fail=[&](const std::string& reason) {
        out.defined=false;out.velocity.clear();out.provenance.clear();
        out.diagnostic=reason;out.scientificHash.clear();return out;
    };
    if(settings.icosphereSubdivisions<0 || settings.icosphereSubdivisions>7 ||
       (settings.maximumIcosphereSubdivisions!=-1 &&
        (settings.maximumIcosphereSubdivisions<settings.icosphereSubdivisions ||
         settings.maximumIcosphereSubdivisions>7)))
        return fail("source motion requires requested subdivisions in [0,7] and maximum subdivisions of -1 or [requested,7]");
    const auto subtype=classifyQuadricSurfaceSubtype(q);
    if(subtype!=QuadricSurfaceSubtype::Ellipsoid &&
       subtype!=QuadricSurfaceSubtype::OneSheetHyperboloid &&
       subtype!=QuadricSurfaceSubtype::TwoSheetHyperboloid)
    {
        const auto tangent=evaluateQuadricStratumTangent(q);
        return fail("RANK_DEFICIENT_STRATUM_DERIVATIVE_UNDEFINED: "+
            std::string(quadricSurfaceSubtypeToken(subtype))+" fixed-stratum tangent dimension="+
            std::to_string(tangent.tangentDimension));
    }
    auto recorded=settings;recorded.retainFixedNormalProvenance=true;
    if(actualSource) {
        if(!settings.rationalSourceCoordinates || settings.algebraicFiniteMap ||
           settings.paperProjectiveHyperboloid || settings.fixedNormalParaboloid ||
           !rationalQuadricConstructionIdentityMatches(q,*actualSource,design,settings) ||
           actualSource->initialIcosphereSubdivisions!=settings.icosphereSubdivisions ||
           actualSource->effectiveIcosphereSubdivisions<settings.icosphereSubdivisions ||
           actualSource->effectiveIcosphereSubdivisions>
               (settings.maximumIcosphereSubdivisions==-1?settings.icosphereSubdivisions:
                   settings.maximumIcosphereSubdivisions))
            return fail("invalid or stale same-operation rational source construction");
        out.triangulation=*actualSource;
    } else out.triangulation=triangulateQuadric(q,design,recorded);
    if(out.triangulation.paperProjectiveHyperboloid)
        return paperProjectiveMotion(q,std::move(out.triangulation),settings.algebraicFiniteMap);
    const auto& t=out.triangulation;
    const auto& vertices=t.closedBoundary.vertices;
    const double radius=t.truncationRadius;
    if (t.type == QuadricSurfaceType::Hyperboloid &&
        t.sourceNormalEdgeByVertex.size() !=
            static_cast<std::size_t>(vertices.rows()))
        return fail("hyperboloid triangulation omitted source normal-edge provenance");
    if (t.type == QuadricSurfaceType::Hyperboloid &&
        t.sourceNormalVertexByVertex.size() !=
            static_cast<std::size_t>(vertices.rows()))
        return fail("hyperboloid triangulation omitted source normal-vertex provenance");
    auto sphereSettings=recorded;
    // If the opt-in cap-topology policy had to refine the normal sphere,
    // rebuild the reference graph at the exact accepted resolution. Reusing
    // the requested coarse graph here would bind Jacobians to a different
    // discrete material map than the one returned by triangulation.
    sphereSettings.icosphereSubdivisions =
        t.effectiveIcosphereSubdivisions >= 0
            ? t.effectiveIcosphereSubdivisions
            : settings.icosphereSubdivisions;
    sphereSettings.maximumIcosphereSubdivisions = -1;
    // Hyperboloid triangulation may double the requested truncation radius
    // until it obtains the two finite boundary loops.  Rebuild the reference
    // normal sphere at that *effective* radius so its inside/crossing graph is
    // the same branch as the returned source mesh.
    sphereSettings.truncationRadius=radius;
    QuadricCoefficients10::Vector sphereParameters;
    sphereParameters<<1,0,0,0,-1,-1,-1,0,0,0;
    const auto sphere=triangulateQuadric(
        QuadricCoefficients10::fromUnnormalized(sphereParameters),design,sphereSettings);
    if(sphere.truncationRadius!=radius)
        return fail("effective truncation radius was not preserved by the reference normal sphere");
    const auto& normals=sphere.fixedNormals;
    QuadricInverseGaussMap map(q);
    const TruncationRootCertifier rootCertifier(q,radius);
    std::vector<QuadricInverseGaussPointResult> points;
    std::vector<bool> inside;
    for(int i=0;i<normals.rows();++i) {
        points.push_back(settings.algebraicFiniteMap
            ? map.algebraicPoint(normals.row(i)) : map.point(normals.row(i)));
        const auto& p=points.back();
        if(p.valid() && t.type==QuadricSurfaceType::Hyperboloid &&
           std::abs(p.point->norm()-radius)<=1e-8*std::max(1.0,radius))
            return fail("original normal at truncation membership transition");
        inside.push_back(p.valid() && p.point->norm()<=radius);
    }
    // The old source-identity recovery scanned every fixed normal for every
    // interface vertex.  That is quadratic in the icosphere size and made
    // levels 6–7 unusable.  Index valid inverse-Gauss points in a conservative
    // spatial grid, then retain the original Euclidean tolerance and ambiguity
    // check on every nearby candidate.  Grid lookup only prunes impossible
    // matches; it does not relax source-vertex identity semantics.
    double identityScale=1.0;
    for(const auto& point:points)
        if(point.valid())identityScale=std::max(identityScale,point.point->norm());
    for(Eigen::Index i=0;i<vertices.rows();++i)
        identityScale=std::max(identityScale,vertices.row(i).norm());
    if(!std::isfinite(identityScale))return fail("nonfinite scale in fixed-normal identity index");
    const double identityCellWidth=1e-12*identityScale;
    using IdentityCell=std::array<std::int64_t,3>;
    const auto identityCell=[&](const Eigen::Vector3d& point,
                                IdentityCell& cell)->bool {
        for(int axis=0;axis<3;++axis) {
            const double scaled=point(axis)/identityCellWidth;
            if(!std::isfinite(scaled) ||
               scaled<=static_cast<double>(std::numeric_limits<std::int64_t>::min()+1) ||
               scaled>=static_cast<double>(std::numeric_limits<std::int64_t>::max()-1))
                return false;
            cell[static_cast<std::size_t>(axis)]=
                static_cast<std::int64_t>(std::floor(scaled));
        }
        return true;
    };
    std::map<IdentityCell,std::vector<int>> fixedNormalIndex;
    for(int j=0;j<static_cast<int>(points.size());++j)if(points[j].valid()) {
        IdentityCell cell{};
        if(!identityCell(*points[j].point,cell))
            return fail("fixed-normal point exceeds spatial identity-index range");
        fixedNormalIndex[cell].push_back(j);
    }
    std::set<std::pair<int,int>> crossing;
    if(t.type==QuadricSurfaceType::Hyperboloid) {
        for(int f=0;f<sphere.closedBoundary.faces.rows();++f)
            for(int k=0;k<3;++k) {
                int a=sphere.closedBoundary.faces(f,k),b=sphere.closedBoundary.faces(f,(k+1)%3);
                if(a>b)std::swap(a,b);
                if(inside[a]!=inside[b])crossing.emplace(a,b);
            }
    }
    const Eigen::Matrix<double,10,10> tangent=Eigen::Matrix<double,10,10>::Identity()-q.values()*q.values().transpose();
    out.velocity.resize(vertices.rows(),SourceVertexJacobian::Zero());
    out.provenance.resize(vertices.rows());
    std::set<int> usedOriginalNormalVertices;
    for(int i=0;i<t.interfaceVertexCount;++i) {
        const Eigen::Vector3d x=vertices.row(i);
        int original=-1;
        if(t.type==QuadricSurfaceType::Hyperboloid) {
            original=t.sourceNormalVertexByVertex[static_cast<std::size_t>(i)];
            if(original>=0) {
                if(original>=normals.rows() ||
                   !usedOriginalNormalVertices.insert(original).second ||
                   t.sourceNormalEdgeByVertex[static_cast<std::size_t>(i)][0]>=0 ||
                   t.sourceNormalEdgeByVertex[static_cast<std::size_t>(i)][1]>=0 ||
                   !points[static_cast<std::size_t>(original)].valid() ||
                   (x-*points[static_cast<std::size_t>(original)].point).norm()>
                       1e-12*std::max(1.0,x.norm()))
                    return fail("invalid native fixed-normal vertex identity");
            }
        } else {
            IdentityCell center{};
            if(!identityCell(x,center))
                return fail("interface point exceeds spatial identity-index range");
            for(std::int64_t dx=-1;dx<=1;++dx)
                for(std::int64_t dy=-1;dy<=1;++dy)
                    for(std::int64_t dz=-1;dz<=1;++dz) {
                        const IdentityCell cell{{center[0]+dx,center[1]+dy,center[2]+dz}};
                        const auto found=fixedNormalIndex.find(cell);
                        if(found==fixedNormalIndex.end())continue;
                        for(const int candidate:found->second)
                            if((*points[static_cast<std::size_t>(candidate)].point-x).norm()<=
                               1e-12*std::max(1.0,x.norm())) {
                                if(original>=0)return fail("ambiguous original fixed-normal identity");
                                original=candidate;
                            }
                    }
        }
        if(original>=0) {
            const auto differential=settings.algebraicFiniteMap
                ? map.algebraicCoefficientDifferential(normals.row(original))
                : map.coefficientDifferential(normals.row(original));
            if(!differential.valid())return fail(differential.diagnostic);
            out.velocity[i]=*differential.coefficientJacobian;
            out.provenance[i]="normal:"+std::to_string(original);
            continue;
        }
        if(t.type!=QuadricSurfaceType::Hyperboloid ||
           std::abs(x.norm()-radius)>1e-10*std::max(1.0,radius))
            return fail("unmatched finite source point or uncertified truncation root");
        const Eigen::Vector3d gradient=q.gradient(x);
        if(!gradient.allFinite() || gradient.norm()<1e-12)
            return fail("singular implicit surface at truncation root");
        const Eigen::Vector3d n=gradient.normalized();
        const auto identity = t.sourceNormalEdgeByVertex[static_cast<std::size_t>(i)];
        if (identity[0] < 0 || identity[1] <= identity[0] ||
            identity[1] >= normals.rows() ||
            crossing.find({identity[0], identity[1]}) == crossing.end() ||
            inside[static_cast<std::size_t>(identity[0])] ==
                inside[static_cast<std::size_t>(identity[1])])
            return fail("truncation vertex lacks a valid native normal-edge identity");
        const Eigen::Vector3d a=normals.row(identity[0]),b=normals.row(identity[1]);
        const Eigen::Vector3d edgeNormal=a.cross(b).normalized();
        const double edgeDot=a.dot(b),denominator=1-edgeDot*edgeDot;
        if (!std::isfinite(denominator) || denominator <= 1e-15 ||
            std::abs(edgeNormal.dot(n)) > 1e-10)
            return fail("native truncation edge does not certify the root normal");
        const double u=(a.dot(n)-edgeDot*b.dot(n))/denominator;
        const double v=(b.dot(n)-edgeDot*a.dot(n))/denominator;
        if (!std::isfinite(u) || !std::isfinite(v) || u<=1e-9 || v<=1e-9)
            return fail("native truncation edge root is not in the edge interior");
        std::string rootCertificate;
        if(!rootCertifier.certify(a,b,&rootCertificate))
            return fail(rootCertificate);
        if(std::abs(q.evaluate(x))>1e-10*std::max(1.0,x.squaredNorm()))
            return fail("native bisection point does not satisfy implicit root equation");
        Eigen::Matrix3d system;
        system.row(0)=gradient.transpose();
        system.row(1)=2*x.transpose();
        system.row(2)=2*edgeNormal.transpose()*q.quadraticMatrix();
        Eigen::Matrix<double,3,10> rhs=Eigen::Matrix<double,3,10>::Zero();
        rhs.row(0)<<-1,-x.x(),-x.y(),-x.z(),-x.x()*x.x(),-x.y()*x.y(),-x.z()*x.z(),
                    -2*x.y()*x.z(),-2*x.x()*x.z(),-2*x.x()*x.y();
        Eigen::Matrix<double,3,10> gradientDerivative=Eigen::Matrix<double,3,10>::Zero();
        gradientDerivative.block<3,3>(0,1)=Eigen::Matrix3d::Identity();
        gradientDerivative(0,4)=2*x.x();gradientDerivative(1,5)=2*x.y();gradientDerivative(2,6)=2*x.z();
        gradientDerivative(1,7)=2*x.z();gradientDerivative(2,7)=2*x.y();
        gradientDerivative(0,8)=2*x.z();gradientDerivative(2,8)=2*x.x();
        gradientDerivative(0,9)=2*x.y();gradientDerivative(1,9)=2*x.x();
        rhs.row(2)=-edgeNormal.transpose()*gradientDerivative;
        rhs=rhs*tangent;
        for(int row=0;row<3;++row) {
            const double scale=system.row(row).norm();
            if(!std::isfinite(scale) || scale<1e-12)return fail("singular truncation constraint");
            system.row(row)/=scale;rhs.row(row)/=scale;
        }
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(system,Eigen::ComputeFullU|Eigen::ComputeFullV);
        if(svd.singularValues()(2)<1e-10*svd.singularValues()(0))
            return fail("ill-conditioned truncation root differential");
        out.velocity[i]=svd.solve(rhs);
        if(!out.velocity[i].allFinite() ||
           (system*out.velocity[i]-rhs).norm()>1e-10*std::max(1.0,rhs.norm()))
            return fail("uncertified truncation root differential");
        out.provenance[i]="edge:"+std::to_string(identity[0])+","+std::to_string(identity[1]);
    }
    for(int i=t.interfaceVertexCount;i<vertices.rows();++i) {
        std::set<int> loop;
        for(int f=t.interfaceFaceCount;f<t.closedBoundary.faces.rows();++f) {
            bool contains=false;
            for(int k=0;k<3;++k)contains|=t.closedBoundary.faces(f,k)==i;
            if(contains)for(int k=0;k<3;++k) {
                const int v=t.closedBoundary.faces(f,k);
                if(v!=i) {
                    if(v<0 || v>=t.interfaceVertexCount)return fail("unsupported cap construction");
                    loop.insert(v);
                }
            }
        }
        if(loop.size()<3)return fail("cap centroid lacks complete loop identity");
        Eigen::Vector3d center=Eigen::Vector3d::Zero();
        std::string key="cap:";
        for(int v:loop) {
            center+=vertices.row(v).transpose();
            out.velocity[i]+=out.velocity[v];key+=std::to_string(v)+",";
        }
        center/=loop.size();out.velocity[i]/=loop.size();
        if((center-vertices.row(i).transpose()).norm()>1e-12*std::max(1.0,center.norm()))
            return fail("cap centroid does not match native construction");
        out.provenance[i]=key;
    }
    std::ostringstream identity;identity.imbue(std::locale::classic());
    identity<<"NativeFullRankSourceMotion-v3-extended-solve-exact-root;"<<q.hash()<<';'<<t.hash<<';'<<std::hexfloat;
    for(int i=0;i<vertices.rows();++i) {
        identity<<out.provenance[i]<<';';
        for(int a=0;a<3;++a)for(int b=0;b<10;++b)identity<<out.velocity[i](a,b)<<';';
    }
    out.defined=true;
    out.diagnostic="regular fixed-normal/truncation-root/cap-centroid construction";
    out.scientificHash=manufacturing::computeContentSha256(identity.str());
    return out;
}

QuadricSourceMotion evaluateQuadricSourceMotion(const QuadricCoefficients10& q,
    const Eigen::MatrixXd& design,const QuadricTriangulationSettings& settings) {
    return sourceMotionImpl(q,design,settings,nullptr);
}
QuadricSourceMotion evaluateQuadricSourceMotion(const QuadricCoefficients10& q,
    const Eigen::MatrixXd& design,const QuadricTriangulationSettings& settings,
    const QuadricTriangulation& actualSource) {
    return sourceMotionImpl(q,design,settings,&actualSource);
}

QuadricRationalSourceMotion evaluateQuadricSourceMotion(
    const QuadricCoefficientAuthority10& authority,const Eigen::MatrixXd& design,
    const QuadricTriangulationSettings& settings,const QuadricTriangulation& source) {
    using Q=CGAL::Gmpq;
    QuadricRationalSourceMotion out;
    const auto fail=[&](const std::string& reason) {
        out.defined=false;out.coefficientJacobian.clear();out.lower.clear();out.upper.clear();
        out.pointProofLower.clear();out.pointProofUpper.clear();out.provenance.clear();
        out.scientificHash.clear();out.diagnostic=reason;return out;
    };
    if(!rationalQuadricConstructionIdentityMatches(authority,source,design,settings))
        return fail("invalid or stale authority same-operation source construction");
    if(source.type!=QuadricSurfaceType::Ellipsoid && source.type!=QuadricSurfaceType::Hyperboloid)
        return fail("authority source motion requires a supported full-rank chart");
    const std::size_t count=source.exactVertices.size();
    if(source.interfaceVertexCount<=0 || source.interfaceVertexCount>int(count) ||
       source.interfaceFaceCount<=0 || source.interfaceFaceCount>source.closedBoundary.faces.rows())
        return fail("invalid authority interface counts");
    const bool hyperboloid=source.type==QuadricSurfaceType::Hyperboloid;
    if(hyperboloid && (source.sourceNormalVertexByVertex.size()!=count ||
                      source.sourceNormalEdgeByVertex.size()!=count))
        return fail("authority hyperboloid omitted same-source feature identities");
    if(!hyperboloid && (source.fixedNormals.rows()!=int(count) ||
        source.interfaceVertexCount!=int(count) || source.interfaceFaceCount!=source.closedBoundary.faces.rows()))
        return fail("authority ellipsoid fixed-normal layout differs from actual source");
    out.triangulation=source;
    out.coefficientJacobian.resize(count);out.lower.resize(count);out.upper.resize(count);
    out.pointProofLower.resize(count);out.pointProofUpper.resize(count);out.provenance.resize(count);
    const TruncationRootCertifier roots(authority,source.truncationRadius);
    const auto samePoint=[&](const QuadricRationalPointResult& point,int vertex) {
        return point.defined && point.coordinates==source.exactVertices[vertex] &&
            point.lower==source.coordinateLower[vertex] && point.upper==source.coordinateUpper[vertex];
    };
    const auto text=[](const Q& value) {
        std::ostringstream stream;stream.imbue(std::locale::classic());stream<<value;return stream.str();
    };
    std::set<int> normals;std::set<std::pair<int,int>> edges;
    try {
        for(int vertex=0;vertex<source.interfaceVertexCount;++vertex) {
            const int normal=hyperboloid ? source.sourceNormalVertexByVertex[vertex] : vertex;
            if(normal>=0) {
                if(normal>=source.fixedNormals.rows() || !normals.insert(normal).second ||
                   (hyperboloid && source.sourceNormalEdgeByVertex[vertex]!=std::array<int,2>{-1,-1}))
                    return fail("authority finite point has invalid fixed-normal provenance");
                const auto differential=QuadricInverseGaussMap::rationalDifferentialFromAuthority(
                    authority,source.fixedNormals.row(normal),settings.compactRationalSourceCoordinates);
                if(!differential.defined)return fail(differential.diagnostic);
                if(!samePoint(differential.point,vertex))
                    return fail("authority finite differential does not match actual rational Boolean operand");
                out.coefficientJacobian[vertex]=differential.coefficientJacobian;
                out.lower[vertex]=differential.lower;out.upper[vertex]=differential.upper;
                out.pointProofLower[vertex]=differential.point.lower;
                out.pointProofUpper[vertex]=differential.point.upper;
                out.provenance[vertex]="normal:"+std::to_string(normal);
            } else {
                const auto edge=source.sourceNormalEdgeByVertex[vertex];
                if(edge[0]<0 || edge[1]<=edge[0] || edge[1]>=source.fixedNormals.rows() ||
                   !edges.emplace(edge[0],edge[1]).second)
                    return fail("authority truncation point has invalid normal-edge provenance");
                const auto differential=roots.rationalDifferential(source.fixedNormals.row(edge[0]),
                    source.fixedNormals.row(edge[1]),settings.compactRationalSourceCoordinates);
                if(!differential.defined)return fail(differential.diagnostic);
                if(!samePoint(differential.point,vertex))
                    return fail("authority root differential does not match actual rational Boolean operand");
                out.coefficientJacobian[vertex]=differential.coefficientJacobian;
                out.lower[vertex]=differential.lower;out.upper[vertex]=differential.upper;
                out.pointProofLower[vertex]=differential.derivativePointLower;
                out.pointProofUpper[vertex]=differential.derivativePointUpper;
                out.provenance[vertex]="edge:"+std::to_string(edge[0])+","+std::to_string(edge[1]);
            }
        }
        for(int center=source.interfaceVertexCount;center<int(count);++center) {
            if(source.sourceNormalVertexByVertex[center]!=-1 ||
               source.sourceNormalEdgeByVertex[center]!=std::array<int,2>{-1,-1})
                return fail("authority cap center has interface feature provenance");
            std::set<int> ring;
            for(int face=source.interfaceFaceCount;face<source.closedBoundary.faces.rows();++face) {
                bool contains=false;for(int corner=0;corner<3;++corner)
                    contains|=source.closedBoundary.faces(face,corner)==center;
                if(contains)for(int corner=0;corner<3;++corner) {
                    const int vertex=source.closedBoundary.faces(face,corner);
                    if(vertex==center)continue;
                    if(vertex<0 || vertex>=source.interfaceVertexCount)
                        return fail("authority cap has unsupported loop construction");
                    ring.insert(vertex);
                }
            }
            if(ring.size()<3)return fail("authority cap lacks a complete same-source loop");
            const Q denominator(static_cast<unsigned long>(ring.size()));
            out.provenance[center]="cap:";
            for(int vertex:ring)out.provenance[center]+=std::to_string(vertex)+",";
            for(int axis=0;axis<3;++axis) {
                Q coordinate(0),low(0),high(0),proofLow(0),proofHigh(0);
                for(int vertex:ring) {
                    coordinate+=Q(source.exactVertices[vertex][axis]);
                    low+=Q(source.coordinateLower[vertex][axis]);high+=Q(source.coordinateUpper[vertex][axis]);
                    proofLow+=Q(out.pointProofLower[vertex][axis]);proofHigh+=Q(out.pointProofUpper[vertex][axis]);
                }
                if(coordinate/denominator!=Q(source.exactVertices[center][axis]) ||
                   low/denominator!=Q(source.coordinateLower[center][axis]) ||
                   high/denominator!=Q(source.coordinateUpper[center][axis]))
                    return fail("authority cap center is not the actual rational loop mean");
                out.pointProofLower[center][axis]=text(proofLow/denominator);
                out.pointProofUpper[center][axis]=text(proofHigh/denominator);
                for(int column=0;column<10;++column) {
                    Q value(0),lower(0),upper(0);
                    for(int vertex:ring) {
                        value+=Q(out.coefficientJacobian[vertex][axis][column]);
                        lower+=Q(out.lower[vertex][axis][column]);upper+=Q(out.upper[vertex][axis][column]);
                    }
                    out.coefficientJacobian[center][axis][column]=text(value/denominator);
                    out.lower[center][axis][column]=text(lower/denominator);
                    out.upper[center][axis][column]=text(upper/denominator);
                }
            }
        }
    } catch(const std::runtime_error& error) {return fail(error.what());}
    std::ostringstream identity;identity.imbue(std::locale::classic());
    identity<<"AuthorityRationalSourceMotion-v1;"<<authority.policy<<';'<<authority.precisionBits<<';'
            <<source.rationalConstructionIdentity<<';';
    for(std::size_t vertex=0;vertex<count;++vertex) {
        identity<<out.provenance[vertex]<<';';
        for(int axis=0;axis<3;++axis) {
            identity<<out.pointProofLower[vertex][axis]<<';'<<out.pointProofUpper[vertex][axis]<<';';
            for(int column=0;column<10;++column)
                identity<<out.coefficientJacobian[vertex][axis][column]<<';'
                        <<out.lower[vertex][axis][column]<<';'<<out.upper[vertex][axis][column]<<';';
        }
    }
    out.defined=true;
    out.scientificHash=manufacturing::computeContentSha256(identity.str());
    out.diagnostic="same-operation rational finite/root/cap differential; fixed chart/topology/radius; not Boolean moment certification";
    return out;
}

QuadricStratumSourceMotion evaluateQuadricStratumSourceMotion(
    const QuadricCoefficients10& q,const Eigen::MatrixXd& design,
    const QuadricTriangulationSettings& settings) {
    QuadricStratumSourceMotion out;
    const auto fail=[&](const std::string& reason) {
        out.defined=false;out.coefficientBasis.resize(0,0);
        out.vertexVelocity.clear();out.provenance.clear();
        out.diagnostic=reason;out.scientificHash.clear();return out;
    };
    const auto tangent=evaluateQuadricStratumTangent(q);
    if(tangent.subtype==QuadricSurfaceSubtype::Plane) {
        if(!tangent.defined || tangent.tangentDimension!=3 ||
           tangent.basis.rows()!=10 || tangent.basis.cols()!=3)
            return fail("plane tangent chart is not three-dimensional");
        const Eigen::Vector3d linear=q.linearVector();
        const double linearNorm=linear.norm();
        if(!std::isfinite(linearNorm) || linearNorm<=1e-12)
            return fail("plane linear coefficient is singular");
        const Eigen::Vector3d normal=linear/linearNorm;
        const auto frame=makeQuadricPlaneFrame(normal);
        if(!frame.chartStable)
            return fail("plane tangent-frame chart boundary: "+frame.diagnostic);
        const Eigen::Vector3d first=frame.first;
        const Eigen::Vector3d second=frame.second;
        const int referenceAxis=frame.chart==0?0:1;
        Eigen::Vector3d reference=Eigen::Vector3d::Zero();
        reference(referenceAxis)=1.0;

        out.triangulation=triangulateQuadric(q,design,settings);
        if(out.triangulation.type!=QuadricSurfaceType::Plane ||
           out.triangulation.interfaceVertexCount!=4 ||
           out.triangulation.closedBoundary.vertices.rows()!=8)
            return fail("unsupported native plane patch/extrusion vertex layout");
        const auto& mesh=out.triangulation.closedBoundary;
        const double radius=out.triangulation.truncationRadius;
        const double constant=q.constant();
        const Eigen::Vector3d origin=-constant*normal/linearNorm;
        const std::array<std::pair<double,double>,4> signs{{
            {-1.0,-1.0},{1.0,-1.0},{1.0,1.0},{-1.0,1.0}}};
        out.coefficientBasis=tangent.basis;
        out.vertexVelocity.resize(8,Eigen::MatrixXd::Zero(3,3));
        out.provenance.resize(8);
        for(int column=0;column<3;++column) {
            const Eigen::Matrix<double,10,1> dq=tangent.basis.col(column);
            if(!tangent.contains(dq,1e-9))
                return fail("plane tangent basis violates its normalized plane constraints");
            const Eigen::Vector3d dLinear=dq.segment<3>(1);
            const double dLinearNorm=normal.dot(dLinear);
            const Eigen::Vector3d dNormal=(dLinear-normal*dLinearNorm)/linearNorm;
            const double dConstant=dq(0);
            const Eigen::Vector3d dOrigin=-dConstant*normal/linearNorm-
                constant*dNormal/linearNorm+constant*normal*dLinearNorm/(linearNorm*linearNorm);

            const Eigen::Vector3d unnormalizedDerivative=dNormal.cross(reference);
            const Eigen::Vector3d dFirst=(unnormalizedDerivative-
                first*first.dot(unnormalizedDerivative))/
                (normal.cross(reference)).norm();
            const Eigen::Vector3d dSecond=dNormal.cross(first)+normal.cross(dFirst);
            const Eigen::Vector3d expectedOrigin=-constant*linear/(linear.squaredNorm());
            if((expectedOrigin-origin).norm()>1e-11*std::max(1.0,origin.norm()))
                return fail("plane origin differential disagrees with native construction");
            for(int vertex=0;vertex<4;++vertex) {
                const auto [signFirst,signSecond]=signs[static_cast<std::size_t>(vertex)];
                const Eigen::Vector3d expected=origin+radius*(signFirst*first+signSecond*second);
                if((expected-mesh.vertices.row(vertex).transpose()).norm()>
                   1e-11*std::max(1.0,expected.norm()))
                    return fail("plane patch vertex disagrees with native construction");
                const Eigen::Vector3d velocity=dOrigin+radius*(signFirst*dFirst+signSecond*dSecond);
                out.vertexVelocity[static_cast<std::size_t>(vertex)].col(column)=velocity;
                out.vertexVelocity[static_cast<std::size_t>(vertex+4)].col(column)=
                    velocity+4.0*radius*dNormal;
            }
        }
        for(int vertex=0;vertex<4;++vertex) {
            out.provenance[static_cast<std::size_t>(vertex)]="plane-corner:"+std::to_string(vertex);
            out.provenance[static_cast<std::size_t>(vertex+4)]="plane-extruded-corner:"+std::to_string(vertex);
        }
        std::ostringstream identity;identity.imbue(std::locale::classic());
        identity<<"NativePlaneStratumSourceMotion-v1;"<<q.hash()<<';'
                <<out.triangulation.hash<<';'<<std::hexfloat;
        for(int row=0;row<out.coefficientBasis.rows();++row)
            for(int column=0;column<out.coefficientBasis.cols();++column)
                identity<<out.coefficientBasis(row,column)<<';';
        for(std::size_t vertex=0;vertex<out.vertexVelocity.size();++vertex) {
            identity<<out.provenance[vertex]<<';';
            for(int axis=0;axis<3;++axis)
                for(int column=0;column<3;++column)
                    identity<<out.vertexVelocity[vertex](axis,column)<<';';
        }
        out.defined=true;
        out.diagnostic="reduced plane stratum source motion; fixed orthogonal-frame branch, square patch and extrusion radius";
        out.scientificHash=manufacturing::computeContentSha256(identity.str());
        return out;
    }
    if(tangent.subtype==QuadricSurfaceSubtype::Paraboloid) {
        if(!tangent.defined || tangent.tangentDimension!=8 ||
           tangent.basis.rows()!=10 || tangent.basis.cols()!=8)
            return fail("paraboloid tangent chart is not eight-dimensional");
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(q.quadraticMatrix());
        if(eigen.info()!=Eigen::Success)return fail("paraboloid eigensystem failed");
        const auto eigenvalues=eigen.eigenvalues();
        const Eigen::Matrix3d eigenvectors=canonicalEigenvectorSigns(eigen.eigenvectors());
        if(!eigenvectorSignPivotsAreSeparated(eigenvectors))
            return fail("paraboloid eigenvector sign-pivot boundary is not differentiable in this chart");
        int nullIndex=0;eigenvalues.cwiseAbs().minCoeff(&nullIndex);
        const double eigenScale=std::max(1.0,eigenvalues.cwiseAbs().maxCoeff());
        for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)
            if(std::abs(eigenvalues(i)-eigenvalues(j))<=1e-10*eigenScale)
                return fail("paraboloid source derivative requires a simple, separated eigenspectrum");
        std::array<int,2> curved{};int cursor=0;
        for(int i=0;i<3;++i)if(i!=nullIndex)curved[cursor++]=i;
        const int uIndex=curved[0],vIndex=curved[1];
        const Eigen::Vector3d u=eigenvectors.col(uIndex),v=eigenvectors.col(vIndex);
        Eigen::Vector3d n=eigenvectors.col(nullIndex);
        const Eigen::Vector3d b=q.linearVector();
        double bn=b.dot(n);
        if(!std::isfinite(bn) || std::abs(bn)<=1e-12*std::max(1.0,b.norm()))
            return fail("rank-two quadric is not a regular paraboloid: null-axis linear coefficient vanishes");
        const double normalSign=bn<0.0?-1.0:1.0;
        n*=normalSign;bn*=normalSign;

        out.triangulation=triangulateQuadric(q,design,settings);
        if(out.triangulation.type!=QuadricSurfaceType::Paraboloid)
            return fail("native construction changed paraboloid type");
        const auto& mesh=out.triangulation.closedBoundary;
        const double radius=out.triangulation.truncationRadius;
        const int resolution=settings.paraboloidGridResolution;
        const int side=resolution+1;
        const bool fixedNormals=settings.fixedNormalParaboloid;
        const int interfaceCount=fixedNormals?out.triangulation.interfaceVertexCount:side*side;
        if(resolution<1 || out.triangulation.interfaceVertexCount!=interfaceCount ||
           mesh.vertices.rows()!=2*interfaceCount)
            return fail("unsupported native paraboloid grid/extruded-cap vertex layout");
        if(fixedNormals) {
            if(!out.triangulation.fixedNormalParaboloid ||
               out.triangulation.sourceNormalVertexByVertex.size()!=static_cast<std::size_t>(mesh.vertices.rows()))
                return fail("missing fixed-normal paraboloid provenance");
            for(int i=0;i<out.triangulation.fixedNormals.rows();++i)
                if(std::abs(out.triangulation.fixedNormals.row(i).dot(n)-0.01)<=1e-8)
                    return fail("fixed-normal paraboloid finite-patch selection boundary");
        }
        const double au=eigenvalues(uIndex),av=eigenvalues(vIndex);
        if(std::abs(au)<=1e-12*eigenScale || std::abs(av)<=1e-12*eigenScale)
            return fail("paraboloid non-null curvature is numerically singular");

        out.coefficientBasis=tangent.basis;
        out.vertexVelocity.resize(static_cast<std::size_t>(mesh.vertices.rows()));
        out.provenance.resize(static_cast<std::size_t>(mesh.vertices.rows()));
        for(auto& velocity:out.vertexVelocity)velocity=Eigen::MatrixXd::Zero(3,8);
        for(int column=0;column<8;++column) {
            const auto dq=tangent.basis.col(column).eval();
            if(!tangent.contains(dq,1e-9))return fail("basis column violates paraboloid constraints");
            const Eigen::Matrix3d dA=quadraticDirection(dq);
            const Eigen::Vector3d db=dq.segment<3>(1);
            const auto eigenvectorDerivative=[&](int index) {
                Eigen::Vector3d derivative=Eigen::Vector3d::Zero();
                for(int j=0;j<3;++j)if(j!=index) {
                    const double denominator=eigenvalues(index)-eigenvalues(j);
                    if(std::abs(denominator)<=1e-10*eigenScale)
                        throw std::runtime_error("paraboloid eigenvector differential is ill-conditioned");
                    derivative+=eigenvectors.col(j)*
                        (eigenvectors.col(j).dot(dA*eigenvectors.col(index))/denominator);
                }
                return derivative;
            };
            Eigen::Vector3d du,dv,dn;
            try {
                du=eigenvectorDerivative(uIndex);
                dv=eigenvectorDerivative(vIndex);
                dn=normalSign*eigenvectorDerivative(nullIndex);
            } catch(const std::exception& error) {return fail(error.what());}
            const double dau=u.dot(dA*u),dav=v.dot(dA*v);
            const double dbu=db.dot(u)+b.dot(du),dbv=db.dot(v)+b.dot(dv);
            const double dbn=db.dot(n)+b.dot(dn);
            const double dc=dq(0);
            for(int index=0;index<interfaceCount;++index) {
                double t=-radius+2.0*radius*(index/side)/resolution;
                double s=-radius+2.0*radius*(index%side)/resolution;
                double ds=0,dt=0;
                if(fixedNormals) {
                    const int id=out.triangulation.sourceNormalVertexByVertex[index];
                    if(id<0 || id>=out.triangulation.fixedNormals.rows())
                        return fail("invalid fixed-normal paraboloid vertex identity");
                    const Eigen::Vector3d normal=out.triangulation.fixedNormals.row(id);
                    const double denominator=normal.dot(n);
                    if(denominator<=0.01)return fail("inadmissible finite paraboloid normal");
                    const double k=bn/denominator;
                    const double dk=dbn/denominator-bn*normal.dot(dn)/(denominator*denominator);
                    s=(k*normal.dot(u)-b.dot(u))/(2*au);
                    t=(k*normal.dot(v)-b.dot(v))/(2*av);
                    ds=(dk*normal.dot(u)+k*normal.dot(du)-dbu)/(2*au)-s*dau/au;
                    dt=(dk*normal.dot(v)+k*normal.dot(dv)-dbv)/(2*av)-t*dav/av;
                }
                const double numerator=q.constant()+b.dot(u)*s+b.dot(v)*t+
                    au*s*s+av*t*t;
                const double w=-numerator/bn;
                double dNumerator=dc+dbu*s+dbv*t+dau*s*s+dav*t*t;
                if(fixedNormals)dNumerator+=(b.dot(u)+2*au*s)*ds+(b.dot(v)+2*av*t)*dt;
                const double dw=-dNumerator/bn+numerator*dbn/(bn*bn);
                const Eigen::Vector3d expected=u*s+v*t+n*w;
                if((expected-mesh.vertices.row(index).transpose()).norm()>
                   1e-10*std::max(1.0,expected.norm()))
                    return fail("paraboloid source point disagrees with native grid construction");
                // Preserve the legacy arithmetic order when this profile is off.
                out.vertexVelocity[static_cast<std::size_t>(index)].col(column)=
                    du*s+dv*t+dn*w+n*dw;
                if(fixedNormals)out.vertexVelocity[static_cast<std::size_t>(index)].col(column)+=u*ds+v*dt;
                out.vertexVelocity[static_cast<std::size_t>(index+interfaceCount)].col(column)=
                    out.vertexVelocity[static_cast<std::size_t>(index)].col(column)+4.0*radius*dn;
            }
        }
        for(int i=0;i<interfaceCount;++i) {
            const std::string id=std::to_string(fixedNormals?out.triangulation.sourceNormalVertexByVertex[i]:i);
            out.provenance[static_cast<std::size_t>(i)]=(fixedNormals?"paraboloid-fixed-normal:":"paraboloid-grid:")+id;
            out.provenance[static_cast<std::size_t>(i+interfaceCount)]=(fixedNormals?"paraboloid-extruded-normal:":"paraboloid-extruded-grid:")+id;
        }
        std::ostringstream identity;identity.imbue(std::locale::classic());
        identity<<"NativeParaboloidStratumSourceMotion-v1;"<<q.hash()<<';'
                <<out.triangulation.hash<<';'<<std::hexfloat;
        for(int i=0;i<out.coefficientBasis.rows();++i)
            for(int j=0;j<out.coefficientBasis.cols();++j)identity<<out.coefficientBasis(i,j)<<';';
        for(std::size_t i=0;i<out.vertexVelocity.size();++i) {
            identity<<out.provenance[i]<<';';
            for(int a=0;a<3;++a)for(int j=0;j<8;++j)
                identity<<out.vertexVelocity[i](a,j)<<';';
        }
        out.defined=true;
        out.diagnostic=fixedNormals?
            "reduced fixed-normal paraboloid source motion; finite patch and fixed selection/extrusion branches; no ideal or cross-type derivative":
            "reduced paraboloid stratum source motion; fixed eigenbasis ordering, parameter grid and extrusion radius";
        out.scientificHash=manufacturing::computeContentSha256(identity.str());
        return out;
    }
    if(tangent.subtype==QuadricSurfaceSubtype::ParabolicCylinder) {
        if(!tangent.defined || tangent.tangentDimension!=6 ||
           tangent.basis.rows()!=10 || tangent.basis.cols()!=6)
            return fail("parabolic-cylinder tangent chart is not six-dimensional");
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(q.quadraticMatrix());
        if(eigen.info()!=Eigen::Success)return fail("parabolic-cylinder eigensystem failed");
        const auto eigenvalues=eigen.eigenvalues();
        const Eigen::Matrix3d eigenvectors=canonicalEigenvectorSigns(eigen.eigenvectors());
        int curvedIndex=0;eigenvalues.cwiseAbs().maxCoeff(&curvedIndex);
        const double eigenScale=std::max(1.0,eigenvalues.cwiseAbs().maxCoeff());
        const double lambda=eigenvalues(curvedIndex);
        if(std::abs(lambda)<=1e-12*eigenScale)
            return fail("parabolic-cylinder curvature is numerically singular");
        for(int other=0;other<3;++other)if(other!=curvedIndex &&
           std::abs(eigenvalues(other))>1e-10*eigenScale)
            return fail("parabolic-cylinder quadratic block left the rank-one stratum");
        const Eigen::Vector3d u=eigenvectors.col(curvedIndex);
        Eigen::Index dominant=0;u.cwiseAbs().maxCoeff(&dominant);
        double secondLargest=0.0;
        for(int axis=0;axis<3;++axis)if(axis!=dominant)
            secondLargest=std::max(secondLargest,std::abs(u(axis)));
        if(std::abs(std::abs(u(dominant))-secondLargest)<=1e-10)
            return fail("parabolic-cylinder principal-axis sign chart is on a pivot boundary");
        const Eigen::Vector3d b=q.linearVector();
        const double bu=b.dot(u);
        const Eigen::Vector3d transverse=b-bu*u;
        const double beta=transverse.norm();
        if(!std::isfinite(beta) || beta<=1e-12*std::max(1.0,b.norm()))
            return fail("parabolic-cylinder transverse linear component is singular");
        const Eigen::Vector3d v=transverse/beta;
        const Eigen::Vector3d w=u.cross(v).normalized();

        out.triangulation=triangulateQuadric(q,design,settings);
        if(out.triangulation.type!=QuadricSurfaceType::ParabolicCylinder)
            return fail("native construction changed parabolic-cylinder type");
        const auto& mesh=out.triangulation.closedBoundary;
        const double radius=out.triangulation.truncationRadius;
        const int resolution=settings.parabolicCylinderGridResolution;
        const int side=resolution+1;
        const int interfaceCount=side*side;
        if(resolution<1 || out.triangulation.interfaceVertexCount!=interfaceCount ||
           mesh.vertices.rows()!=2*interfaceCount)
            return fail("unsupported native parabolic-cylinder grid/extruded-cap layout");

        out.coefficientBasis=tangent.basis;
        out.vertexVelocity.resize(static_cast<std::size_t>(mesh.vertices.rows()));
        out.provenance.resize(static_cast<std::size_t>(mesh.vertices.rows()));
        for(auto& velocity:out.vertexVelocity)velocity=Eigen::MatrixXd::Zero(3,6);
        for(int column=0;column<6;++column) {
            const auto dq=tangent.basis.col(column).eval();
            if(!tangent.contains(dq,1e-9))
                return fail("basis column violates parabolic-cylinder rank-one constraints");
            const Eigen::Matrix3d dA=quadraticDirection(dq);
            const Eigen::Vector3d db=dq.segment<3>(1);
            const Eigen::Vector3d du=(Eigen::Matrix3d::Identity()-u*u.transpose())*
                dA*u/lambda;
            const double dlambda=u.dot(dA*u);
            const double dbu=db.dot(u)+b.dot(du);
            const Eigen::Vector3d dTransverse=db-dbu*u-bu*du;
            const double dbeta=v.dot(dTransverse);
            const Eigen::Vector3d dv=(dTransverse-v*dbeta)/beta;
            const Eigen::Vector3d dw=du.cross(v)+u.cross(dv);
            for(int row=0;row<=resolution;++row)for(int col=0;col<=resolution;++col) {
                const int index=row*side+col;
                const double r=-radius+2.0*radius*row/resolution;
                const double z=-radius+2.0*radius*col/resolution;
                const double numerator=q.constant()+bu*r+lambda*r*r;
                const double s=-numerator/beta;
                const double dNumerator=dq(0)+dbu*r+dlambda*r*r;
                const double ds=-dNumerator/beta+
                    numerator*dbeta/(beta*beta);
                const Eigen::Vector3d expected=u*r+v*s+w*z;
                if((expected-mesh.vertices.row(index).transpose()).norm()>
                   1e-10*std::max(1.0,expected.norm()))
                    return fail("parabolic-cylinder source point disagrees with native grid construction");
                out.vertexVelocity[static_cast<std::size_t>(index)].col(column)=
                    du*r+dv*s+v*ds+dw*z;
                out.vertexVelocity[static_cast<std::size_t>(index+interfaceCount)].col(column)=
                    out.vertexVelocity[static_cast<std::size_t>(index)].col(column)+4.0*radius*dv;
            }
        }
        for(int i=0;i<interfaceCount;++i) {
            out.provenance[static_cast<std::size_t>(i)]="parabolic-cylinder-grid:"+std::to_string(i);
            out.provenance[static_cast<std::size_t>(i+interfaceCount)]="parabolic-cylinder-extruded:"+std::to_string(i);
        }
        std::ostringstream identity;identity.imbue(std::locale::classic());
        identity<<"NativeParabolicCylinderStratumSourceMotion-v1;"<<q.hash()<<';'
                <<out.triangulation.hash<<';'<<std::hexfloat;
        for(int i=0;i<out.coefficientBasis.rows();++i)
            for(int j=0;j<out.coefficientBasis.cols();++j)identity<<out.coefficientBasis(i,j)<<';';
        for(std::size_t i=0;i<out.vertexVelocity.size();++i) {
            identity<<out.provenance[i]<<';';
            for(int a=0;a<3;++a)for(int j=0;j<6;++j)
                identity<<out.vertexVelocity[i](a,j)<<';';
        }
        out.defined=true;
        out.diagnostic="reduced parabolic-cylinder source motion; fixed principal-axis chart, parabolic graph grid and transverse closure branch";
        out.scientificHash=manufacturing::computeContentSha256(identity.str());
        return out;
    }
    if(tangent.subtype==QuadricSurfaceSubtype::HyperbolicCylinder) {
        if(!tangent.defined || tangent.tangentDimension!=7 ||
           tangent.basis.rows()!=10 || tangent.basis.cols()!=7)
            return fail("hyperbolic-cylinder tangent chart is not seven-dimensional");
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(q.quadraticMatrix());
        if(eigen.info()!=Eigen::Success)return fail("hyperbolic-cylinder eigensystem failed");
        const auto values=eigen.eigenvalues();
        const Eigen::Matrix3d vectors=canonicalEigenvectorSigns(eigen.eigenvectors());
        if(!eigenvectorSignPivotsAreSeparated(vectors))
            return fail("hyperbolic-cylinder eigenvector sign-pivot boundary is not differentiable in this chart");
        int nullIndex=0;values.cwiseAbs().minCoeff(&nullIndex);
        const double eigenScale=std::max(1.0,values.cwiseAbs().maxCoeff());
        for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)
            if(std::abs(values(i)-values(j))<=1e-10*eigenScale)
                return fail("hyperbolic-cylinder source derivative requires a simple, separated eigenspectrum");
        std::array<int,2> nonzero{};int cursor=0;
        for(int i=0;i<3;++i)if(i!=nullIndex)nonzero[cursor++]=i;
        int longitudinalPosition=-1;
        std::array<double,2> centers{};
        for(int position=0;position<2;++position) {
            const int index=nonzero[position];
            const double component=q.linearVector().dot(vectors.col(index));
            centers[position]=-component/(2.0*values(index));
        }
        const double shifted=q.constant()-
            std::pow(q.linearVector().dot(vectors.col(nonzero[0])),2)/(4.0*values(nonzero[0]))-
            std::pow(q.linearVector().dot(vectors.col(nonzero[1])),2)/(4.0*values(nonzero[1]));
        for(int position=0;position<2;++position)
            if(-shifted/values(nonzero[position])>0.0)longitudinalPosition=position;
        if(longitudinalPosition<0 || std::abs(shifted)<=1e-12*eigenScale)
            return fail("hyperbolic-cylinder branch scales are singular");
        const int transversePosition=1-longitudinalPosition;
        const int longitudinalIndex=nonzero[longitudinalPosition];
        const int transverseIndex=nonzero[transversePosition];
        const Eigen::Vector3d longitudinal=vectors.col(longitudinalIndex);
        const Eigen::Vector3d transverse=vectors.col(transverseIndex);
        const Eigen::Vector3d nullAxis=vectors.col(nullIndex);
        const double longitudinalCenter=centers[longitudinalPosition];
        const double transverseCenter=centers[transversePosition];
        const double longitudinalEigenvalue=values(longitudinalIndex);
        const double transverseEigenvalue=values(transverseIndex);
        const double longitudinalScaleSquared=-shifted/longitudinalEigenvalue;
        const double transverseScaleSquared=shifted/transverseEigenvalue;
        if(!(longitudinalScaleSquared>0.0) || !(transverseScaleSquared>0.0))
            return fail("hyperbolic-cylinder real branch scales are invalid");
        const double longitudinalScale=std::sqrt(longitudinalScaleSquared);
        const double transverseScale=std::sqrt(transverseScaleSquared);

        out.triangulation=triangulateQuadric(q,design,settings);
        if(out.triangulation.type!=QuadricSurfaceType::HyperbolicCylinder)
            return fail("native construction changed hyperbolic-cylinder type");
        const auto& mesh=out.triangulation.closedBoundary;
        const double radius=out.triangulation.truncationRadius;
        const int resolution=settings.hyperbolicCylinderGridResolution;
        const int side=resolution+1;
        const std::size_t perBranch=static_cast<std::size_t>(side)*static_cast<std::size_t>(side);
        const int expectedInterfaceVertices=static_cast<int>(2*perBranch);
        if(out.triangulation.interfaceVertexCount!=expectedInterfaceVertices ||
           mesh.vertices.rows()!=expectedInterfaceVertices+2 ||
           out.triangulation.boundaryLoopCount!=2)
            return fail("unsupported native hyperbolic-cylinder grid/cap layout");

        out.coefficientBasis=tangent.basis;
        out.vertexVelocity.resize(static_cast<std::size_t>(mesh.vertices.rows()));
        out.provenance.resize(static_cast<std::size_t>(mesh.vertices.rows()));
        for(auto& velocity:out.vertexVelocity)velocity=Eigen::MatrixXd::Zero(3,7);
        for(int column=0;column<7;++column) {
            const auto dq=tangent.basis.col(column).eval();
            if(!tangent.contains(dq,1e-9))
                return fail("basis column violates hyperbolic-cylinder constraints");
            const Eigen::Matrix3d dA=quadraticDirection(dq);
            const Eigen::Vector3d db=dq.segment<3>(1);
            std::array<Eigen::Vector3d,3> de;
            std::array<double,3> dlambda{};
            for(int i=0;i<3;++i) {
                de[i].setZero();
                dlambda[i]=vectors.col(i).dot(dA*vectors.col(i));
                for(int j=0;j<3;++j)if(j!=i) {
                    const double denominator=values(i)-values(j);
                    if(std::abs(denominator)<=1e-10*eigenScale)
                        return fail("hyperbolic-cylinder eigenvector differential is ill-conditioned");
                    de[i]+=vectors.col(j)*
                        (vectors.col(j).dot(dA*vectors.col(i))/denominator);
                }
            }
            const double bl=q.linearVector().dot(longitudinal);
            const double bt=q.linearVector().dot(transverse);
            const double dbl=db.dot(longitudinal)+q.linearVector().dot(de[longitudinalIndex]);
            const double dbt=db.dot(transverse)+q.linearVector().dot(de[transverseIndex]);
            const double dl=dlambda[longitudinalIndex];
            const double dt=dlambda[transverseIndex];
            const double dLongitudinalCenter=-dbl/(2.0*longitudinalEigenvalue)+
                bl*dl/(2.0*longitudinalEigenvalue*longitudinalEigenvalue);
            const double dTransverseCenter=-dbt/(2.0*transverseEigenvalue)+
                bt*dt/(2.0*transverseEigenvalue*transverseEigenvalue);
            const double dShift=dq(0)-bl*dbl/(2.0*longitudinalEigenvalue)+
                bl*bl*dl/(4.0*longitudinalEigenvalue*longitudinalEigenvalue)-
                bt*dbt/(2.0*transverseEigenvalue)+
                bt*bt*dt/(4.0*transverseEigenvalue*transverseEigenvalue);
            const double dLongitudinalScale=0.5/longitudinalScale*
                (-dShift/longitudinalEigenvalue+
                 shifted*dl/(longitudinalEigenvalue*longitudinalEigenvalue));
            const double dTransverseScale=0.5/transverseScale*
                (dShift/transverseEigenvalue-
                 shifted*dt/(transverseEigenvalue*transverseEigenvalue));
            for(int branch=0;branch<2;++branch) {
                const double branchSign=branch==0?-1.0:1.0;
                for(int row=0;row<=resolution;++row) {
                    const double y=-2.0*radius+4.0*radius*row/resolution;
                    const double g=std::sqrt(1.0+(y/transverseScale)*(y/transverseScale));
                    const double longitudinalOffset=branchSign*longitudinalScale*g;
                    const double dLongitudinalOffset=branchSign*(dLongitudinalScale*g-
                        longitudinalScale*y*y*dTransverseScale/
                        (transverseScale*transverseScale*transverseScale*g));
                    for(int col=0;col<=resolution;++col) {
                        const double z=-2.0*radius+4.0*radius*col/resolution;
                        const int vertex=branch*static_cast<int>(perBranch)+row*side+col;
                        const Eigen::Vector3d expected=
                            longitudinal*(longitudinalCenter+longitudinalOffset)+
                            transverse*(transverseCenter+y)+nullAxis*z;
                        if((expected-mesh.vertices.row(vertex).transpose()).norm()>
                           1e-10*std::max(1.0,expected.norm()))
                            return fail("hyperbolic-cylinder source point disagrees with native grid construction");
                        out.vertexVelocity[static_cast<std::size_t>(vertex)].col(column)=
                            de[longitudinalIndex]*(longitudinalCenter+longitudinalOffset)+
                            longitudinal*(dLongitudinalCenter+dLongitudinalOffset)+
                            de[transverseIndex]*(transverseCenter+y)+
                            transverse*dTransverseCenter+de[nullIndex]*z;
                        out.provenance[static_cast<std::size_t>(vertex)]=
                            "hyperbolic-cylinder-branch:"+std::to_string(branch)+
                            ":grid:"+std::to_string(row)+":"+std::to_string(col);
                    }
                }
            }
            for(int cap=0;cap<2;++cap) {
                const int capVertex=expectedInterfaceVertices+cap;
                std::set<int> loopVertices;
                for(int face=out.triangulation.interfaceFaceCount;
                    face<mesh.faces.rows();++face) {
                    bool contains=false;
                    for(int corner=0;corner<3;++corner)
                        contains|=mesh.faces(face,corner)==capVertex;
                    if(contains)for(int corner=0;corner<3;++corner) {
                        const int vertex=mesh.faces(face,corner);
                        if(vertex!=capVertex)loopVertices.insert(vertex);
                    }
                }
                if(loopVertices.size()!=static_cast<std::size_t>(4*resolution))
                    return fail("hyperbolic-cylinder cap does not identify one complete boundary loop");
                Eigen::Vector3d center=Eigen::Vector3d::Zero();
                Eigen::Vector3d velocity=Eigen::Vector3d::Zero();
                for(const int vertex:loopVertices) {
                    center+=mesh.vertices.row(vertex).transpose();
                    velocity+=out.vertexVelocity[static_cast<std::size_t>(vertex)].col(column);
                }
                center/=static_cast<double>(loopVertices.size());
                velocity/=static_cast<double>(loopVertices.size());
                const double branchCoordinate=(center-longitudinal*longitudinalCenter-
                    transverse*transverseCenter).dot(longitudinal);
                const double branchSign=branchCoordinate<0.0?-1.0:1.0;
                const Eigen::Vector3d expected=center+branchSign*4.0*radius*longitudinal;
                if((expected-mesh.vertices.row(capVertex).transpose()).norm()>
                   1e-10*std::max(1.0,expected.norm()))
                    return fail("hyperbolic-cylinder cap vertex disagrees with native centroid construction");
                out.vertexVelocity[static_cast<std::size_t>(capVertex)].col(column)=
                    velocity+branchSign*4.0*radius*de[longitudinalIndex];
                out.provenance[static_cast<std::size_t>(capVertex)]=
                    "hyperbolic-cylinder-cap-centroid:"+std::to_string(cap);
            }
        }
        std::ostringstream identity;identity.imbue(std::locale::classic());
        identity<<"NativeHyperbolicCylinderStratumSourceMotion-v1;"<<q.hash()<<';'
                <<out.triangulation.hash<<';'<<std::hexfloat;
        for(int i=0;i<out.coefficientBasis.rows();++i)
            for(int j=0;j<out.coefficientBasis.cols();++j)
                identity<<out.coefficientBasis(i,j)<<';';
        for(std::size_t i=0;i<out.vertexVelocity.size();++i) {
            identity<<out.provenance[i]<<';';
            for(int a=0;a<3;++a)for(int j=0;j<7;++j)
                identity<<out.vertexVelocity[i](a,j)<<';';
        }
        out.defined=true;
        out.diagnostic="reduced hyperbolic-cylinder source motion; fixed eigenbasis ordering, rectangular grid and exterior cap branches";
        out.scientificHash=manufacturing::computeContentSha256(identity.str());
        return out;
    }
    if(tangent.subtype!=QuadricSurfaceSubtype::EllipticCylinder)
        return fail("RANK_DEFICIENT_STRATUM_UNSUPPORTED:"+
                    std::string(quadricSurfaceSubtypeToken(tangent.subtype)));
    if(!tangent.defined || tangent.tangentDimension!=7 ||
       tangent.basis.rows()!=10 || tangent.basis.cols()!=7)
        return fail("elliptic-cylinder tangent chart is not seven-dimensional");

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(q.quadraticMatrix());
    if(eigen.info()!=Eigen::Success)return fail("cylinder eigensystem failed");
    const auto values=eigen.eigenvalues();
    const Eigen::Matrix3d eigenvectors=canonicalEigenvectorSigns(eigen.eigenvectors());
    if(!eigenvectorSignPivotsAreSeparated(eigenvectors))
        return fail("cylinder eigenvector sign-pivot boundary is not differentiable in this chart");
    int nullIndex=0;values.cwiseAbs().minCoeff(&nullIndex);
    const double eigenScale=std::max(1.0,values.cwiseAbs().maxCoeff());
    for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)
        if(std::abs(values(i)-values(j))<=1e-10*eigenScale)
            return fail("cylinder source derivative requires a simple, separated eigenspectrum");

    out.triangulation=triangulateQuadric(q,design,settings);
    if(out.triangulation.type!=QuadricSurfaceType::EllipticCylinder)
        return fail("native construction changed elliptic-cylinder type");
    const auto& mesh=out.triangulation.closedBoundary;
    const double radius=out.triangulation.truncationRadius;
    const Eigen::Vector3d b=q.linearVector();
    const double constant=q.constant();
    std::array<int,2> curved{};int cursor=0;
    for(int i=0;i<3;++i)if(i!=nullIndex)curved[cursor++]=i;
    const Eigen::Vector3d n=eigenvectors.col(nullIndex);
    const Eigen::Vector3d u=eigenvectors.col(curved[0]);
    const Eigen::Vector3d v=eigenvectors.col(curved[1]);
    const double au=values(curved[0]),av=values(curved[1]);
    if(std::abs(au)<=1e-12*eigenScale||std::abs(av)<=1e-12*eigenScale)
        return fail("cylinder cross-section curvature is rank deficient");
    const double bu=b.dot(u),bv=b.dot(v);
    const double centerU=-bu/(2*au),centerV=-bv/(2*av);
    const double shifted=constant-bu*bu/(4*au)-bv*bv/(4*av);
    const double r2u=-shifted/au,r2v=-shifted/av;
    if(!(r2u>0&&r2v>0)||!std::isfinite(r2u)||!std::isfinite(r2v))
        return fail("elliptic-cylinder cross-section is not real and bounded");
    const double ru=std::sqrt(r2u),rv=std::sqrt(r2v);
    const int resolution=settings.cylinderAngularResolution;
    if(mesh.vertices.rows()!=2*resolution+2 ||
       out.triangulation.interfaceVertexCount!=2*resolution)
        return fail("unsupported native cylinder vertex layout");

    out.coefficientBasis=tangent.basis;
    out.vertexVelocity.resize(static_cast<std::size_t>(mesh.vertices.rows()));
    out.provenance.resize(static_cast<std::size_t>(mesh.vertices.rows()));
    for(int vertex=0;vertex<mesh.vertices.rows();++vertex) {
        out.vertexVelocity[static_cast<std::size_t>(vertex)]=
            Eigen::MatrixXd::Zero(3,tangent.tangentDimension);
    }
    for(int column=0;column<tangent.tangentDimension;++column) {
        const auto dq=tangent.basis.col(column).eval();
        if(!tangent.contains(dq,1e-9))return fail("basis column violates cylinder constraints");
        const Eigen::Matrix3d dA=quadraticDirection(dq);
        const Eigen::Vector3d db=dq.segment<3>(1);
        std::array<Eigen::Vector3d,3> de;
        for(int i=0;i<3;++i) {
            de[i].setZero();
            for(int j=0;j<3;++j)if(j!=i) {
                const double denominator=values(i)-values(j);
                if(std::abs(denominator)<=1e-10*eigenScale)
                    return fail("eigenvector differential is ill-conditioned");
                de[i]+=eigenvectors.col(j)*
                    (eigenvectors.col(j).dot(dA*eigenvectors.col(i))/denominator);
            }
        }
        const double dau=eigenvectors.col(curved[0]).dot(
            dA*eigenvectors.col(curved[0]));
        const double dav=eigenvectors.col(curved[1]).dot(
            dA*eigenvectors.col(curved[1]));
        const double dbu=db.dot(u)+b.dot(de[curved[0]]);
        const double dbv=db.dot(v)+b.dot(de[curved[1]]);
        const double dcenterU=-dbu/(2*au)+bu*dau/(2*au*au);
        const double dcenterV=-dbv/(2*av)+bv*dav/(2*av*av);
        const double dshift=dq(0)-bu*dbu/(2*au)+bu*bu*dau/(4*au*au)
                                  -bv*dbv/(2*av)+bv*bv*dav/(4*av*av);
        const double dru=.5/ru*(-dshift/au+shifted*dau/(au*au));
        const double drv=.5/rv*(-dshift/av+shifted*dav/(av*av));
        for(int layer=0;layer<2;++layer)for(int segment=0;segment<resolution;++segment) {
            const int index=layer*resolution+segment;
            const double angle=2*std::acos(-1.0)*(segment+.5)/resolution;
            const double c=std::cos(angle),s=std::sin(angle);
            const double height=layer==0?-2*radius:2*radius;
            const Eigen::Vector3d derivative=
                de[curved[0]]*(centerU+ru*c)+u*(dcenterU+dru*c)+
                de[curved[1]]*(centerV+rv*s)+v*(dcenterV+drv*s)+de[nullIndex]*height;
            const Eigen::Vector3d expected=u*(centerU+ru*c)+v*(centerV+rv*s)+n*height;
            if((expected-mesh.vertices.row(index).transpose()).norm()>
               1e-11*std::max(1.0,expected.norm()))
                return fail("cylinder source point disagrees with native construction");
            out.vertexVelocity[static_cast<std::size_t>(index)].col(column)=derivative;
        }
        const int lower=2*resolution,upper=lower+1;
        out.vertexVelocity[static_cast<std::size_t>(lower)].col(column)=
            de[curved[0]]*centerU+u*dcenterU+de[curved[1]]*centerV+v*dcenterV-de[nullIndex]*(2*radius);
        out.vertexVelocity[static_cast<std::size_t>(upper)].col(column)=
            de[curved[0]]*centerU+u*dcenterU+de[curved[1]]*centerV+v*dcenterV+de[nullIndex]*(2*radius);
    }
    for(int i=0;i<mesh.vertices.rows();++i)
        out.provenance[static_cast<std::size_t>(i)]=
            i<2*resolution?"cylinder-ring:"+std::to_string(i):"cylinder-cap-center:"+std::to_string(i-2*resolution);

    std::ostringstream identity;identity.imbue(std::locale::classic());
    identity<<"NativeEllipticCylinderStratumSourceMotion-v1;"<<q.hash()<<';'
            <<out.triangulation.hash<<';'<<std::hexfloat;
    for(int i=0;i<out.coefficientBasis.rows();++i)
        for(int j=0;j<out.coefficientBasis.cols();++j)identity<<out.coefficientBasis(i,j)<<';';
    for(std::size_t i=0;i<out.vertexVelocity.size();++i) {
        identity<<out.provenance[i]<<';';
        for(int a=0;a<3;++a)for(int j=0;j<tangent.tangentDimension;++j)
            identity<<out.vertexVelocity[i](a,j)<<';';
    }
    out.defined=true;
    out.diagnostic="reduced elliptic-cylinder stratum source motion; fixed triangulation/eigenvalue ordering and truncation radius branches";
    out.scientificHash=manufacturing::computeContentSha256(identity.str());
    return out;
}
}
