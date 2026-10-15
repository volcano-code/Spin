#include "spin/CanonicalFullRankHollowing.h"
#include <algorithm>
#include <array>
#include <numeric>
#include <tuple>

namespace spin {
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity, const TriangleMesh& material,
    const Eigen::VectorXi& births, double density, bool self, bool sourceAudit, bool materialAudit,
    const ExactBooleanVertexCoordinates& exactMaterial,
    const ExactBooleanVertexCoordinates& exactSource,
    const ExactBooleanVertexCoordinates& exactDesign) {
    return differentiateCanonicalPolyhedralMaterialMoments(design,source,velocity,material,births,
        density,self,sourceAudit,materialAudit,exactMaterial,exactSource,exactDesign,{});
}
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity, const TriangleMesh& material,
    const Eigen::VectorXi& births, double density, bool self, bool sourceAudit, bool materialAudit,
    const ExactBooleanVertexCoordinates& exactMaterial,
    const ExactBooleanVertexCoordinates& exactSource,
    const ExactBooleanVertexCoordinates& exactDesign,
    const std::function<void(const std::string&)>& diagnosticObserver) {
    // The producer already orders by rational coordinates and oriented facets.
    // A second binary64 ordering here would lose that authoritative identity.
    return differentiatePolyhedralMaterialMoments(design,source,velocity,material,births,
        density,self,sourceAudit,materialAudit,exactMaterial,exactSource,exactDesign,diagnosticObserver);
}
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density) {
    return differentiateCanonicalPolyhedralMaterialMoments(design,source,velocity,material,birthFaces,density,false);
}
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection) {
    return differentiateCanonicalPolyhedralMaterialMoments(design,source,velocity,
        material,birthFaces,density,exactBinary64SelfIntersection,false);
}
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary) {
    return differentiateCanonicalPolyhedralMaterialMoments(design,source,velocity,material,
        birthFaces,density,exactBinary64SelfIntersection,exactBinary64SourceBoundary,false);
}
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary) {
    return differentiateCanonicalPolyhedralMaterialMoments(design,source,velocity,material,
        birthFaces,density,exactBinary64SelfIntersection,exactBinary64SourceBoundary,
        exactBinary64MaterialBoundary,{});
}
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary,const ExactBooleanVertexCoordinates& exactVertices) {
    // Check before indexing/sorting. The original differential performs all
    // remaining strict boundary, source-motion and provenance audits afterward.
    if(material.vertices.cols()!=3 || material.faces.cols()!=3 ||
       material.vertices.rows()<4 || material.faces.rows()<4 ||
       !material.vertices.allFinite() || birthFaces.size()!=material.faces.rows() ||
       (material.faces.array()<0).any() ||
       (material.faces.array()>=material.vertices.rows()).any() ||
       (!exactVertices.empty() && exactVertices.size()!=static_cast<std::size_t>(material.vertices.rows()))) {
        QuadricDiscreteMomentDifferential out;
        out.diagnostic="invalid canonical material mesh or provenance shape";
        return out;
    }
    TriangleMesh ordered;
    ordered.vertices.resizeLike(material.vertices);
    ordered.faces.resizeLike(material.faces);
    std::vector<int> order(material.vertices.rows()), remap(order.size());
    std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](int a,int b){
        for(int c=0;c<3;++c) {
            if(material.vertices(a,c)<material.vertices(b,c))return true;
            if(material.vertices(b,c)<material.vertices(a,c))return false;
        }
        return a<b;
    });
    ExactBooleanVertexCoordinates orderedExact(exactVertices.size());
    for(std::size_t i=0;i<order.size();++i) {
        remap[order[i]]=static_cast<int>(i);
        for(int c=0;c<3;++c) {
            const double value=material.vertices(order[i],c);
            ordered.vertices(i,c)=value==0.0?0.0:value;
        }
        if(!exactVertices.empty())orderedExact[i]=exactVertices[order[i]];
    }
    struct Face { std::array<int,3> vertices; int birth; };
    std::vector<Face> faces;
    faces.reserve(material.faces.rows());
    for(int f=0;f<material.faces.rows();++f) {
        std::array<int,3> ids={{remap[material.faces(f,0)],
            remap[material.faces(f,1)],remap[material.faces(f,2)]}};
        const int first=static_cast<int>(std::min_element(ids.begin(),ids.end())-ids.begin());
        faces.push_back({{{ids[first],ids[(first+1)%3],ids[(first+2)%3]}},birthFaces(f)});
    }
    std::sort(faces.begin(),faces.end(),[](const Face& a,const Face& b){
        return std::tie(a.vertices,a.birth)<std::tie(b.vertices,b.birth);
    });
    Eigen::VectorXi provenance(birthFaces.size());
    for(std::size_t f=0;f<faces.size();++f) {
        for(int c=0;c<3;++c)ordered.faces(f,c)=faces[f].vertices[c];
        provenance(f)=faces[f].birth;
    }
    return differentiatePolyhedralMaterialMoments(
        design,source,velocity,ordered,provenance,density,exactBinary64SelfIntersection,
        exactBinary64SourceBoundary,exactBinary64MaterialBoundary,orderedExact);
}
}
