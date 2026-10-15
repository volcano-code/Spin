#include "verification/ExactBoundaryAudit.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

template<class T> void emit(std::ostream& out,const T& value) { out<<value; }
void emit(std::ostream& out,const std::string& value) { out<<std::quoted(value); }
template<class T,std::size_t N> void emit(std::ostream& out,const std::array<T,N>& values) {
    out<<'['; bool first=true;
    for(const auto& value:values) { if(!first) out<<','; first=false; emit(out,value); }
    out<<']';
}
template<class T> void emit(std::ostream& out,const std::vector<T>& values) {
    out<<'['; bool first=true;
    for(const auto& value:values) { if(!first) out<<','; first=false; emit(out,value); }
    out<<']';
}

int main(int argc,char** argv) {
    try {
        if(argc!=3) throw std::invalid_argument("boundary input and report output required");
        const auto begin=std::chrono::steady_clock::now();
        std::ifstream input(argv[1]);
        std::size_t nv,nf;
        if(!(input>>nv>>nf) || !nv || !nf || nv>2147483647 || nf>1073741823)
            throw std::invalid_argument("invalid boundary sizes");
        std::vector<std::array<std::string,3>> vertices(nv);
        for(auto& row:vertices) for(auto& value:row)
            if(!(input>>value)) throw std::invalid_argument("truncated rational vertices");
        Eigen::MatrixXi faces(nf,3);
        for(int i=0;i<faces.rows();++i) for(int k=0;k<3;++k)
            if(!(input>>faces(i,k))) throw std::invalid_argument("truncated faces");
        std::string excess;
        if(input>>excess) throw std::invalid_argument("unexpected trailing boundary data");
        const auto loaded=std::chrono::steady_clock::now();
        const auto checked=verification::auditExactTriangleBoundaryRational(
            vertices,faces,verification::ExactBoundaryOrientation::AnyNonzero);
        const auto completed=std::chrono::steady_clock::now();
        const auto& r=checked.facts();
        std::ofstream out(argv[2],std::ios::out);
        if(!out) throw std::runtime_error("cannot open report");
        out<<std::boolalpha<<std::setprecision(17);
        out<<"{\"scope\":\"EXACT_BOUNDARY_AUDIT_ONLY_NO_MATERIAL_ADMISSION\",\"load_seconds\":"
           <<std::chrono::duration<double>(loaded-begin).count()<<",\"audit_seconds\":"
           <<std::chrono::duration<double>(completed-loaded).count()<<",\"facts\":{";
        bool first=true;
        auto field=[&](const char* name,const auto& value) {
            if(!first) out<<','; first=false;
            out<<std::quoted(name)<<':'; emit(out,value);
        };
        field("algorithm",r.algorithm); field("backend",r.backend);
        field("inputDigest",r.inputDigest); field("orientation",r.orientation);
        field("intersectionError",r.intersectionError); field("boundaryValid",r.boundaryValid);
        field("closed",r.closed); field("edgeManifold",r.edgeManifold);
        field("vertexManifold",r.vertexManifold); field("consistentlyOriented",r.consistentlyOriented);
        field("componentOrientationsValid",r.componentOrientationsValid);
        field("intersectionCompleted",r.intersectionCompleted); field("intersectionFree",r.intersectionFree);
        field("exactPairTests",r.exactPairTests); field("unreferencedVertices",r.unreferencedVertices);
        field("degenerateFaces",r.degenerateFaces); field("nonManifoldVertices",r.nonManifoldVertices);
        field("boundaryEdges",r.boundaryEdges); field("nonManifoldEdges",r.nonManifoldEdges);
        field("inconsistentEdges",r.inconsistentEdges); field("intersectingPairs",r.intersectingPairs);
        field("faceComponents",r.faceComponents); field("componentVolumeSigns",r.componentVolumeSigns);
        field("componentSignedVolumesExact",r.componentSignedVolumesExact);
        out<<"}}\n";
        if(!out) throw std::runtime_error("report write failed");
        return r.intersectionCompleted?0:1;
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
