// Controller-only export of actual public API values for an independent Python
// reintegration. Not a solver/demo or publication path; stdout is cage geometry.
#include "spin2014/SourceSurfaceCage.h"
#include <iomanip>
#include <iostream>

template<class Matrix> void array(const Matrix& a) {
    std::cout << '[';
    for(Eigen::Index i=0;i<a.rows();++i) {
        if(i) std::cout << ',';
        std::cout << '[';
        for(Eigen::Index j=0;j<a.cols();++j) {
            if(j) std::cout << ',';
            std::cout << a(i,j);
        }
        std::cout << ']';
    }
    std::cout << ']';
}
int main() {
    try {
        Eigen::MatrixXd v(8,3);
        for(int i=0;i<8;++i) for(int axis=0;axis<3;++axis)
            v(i,axis)=((i>>axis)&1)?.75:.25;
        Eigen::MatrixXi f(12,3);
        f << 0,2,1, 1,2,3, 4,5,6, 5,7,6,
             0,4,2, 2,4,6, 1,3,5, 3,7,5,
             0,1,4, 1,5,4, 2,6,3, 3,6,7;
        spin2014::SourceSurfaceCageConfig config;
        config.cageLevel=3;
        const spin2014::SourceSurfaceCage cage(v,f,config);
        std::cout << std::setprecision(17)
                  << "{\"scope\":\"SYNTHETIC_CAGE_GEOMETRY_ONLY_NOT_MATERIAL\","
                  << "\"source_sha256\":\"" << cage.sourceSha256()
                  << "\",\"config_sha256\":\"" << cage.configSha256()
                  << "\",\"output_sha256\":\"" << cage.outputSha256()
                  << "\",\"selected_count\":" << cage.selectedCells().size()
                  << ",\"handle_count\":" << cage.handleVertexIndices().size()
                  << ",\"source_nodes\":" << cage.sourceNodeCount()
                  << ",\"vertices\":";
        array(cage.boundaryVertices()); std::cout << ",\"faces\":";
        array(cage.boundaryFaces()); std::cout << "}\n";
        return std::cout.good()?0:1;
    } catch(const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
