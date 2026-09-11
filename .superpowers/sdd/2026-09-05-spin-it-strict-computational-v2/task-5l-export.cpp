#include "spin2014/SourceCageFem.h"
#include "source_surface_test_utils.h"
#include <iomanip>
#include <iostream>
template<class Matrix> void matrix(const Matrix& m) {
    std::cout<<'[';
    for(int i=0;i<m.rows();++i) {
        if(i)std::cout<<',';std::cout<<'[';
        for(int j=0;j<m.cols();++j) {if(j)std::cout<<',';std::cout<<m(i,j);}
        std::cout<<']';
    }
    std::cout<<']';
}
int main(int argc,char** argv) {
    try {
        if(argc!=2 || (std::string(argv[1])!="voronoi" && std::string(argv[1])!="barycentric"))
            throw std::invalid_argument("mass rule required");
        auto mesh=surface_test::cube(.3,.7);
        const Eigen::Vector3d origin(-.3,-.2,-.1);
        mesh.v=(1.7*mesh.v).rowwise()+origin.transpose();
        spin2014::SourceSurfaceCageConfig cc;cc.cageLevel=2;cc.rootOrigin=origin;cc.rootSide=1.7;
        spin2014::SourceCageGridConfig gc;gc.fineLevel=3;
        const spin2014::SourceCageGrid grid(spin2014::SourceSurfaceCage(mesh.v,mesh.f,cc),gc);
        spin2014::SourceCageFemConfig fc;
        if(std::string(argv[1])=="barycentric")fc.massRule=spin2014::SourceCageFemMassRule::BARYCENTRIC_TET;
        const spin2014::SourceCageFem fem(grid,fc);
        std::cout<<std::setprecision(17)<<"{\"scope\":\"SYNTHETIC_BOUNDARY_ONLY\",\"rule\":\""<<argv[1]
                 <<"\",\"fem_hash\":\""<<fem.outputSha256()<<"\",\"grid_hash\":\""<<grid.outputSha256()
                 <<"\",\"nodes\":";matrix(grid.nodePositions());
        std::cout<<",\"keys\":[";bool first=true;
        for(const auto& key:grid.nodeKeys()) {if(!first)std::cout<<',';first=false;std::cout<<'['<<key[0]<<','<<key[1]<<','<<key[2]<<']';}
        std::cout<<"],\"handles\":[";first=true;
        for(auto n:grid.handleNodeIndices()) {if(!first)std::cout<<',';first=false;std::cout<<n;}
        std::cout<<"],\"boundary\":[";first=true;
        for(const auto& row:fem.boundaryConstraints()) {
            if(!first)std::cout<<',';first=false;std::cout<<"{\"node\":"<<row.nodeIndex<<",\"entries\":[";bool inner=true;
            for(const auto& e:row.entries) {if(!inner)std::cout<<',';inner=false;std::cout<<'['<<e.handleIndex<<','<<e.weight<<']';}
            std::cout<<"]}";
        }
        std::cout<<"]}\n";return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
