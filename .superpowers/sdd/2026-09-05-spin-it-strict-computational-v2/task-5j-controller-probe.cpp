// Export actual public native values for the independent rational checker.
// Diagnostic only: this does not solve BBW or grant paper-case acceptance.
#include "spin2014/SourceCageFem.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

template<class Range> void rowJson(const Range& r) {
    std::cout << '['; bool first=true;
    for(const auto& x:r) {if(!first) std::cout << ',';first=false;std::cout << x;}
    std::cout << ']';
}
template<class Matrix> void matrixJson(const Matrix& m) {
    std::cout << '[';
    for(Eigen::Index i=0;i<m.rows();++i) {
        if(i) std::cout << ',';
        std::cout << '[';
        for(Eigen::Index j=0;j<m.cols();++j) {if(j) std::cout << ',';std::cout << m(i,j);}
        std::cout << ']';
    }
    std::cout << ']';
}
void vectorJson(const Eigen::VectorXd& v) {
    std::cout << '[';
    for(Eigen::Index i=0;i<v.size();++i) {if(i) std::cout << ',';std::cout << v(i);}
    std::cout << ']';
}
int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::invalid_argument("explicit mass-rule argument required");
        spin2014::SourceCageFemConfig config;
        const std::string rule=argv[1];
        if(rule=="BARYCENTRIC_TET") config.massRule=spin2014::SourceCageFemMassRule::BARYCENTRIC_TET;
        else if(rule=="RESTRICTED_GRID_VORONOI") config.massRule=spin2014::SourceCageFemMassRule::RESTRICTED_GRID_VORONOI;
        else throw std::invalid_argument("unknown mass rule");
        Eigen::MatrixXd v(8,3);
        for(int i=0;i<8;++i) for(int a=0;a<3;++a) v(i,a)=((i>>a)&1)?.4:.3;
        Eigen::MatrixXi f(12,3);
        f << 0,2,1,1,2,3,4,5,6,5,7,6,0,4,2,2,4,6,1,3,5,3,7,5,0,1,4,1,5,4,2,6,3,3,6,7;
        spin2014::SourceSurfaceCageConfig coarse;coarse.cageLevel=2;
        spin2014::SourceCageGridConfig fine;fine.fineLevel=3;
        const spin2014::SourceCageGrid grid(spin2014::SourceSurfaceCage(v,f,coarse),fine);
        const spin2014::SourceCageFem fem(grid,config);
        Eigen::MatrixXd field(grid.nodePositions().rows(),1);
        for(Eigen::Index i=0;i<field.rows();++i) field(i,0)=((7*i)%11-5)/8.;
        const auto action=fem.applyQ(field);
        // Dense diagnostic export is deliberately bounded to exactly27 nodes.
        if(field.rows()!=27) throw std::runtime_error("unexpected diagnostic grid size");
        const Eigen::MatrixXd stiffness(fem.stiffness());
        std::cout << std::setprecision(17)
                  << "{\"scope\":\"SYNTHETIC_FEM_NOT_BBW_SOLVE\",\"mass_rule\":\"" << rule
                  << "\",\"config_sha256\":\"" << fem.configSha256()
                  << "\",\"output_sha256\":\"" << fem.outputSha256() << "\",\"nodes\":";
        matrixJson(grid.nodePositions());std::cout << ",\"tetrahedra\":[";
        bool first=true;
        for(const auto& t:fem.tetrahedra()) {if(!first) std::cout << ',';first=false;rowJson(t);}
        std::cout << "],\"mass\":";vectorJson(fem.lumpedMass());
        std::cout << ",\"stiffness\":";matrixJson(stiffness);
        std::cout << ",\"field\":";vectorJson(field.col(0));
        std::cout << ",\"action\":";vectorJson(action.col(0));
        std::cout << ",\"energy\":" << fem.energy(field) << ",\"handle_nodes\":";
        rowJson(grid.handleNodeIndices());std::cout << ",\"boundary\":[";first=true;
        for(const auto& boundary:fem.boundaryConstraints()) {
            if(!first) std::cout << ',';first=false;
            Eigen::VectorXd weights=Eigen::VectorXd::Zero(grid.handleNodeIndices().size());
            for(const auto& entry:boundary.entries) weights(entry.handleIndex)=entry.weight;
            std::cout << "{\"node\":" << boundary.nodeIndex << ",\"weights\":";
            vectorJson(weights);std::cout << '}';
        }
        std::cout << "]}\n";
        return std::cout.good()?0:1;
    } catch(const std::exception& e) {std::cerr << e.what() << '\n';return 1;}
}
