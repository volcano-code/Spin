// Diagnostic actual API export for independent Python checks, not BBW evidence.
#include "spin2014/SourceCageGrid.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

template<class Matrix> void matrixJson(const Matrix& m) {
    std::cout << '[';
    for(Eigen::Index i=0;i<m.rows();++i) {
        if(i) std::cout << ',';
        std::cout << '[';
        for(Eigen::Index j=0;j<m.cols();++j) {
            if(j) std::cout << ',';
            std::cout << m(i,j);
        }
        std::cout << ']';
    }
    std::cout << ']';
}
template<class Range> void rowJson(const Range& row) {
    std::cout << '['; bool first=true;
    for(const auto& x:row) {if(!first) std::cout << ','; first=false; std::cout << x;}
    std::cout << ']';
}
template<class Range> void rowsJson(const Range& rows) {
    std::cout << '['; bool first=true;
    for(const auto& row:rows) {if(!first) std::cout << ','; first=false; rowJson(row);}
    std::cout << ']';
}
int main() {
    try {
        Eigen::MatrixXd v(8,3);
        for(int i=0;i<8;++i) for(int a=0;a<3;++a) v(i,a)=((i>>a)&1)?.75:.25;
        Eigen::MatrixXi f(12,3);
        f << 0,2,1,1,2,3,4,5,6,5,7,6,0,4,2,2,4,6,1,3,5,3,7,5,0,1,4,1,5,4,2,6,3,3,6,7;
        spin2014::SourceSurfaceCageConfig coarse;
        coarse.cageLevel=3;
        spin2014::SourceCageGridConfig fine;
        fine.fineLevel=4;
        const spin2014::SourceCageGrid grid(spin2014::SourceSurfaceCage(v,f,coarse),fine);
        Eigen::MatrixXd fields(grid.nodePositions().rows(),9);
        for(Eigen::Index i=0;i<fields.rows();++i) {
            const double x=grid.nodePositions()(i,0),y=grid.nodePositions()(i,1),z=grid.nodePositions()(i,2);
            fields.row(i) << 1,x,y,z,x*y,x*z,y*z,x*y*z,2;
        }
        Eigen::MatrixXd queries(6,3);
        queries << .3,.45,.6, .5,.45,.6, .5,.5,.5, .28125,.34375,.40625, .125,.125,.125, .875,.875,.875;
        const auto values=grid.interpolate(fields,queries);
        bool outsideRejected=false;
        try {grid.embed(Eigen::Vector3d(std::nextafter(.125,-std::numeric_limits<double>::infinity()),.5,.5));}
        catch(const std::invalid_argument&) {outsideRejected=true;}
        catch(const std::out_of_range&) {outsideRejected=true;}
        std::cout << std::setprecision(17) << "{\"scope\":\"SYNTHETIC_GRID_INTERPOLATION_NOT_BBW\",\"output_sha256\":\""
                  << grid.outputSha256() << "\",\"nodes\":";
        matrixJson(grid.nodePositions()); std::cout << ",\"cells\":"; rowsJson(grid.cellNodes());
        std::cout << ",\"edges\":"; rowsJson(grid.edges());
        std::cout << ",\"handle_nodes\":"; rowJson(grid.handleNodeIndices());
        std::cout << ",\"outside_rejected\":" << (outsideRejected?"true":"false") << ",\"queries\":[";
        for(Eigen::Index i=0;i<queries.rows();++i) {
            if(i) std::cout << ',';
            const auto embed=grid.embed(queries.row(i).transpose());
            std::cout << "{\"point\":[" << queries(i,0) << ',' << queries(i,1) << ',' << queries(i,2) << "],\"cell_key\":";
            rowJson(embed.cellKey); std::cout << ",\"values\":[";
            for(Eigen::Index j=0;j<values.cols();++j) {if(j) std::cout << ','; std::cout << values(i,j);}
            std::cout << "]}";
        }
        std::cout << "]}\n";
        return std::cout.good()?0:1;
    } catch(const std::exception& e) {std::cerr << e.what() << '\n';return 1;}
}
