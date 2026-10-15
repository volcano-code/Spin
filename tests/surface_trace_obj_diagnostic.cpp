#include "surface_trace_independent_oracle.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::runtime_error("requires triangle OBJ path");
        std::ifstream input(argv[1]);
        if(!input) throw std::runtime_error("missing input");
        std::vector<Eigen::Vector3d> v;
        std::vector<Eigen::Vector3i> f;
        std::string line;
        while(std::getline(input,line)) {
            std::istringstream row(line); std::string type,extra; row>>type;
            if(type.empty() || type[0]=='#') continue;
            if(type=="v") {
                Eigen::Vector3d p;
                if(!(row>>p[0]>>p[1]>>p[2]) || (row>>extra) || !p.allFinite())
                    throw std::runtime_error("invalid vertex");
                v.push_back(p);
            } else if(type=="f") {
                Eigen::Vector3i face;
                if(!(row>>face[0]>>face[1]>>face[2]) || (row>>extra) || face.minCoeff()<1)
                    throw std::runtime_error("requires positive integer triangular indices");
                f.push_back((face.array()-1).matrix());
            } else throw std::runtime_error("unsupported OBJ record");
        }
        if(v.empty() || f.empty()) throw std::runtime_error("empty input");
        Eigen::MatrixXd vertices(v.size(),3); Eigen::MatrixXi faces(f.size(),3);
        for(std::size_t i=0;i<v.size();++i) vertices.row(i)=v[i];
        for(std::size_t i=0;i<f.size();++i) faces.row(i)=f[i];
        // Normalize residual units only: rounding transformed coordinates can
        // collapse valid near-coincident vertices. Keep the original bytes.
        const Eigen::RowVector3d low=vertices.colwise().minCoeff(), high=vertices.colwise().maxCoeff();
        const Eigen::RowVector3d center=Eigen::RowVector3d::Zero();
        const double length=(high-low).maxCoeff();
        if(!std::isfinite(length) || length<=0) throw std::runtime_error("invalid input scale");
        surface_trace_oracle::verifyCubature();
        const auto errors=surface_trace_oracle::verifyPrefixes(vertices,faces,1.,length);
        std::cout<<std::setprecision(17)<<"{\"vertices\":"<<v.size()<<",\"faces\":"<<errors.faces
            <<",\"max_moment_error\":"<<errors.moment<<",\"max_derivative_error\":"<<errors.derivative
            <<",\"input_scale\":"<<length<<",\"normalization_origin\":["<<center[0]<<","<<center[1]<<","<<center[2]
            <<"],\"density_in_normalized_check\":1,\"scope\":\"RECONSTRUCTED_INPUT_PREFIX_NUMERICS_ONLY\"}\n";
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
