#include <cstdint>
#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
    try {
        if(argc!=5) throw std::runtime_error("operation inputA inputB output");
        Eigen::MatrixXd a,b,v; Eigen::MatrixXi fa,fb,f; Eigen::VectorXi j;
        if(!igl::read_triangle_mesh(argv[2],a,fa)||!igl::read_triangle_mesh(argv[3],b,fb))
            throw std::runtime_error("invalid input file");
        std::string op=argv[1];
        igl::MeshBooleanType type;
        if(op=="union") type=igl::MESH_BOOLEAN_TYPE_UNION;
        else if(op=="difference") type=igl::MESH_BOOLEAN_TYPE_MINUS;
        else if(op=="intersect") type=igl::MESH_BOOLEAN_TYPE_INTERSECT;
        else throw std::runtime_error("unsupported operation");
        if(!igl::copyleft::cgal::mesh_boolean(a,fa,b,fb,type,v,f,j))
            throw std::runtime_error("CGAL Boolean failed; no fallback");
        std::ofstream out(argv[4]); out<<std::setprecision(17);
        for(int i=0;i<v.rows();++i) out<<"v "<<v(i,0)<<" "<<v(i,1)<<" "<<v(i,2)<<"\n";
        for(int i=0;i<f.rows();++i) out<<"f "<<f(i,0)+1<<" "<<f(i,1)+1<<" "<<f(i,2)+1<<"\n";
        if(!out) throw std::runtime_error("write failed");
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
