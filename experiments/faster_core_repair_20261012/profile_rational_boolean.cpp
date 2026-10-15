// Diagnostic only: same Epeck operands and libigl operation as the native wrapper.
#include "verification/ExactRationalCoordinates.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <locale>
#include <stdexcept>
using Matrix=Eigen::Matrix<CGAL::Epeck::FT,Eigen::Dynamic,Eigen::Dynamic>;
static const auto start=std::chrono::steady_clock::now();
void stage(const char* name) {
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"{\"stage\":\""<<name<<"\",\"seconds\":"<<seconds<<"}\n"<<std::flush;
}
void read(std::istream& input,Matrix& vertices,Eigen::MatrixXi& faces) {
    int n,m;if(!(input>>n>>m)||n<=0||m<=0||n>2000000||m>4000000)
        throw std::invalid_argument("invalid operand dimensions");
    vertices.resize(n,3);faces.resize(m,3);
    for(int i=0;i<n;++i)for(int j=0;j<3;++j) {
        std::string value;if(!(input>>value))throw std::invalid_argument("missing coordinate");
        vertices(i,j)=CGAL::Epeck::FT(verification::parseExactRationalCoordinate(value));
    }
    for(int i=0;i<m;++i)for(int j=0;j<3;++j)
        if(!(input>>faces(i,j))||faces(i,j)<0||faces(i,j)>=n)
            throw std::invalid_argument("invalid face index");
}
int main(int argc,char** argv) {
    if(argc!=3)throw std::invalid_argument("operands and output paths required");
    std::ifstream input(argv[1]);Matrix a,b,result;Eigen::MatrixXi fa,fb,faces;Eigen::VectorXi birth;
    stage("parse-operands-begin");read(input,a,fa);read(input,b,fb);stage("parse-operands-end");
    stage("mesh-boolean-begin");
    const bool succeeded=igl::copyleft::cgal::mesh_boolean(a,fa,b,fb,igl::MESH_BOOLEAN_TYPE_MINUS,result,faces,birth);
    stage("mesh-boolean-end");
    if(!succeeded)throw std::runtime_error("actual libigl Boolean failed");
    stage("serialize-result-begin");std::ofstream output(argv[2]);output.imbue(std::locale::classic());
    output<<"{\"succeeded\":true,\"vertices\":[";
    for(int i=0;i<result.rows();++i) {
        if(i)output<<',';output<<'[';
        for(int j=0;j<3;++j){if(j)output<<',';output<<'\"'<<CGAL::exact(result(i,j))<<'\"';}
        output<<']';
    }
    output<<"],\"faces\":[";
    for(int i=0;i<faces.rows();++i){if(i)output<<',';output<<'['<<faces(i,0)<<','<<faces(i,1)<<','<<faces(i,2)<<']';}
    output<<"],\"birth_faces\":[";
    for(int i=0;i<birth.size();++i){if(i)output<<',';output<<birth(i);}
    output<<"]}\n";output.close();if(!output)throw std::runtime_error("output persistence failed");
    stage("serialize-result-end");
}
