// Input construction candidate only; never grants mesh or optimization admission.
#include "verification/ExactRationalCoordinates.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <fstream>
#include <iostream>
#include <locale>
#include <stdexcept>
using Matrix=Eigen::Matrix<CGAL::Epeck::FT,Eigen::Dynamic,Eigen::Dynamic>;
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
    try {
        if(argc!=4)throw std::invalid_argument("operands output operation required");
        const std::string operation(argv[3]);igl::MeshBooleanType type;
        if(operation=="union")type=igl::MESH_BOOLEAN_TYPE_UNION;
        else if(operation=="difference")type=igl::MESH_BOOLEAN_TYPE_MINUS;
        else if(operation=="intersection")type=igl::MESH_BOOLEAN_TYPE_INTERSECT;
        else throw std::invalid_argument("unsupported operation");
        std::ifstream input(argv[1]);Matrix a,b,result;Eigen::MatrixXi fa,fb,faces;Eigen::VectorXi birth;
        read(input,a,fa);read(input,b,fb);std::string extra;
        if(input>>extra)throw std::invalid_argument("trailing operand data");
        if(!igl::copyleft::cgal::mesh_boolean(a,fa,b,fb,type,result,faces,birth))
            throw std::runtime_error("actual libigl Boolean failed");
        std::ofstream output(argv[2]);output.imbue(std::locale::classic());
        output<<"{\"succeeded\":true,\"operation\":\""<<operation<<"\",\"vertices\":[";
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
        return 0;
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 2; }
}
