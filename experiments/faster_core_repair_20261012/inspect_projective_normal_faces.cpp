// Actual frozen sphere and inverse-Gauss roots, source diagnostics only.
#include "spin/QuadricInverseGaussMap.h"
#include "spin/detail/FixedNormalSphere.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
    if(argc!=11)throw std::invalid_argument("supply ten saved coefficients");
    spin::QuadricCoefficients10::Vector values;
    for(int i=0;i<10;++i)values(i)=std::stod(argv[i+1]);
    const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
    const spin::QuadricInverseGaussMap map(q);
    const auto sphere=spin::detail::fixedNormalSphere(2);
    std::cout<<std::setprecision(17)<<"{\"multiple_faces\":[";
    bool first=true;
    for(std::size_t f=0;f<sphere.faces.size();++f) {
        const auto& face=sphere.faces[f];
        std::array<spin::QuadricIdealEdgeResult,3> edges;
        int roots=0;bool single=true;
        for(int k=0;k<3;++k) {
            edges[k]=map.idealEdgeDifferential(sphere.vertices[face[k]],sphere.vertices[face[(k+1)%3]]);
            if(!edges[k].defined)throw std::runtime_error(edges[k].diagnostic);
            roots+=int(edges[k].points.size());single&=edges[k].points.size()<=1;
        }
        if(roots==0 || (roots==2 && single))continue;
        if(!first)std::cout<<',';first=false;
        std::cout<<"{\"face\":"<<f<<",\"corners\":[";
        for(int k=0;k<3;++k) {
            if(k)std::cout<<',';
            const auto& n=sphere.vertices[face[k]];
            std::cout<<"{\"id\":"<<face[k]<<",\"finite\":"
                <<(map.point(n).valid()?"true":"false")<<",\"normal\":["
                <<n.x()<<','<<n.y()<<','<<n.z()<<"]}";
        }
        std::cout<<"],\"edges\":[";
        for(int k=0;k<3;++k) {
            if(k)std::cout<<',';
            std::cout<<"{\"from\":"<<face[k]<<",\"to\":"<<face[(k+1)%3]<<",\"roots\":[";
            for(std::size_t j=0;j<edges[k].points.size();++j) {
                if(j)std::cout<<',';
                const auto& p=edges[k].points[j];
                std::cout<<"{\"parameter\":"<<p.parameter<<",\"direction\":["
                    <<p.direction.x()<<','<<p.direction.y()<<','<<p.direction.z()<<"]}";
            }
            std::cout<<"]}";
        }
        std::cout<<"]}";
    }
    std::cout<<"],\"mapping_evaluations\":0,\"optimization_admission\":false}\n";
}
