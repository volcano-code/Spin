// Read-only boundary arithmetic reproducer; no FEM/BBW completion authority.
#include "spin2014/SourceCageFem.h"
#include "source_surface_test_utils.h"
#include <iomanip>
#include <iostream>
using namespace spin2014;
int main() {
    auto mesh=surface_test::cube(.3,.7);
    const Eigen::Vector3d origin(-.3,-.2,-.1);
    mesh.v=(1.7*mesh.v).rowwise()+origin.transpose();
    SourceSurfaceCageConfig c;c.cageLevel=2;c.rootOrigin=origin;c.rootSide=1.7;
    SourceCageGridConfig gc;gc.fineLevel=3;
    const SourceCageGrid grid(SourceSurfaceCage(mesh.v,mesh.f,c),gc);
    using V=Eigen::Matrix<long double,3,1>;
    std::cout << std::setprecision(25);
    int disagreement=0;
    for(int t=0;t<grid.cage().boundaryFaces().rows();++t) {
        V p[3]; auto lo=grid.cage().boundaryVertexKeys()[grid.cage().boundaryFaces()(t,0)],hi=lo;
        for(int j=0;j<3;++j) {
            auto vertex=grid.cage().boundaryFaces()(t,j);
            p[j]=grid.cage().boundaryVertices().row(vertex).transpose().cast<long double>();
            for(int a=0;a<3;++a) {lo[a]=std::min(lo[a],grid.cage().boundaryVertexKeys()[vertex][a]);hi[a]=std::max(hi[a],grid.cage().boundaryVertexKeys()[vertex][a]);}
        }
        int axis=0;while(lo[axis]!=hi[axis])++axis;
        const int u=(axis+1)%3,v=(axis+2)%3;
        const auto cross=[&](const V& x,const V& y){return x[u]*y[v]-x[v]*y[u];};
        const long double det=cross(p[1]-p[0],p[2]-p[0]);
        for(int n=0;n<grid.nodePositions().rows();++n) {
            bool inside=true;for(int a=0;a<3;++a)inside&=grid.nodeKeys()[n][a]>=2*lo[a] && grid.nodeKeys()[n][a]<=2*hi[a];
            if(!inside)continue;
            const V x=grid.nodePositions().row(n).transpose().cast<long double>();
            const long double w1=cross(x-p[0],p[2]-p[0])/det,w2=cross(p[1]-p[0],x-p[0])/det;
            const long double w0=1-w1-w2;
            const long double e0=cross(p[1]-x,p[2]-x)/det,e1=cross(p[2]-x,p[0]-x)/det,e2=cross(p[0]-x,p[1]-x)/det;
            if((w0<0 || w1<0 || w2<0) && e0>=0 && e1>=0 && e2>=0) {
                ++disagreement;
                std::cout<<"node="<<n<<" triangle="<<t<<" point="<<x.transpose()<<" old="<<w0<<','<<w1<<','<<w2<<" oriented="<<e0<<','<<e1<<','<<e2<<'\n';
            }
        }
    }
    std::cout << "false_negative_triangle_tests=" << disagreement << '\n';
    return disagreement>0?0:1;
}
