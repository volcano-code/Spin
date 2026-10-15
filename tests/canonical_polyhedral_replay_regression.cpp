#include "spin/PolyhedralMomentDifferential.h"
#include <Eigen/Geometry>
#include <iostream>
#include <iomanip>
#include <stdexcept>

namespace {
spin::TriangleMesh tetra(double s) {
 spin::TriangleMesh m;m.vertices.resize(4,3);m.faces.resize(4,3);
 m.vertices<<.13,.17,.19,1.13,.23,.31,.29,1.31,.37,.41,.43,1.47;
 m.vertices*=s;m.faces<<0,2,1,0,1,3,0,3,2,1,2,3;return m;
}
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
}
int main(){try {
 std::cout<<std::setprecision(17);
 auto source=tetra(1),design=tetra(8);design.vertices.rowwise()-=Eigen::RowVector3d(3,3,3);
 std::vector<spin::SourceVertexJacobian> velocity(4);
 for(int v=0;v<4;++v)for(int a=0;a<3;++a)for(int j=0;j<10;++j)
  velocity[v](a,j)=std::sin(.31*(v+1)*(a+2)*(j+1));
 Eigen::VectorXi births(4);births<<4,5,6,7;
 auto reference=spin::differentiatePolyhedralMaterialMoments(design,source,velocity,source,births,1.7);
 check(reference.defined,reference.diagnostic.c_str());
 for(int shift=1;shift<=2;++shift) {
  auto permuted=source;for(int f=0;f<4;++f)for(int a=0;a<3;++a)permuted.faces(f,a)=source.faces(f,(a+shift)%3);
  auto value=spin::differentiatePolyhedralMaterialMoments(design,source,velocity,permuted,births,1.7);
  check(value.defined,value.diagnostic.c_str());
  std::cout<<"cyclic_shift="<<shift<<" J_difference="<<(value.jacobian-reference.jacobian).norm()<<'\n';
  check(value.jacobian==reference.jacobian,"cyclic orientation-preserving triangle permutation changed exact Jacobian");
 }
 // Reverse both row order and vertex numbering, carrying birth labels. The
 // surface and orientation are unchanged, so this must be exactly replayable.
 auto permuted=source;Eigen::VectorXi permutedBirths(4);
 for(int v=0;v<4;++v)permuted.vertices.row(v)=source.vertices.row(3-v);
 for(int f=0;f<4;++f) {
  permutedBirths(f)=births(3-f);
  for(int a=0;a<3;++a)permuted.faces(f,a)=3-source.faces(3-f,(a+1)%3);
 }
 const auto relabeled=spin::differentiatePolyhedralMaterialMoments(design,source,velocity,permuted,permutedBirths,1.7);
 check(relabeled.defined && relabeled.jacobian==reference.jacobian,"material relabeling changed Jacobian or lost birth alignment");
 permuted.faces.col(0).swap(permuted.faces.col(2));
 check(!spin::differentiatePolyhedralMaterialMoments(design,source,velocity,permuted,permutedBirths,1.7).defined,"canonicalization reversed an inward mesh to force acceptance");
 for(int variant=0;variant<5;++variant) {
  auto bad=source;auto badBirths=births;
  if(variant==0)bad.faces(0,0)=-1;
  if(variant==1)bad.faces(0,0)=4;
  if(variant==2)bad.vertices(0,0)=std::numeric_limits<double>::quiet_NaN();
  if(variant==3)badBirths.resize(3);
  if(variant==4)bad.faces.row(0)=bad.faces.row(1);
  check(!spin::differentiatePolyhedralMaterialMoments(design,source,velocity,bad,badBirths,1.7).defined,"invalid material/provenance was silently repaired or indexed");
 }
 std::cout<<"[PASS] canonical polyhedral replay\n";return 0;
}catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<'\n';return 1;}}
