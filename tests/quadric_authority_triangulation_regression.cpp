#include "spin/QuadricTriangulator.h"
#include "spin/QuadricCoefficientAuthority.h"
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>
#include <set>

int main() {
 try {
  using A=spin::QuadricCoefficientAuthority10;
  A::RationalVector q; q.fill("0/1");q[0]="-1/2";q[4]=q[5]=q[6]="1/2";
  A::RationalVector direction;direction.fill("0/1");direction[0]="1/1";
  const auto initial=A::fromNormalizedGeometry(q);
  const auto updated=initial.geodesicUpdate(direction,"1/10000000000000000000000000");
  Eigen::MatrixXd design(8,3);
  for(int i=0;i<8;++i)for(int a=0;a<3;++a)design(i,a)=(i&(1<<a))?.1:-.1;
  spin::QuadricTriangulationSettings settings;
  settings.icosphereSubdivisions=0;settings.maximumIcosphereSubdivisions=-1;
  settings.rationalSourceCoordinates=true;settings.compactRationalSourceCoordinates=true;
  const auto mesh=spin::triangulateQuadric(updated,design,settings);
  if(mesh.exactVertices.size()!=12 || mesh.closedBoundary.faces.rows()!=20)
   throw std::runtime_error("authority ellipsoid lost fixed normal connectivity");
  for(const auto& p:mesh.exactVertices) {
   mpq_class squared=0;for(const auto& s:p)squared+=mpq_class(s)*mpq_class(s);
   if(squared>=1)throw std::runtime_error("complete source used double coefficients");
  }
  if(!spin::rationalQuadricConstructionIdentityMatches(updated,mesh,design,settings))
   throw std::runtime_error("authority construction identity did not verify");
  auto explicitFixed=settings;
  explicitFixed.maximumIcosphereSubdivisions=explicitFixed.icosphereSubdivisions;
  const auto fixedMesh=spin::triangulateQuadric(updated,design,explicitFixed);
  if(fixedMesh.exactVertices!=mesh.exactVertices ||
     fixedMesh.closedBoundary.faces!=mesh.closedBoundary.faces ||
     !spin::rationalQuadricConstructionIdentityMatches(updated,fixedMesh,design,explicitFixed))
   throw std::runtime_error("explicit equal maximum did not preserve fixed authority construction");
  auto adaptive=settings;adaptive.maximumIcosphereSubdivisions=settings.icosphereSubdivisions+1;
  bool adaptiveRejected=false;
  try {spin::triangulateQuadric(updated,design,adaptive);}
  catch(const std::invalid_argument&){adaptiveRejected=true;}
  if(!adaptiveRejected)throw std::runtime_error("authority source silently allowed adaptive normal refinement");
  if(spin::rationalQuadricConstructionIdentityMatches(initial,mesh,design,settings))
   throw std::runtime_error("old coefficients reused the new source identity");
  auto wrongRecipe=settings;wrongRecipe.icosphereSubdivisions=1;
  if(spin::rationalQuadricConstructionIdentityMatches(updated,mesh,design,wrongRecipe))
   throw std::runtime_error("different fixed normal level reused source identity");
  wrongRecipe=settings;wrongRecipe.rationalSourceCoordinates=false;
  if(spin::rationalQuadricConstructionIdentityMatches(updated,mesh,design,wrongRecipe))
   throw std::runtime_error("double source policy reused authority construction");
  auto changed=mesh;changed.exactVertices[0][0]="0/1";
  if(spin::rationalQuadricConstructionIdentityMatches(updated,changed,design,settings))
   throw std::runtime_error("changed rational boundary retained its identity");
  for(const auto& constant:{"-1/2","1/2"}) {
   auto hq=q;hq[0]=constant;hq[6]="-1/2";
   const auto hyper=A::fromNormalizedGeometry(hq).geodesicUpdate(
       direction,"1/10000000000000000000000000");
   auto hs=settings;hs.icosphereSubdivisions=1;hs.truncationRadius=8;
   const auto hm=spin::triangulateQuadric(hyper,design,hs);
   if(hm.type!=spin::QuadricSurfaceType::Hyperboloid || hm.capFaceCount<=0 ||
      !spin::rationalQuadricConstructionIdentityMatches(hyper,hm,design,hs))
    throw std::runtime_error("authority hyperboloid classification/caps/identity failed");
   for(int center=hm.interfaceVertexCount;center<hm.closedBoundary.vertices.rows();++center) {
    std::set<int> ring;
    for(int f=hm.interfaceFaceCount;f<hm.closedBoundary.faces.rows();++f)
     for(int k=0;k<3;++k)if(hm.closedBoundary.faces(f,k)==center)
      for(int j=0;j<3;++j)if(hm.closedBoundary.faces(f,j)!=center)
       ring.insert(hm.closedBoundary.faces(f,j));
    if(ring.empty())throw std::runtime_error("authority cap has no ring");
    for(int a=0;a<3;++a) {
     mpq_class sum=0;for(int v:ring)sum+=mpq_class(hm.exactVertices[v][a]);
     if(sum/static_cast<unsigned long>(ring.size())!=mpq_class(hm.exactVertices[center][a]))
      throw std::runtime_error("authority cap center did not use its exact source ring");
    }
   }
  }
  std::cout<<"PASS: authority triangulation, fixed normals, exact audit and identity\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
