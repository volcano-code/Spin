#include "spin/QuadricMomentEvaluator.h"
#include "physics/MassProperties.h"
#include <gmpxx.h>
#include <map>
#include <iostream>
#include <stdexcept>

static spin::TriangleMesh cube(double extent) {
 spin::TriangleMesh m;m.vertices.resize(8,3);
 m.vertices<<-extent,-extent,-extent,extent,-extent,-extent,extent,extent,-extent,-extent,extent,-extent,
             -extent,-extent,extent,extent,-extent,extent,extent,extent,extent,-extent,extent,extent;
 m.faces.resize(12,3);m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
 return m;
}

// Independent checker: divergence field (x^(a+1)y^b z^c/(a+1),0,0),
// followed by exact barycentric monomial integration on every triangle.
// Does not call the production signed-tetrahedron mass-moment function.
static std::array<mpq_class,10> checkMoments(const spin::QuadricMomentEvaluation& e,unsigned density) {
 using Powers=std::array<int,3>;
 using Polynomial=std::map<Powers,mpq_class>;
 const std::array<Powers,10> powers{{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},
     {{2,0,0}},{{0,2,0}},{{0,0,2}},{{0,1,1}},{{1,0,1}},{{1,1,0}}}};
 const auto factorial=[](int n){unsigned long r=1;for(int i=2;i<=n;++i)r*=i;return r;};
 std::array<mpq_class,10> result{};
 for(int f=0;f<e.materialBoundaryTarget.faces.rows();++f) {
  std::array<std::array<mpq_class,3>,3> p;
  for(int k=0;k<3;++k)for(int a=0;a<3;++a)
   p[k][a]=mpq_class(e.materialExactVerticesTarget.at(e.materialBoundaryTarget.faces(f,k))[a]);
  const mpq_class normalX=(p[1][1]-p[0][1])*(p[2][2]-p[0][2])-
                           (p[1][2]-p[0][2])*(p[2][1]-p[0][1]);
  for(int m=0;m<10;++m) {
   Polynomial polynomial;polynomial[Powers{{0,0,0}}]=1;
   for(int axis=0;axis<3;++axis)for(int power=0;power<powers[m][axis]+(axis==0);++power) {
    Polynomial next;
    for(const auto& term:polynomial)for(int corner=0;corner<3;++corner) {
     auto exponent=term.first;++exponent[corner];
     next[exponent]+=term.second*p[corner][axis];
    }
    polynomial=std::move(next);
   }
   mpq_class integral=0;
   for(const auto& term:polynomial) {
    const auto& a=term.first;
    integral+=term.second*mpq_class(factorial(a[0])*factorial(a[1])*factorial(a[2]),
                                  factorial(a[0]+a[1]+a[2]+2));
   }
   result[m]+=normalX*integral*density/(powers[m][0]+1);
  }
 }
 return result;
}

int main() {
 try {
  using A=spin::QuadricCoefficientAuthority10;
  A::RationalVector q;q.fill("0/1");q[0]="-1/2";q[4]=q[5]=q[6]="1/2";
  A::RationalVector direction;direction.fill("0/1");direction[0]="1/1";
  const auto initial=A::fromNormalizedGeometry(q);
  const auto updated=initial.geodesicUpdate(direction,"1/10000000000000000000000000");
  spin::SpinDesignRequest request;request.outerSurface=cube(3);request.designDomain=cube(2);
  const auto outer=computeSurfaceMassMoments(request.outerSurface.vertices,request.outerSurface.faces,1);
  const auto inner=computeSurfaceMassMoments(request.designDomain.vertices,request.designDomain.faces,1);
  request.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(outer)-physics::RawMassMoments10::fromMassMoments(inner);
  spin::QuadricMomentEvaluatorSettings settings;
  settings.triangulation.icosphereSubdivisions=0;settings.triangulation.rationalSourceCoordinates=true;
  settings.triangulation.compactRationalSourceCoordinates=true;
  settings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
  settings.exactBinary64SourceBoundary=settings.exactBinary64MaterialBoundary=true;
  const auto base=spin::evaluateQuadricMaterialMoments(request,initial,settings);
  const auto changed=spin::evaluateQuadricMaterialMoments(request,updated,settings);
  if(mpq_class(changed.variableMomentsExactTarget.at(0))<=mpq_class(base.variableMomentsExactTarget.at(0)))
   throw std::runtime_error("actual Boolean mass lost the authoritative update");
  if(changed.derivativeDefined || changed.derivativeDiagnostic.empty())
   throw std::runtime_error("unverified authority Jacobian was granted");
  if(changed.derivativeScientificHash==base.derivativeScientificHash ||
     changed.quadricTriangulation.rationalConstructionIdentity==base.quadricTriangulation.rationalConstructionIdentity)
   throw std::runtime_error("material identity lost coefficient authority");
  const auto checked=checkMoments(changed,1);
  for(int i=0;i<10;++i)if(checked[i]!=mpq_class(changed.variableMomentsExactTarget.at(i)))
   throw std::runtime_error("independent divergence integration disagrees with actual material");
  request.materialDensity=2;
  const auto dense=spin::evaluateQuadricMaterialMoments(request,updated,settings);
  const auto denseChecked=checkMoments(dense,2);
  for(int i=0;i<10;++i)if(denseChecked[i]!=2*checked[i] ||
                           denseChecked[i]!=mpq_class(dense.variableMomentsExactTarget.at(i)))
   throw std::runtime_error("material density was lost or multiplied twice");
  std::cout<<"PASS: authority Boolean, exact ten moments, independent divergence and density\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
