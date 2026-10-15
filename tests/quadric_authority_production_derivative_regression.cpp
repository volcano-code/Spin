#include "spin/QuadricMomentEvaluator.h"
#include "spin/PolyhedralMomentDifferential.h"
#include "physics/MassProperties.h"
#include <gmpxx.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

using A=spin::QuadricCoefficientAuthority10;
using Q=mpq_class;
static void require(bool c,const std::string& why){if(!c)throw std::runtime_error(why);}
static spin::TriangleMesh cube(double r) {
 spin::TriangleMesh m;m.vertices.resize(8,3);m.faces.resize(12,3);
 m.vertices<<-r,-r,-r,r,-r,-r,r,r,-r,-r,r,-r,-r,-r,r,r,-r,r,r,r,r,-r,r,r;
 m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
 return m;
}
// Independent divergence field in the x direction and cubic triangle rule.
// This does not use the production oriented-tetrahedron moment integrator.
static std::array<Q,10> independentMoments(const spin::QuadricMomentEvaluation& e) {
 const std::array<std::array<Q,3>,4> bary{{{{Q(1,3),Q(1,3),Q(1,3)}},
     {{Q(3,5),Q(1,5),Q(1,5)}},{{Q(1,5),Q(3,5),Q(1,5)}},{{Q(1,5),Q(1,5),Q(3,5)}}}};
 // mpq_class numerator/denominator construction assumes reduced inputs.
 // -9/16 is the same cubic-rule weight as -27/48, in canonical form.
 const std::array<Q,4> weights{{Q(-9,16),Q(25,48),Q(25,48),Q(25,48)}};
 std::array<Q,10> moments{};
 for(int f=0;f<e.materialBoundaryTarget.faces.rows();++f) {
  std::array<std::array<Q,3>,3> p;
  for(int k=0;k<3;++k)for(int a=0;a<3;++a)
   p[k][a]=Q(e.materialExactVerticesTarget[e.materialBoundaryTarget.faces(f,k)][a]);
  const Q areaX=((p[1][1]-p[0][1])*(p[2][2]-p[0][2])-
                 (p[1][2]-p[0][2])*(p[2][1]-p[0][1]))/2;
  for(int s=0;s<4;++s) {
   std::array<Q,3> position{};
   for(int a=0;a<3;++a)for(int k=0;k<3;++k)position[a]+=bary[s][k]*p[k][a];
   const Q x=position[0],y=position[1],z=position[2];
   const std::array<Q,10> field{{x,x*x/2,x*y,x*z,x*x*x/3,x*y*y,x*z*z,x*y*z,x*x*z/2,x*x*y/2}};
   for(int m=0;m<10;++m)moments[m]+=areaX*weights[s]*field[m];
  }
 }
 for(int m=0;m<10;++m)if(moments[m]!=Q(e.variableMomentsExactTarget[m])) {
   const Q difference=moments[m]-Q(e.variableMomentsExactTarget[m]);
   std::ostringstream diagnostic;diagnostic<<"independent moment "<<m<<" disagrees; approximate difference="
       <<difference.get_d()<<"; faces="<<e.materialBoundaryTarget.faces.rows()
       <<"; sign="<<mpq_sgn(difference.get_mpq_t())<<"; stored="<<e.variableMomentsExactTarget[m].substr(0,60)
       <<"; checker="<<moments[m].get_str().substr(0,60);
   throw std::runtime_error(diagnostic.str());
 }
 return moments;
}
int main() {
 try {
  spin::SpinDesignRequest request;request.outerSurface=cube(3);request.designDomain=cube(2);
  request.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(
      request.outerSurface.vertices,request.outerSurface.faces,1))-
      physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(request.designDomain.vertices,request.designDomain.faces,1));
  spin::QuadricMomentEvaluatorSettings settings;settings.triangulation.icosphereSubdivisions=0;
  settings.triangulation.rationalSourceCoordinates=true;settings.triangulation.compactRationalSourceCoordinates=true;
  settings.triangulation.truncationRadius=8;settings.triangulation.maximumRadiusDoublings=0;
  settings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
  settings.exactBinary64SourceBoundary=settings.exactBinary64MaterialBoundary=true;
  A::RationalVector q;q.fill("0/1");q[0]="-1/2";q[4]=q[5]=q[6]="1/2";
  const auto initial=A::fromNormalizedGeometry(q);
  A::RationalVector d;d.fill("0/1");d[0]="1/1";
  const auto changed=initial.geodesicUpdate(d,"1/10000000000000000000000000");
  spin::QuadricAuthorityMomentAttachment baseAttachment,changedAttachment;
  const auto base=spin::evaluateQuadricMaterialMoments(request,initial,settings,baseAttachment);
  const auto value=spin::evaluateQuadricMaterialMoments(request,changed,settings,changedAttachment);
  require(base.derivativeDefined && value.derivativeDefined && baseAttachment.defined && changedAttachment.defined,
      "same-operation authority derivative unavailable: "+value.derivativeDiagnostic);
  require(changedAttachment.coefficientAuthority==changed.geometry(),"production attachment lost coefficients");
  require(changedAttachment.differential.coefficientJacobian[0][0]!=baseAttachment.differential.coefficientJacobian[0][0],
      "production Jacobian rounded away authoritative update");
  require(changedAttachment.scientificHash!=baseAttachment.scientificHash,"authority derivative identity not bound");
  std::cout<<"ellipsoid authority derivative hash="<<changedAttachment.scientificHash<<'\n';
  // The explicit historical authority value-only entry remains fail-closed.
  require(!spin::evaluateQuadricMaterialMoments(request,changed,settings).derivativeDefined,
      "value-only interface silently became Newton-ready");
  // Independent material integration at actual plus/minus evaluations. This
  // local one-direction regression does not replace the 81-point protocol.
  double maxError=0;
  for(const std::string step:{"1/10000","1/20000","1/40000"}) {
   const Q h(step);auto negative=d;negative[0]="-1/1";
   const auto plus=spin::evaluateQuadricMaterialMoments(request,initial.geodesicUpdate(d,step),settings);
   const auto minus=spin::evaluateQuadricMaterialMoments(request,initial.geodesicUpdate(negative,step),settings);
   const auto plusMoments=independentMoments(plus),minusMoments=independentMoments(minus);
   for(int m=0;m<10;++m) {
    // Use the actual normalized coefficient displacement, not a fabricated
    // unit coordinate direction or an assumed normalization derivative.
    const auto p=initial.geodesicUpdate(d,step),n=initial.geodesicUpdate(negative,step);
    Q predicted=0;
    for(int j=0;j<10;++j)predicted+=Q(baseAttachment.differential.coefficientJacobian[m][j])*
        (Q(p.geometry()[j])-Q(n.geometry()[j]))/(2*h);
    const Q fd=(plusMoments[m]-minusMoments[m])/(2*h);
    const double error=std::abs(Q(fd-predicted).get_d())/std::max(1.,std::abs(predicted.get_d()));
    maxError=std::max(maxError,error);require(error<=1e-6,"actual authority material difference disagrees with Jacobian");
   }
  }
  std::cout<<"local production maximum difference="<<maxError<<'\n';
  q[6]="-1/2";const auto hyper=A::fromNormalizedGeometry(q);
  spin::QuadricAuthorityMomentAttachment hyperAttachment;
  const auto hyperValue=spin::evaluateQuadricMaterialMoments(request,hyper,settings,hyperAttachment);
  require(hyperValue.derivativeDefined && hyperAttachment.defined,
      "actual hyperboloid root/cap/Boolean differential unavailable: "+hyperValue.derivativeDiagnostic);
  require(hyperValue.constructedIntersectionVertexCount>0,"hyperboloid fixture missed intersection transport");
  independentMoments(hyperValue);
  std::cout<<"hyperboloid authority derivative hash="<<hyperAttachment.scientificHash<<'\n';
  std::cout<<"PASS: authority actual material, exact Jacobian and local finite differences\n";
 } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
