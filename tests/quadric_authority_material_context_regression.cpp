#include "spin/QuadricMomentEvaluator.h"
#include "spin2014/QuadricAuthorityNewtonState.h"
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>

static void require(bool c,const char* message){if(!c)throw std::runtime_error(message);}
static spin::TriangleMesh cube(double r) {
 spin::TriangleMesh m;m.vertices.resize(8,3);m.faces.resize(12,3);
 m.vertices<<-r,-r,-r,r,-r,-r,r,r,-r,-r,r,-r,-r,-r,r,r,-r,r,r,r,r,-r,r,r;
 m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
 return m;
}
int main() {
 try {
  spin::SpinDesignRequest request;request.outerSurface=cube(3);request.designDomain=cube(2);
  request.materialDensity=3;request.mode=spin::SpinMode::SpinningTop;
  const auto context=spin::makeQuadricAuthorityMaterialContext(request,2);
  require(context.fixedMomentsTarget()[0]=="304/1","fixed shell must use its own density once");
  for(unsigned i=4;i<7;++i)require(context.fixedMomentsTarget()[i]=="3376/3","fixed shell second moments wrong");
  for(unsigned i:{1U,2U,3U,7U,8U,9U})require(context.fixedMomentsTarget()[i]=="0/1","symmetric shell moments wrong");
  require(context.request().fixedShellMoments.mass()==304,"compatibility request did not derive actual shell");
  auto changedRequest=request;changedRequest.fixedShellMoments=context.request().fixedShellMoments.scaledDensity(99);
  const auto numbersIgnored=spin::makeQuadricAuthorityMaterialContext(changedRequest,2);
  require(numbersIgnored.identity()==context.identity(),"unowned caller moments influenced geometric context");
  const auto changedDensity=spin::makeQuadricAuthorityMaterialContext(request,4);
  require(changedDensity.identity()!=context.identity() && changedDensity.fixedMomentsTarget()[0]=="608/1",
      "changed actual density reused old material identity");
  request.targetFrame.origin=Eigen::Vector3d(.25,-.5,.125);
  request.targetFrame.axes<<0,-1,0,1,0,0,0,0,1;
  const auto transformed=spin::makeQuadricAuthorityMaterialContext(request,2);
  require(transformed.fixedMomentsTarget()[1]=="152/1" &&
          transformed.fixedMomentsTarget()[2]=="76/1" &&
          transformed.fixedMomentsTarget()[3]=="-38/1","fixed moments not integrated in actual target coordinates");
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<4;++j) {
    mpq_class composed=0;
    for(unsigned k=0;k<3;++k)composed+=mpq_class(transformed.targetToModel()[i][k])*
        mpq_class(transformed.modelToTarget()[k][j]);
    if(j==3)composed+=mpq_class(transformed.targetToModel()[i][3]);
    require(composed==(i==j?1:0),"affine inverse not bound to exact forward map");
  }
  spin::QuadricMomentEvaluatorSettings settings;settings.triangulation.icosphereSubdivisions=0;
  settings.triangulation.rationalSourceCoordinates=true;settings.triangulation.compactRationalSourceCoordinates=true;
  settings.triangulation.truncationRadius=8;settings.triangulation.maximumRadiusDoublings=0;
  settings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
  settings.exactBinary64SourceBoundary=settings.exactBinary64MaterialBoundary=true;
  spin::QuadricCoefficientAuthority10::RationalVector q;q.fill("0/1");q[0]="-1/2";q[4]=q[5]=q[6]="1/2";
  const auto authority=spin::QuadricCoefficientAuthority10::fromNormalizedGeometry(q);
  const auto evaluated=spin::evaluateQuadricMaterialMoments(context,authority,settings);
  require(evaluated.derivative.defined && evaluated.material.derivativeDefined,"actual authority material derivative missing");
  for(unsigned i=0;i<10;++i) {
    mpq_class expected=mpq_class(context.fixedMomentsTarget()[i])+mpq_class(evaluated.material.variableMomentsExactTarget[i]);
    expected.canonicalize();
    require(mpq_class(evaluated.totalMomentsTarget[i])==expected,"total moments lost geometry-bound contributions");
  }
  spin2014::QuadricObjectiveSettings objective;objective.objective=spin2014::Spin2014Objective::TopMassAndLowering;
  const auto state=spin2014::linearizeSpin2014AuthorityNewtonState(objective,authority,evaluated.totalMomentsTarget,
      evaluated.derivative.differential.coefficientJacobian,evaluated.scientificHash);
  require(state.evaluationMapIdentity==evaluated.scientificHash && state.matrix.size()==15,
      "real geometry/moments/Jacobian did not reach authority Eq.14");
  bool rejected=false;
  try {(void)spin::makeQuadricAuthorityMaterialContext(request,-1);}catch(const std::exception&){rejected=true;}
  require(rejected,"negative fixed density accepted");
  auto outside=request;outside.designDomain=cube(4);
  rejected=false;
  try {(void)spin::makeQuadricAuthorityMaterialContext(outside,2);}catch(const std::exception&){rejected=true;}
  require(rejected,"design outside outer domain accepted");
  auto reversed=request;
  for(int f=0;f<reversed.designDomain.faces.rows();++f)
      std::swap(reversed.designDomain.faces(f,1),reversed.designDomain.faces(f,2));
  rejected=false;
  try {(void)spin::makeQuadricAuthorityMaterialContext(reversed,2);}catch(const std::exception&){rejected=true;}
  require(rejected,"silently reversed invalid input winding");
  std::cout<<"geometry-owned fixed shell, two densities, affine inverse and actual material Eq.14 passed\n";
  std::cout<<"context hash="<<context.identity()<<"; state hash="<<state.scientificHash<<'\n';
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
