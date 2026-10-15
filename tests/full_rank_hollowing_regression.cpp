#include "spin/FullRankHollowingSolver.h"
#include "spin/DualMonomialCoordinates.h"
#include "spin/QuadricSourceMotion.h"
#include <Eigen/Geometry>
#include <array>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactRationalCoordinates.h"
#endif
namespace {
using V=spin::QuadricCoefficients10::Vector;
void require(bool b,const std::string& s) { if(!b) throw std::runtime_error(s); }
spin::TriangleMesh box(double a,double b) {
 spin::TriangleMesh m;m.vertices.resize(8,3);m.faces.resize(12,3);
 m.vertices<<a,a,a,b,a,a,b,b,a,a,b,a,a,a,b,b,a,b,b,b,b,a,b,b;
 m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;return m;
}
spin::SpinDesignRequest request() {
 spin::SpinDesignRequest r;r.designDomain=box(-1,1);r.outerSurface=box(-2,2);
 r.minimumWallThickness=.2;r.materialDensity=1.7;
 r.settings.equationProfile=spin::SpinEquationProfile::DualMonomialRawV1;
 return r;
}
spin::QuadricMomentEvaluatorSettings settings() {
 spin::QuadricMomentEvaluatorSettings s;s.triangulation.icosphereSubdivisions=2;
 s.triangulation.truncationRadius=8;return s;
}
spin::QuadricCoefficients10 fixture() { V p;p<<-.31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03;
 return spin::QuadricCoefficients10::fromUnnormalized(p); }
void adapter() {
 auto r=request();auto es=settings();const auto q=fixture();
 r.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.outerSurface.vertices,r.outerSurface.faces,r.materialDensity))-
 physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.designDomain.vertices,r.designDomain.faces,r.materialDensity));
 const auto reference=spin::evaluateFullRankMaterialMoments(r,q,es);
 auto discreteSettings=es;
 discreteSettings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
 const auto unified=spin::evaluateQuadricMaterialMoments(r,q,discreteSettings);
 require(reference.coefficientJacobian.has_value(),"fixture derivative unavailable");
 require(reference.effectiveDerivativeMethod=="discrete-fullrank-quadric-boolean-transport-v2",
         "strict full-rank path did not bind its Boolean moments to the discrete derivative");
 require(unified.derivativeDefined && reference.geometryHash==unified.geometryHash &&
         spin::sameOrientedTriangleBoundary(reference.materialBoundaryTarget,unified.materialBoundaryTarget),
         "strict adapter geometry is not the actual value-side Boolean output");
 require(reference.coefficientJacobian->allFinite() && !reference.sourceScientificHash.empty(),
         "strict adapter did not return a finite derivative/source identity");
 require((reference.coefficientJacobian.value()-unified.momentJacobian).norm()<=
         1e-10*std::max(1.0,unified.momentJacobian.norm()),
         "strict derivative disagrees with the same value-side Boolean map");
 require(reference.topologyHash==unified.provenanceTopologyHash,
         "strict adapter topology identity differs from the value-side Boolean provenance");
 auto evaluator=spin::makeFullRankMomentEvaluator(r,es);const auto value=evaluator(q);
 require(value.derivativeDefined,"adapter lost available native Jacobian");
 require(value.rawMoments==reference.totalMomentsTarget.values(),"adapter omitted native total/fixed-shell moments");
 require(value.momentJacobian==*reference.coefficientJacobian,"adapter changed native-p Jacobian");
 require(value.geometryHash==reference.geometryHash,"adapter lost geometry identity");
 if (!(value.topologyHash.size()==64 && value.topologyHash==reference.topologyHash)) {
     std::cerr << "[DETAIL] adapter topology value='" << value.topologyHash
               << "' reference='" << reference.topologyHash << "'\n";
     require(false,"adapter lost provenance-aware Boolean topology identity");
 }
 r.materialDensity=9;es.triangulation.truncationRadius=20;
 const auto captured=evaluator(q);
 require(captured.rawMoments==value.rawMoments && captured.momentJacobian==value.momentJacobian && captured.geometryHash==value.geometryHash,"caller mutation leaked into captured map");
 require(value.residualEvaluationSeconds>0 && value.derivativeEvaluationSeconds==0,"coupled evaluation timing must be counted once");
}
void chain() {
 const auto r=request();const auto es=settings();const auto q=fixture();
 const auto value=spin::makeFullRankMomentEvaluator(r,es)(q);
 const auto s=spin::DualMonomialCoordinates10::fromGeometry(q);
 const Eigen::Matrix<double,10,10> predicted=value.momentJacobian*s.geometryJacobian();
 double previous=1;
 for(double h:{1e-4,3e-5,1e-5}) {
  Eigen::Matrix<double,10,10> fd;
  for(int j=0;j<10;++j) {
   V plus=s.values(),minus=s.values();plus(j)+=h;minus(j)-=h;
   const auto a=spin::evaluateQuadricMaterialMoments(r,spin::DualMonomialCoordinates10::fromUnnormalized(plus).geometry(),es);
   const auto b=spin::evaluateQuadricMaterialMoments(r,spin::DualMonomialCoordinates10::fromUnnormalized(minus).geometry(),es);
   fd.col(j)=(a.totalMomentsTarget.values()-b.totalMomentsTarget.values())/(2*h);
  }
  const double error=(fd-predicted).norm()/predicted.norm();
  std::cout<<"dual-chain h="<<h<<" relative="<<error<<'\n';
  if(h==1e-5)require(error<=1e-6,"double D scaling or missing normalization in adapter composition");
  if(h<1e-4 && previous>1e-8) require(error<.3*previous,"dual finite difference failed to converge");
  previous=error;
 }
}
template<class F> void throwsInvalid(F f) {bool caught=false;try{f();}catch(const std::invalid_argument&){caught=true;}require(caught,"invalid request/profile/settings did not throw at boundary");}
void negative() {
 auto r=request();const auto es=settings();
 auto invalid=r;invalid.settings.equationProfile=spin::SpinEquationProfile::LegacyNativeNormalizedV1;
 throwsInvalid([&]{spin::makeFullRankMomentEvaluator(invalid,es);});
 throwsInvalid([&]{spin::solveFullRankQuadricHollowing(invalid,es,{fixture()});});
 invalid=r;invalid.settings.equationTolerance=0;
 throwsInvalid([&]{spin::makeFullRankMomentEvaluator(invalid,es);});
 invalid=r;invalid.materialDensity=-1;
 throwsInvalid([&]{spin::solveFullRankQuadricHollowing(invalid,es,{fixture()});});
 invalid=r;invalid.staticObjective=static_cast<spin::StaticObjective>(999);
 throwsInvalid([&]{spin::makeFullRankMomentEvaluator(invalid,es);});
 auto bad=es;bad.derivativeMethod=static_cast<spin::QuadricMomentDerivativeMethod>(999);
 throwsInvalid([&]{spin::makeFullRankMomentEvaluator(r,bad);});
 V p;p<<.2,0,0,-1,0,0,0,0,0,0;
 const auto plane=spin::QuadricCoefficients10::fromUnnormalized(p);
 // A genuine rank-two cylinder replaces the formerly misclassified regular
 // inactive-centroid fixture; the latter is now tested in selector-fallback.
 p<<-.4,0,0,0,1,1,0,0,0,0;
 r.settings.maximumRestarts=1;r.settings.maximumFunctionEvaluationsPerRestart=2;
 for(const auto& q:{plane,spin::QuadricCoefficients10::fromUnnormalized(p)}) {
  const auto map=spin::makeFullRankMomentEvaluator(r,es)(q);
  require(!map.derivativeDefined && !map.derivativeDiagnostic.empty() && map.rawMoments(0)>0,"unsupported derivative fabricated or native values lost");
  const auto solved=spin::solveFullRankQuadricHollowing(r,es,{q});
  require(solved.numerical.status==spin::SpinSolverStatus::GeometryEvaluationFailed,"undefined derivative failure semantics changed");
  require(!solved.finalGeometryVerified && !solved.finalEvaluation,"failed solve returned final material");
 }
}
void solver() {
 auto r=request();const auto es=settings();r.mode=spin::SpinMode::SpinningTop;r.centerOfMassWeight=.1;
 r.settings.maximumIterations=1;r.settings.maximumRestarts=1;r.settings.maximumFunctionEvaluationsPerRestart=3;
 const auto out=spin::solveFullRankQuadricHollowing(r,es,{fixture()});
 require(out.numerical.functionEvaluations>0,"solver never evaluated real geometry");
 require(out.numerical.status!=spin::SpinSolverStatus::Infeasible && out.numerical.status!=spin::SpinSolverStatus::Converged,"small-budget probe was mislabeled");
 require(!out.finalGeometryVerified && !out.finalEvaluation && out.finalVerificationSeconds==0,"failed numerical search exposed certified material");
 require(out.numerical.functionEvaluations<=3 && out.numerical.evaluationTimings.size()==static_cast<std::size_t>(out.numerical.functionEvaluations),"evaluation attempt accounting lost");
 require(!out.numerical.trace.empty(),"real fixture had no accepted initial state");
 for(const auto& t:out.numerical.trace)require(!t.materialGeometryHash.empty() && t.equationProfile==r.settings.equationProfile,"accepted trace missing geometry/profile");
 for(const auto& e:out.numerical.events)if(e.rawEquationResidual)require(e.materialGeometryHash && !e.materialGeometryHash->empty() && e.equationProfile==r.settings.equationProfile,"evaluated event missing geometry/profile");
 require(out.totalWallSeconds>=out.numerical.evaluatorWallSeconds && out.numerical.residualEvaluationSeconds>0,"wall accounting missing work");
 const auto again=spin::solveFullRankQuadricHollowing(r,es,{fixture()});
 require(!out.inputHash.empty() && out.inputHash==again.inputHash && out.outputHash==again.outputHash,"scientific hashes unstable");
 auto changed=r;changed.materialDensity*=2;
 require(spin::solveFullRankQuadricHollowing(changed,es,{fixture()}).inputHash!=out.inputHash,"density not bound");
 changed=r;changed.settings.deterministicSeed++;
 require(spin::solveFullRankQuadricHollowing(changed,es,{fixture()}).inputHash!=out.inputHash,"seed not bound");
 auto modified=es;modified.triangulation.retainFixedNormalProvenance=true;
 require(spin::solveFullRankQuadricHollowing(r,modified,{fixture()}).inputHash!=out.inputHash,"evaluator metadata not bound");
 const auto defaults=spin::deterministicInitialQuadrics(r.mode,r.settings);
 require(spin::solveFullRankQuadricHollowing(r,es).inputHash!=spin::solveFullRankQuadricHollowing(r,es,defaults).inputHash,"default restart marker missing");
 std::cout<<"budget status="<<spin::spinSolverStatusToken(out.numerical.status)<<" evaluations="<<out.numerical.functionEvaluations<<" raw="<<*out.numerical.rawEquationResidual<<'\n';
}
void integration() {
 const auto b=spin::integrateFullRankBoundaryMoments(box(-1,1),1.7);
 V expected=V::Zero();expected(0)=13.6;expected.segment<3>(4).setConstant(13.6/3);
 require((b.values()-expected).norm()<1e-13,"independent density-weighted box moments wrong");
 spin::TriangleMesh t;t.vertices.resize(4,3);t.vertices<<0,0,0,1,0,0,0,1,0,0,0,1;
 t.faces.resize(4,3);t.faces<<0,2,1,0,1,3,0,3,2,1,2,3;
 const auto raw=spin::integrateFullRankBoundaryMoments(t,2);
 expected<<1./3,1./12,1./12,1./12,1./30,1./30,1./30,1./60,1./60,1./60;
 require((raw.values()-expected).norm()<1e-14,"independent analytic tetra polynomial moments wrong");
 const Eigen::Matrix3d rotation=Eigen::AngleAxisd(.43,Eigen::Vector3d(1,2,-1).normalized()).toRotationMatrix();
 const Eigen::Vector3d translation(.3,-.2,.1);
 for(int i=0;i<t.vertices.rows();++i)t.vertices.row(i)=rotation*t.vertices.row(i).transpose()+translation;
 require((spin::integrateFullRankBoundaryMoments(t,2).values()-raw.transformed(rotation,translation).values()).norm()<1e-14,"signed-tetra affine covariance wrong");
 t.faces.col(0).swap(t.faces.col(2));bool caught=false;
 try{spin::integrateFullRankBoundaryMoments(t,2);}catch(const std::runtime_error&){caught=true;}
 require(caught,"inward material accepted by strict independent audit");
 throwsInvalid([&]{spin::integrateFullRankBoundaryMoments(box(-1,1),0);});
}
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
void rationalIntegration() {
 // Break caught: independently validating the authority through its rounded
 // display, or integrating density twice. Distinct exact x points round to
 // the SAME binary64 point. This is a geometry unit test, not a solve.
 using Q=CGAL::Gmpq;
 spin::ExactBooleanVertexCoordinates points={
  {"1/1","0/1","0/1"},
  {"1208925819614629174706177/1208925819614629174706176","0/1","0/1"},
  {"1/1","1/1","0/1"},{"1/1","0/1","1/1"}};
 Eigen::MatrixXi faces(4,3);faces<<0,2,1,0,1,3,0,3,2,1,2,3;
 const auto actual=spin::integrateFullRankBoundaryMoments(points,faces,2.5);
 const Q epsilon=verification::parseExactRationalCoordinate("1/1208925819614629174706176");
 const Q mass=Q(5)*epsilon/12;
 const std::array<Q,10> expected{{mass,mass*(Q(1)+epsilon/4),mass/4,mass/4,
  mass*(Q(1)+epsilon/2+epsilon*epsilon/10),mass/10,mass/10,mass/20,
  mass*(Q(1)/4+epsilon/20),mass*(Q(1)/4+epsilon/20)}};
 require(actual.size()==10,"rational independent integration lost mass moment components");
 for(int i=0;i<10;++i)require(verification::parseExactRationalCoordinate(actual[i])==expected[i],
                            "independent rational tetrahedron/density moment differs from hand integral");
 spin::TriangleMesh display;display.vertices.resize(4,3);
 display.vertices<<1,0,0,1,0,0,1,1,0,1,0,1;display.faces=faces;
 bool caught=false;try{spin::integrateFullRankBoundaryMoments(display,2.5);}
 catch(const std::runtime_error&){caught=true;}
 require(caught,"invalid floating derivative was blessed by authority unit test");
 // Signed nested shells are actual material, not an overlapping component
 // sum. Their independently derived volume is 64-8; Qxx=256/3-8/3.
 auto outer=box(-2,2),inner=box(-1,1);
 spin::ExactBooleanVertexCoordinates shells;
 for(const auto* mesh:{&outer,&inner})for(int v=0;v<mesh->vertices.rows();++v)
  shells.push_back({std::to_string(int(mesh->vertices(v,0)))+"/1",
                    std::to_string(int(mesh->vertices(v,1)))+"/1",
                    std::to_string(int(mesh->vertices(v,2)))+"/1"});
 Eigen::MatrixXi shellFaces(24,3);shellFaces.topRows(12)=outer.faces;
 for(int f=0;f<12;++f)for(int a=0;a<3;++a)shellFaces(12+f,a)=8+inner.faces(f,2-a);
 const auto shell=spin::integrateFullRankBoundaryMoments(shells,shellFaces,4.);
 const std::array<Q,10> shellExpected{{Q(224),Q(0),Q(0),Q(0),Q(992)/3,Q(992)/3,Q(992)/3,Q(0),Q(0),Q(0)}};
 for(int i=0;i<10;++i)require(verification::parseExactRationalCoordinate(shell[i])==shellExpected[i],
                            "signed shell integration omitted cavity or density");
 auto bad=faces;bad(0,0)=bad(0,1);caught=false;
 try{spin::integrateFullRankBoundaryMoments(points,bad,2.5);}
 catch(const std::exception&){caught=true;}
 require(caught,"invalid exact face accepted as authoritative material");
 throwsInvalid([&]{spin::integrateFullRankBoundaryMoments(points,faces,0.);});
}
#endif
void modes() {
 const auto check=[](spin::SpinDesignRequest r,const spin::SpinMomentProblem& expected) {
  r.settings.maximumRestarts=1;r.settings.maximumFunctionEvaluationsPerRestart=1;
  const auto map=spin::evaluateFullRankMaterialMoments(r,fixture(),settings());
  const auto reference=expected.evaluate(map.totalMomentsTarget.values());
  const auto out=spin::solveFullRankQuadricHollowing(r,settings(),{fixture()});
  const auto initial=std::find_if(out.numerical.events.begin(),out.numerical.events.end(),
    [](const spin::SpinSolverEvent& event) {
      return event.evaluationIndex==1 && event.objective.has_value();
    });
  require(initial!=out.numerical.events.end() && initial->objective &&
          initial->physicalRelativeResidual && initial->multiplierMu,
          "mode dispatch lost the actual initial evaluation");
  require(*initial->objective==reference.objective,
          "request mode/objective dispatched to wrong moment problem");
  require(std::abs(*initial->physicalRelativeResidual-
          reference.constraints.norm()/std::max(1.,map.totalMomentsTarget.values().norm()))<1e-15,
          "mode dispatch changed physical constraints");
  if (out.numerical.trace.empty())
    require(initial->reason=="initial_multiplier_nonpositive" &&
            *initial->multiplierMu<=0.0 && !out.finalEvaluation,
            "strict initializer rejected valid geometry for a reason other than nonpositive multiplier");
 };
 auto r=request();check(r,spin::SpinMomentProblem::staticBalance(spin::StaticObjective::MaximumMass));
 r.staticObjective=spin::StaticObjective::LowerCenterOfMass;check(r,spin::SpinMomentProblem::staticBalance(spin::StaticObjective::LowerCenterOfMass));
 r.staticObjective=spin::StaticObjective::ParetoTradeoff;r.staticParetoWeight=.4;check(r,spin::SpinMomentProblem::staticBalancePareto(.4));
 r=request();r.mode=spin::SpinMode::FidgetSpinner;check(r,spin::SpinMomentProblem::fidgetSpinner());
 r.mode=spin::SpinMode::YoYo;check(r,spin::SpinMomentProblem::yoYo());
 r.mode=spin::SpinMode::SpinningTop;r.centerOfMassWeight=.1;check(r,spin::SpinMomentProblem::spinningTop(.1));
 r=request();r.mode=spin::SpinMode::BuoyantStability;r.buoyancy.liquidDensity=1;r.buoyancy.displacedLiquidVolume=2.3;
 r.buoyancy.centerOfBuoyancyTarget<<.1,-.2,0;r.buoyancy.waterlineArea=2;r.buoyancy.waterlineCentralSecondMoment=Eigen::Matrix2d::Identity();
 check(r,spin::SpinMomentProblem::buoyantStability(2.3,Eigen::Vector2d(.1,-.2)));
}
void staticLowerTrace() {
 auto r=request();r.mode=spin::SpinMode::StaticBalance;
 r.staticObjective=spin::StaticObjective::LowerCenterOfMass;
 r.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(
     computeSurfaceMassMoments(r.outerSurface.vertices,r.outerSurface.faces,r.materialDensity))-
     physics::RawMassMoments10::fromMassMoments(
     computeSurfaceMassMoments(r.designDomain.vertices,r.designDomain.faces,r.materialDensity));
 r.settings.maximumRestarts=1;r.settings.maximumFunctionEvaluationsPerRestart=60;
 const auto out=spin::solveFullRankQuadricHollowing(r,settings());
 std::cout<<"static-lowering diagnostic status="<<spin::spinSolverStatusToken(out.numerical.status)
     <<" evaluations="<<out.numerical.functionEvaluations<<" residual="<<out.numerical.normalizedEquationResidual
     <<" final="<<out.finalGeometryVerified<<" diagnostic="<<out.numerical.diagnostic<<'\n';
 const auto begin=out.numerical.events.size()>12?out.numerical.events.size()-12:0;
 for(std::size_t i=begin;i<out.numerical.events.size();++i) {
  const auto& e=out.numerical.events[i];
  std::cout<<"event index="<<e.evaluationIndex<<" restart="<<e.restart<<" iteration="<<e.iteration
      <<" type="<<e.type<<" outcome="<<e.outcome<<" reason="<<e.reason
      <<" residual="<<e.normalizedEquationResidual.value_or(-1.)
      <<" q=";
  if(e.quadricCoefficients)for(const double value:*e.quadricCoefficients)std::cout<<value<<',';
  else std::cout<<"missing";
  std::cout<<" detail="<<e.detail<<'\n';
 }
 require(out.numerical.status==spin::SpinSolverStatus::ResourceLimitReached&&
     out.numerical.functionEvaluations==60&&!out.finalGeometryVerified,
     "diagnostic replay changed the frozen lower-center-of-mass budget outcome");
}
void nearParabolicTruncationProbe() {
 V p;
 p<<-0.03229148714369761,-3.5981802826533382e-05,0.00013923083159621887,
     0.99947721704275616,-0.00011842626037807066,0.0010909647938445234,
     0.00047307981633996752,-0.00023316715451590989,0.0005345316953516031,
     -0.00087387989594205962;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
 const auto domain=box(-1,1);
 for(int doublings:{8,16,24,32}) {
  auto es=settings();es.triangulation.icosphereSubdivisions=2;
  es.triangulation.maximumRadiusDoublings=doublings;
  try {
   const auto tri=spin::triangulateQuadric(q,domain.vertices,es.triangulation);
   std::cout<<"near-parabolic doublings="<<doublings<<" status=OK loops="<<tri.boundaryLoopCount
       <<" radius="<<tri.truncationRadius<<" faces="<<tri.closedBoundary.faces.rows()<<'\n';
   return;
  } catch(const std::exception& error) {
   std::cout<<"near-parabolic doublings="<<doublings<<" status=FAIL reason="<<error.what()<<'\n';
  }
 }
}
void offTargetYoyo() {
 auto r=request();auto es=settings();r.mode=spin::SpinMode::YoYo;
 r.settings.maximumRestarts=1;
 r.settings.maximumFunctionEvaluationsPerRestart=60;
 V p;p<<0,0,0,1,1,1,-1,0,0,0;
 const auto target=spin::QuadricCoefficients10::fromUnnormalized(p);
 const auto targetEvaluation=spin::evaluateFullRankMaterialMoments(r,target,es);
 MassMoments desired;
 desired.mass=20.0;
 desired.firstMoment=Eigen::Vector3d::Zero();
 desired.secondMoment=Eigen::Vector3d(100.0,100.0,200.0).asDiagonal();
 r.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(desired)-
     targetEvaluation.variableMomentsTarget;
 p=target.values();p(0)+=.015;p(6)-=.025;
 const auto seed=spin::QuadricCoefficients10::fromUnnormalized(p);
 require(seed.hash()!=target.hash(),"off-target seed equals the known solution");
 const auto result=spin::solveFullRankQuadricHollowing(r,es,{seed});
 const auto repeated=spin::solveFullRankQuadricHollowing(r,es,{seed});
 std::cout<<"off-target-yoyo status="
          <<spin::spinSolverStatusToken(result.numerical.status)
          <<" residual="<<result.numerical.normalizedEquationResidual
          <<" trace="<<result.numerical.trace.size()
          <<" final_verified="<<result.finalGeometryVerified
          <<" output_hash="<<result.outputHash
          <<" trace_hash="<<result.numerical.traceHash
          <<" diagnostic="<<result.verificationDiagnostic<<'\n';
 bool sawTopologyRelinearization=false;
 for(const auto& e:result.numerical.events) sawTopologyRelinearization=sawTopologyRelinearization||
   (e.type=="evaluation" && e.outcome=="accepted" && e.reason=="topology_transition_relinearized" &&
    e.detail.find("fresh derivative")!=std::string::npos);
 require(result.numerical.status==spin::SpinSolverStatus::Converged &&
             sawTopologyRelinearization && result.finalGeometryVerified && result.finalEvaluation,
         "off-target yoyo did not record a residual-decreasing topology relinearization");
 require(repeated.outputHash==result.outputHash &&
             repeated.numerical.eventHash==result.numerical.eventHash,
         "off-target yoyo topology-transition replay changed the scientific result");
}
void audit(bool strictReplay) {
 auto r=request();auto es=settings();r.designDomain=box(-.8,.8);r.outerSurface=box(-1,1);
 r.designDomain.vertices.col(2).array()+=2;r.outerSurface.vertices.col(2).array()+=2;
 r.minimumWallThickness=.15;r.materialDensity=1;r.mode=spin::SpinMode::SpinningTop;r.centerOfMassWeight=1;
 r.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.outerSurface.vertices,r.outerSurface.faces,1))-
 physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.designDomain.vertices,r.designDomain.faces,1));
 r.settings.maximumRestarts=1;r.settings.maximumFunctionEvaluationsPerRestart=20;es.triangulation.icosphereSubdivisions=1;
 const auto seed=fixture();const Eigen::Vector3d c(0,0,2);V p=seed.values();
 p(0)=seed.constant()-seed.linearVector().dot(c)+c.dot(seed.quadraticMatrix()*c);
 p.segment<3>(1)=seed.linearVector()-2*seed.quadraticMatrix()*c;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
 if(strictReplay) {
  const auto mapA=spin::makeFullRankMomentEvaluator(r,es)(q);
  const auto mapB=spin::makeFullRankMomentEvaluator(r,es)(q);
  std::cout<<"same-q replay geometry="<<(mapA.geometryHash==mapB.geometryHash)<<" moment-diff="<<(mapA.rawMoments-mapB.rawMoments).norm()<<" J-diff="<<(mapA.momentJacobian-mapB.momentJacobian).norm()<<'\n';
 }
 const auto solved=spin::solveFullRankQuadricHollowing(r,es,{q});
 bool sawTopologyRelinearization=false;
 for(const auto& e:solved.numerical.events) sawTopologyRelinearization=sawTopologyRelinearization||
   (e.type=="evaluation" && e.outcome=="accepted" && e.reason=="topology_transition_relinearized" &&
    e.detail.find("fresh derivative")!=std::string::npos);
 require(sawTopologyRelinearization,
         "spinning-top fixture did not record its topology relinearization");
 if(solved.numerical.status==spin::SpinSolverStatus::Converged)
  require(solved.finalGeometryVerified && solved.finalEvaluation,
          "converged topology-relinearized fixture lacks final geometry verification");
 else
  require(!solved.finalGeometryVerified && !solved.finalEvaluation,
          "nonconverged topology-relinearized fixture exposed final geometry");
 if(!strictReplay)return;
 const auto repeated=spin::solveFullRankQuadricHollowing(r,es,{q});
 if(repeated.outputHash!=solved.outputHash) {
  std::cout<<"hash-debug input="<<(repeated.inputHash==solved.inputHash)<<" event="<<(repeated.numerical.eventHash==solved.numerical.eventHash)<<" trace="<<(repeated.numerical.traceHash==solved.numerical.traceHash)<<" final="<<repeated.finalGeometryVerified<<" diagnostic="<<repeated.verificationDiagnostic<<'\n';
  if(repeated.finalEvaluation) std::cout<<"hash-debug geometry="<<(repeated.finalEvaluation->geometryHash==solved.finalEvaluation->geometryHash)<<" scientific="<<(repeated.finalEvaluation->scientificHash==solved.finalEvaluation->scientificHash)<<" raw="<<*repeated.independentlyRecomputedRawResidual<<" original="<<*solved.independentlyRecomputedRawResidual<<" moment-diff="<<(repeated.independentTotalMoments->values()-solved.independentTotalMoments->values()).norm()<<'\n';
  for(std::size_t i=0;i<std::min(repeated.numerical.trace.size(),solved.numerical.trace.size());++i) {
   const auto& a=solved.numerical.trace[i];const auto& b=repeated.numerical.trace[i];
   std::cout<<"trace-replay i="<<i<<" q="<<(a.quadricHash==b.quadricHash)<<" geometry="<<(a.materialGeometryHash==b.materialGeometryHash)<<" raw-a="<<*a.rawEquationResidual<<" raw-b="<<*b.rawEquationResidual<<'\n';
  }
 }
 require(repeated.finalGeometryVerified==solved.finalGeometryVerified && repeated.outputHash==solved.outputHash &&
             repeated.numerical.eventHash==solved.numerical.eventHash,
         "topology-transition output scientific identity unstable");
}
void strictFallbackSelector() {
 // Independent literal selector predicates. A thin rectangular solid with
 // q=x has center x=0; only its outward +x faces are interface witnesses.
 // These test the selector certificate ONLY, not a quadric/source/J seal.
 V linear=V::Zero();linear(1)=1;
 const auto plane=spin::QuadricCoefficients10::fromUnnormalized(linear);
 spin::QuadricTriangulation witness;
 witness.closedBoundary=box(-1,1);
 const double threshold=1e-10;
 for(int i=0;i<8;++i)
  witness.closedBoundary.vertices(i,0)=witness.closedBoundary.vertices(i,0)*threshold;
 witness.faceIsCap=Eigen::VectorXi::Ones(12);
 witness.faceIsCap(10)=witness.faceIsCap(11)=0;
 require(spin::quadricSelectorHasStrictInactiveCentroidMargin(plane,witness,false),
     "exact strict inactive selector witness rejected");
 require(!spin::quadricSelectorHasStrictInactiveCentroidMargin(plane,witness,true),
     "selector certificate accepted the wrong material side");
 auto boundary=witness;
 // Exact center x=threshold (x endpoints 0 and 2*threshold) MUST fail.
 for(int i=0;i<8;++i)
  boundary.closedBoundary.vertices(i,0)=witness.closedBoundary.vertices(i,0)<0 ? 0 : 2*threshold;
 require(!spin::quadricSelectorHasStrictInactiveCentroidMargin(plane,boundary,false),
     "selector availability equality granted a strict margin");
 auto mixed=witness;mixed.faceIsCap(8)=mixed.faceIsCap(9)=0;
 require(!spin::quadricSelectorHasStrictInactiveCentroidMargin(plane,mixed,false),
     "mixed interface signs granted selector stability");
 auto capOnly=witness;capOnly.faceIsCap.setOnes();
 require(!spin::quadricSelectorHasStrictInactiveCentroidMargin(plane,capOnly,false),
     "missing interface witnesses granted selector stability");
 auto degenerate=witness;
 degenerate.closedBoundary.faces(10,1)=degenerate.closedBoundary.faces(10,0);
 require(!spin::quadricSelectorHasStrictInactiveCentroidMargin(plane,degenerate,false),
     "zero-area interface granted selector stability");
 // Regression for a false sign-transition rejection: the actual selector
 // uses the strictly available oriented-face fallback, not sign(q(centroid)).
 auto r=request();auto es=settings();
 es.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
 Eigen::Matrix3d a;a<<.83,.03,-.04,.03,1.17,.06,-.04,.06,-.72;
 const Eigen::Vector3d center(.1,.2,.37717397599018909);
 V p;p(0)=center.dot(a*center)+.31;p.segment<3>(1)=-2*a*center;
 p.segment<3>(4)=a.diagonal();p(7)=a(1,2);p(8)=a(0,2);p(9)=a(0,1);
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
 const auto source=spin::evaluateQuadricSourceMotion(q,r.designDomain.vertices,es.triangulation);
 require(source.defined,"strict fallback fixture lacks regular native source motion");
 Eigen::Vector3d centroid;
 require(computeSurfaceMassMoments(source.triangulation.closedBoundary.vertices,
     source.triangulation.closedBoundary.faces,1).tryCenterOfMass(centroid,1e-15),
     "strict fallback fixture has no positive source mass");
 require(std::abs(q.evaluate(centroid))<1e-11,"fixture not inside inactive-centroid band");
 const auto native=spin::evaluateQuadricMaterialMoments(r,q,es);
 require(native.derivativeDefined && !native.quadricSolidIsPositiveRegion,
     "actual local map does not supply the expected negative-region derivative");
 const auto base=spin::evaluateFullRankMaterialMoments(r,q,es);
 require(base.coefficientJacobian.has_value(),
     "strict inactive-centroid fallback incorrectly withholds J: "+base.diagnostic);
 require(base.geometryHash==native.geometryHash &&
     base.totalMomentsTarget.values()==native.totalMomentsTarget.values(),
     "selector admission changed actual values or Boolean geometry");
 double previous=1;
 for(double h:{1e-4,3e-5,1e-5}) {
  Eigen::Matrix<double,10,10> fd;
  for(int col=0;col<10;++col){
   V plus=q.values(),minus=q.values();plus(col)+=h;minus(col)-=h;
   const auto vp=spin::evaluateQuadricMaterialMoments(r,
       spin::QuadricCoefficients10::fromUnnormalized(plus),es);
   const auto vm=spin::evaluateQuadricMaterialMoments(r,
       spin::QuadricCoefficients10::fromUnnormalized(minus),es);
   require(vp.quadricSolidIsPositiveRegion==native.quadricSolidIsPositiveRegion &&
       vm.quadricSolidIsPositiveRegion==native.quadricSolidIsPositiveRegion &&
       vp.provenanceTopologyHash==native.provenanceTopologyHash &&
       vm.provenanceTopologyHash==native.provenanceTopologyHash,
       "strict fallback differential crossed a material/source-feature branch");
   fd.col(col)=(vp.totalMomentsTarget.values()-vm.totalMomentsTarget.values())/(2*h);
  }
  const double error=(fd-*base.coefficientJacobian).norm()/base.coefficientJacobian->norm();
  std::cout<<"strict-selector-fallback h="<<h<<" relative="<<error<<'\n';
  require(error<=1e-6,"strict fallback actual Boolean derivative failed");
  if(h<1e-4)require(error<.3*previous,"strict fallback differences do not converge");
  previous=error;
 }
 const auto repeated=spin::evaluateFullRankMaterialMoments(r,q,es);
 require(repeated.coefficientJacobian && repeated.scientificHash==base.scientificHash,
     "strict fallback evaluation identity is nondeterministic");
}
void boundaryIdentity() {
 auto original=box(-1,1);
 auto reindexed=original;
 std::array<int,8> remap{};
 reindexed.vertices.resize(original.vertices.rows(),3);
 for(int old=0;old<original.vertices.rows();++old) {
  const int next=static_cast<int>(original.vertices.rows())-1-old;
  remap[static_cast<std::size_t>(old)]=next;
  reindexed.vertices.row(next)=original.vertices.row(old);
 }
 reindexed.faces.resize(original.faces.rows(),3);
 for(int old=0;old<original.faces.rows();++old) {
  const int next=static_cast<int>(original.faces.rows())-1-old;
  for(int corner=0;corner<3;++corner)
   reindexed.faces(next,corner)=remap[static_cast<std::size_t>(original.faces(old,(corner+1)%3))];
 }
 require(spin::sameOrientedTriangleBoundary(original,reindexed),
         "exact oriented boundary identity depends on mesh storage order");

 auto alternateDiagonal=original;
 alternateDiagonal.faces.row(0)<<0,3,1;
 alternateDiagonal.faces.row(1)<<1,3,2;
 const auto originalMoments=computeSurfaceMassMoments(original.vertices,original.faces,1.0);
 const auto alternateMoments=computeSurfaceMassMoments(alternateDiagonal.vertices,alternateDiagonal.faces,1.0);
 require((physics::RawMassMoments10::fromMassMoments(originalMoments).values()-
          physics::RawMassMoments10::fromMassMoments(alternateMoments).values()).norm()<1e-12,
         "coplanar diagonal control unexpectedly changed its mass moments");
 require(!spin::sameOrientedTriangleBoundary(original,alternateDiagonal),
         "boundary identity accepted a different triangulation with matching moments");

 auto perturbed=original;
 perturbed.vertices(0,0)+=1e-14;
 const auto perturbedMoments=computeSurfaceMassMoments(perturbed.vertices,perturbed.faces,1.0);
 require((physics::RawMassMoments10::fromMassMoments(originalMoments).values()-
          physics::RawMassMoments10::fromMassMoments(perturbedMoments).values()).norm()<1e-10,
         "coordinate perturbation is not below the former moment-only comparison threshold");
 require(!spin::sameOrientedTriangleBoundary(original,perturbed),
         "boundary identity accepted distinct binary64 geometry with nearly equal moments");

 auto reversed=original;
 std::swap(reversed.faces(0,1),reversed.faces(0,2));
 require(!spin::sameOrientedTriangleBoundary(original,reversed),
         "boundary identity ignored triangle orientation");
 std::cout<<"[PASS] exact oriented boundary identity ignores storage order but rejects map changes\n";
}
}
int main(int argc,char** argv) {try {
 std::cout<<std::setprecision(17);const std::string mode=argc>1?argv[1]:"adapter";
 if(mode=="selector-fallback") {
  strictFallbackSelector();std::cout<<"[PASS] "<<mode<<'\n';return 0;
 }
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
 if(mode=="rational-integration") {
  rationalIntegration();std::cout<<"[PASS] "<<mode<<'\n';return 0;
 }
#endif
 if(mode=="adapter")adapter();else if(mode=="chain")chain();else if(mode=="negative")negative();else if(mode=="solver")solver();else if(mode=="integration")integration();else if(mode=="modes")modes();else if(mode=="offtarget-yoyo")offTargetYoyo();else if(mode=="static-lowering-trace")staticLowerTrace();else if(mode=="near-parabolic-truncation")nearParabolicTruncationProbe();else if(mode=="audit")audit(false);else if(mode=="replay")audit(true);else if(mode=="boundary-identity")boundaryIdentity();else throw std::runtime_error("test not implemented: "+mode);
 std::cout<<"[PASS] "<<mode<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<'\n';return 1;}}
