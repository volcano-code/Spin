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
  spin::QuadricMomentEvaluatorSettings evaluation;
  evaluation.triangulation.icosphereSubdivisions=0;
  evaluation.triangulation.rationalSourceCoordinates=true;
  evaluation.triangulation.compactRationalSourceCoordinates=true;
  evaluation.triangulation.truncationRadius=8;evaluation.triangulation.maximumRadiusDoublings=0;
  evaluation.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
  evaluation.exactBinary64SourceBoundary=evaluation.exactBinary64MaterialBoundary=true;
  using A=spin::QuadricCoefficientAuthority10;
  A::RationalVector q;q.fill("0/1");q[0]="-1/2";q[1]="1/10";
  q[4]="1/2";q[5]="13/20";q[6]="2/5";
  const auto coefficient=A::fromUnnormalizedGeometry(q);
  spin2014::QuadricObjectiveSettings objective;
  objective.objective=spin2014::Spin2014Objective::TopMassAndLowering;
  const auto base=spin2014::evaluateSpin2014AuthorityNewtonPoint(context,objective,coefficient,evaluation);
  auto changedEvaluation=evaluation;
  changedEvaluation.minimumInterfaceGradient*=2;
  A::RationalVector recipeDirection;recipeDirection.fill("0/1");recipeDirection[4]="1/1";
  bool recipeRejected=false;
  try{(void)spin2014::evaluateSpin2014AuthorityNewtonTrial(context,objective,base,
      recipeDirection,"1/1000",changedEvaluation);}
  catch(const std::invalid_argument&){recipeRejected=true;}
  require(recipeRejected,"trial reused the base Jacobian under a changed evaluator recipe");
  require(mpq_class(base.state.multiplierMu)>0,"frozen attempt fixture has inadmissible initial multiplier");
  const auto oneCall=spin2014::runSpin2014AuthorityNewtonAttempt(context,objective,coefficient,evaluation,
      "1/1000000000000",1);
  require(oneCall.productionEvaluations==1 && oneCall.status=="evaluation_budget_exhausted" &&
      oneCall.trials.empty() && oneCall.acceptedUpdates==0 && !oneCall.stationaryCandidate.has_value(),
      "attempt did not count initial actual construction or bypassed one-call budget");
  require(oneCall.initialTotalMoments.has_value() &&
      *oneCall.initialTotalMoments==base.material.totalMomentsTarget &&
      oneCall.initialObjective==base.state.objective.objective &&
      oneCall.initialMaterialIdentity==base.material.scientificHash,
      "attempt omitted actual initial material/objective evidence");
  bool invalidBudget=false;
  try{(void)spin2014::runSpin2014AuthorityNewtonAttempt(context,objective,coefficient,evaluation,
      "1/1000000000000",61);}catch(const std::invalid_argument&){invalidBudget=true;}
  require(invalidBudget,"attempt accepted more than 60 evaluations");
  auto restartRequest=request;restartRequest.settings.maximumRestarts=2;
  const auto restartContext=spin::makeQuadricAuthorityMaterialContext(restartRequest,2);
  const auto restarts=spin2014::runSpin2014AuthorityNewtonRestarts(restartContext,objective,evaluation,
      "1/1000000000000",1);
  require(restarts.attempts.size()==2 && restarts.productionEvaluations==2 &&
      !restarts.stationaryAttemptIndex.has_value(),"restart ledger dropped a failed seeded attempt");
  require(restarts.attempts[0].initialGeometryAuthority!=restarts.attempts[1].initialGeometryAuthority,
      "restart reused the same seed start");
  A::RationalVector direction;direction.fill("0/1");direction[4]="1/1";
  const auto trial=spin2014::evaluateSpin2014AuthorityNewtonTrial(context,objective,base,direction,"1/1000",evaluation);
  require(trial.point.state.geometryAuthority!=base.state.geometryAuthority,"trial did not update authority coefficients");
  require(trial.point.material.material.geometryHash!=base.material.material.geometryHash &&
          trial.point.material.totalMomentsTarget!=base.material.totalMomentsTarget,"trial reused base Boolean/moments");
  require(trial.point.state.evaluationMapIdentity==trial.point.material.scientificHash &&
          trial.point.state.evaluationMapIdentity!=base.state.evaluationMapIdentity,"trial reused old Jacobian identity");
  const auto& m=trial.point.material.totalMomentsTarget;
  const double mass=mpq_class(m[0]).get_d(),x=mpq_class(m[1]).get_d(),y=mpq_class(m[2]).get_d(),z=mpq_class(m[3]).get_d();
  const double xx=mpq_class(m[4]).get_d(),yy=mpq_class(m[5]).get_d(),zz=mpq_class(m[6]).get_d();
  const double xy=mpq_class(m[9]).get_d();
  const double a=yy+zz-(y*y+z*z)/mass,b=xx+zz-(x*x+z*z)/mass;
  const double c=xx+yy-(x*x+y*y)/mass,d=-xy+x*y/mass;
  const double independent=(a*a+b*b+2*d*d)/(c*c)+z*z;
  require(std::abs(mpq_class(trial.point.state.objective.objective).get_d()-independent)<1e-12,
      "actual nonlinear objective did not match independent final moments");
  const bool decrease=mpq_class(trial.point.state.rawEquationResidual)<mpq_class(base.state.rawEquationResidual);
  const bool badStationary=mpq_class(trial.point.state.rawEquationResidual)<=mpq_class(1,1000000000000L) &&
      mpq_class(trial.point.state.multiplierMu)<=0;
  require(trial.accepted==(decrease && !badStationary),"trial acceptance not tied to actual nonlinear raw residual");
  const auto newtonDirection=spin2014::solveSpin2014AuthorityNewtonDirection(objective,base.state,"1/1000000000000");
  require(newtonDirection.defined,newtonDirection.diagnostic.c_str());
  const auto advance=spin2014::advanceSpin2014AuthorityNewtonPoint(context,objective,base,evaluation,
      "1/1000000000000",60);
  require(advance.accepted && advance.finalTrial.has_value(),"actual protected Newton backtracking did not reduce nonlinear residual");
  require(advance.trials.size()>1 && advance.trials.front().reason=="geometry_or_derivative_failure",
      "frozen-cap geometry failure was lost from backtracking ledger");
  require(advance.trials.front().diagnostic.find("cap triangle")!=std::string::npos,
      "actual source construction exception lost");
  require(advance.trials.back().accepted && advance.productionEvaluations==advance.trials.size(),
      "actual failed/successful calls not accounted");
  const auto exhausted=spin2014::advanceSpin2014AuthorityNewtonPoint(context,objective,base,evaluation,
      "1/1000000000000",1);
  require(!exhausted.accepted && exhausted.productionEvaluations==1 && exhausted.status=="evaluation_budget_exhausted",
      "failed first trial bypassed evaluation budget");
  auto stale=base;stale.state.evaluationMapIdentity="old-map";
  const auto rejectedValid=spin2014::advanceSpin2014AuthorityNewtonPoint(context,objective,base,evaluation,
      "1/1000000000000",2);
  require(!rejectedValid.accepted && rejectedValid.productionEvaluations==2 &&
      rejectedValid.trials.back().reason=="no_residual_decrease",
      "frozen second trial did not expose the valid-but-rejected material branch");
  require(!rejectedValid.finalTrial.has_value(),
      "rejected material escaped as a committable final trial");
  bool rejected=false;
  try{(void)spin2014::evaluateSpin2014AuthorityNewtonTrial(context,objective,stale,direction,"1/1000",evaluation);}
  catch(const std::exception&){rejected=true;}
  require(rejected,"stale material/state attachment accepted");
  rejected=false;
  try{(void)spin2014::evaluateSpin2014AuthorityNewtonTrial(context,objective,base,direction,"3/4",evaluation);}
  catch(const std::exception&){rejected=true;}
  require(rejected,"trial bypassed frozen maximum geodesic step");
  std::cout<<"actual authority trial Boolean, moments, objective, fresh Eq.14 and nonlinear decision passed\n";
  std::cout<<"accepted="<<trial.accepted<<"; reason="<<trial.reason<<"; trial hash="<<trial.scientificHash<<'\n';
  std::cout<<"actual Newton backtracking trials="<<advance.trials.size()<<"; hash="<<advance.scientificHash<<'\n';
  for(const auto& record:advance.trials)
      std::cout<<"arc="<<record.arcLength<<"; reason="<<record.reason<<"; accepted="<<record.accepted<<'\n';
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
