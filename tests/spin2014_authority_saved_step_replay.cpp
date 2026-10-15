/** Frozen diagnostic replay only. Does not initialize a new optimization. */
#include "spin2014/QuadricAuthorityNewtonState.h"
#include "spin/DualMonomialCoordinates.h"
#include "ReplayResourceAdmission.h"
#include <nlohmann/json.hpp>
#include <gmpxx.h>
#include <fstream>
#include <iostream>
#include <stdexcept>

using Json=nlohmann::json;
static spin::QuadricCoefficientAuthority10::RationalVector vector10(const Json& values) {
 if(values.size()!=10)throw std::invalid_argument("incomplete frozen coefficient vector");
 spin::QuadricCoefficientAuthority10::RationalVector result;
 for(unsigned i=0;i<10;++i)result[i]=values.at(i).get<std::string>();
 return result;
}
static spin::TriangleMesh mesh(const Json& request,const char* vertexKey,const char* faceKey) {
 spin::TriangleMesh result;const auto& v=request.at(vertexKey);const auto& f=request.at(faceKey);
 result.vertices.resize(v.size(),3);result.faces.resize(f.size(),3);
 for(unsigned i=0;i<v.size();++i)for(unsigned j=0;j<3;++j)result.vertices(i,j)=v.at(i).at(j).get<double>();
 for(unsigned i=0;i<f.size();++i)for(unsigned j=0;j<3;++j)result.faces(i,j)=f.at(i).at(j).get<int>();
 return result;
}
int main(int argc,char** argv) {
 try {
  if(argc!=3 || (std::string(argv[2])!="--inspect" && std::string(argv[2])!="--actual" &&
                 std::string(argv[2])!="--actual-current-direction"))
      throw std::invalid_argument("usage: replay INPUT --inspect|--actual|--actual-current-direction");
  const bool currentDirection=std::string(argv[2])=="--actual-current-direction";
  if(std::string(argv[2])!="--inspect")replay::requireCurrentResourceAdmission();
  std::ifstream source(argv[1]);if(!source)throw std::invalid_argument("missing frozen fixture");
  Json input;source>>input;
  if(input.at("schema")!="spin2014-authority-saved-step-input-v1" || input.at("scientific_validation")!="NOT_GRANTED")
      throw std::invalid_argument("wrong diagnostic fixture scope");
  const auto& e=input.at("evaluator");spin::QuadricMomentEvaluatorSettings settings;
  if(e.at("derivative_method")!="DISCRETE_FULL_RANK_QUADRIC" || e.at("cap_radius_policy")!="REJECT_UNCERTIFIED" ||
      e.at("fixed_normal_paraboloid").get<bool>())throw std::invalid_argument("unsupported frozen diagnostic profile");
  settings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
  settings.minimumInterfaceGradient=e.at("minimum_interface_gradient").get<double>();
  settings.triangulation.icosphereSubdivisions=e.at("icosphere_subdivisions").get<int>();
  settings.triangulation.maximumIcosphereSubdivisions=e.at("maximum_icosphere_subdivisions").get<int>();
  settings.triangulation.maximumRadiusDoublings=e.at("maximum_radius_doublings").get<int>();
  settings.triangulation.truncationRadius=e.at("truncation_radius").get<double>();
  settings.triangulation.paraboloidGridResolution=e.at("paraboloid_grid_resolution").get<int>();
  settings.triangulation.cylinderAngularResolution=e.at("cylinder_angular_resolution").get<int>();
  settings.triangulation.hyperbolicCylinderGridResolution=e.at("hyperbolic_cylinder_grid_resolution").get<int>();
  settings.triangulation.parabolicCylinderGridResolution=e.at("parabolic_cylinder_grid_resolution").get<int>();
  settings.triangulation.retainFixedNormalProvenance=e.at("retain_fixed_normal_provenance").get<bool>();
  settings.triangulation.paperProjectiveHyperboloid=e.at("paper_projective_hyperboloid").get<bool>();
  settings.triangulation.algebraicFiniteMap=e.at("algebraic_finite_map").get<bool>();
  settings.triangulation.rationalSourceCoordinates=e.at("rational_source_coordinates").get<bool>();
  settings.triangulation.compactRationalSourceCoordinates=e.at("compact_rational_source_coordinates").get<bool>();
  settings.exactBinary64SelfIntersection=e.at("exact_binary64_self_intersection").get<bool>();
  settings.exactBinary64SourceBoundary=e.at("exact_binary64_source_boundary").get<bool>();
  settings.exactBinary64MaterialBoundary=e.at("exact_binary64_material_boundary").get<bool>();
  const auto q=vector10(input.at("source_coefficients"));
  const auto recorded=vector10(input.at("source_rejected_coefficients"));
  const auto& equation=input.at("source_linear_system");
  Eigen::MatrixXd matrix(15,15);Eigen::VectorXd rhs(15);
  for(int i=0;i<15;++i){rhs(i)=mpq_class(equation.at("right_hand_side").at(i).get<std::string>()).get_d();
    for(int j=0;j<15;++j)matrix(i,j)=mpq_class(equation.at("matrix").at(i).at(j).get<std::string>()).get_d();}
  const auto solved=spin::solveSpinNewtonLinearSystem(matrix,rhs);
  if(!solved.solved)throw std::runtime_error("archived protected equation cannot be restored: "+solved.diagnostic);
  spin::QuadricCoefficients10::Vector display;
  for(int i=0;i<10;++i)display(i)=mpq_class(q[i]).get_d();
  const auto legacy=spin::QuadricCoefficients10::fromNormalized(display);
  const auto coordinates=spin::DualMonomialCoordinates10::fromGeometry(legacy);
  const auto direction=coordinates.projectToTangent(solved.solution.head<10>());
  const auto arc=input.at("source_arc_length").get<std::string>();
  const auto rebuilt=coordinates.geodesicUpdate(direction,mpq_class(arc).get_d()).geometry();
  bool equal=true;double maximumDifference=0;
  Json directionTape=Json::array();
  spin::QuadricCoefficientAuthority10::RationalVector authorityDirection;
  for(int i=0;i<10;++i) {
    const double old=mpq_class(recorded[i]).get_d();
    equal=equal && rebuilt.values()(i)==old;
    maximumDifference=std::max(maximumDifference,std::abs(rebuilt.values()(i)-old));
    mpq_class value(direction(i));value.canonicalize();
    authorityDirection[i]=value.get_num().get_str()+"/"+value.get_den().get_str();
    directionTape.push_back(authorityDirection[i]);
  }
  Json result={{"scope","FROZEN_DIAGNOSTIC_REPLAY_NOT_NEW_OPTIMIZATION"},
      {"old_trial_coefficients_equal",equal},{"old_trial_maximum_difference",maximumDifference},
      {"source_direction",directionTape},{"source_arc_length",arc},{"actual_material_evaluations",0},
      {"scientific_validation","NOT_GRANTED"}};
  if(std::string(argv[2])=="--inspect"){std::cout<<result.dump()<<'\n';return equal?0:2;}
  if(!equal)throw std::runtime_error("old trial direction reconstruction mismatch; actual replay not authorized");
  const auto& frozen=input.at("frozen_input");const auto& r=frozen.at("request");
  spin::SpinDesignRequest request;request.mode=spin::SpinMode::SpinningTop;
  request.outerSurface=mesh(r,"outer_vertices","outer_faces");request.designDomain=mesh(r,"design_vertices","design_faces");
  request.materialDensity=r.at("material_density").get<double>();
  for(int i=0;i<3;++i){request.targetFrame.origin(i)=r.at("target_origin").at(i).get<double>();
    for(int j=0;j<3;++j)request.targetFrame.axes(i,j)=r.at("target_axes").at(i).at(j).get<double>();}
  const auto context=spin::makeQuadricAuthorityMaterialContext(request,mpq_class(input.at("fixed_material_density").get<std::string>()).get_d());
  for(unsigned i=0;i<10;++i)if(context.fixedMomentsTarget()[i]!=input.at("independently_integrated_shell").at(i).get<std::string>())
      throw std::runtime_error("geometry-owned fixed shell differs from independently checked source");
  spin2014::QuadricObjectiveSettings objective;const auto& o=frozen.at("objective");
  if(o.at("objective")!="TOP_INERTIA_ONLY")throw std::invalid_argument("wrong frozen objective");
  objective.objective=spin2014::Spin2014Objective::TopInertiaOnly;
  objective.inertiaWeight=o.at("inertia_weight").get<double>();objective.centerWeight=o.at("center_weight").get<double>();
  objective.declaredLengthUnit=o.at("declared_length_unit").get<std::string>();objective.declaredMassUnit=o.at("declared_mass_unit").get<std::string>();
  objective.solverMassUnitInDeclaredUnits=o.at("solver_mass_unit_in_declared_units").get<double>();
  objective.solverLengthUnitInDeclaredUnits=o.at("solver_length_unit_in_declared_units").get<double>();
  const auto authority=spin::QuadricCoefficientAuthority10::fromUnnormalizedGeometry(q);
  const auto base=spin2014::evaluateSpin2014AuthorityNewtonPoint(context,objective,authority,settings);
  auto trialDirection=authorityDirection;
  auto trialArc=arc;
  if(currentDirection) {
    const auto production=spin2014::solveSpin2014AuthorityNewtonDirection(objective,base.state,"1/1000000000000");
    if(!production.defined)throw std::runtime_error("current production direction unavailable: "+production.diagnostic);
    trialDirection=production.dualDirection;
    trialArc=production.initialArcLength;
    result["scope"]="SUPPLEMENTARY_CURRENT_PRODUCTION_DIRECTION_NOT_OLD_GATE";
    result["production_direction"]=trialDirection;
    result["production_arc_length"]=trialArc;
    result["production_direction_identity"]=production.scientificHash;
    result["original_linear_residual_upper"]=production.linearSolve.residualInfinityUpper;
    result["original_old_direction_gate"]="FAILED_UNCHANGED";
  }
  const auto trial=spin2014::evaluateSpin2014AuthorityNewtonTrial(context,objective,base,trialDirection,trialArc,settings);
  result["actual_material_evaluations"]=2;result["accepted_actual_material_trial"]=trial.accepted;
  result["raw_before"]=base.state.rawEquationResidual;result["raw_after"]=trial.point.state.rawEquationResidual;
  result["base_identity"]=base.state.scientificHash;result["trial_identity"]=trial.scientificHash;
  std::cout<<result.dump()<<'\n';return trial.accepted?0:3;
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
