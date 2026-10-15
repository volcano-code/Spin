#include "spin/QuadricSourceMotion.h"
#include "spin/QuadricInverseGaussMap.h"
#include "spin/QuadricStratumTangent.h"
#include "spin/QuadricMomentEvaluator.h"
#include "spin/TruncationRootCertificate.h"
#include "spin/ExactSturmTarski.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include "manufacturing/ManufacturingGeometry.h"
#include <algorithm>
#include <array>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool b,const std::string& s) {if(!b)throw std::runtime_error(s);}
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
void checkReuseRationalSource() {
 // Break caught: derivative construction reruns the source map, or accepts a
 // stale/mutated Boolean operand instead of the actual authority construction.
 Eigen::MatrixXd domain(2,3);domain<<-1,-1,-1,1,1,1;
 spin::QuadricTriangulationSettings settings;
 settings.icosphereSubdivisions=1;settings.truncationRadius=8;
 settings.rationalSourceCoordinates=true;
 spin::QuadricCoefficients10::Vector values;
 for(bool hyperboloid:{false,true}) {
  if(hyperboloid)values<<-.2,.07,-.04,.03,1.2,.9,-1.1,.02,.03,-.04;
  else values<<-.25,0,0,0,1,1,1,0,0,0;
  const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
  const auto actual=spin::triangulateQuadric(q,domain,settings);
  const auto reused=spin::evaluateQuadricSourceMotion(q,domain,settings,actual);
  require(reused.defined,"same-operation rational source motion unavailable: "+reused.diagnostic);
  const auto rebuilt=spin::evaluateQuadricSourceMotion(q,domain,settings);
  require(rebuilt.defined && reused.scientificHash==rebuilt.scientificHash &&
      reused.provenance==rebuilt.provenance &&
      reused.triangulation.exactVertices==actual.exactVertices,
      "reuse changed scientific output or lost actual rational source");
  require(reused.velocity.size()==rebuilt.velocity.size(),"reuse lost source motion rows");
  for(std::size_t v=0;v<reused.velocity.size();++v)
   require(reused.velocity[v].isApprox(rebuilt.velocity[v],0),"reuse changed source Jacobian");
  const auto reject=[&](const spin::QuadricTriangulation& bad) {
   const auto result=spin::evaluateQuadricSourceMotion(q,domain,settings,bad);
   require(!result.defined && result.scientificHash.empty(),"stale authority construction was accepted");
  };
  auto bad=actual;bad.closedBoundary.vertices(0,0)+=.01;reject(bad);
  bad=actual;bad.exactVertices[0][0]="999/1";reject(bad);
  bad=actual;bad.coordinateLower[0][0]="-999/1";reject(bad);
  bad=actual;std::swap(bad.closedBoundary.faces(0,1),bad.closedBoundary.faces(0,2));reject(bad);
  bad=actual;bad.fixedNormals(0,0)+=.01;reject(bad);
  if(hyperboloid) {
   bad=actual;bad.sourceNormalEdgeByVertex[0]={{0,1}};reject(bad);
  }
  auto wrong=q.values();wrong(1)+=.01;
  require(!spin::evaluateQuadricSourceMotion(
      spin::QuadricCoefficients10::fromUnnormalized(wrong),domain,settings,actual).defined,
      "authority source from different coefficients was accepted");
  auto different=settings;different.icosphereSubdivisions=2;
  require(!spin::evaluateQuadricSourceMotion(q,domain,different,actual).defined,
      "authority source from different normal resolution was accepted");
  different=settings;different.truncationRadius=16;
  require(!spin::evaluateQuadricSourceMotion(q,domain,different,actual).defined,
      "authority source outside the requested radius budget was accepted");
 }
 std::cout<<"[PASS] same-operation rational source motion, bitwise scientific equivalence and authority mismatch guards\n";
}
#endif
void checkLargeRadiusTruncationAccuracy() {
 // Captured coefficients, using the already-supported level-2 map. Large
 // truncation radii magnify the radius-qualified seed's coordinate noise.
 spin::QuadricCoefficients10::Vector values;
 values<<.4091598840032713,-.36688397145905616,-.019533580749283127,
   -.3602578587235418,.14803212299816212,-.34789822850316837,
   -.25957403943858653,-.3788370952579597,.3838860742170953,.2580806535816809;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
 Eigen::MatrixXd domain(2,3);domain<<-.25,-.25,-.25,.25,.25,.25;
 spin::QuadricTriangulationSettings settings;
 settings.icosphereSubdivisions=2;
 settings.maximumIcosphereSubdivisions=2;
 settings.truncationRadius=221.70250336881628;
 settings.maximumRadiusDoublings=0;
 const auto base=spin::evaluateQuadricSourceMotion(q,domain,settings);
 require(base.defined,"large-radius source motion unavailable: "+base.diagnostic);
 for(double h:{3e-6,1e-6,3e-7}) {
  for(int j=0;j<10;++j) {
   double error2=0,scale2=0;
   auto plus=q.values(),minus=q.values();plus(j)+=h;minus(j)-=h;
   const auto p=spin::triangulateQuadric(spin::QuadricCoefficients10::fromUnnormalized(plus),domain,settings);
   const auto m=spin::triangulateQuadric(spin::QuadricCoefficients10::fromUnnormalized(minus),domain,settings);
   require(p.closedBoundary.faces==base.triangulation.closedBoundary.faces &&
           m.closedBoundary.faces==base.triangulation.closedBoundary.faces,
           "large-radius FD changed source topology");
   const Eigen::MatrixXd fd=(p.closedBoundary.vertices-m.closedBoundary.vertices)/(2*h);
   for(int i=0;i<fd.rows();++i) {
    error2+=(fd.row(i).transpose()-base.velocity[i].col(j)).squaredNorm();
    scale2+=base.velocity[i].col(j).squaredNorm();
   }
   const double relative=std::sqrt(error2/scale2);
   std::cout<<"large-radius source h="<<h<<" column="<<j<<" relative="<<relative<<'\n';
   require(relative<=1e-6,"large-radius truncation coordinate noise defeats the actual source derivative");
  }
 }
}
void checkRepresentedCoordinateSideWitness() {
 // Plane z and this represented triangle have exact cross product
 // (0,-eps,+eps), although the rounded double subtraction loses eps.
 spin::QuadricTriangulation triangulation;
 triangulation.closedBoundary.vertices.resize(3,3);
 triangulation.closedBoundary.vertices<<1,1,1,0,0,0,1e-17,0,0;
 triangulation.closedBoundary.faces.resize(1,3);
 triangulation.closedBoundary.faces<<0,1,2;
 triangulation.faceIsCap.resize(1);
 triangulation.faceIsCap<<0;
 spin::QuadricCoefficients10::Vector coefficients;
 coefficients<<0,0,0,1,0,0,0,0,0,0;
 const auto plane=spin::QuadricCoefficients10::fromUnnormalized(coefficients);
 require(!spin::quadricSolidInteriorIsPositiveRegion(plane,triangulation),
         "represented nonzero area must preserve the negative interior sign");
 std::swap(triangulation.closedBoundary.faces(0,1),triangulation.closedBoundary.faces(0,2));
 require(spin::quadricSolidInteriorIsPositiveRegion(plane,triangulation),
         "reversing the represented interface must reverse the selected side");
 triangulation.closedBoundary.faces<<0,1,2;
 const auto rejected=[&](const spin::QuadricCoefficients10& q,
                         const spin::QuadricTriangulation& t) {
   bool failed=false;
   try { (void)spin::quadricSolidInteriorIsPositiveRegion(q,t); }
   catch(const std::runtime_error&) { failed=true; }
   require(failed,"exact fallback accepted genuinely invalid or ambiguous interface");
 };
 auto collinear=triangulation;
 collinear.closedBoundary.vertices.row(2)<<1e-17,1e-17,1e-17;
 rejected(plane,collinear);
 coefficients<<1,0,0,0,0,0,0,0,0,0;
 rejected(spin::QuadricCoefficients10::fromUnnormalized(coefficients),triangulation);
 coefficients<<0,1,0,0,0,0,0,0,0,0;
 rejected(spin::QuadricCoefficients10::fromUnnormalized(coefficients),triangulation);
 auto mixed=triangulation;
 mixed.closedBoundary.faces.resize(2,3);
 mixed.closedBoundary.faces<<0,1,2,0,2,1;
 mixed.faceIsCap.resize(2); mixed.faceIsCap<<0,0;
 rejected(plane,mixed);
 // A well-conditioned +z normal with a nearly tangential plane gradient
 // checks both sides of the unchanged 1e-10 cosine ambiguity threshold.
 auto unitFace=triangulation;
 unitFace.closedBoundary.vertices<<0,0,0,1,0,0,0,1,0;
 coefficients<<0,1,0,0.5e-10,0,0,0,0,0,0;
 rejected(spin::QuadricCoefficients10::fromUnnormalized(coefficients),unitFace);
 coefficients<<0,1,0,2e-10,0,0,0,0,0,0;
 require(!spin::quadricSolidInteriorIsPositiveRegion(
             spin::QuadricCoefficients10::fromUnnormalized(coefficients),unitFace),
         "unchanged angle threshold rejected a non-ambiguous interface");
 std::cout<<"[PASS] exact represented-coordinate side witness and protected rejection guards\n";
}
// Captured rejected full Newton trial (evaluation 61) of the frozen centered
// axis seed-1 request. A nearly ideal normal produced a tiny interior root
// seed; it must not prevent constructing the certified sphere crossing.
void checkCapturedFullStepTruncationSourceMotion() {
 spin::QuadricCoefficients10::Vector coefficients;
 coefficients<<-1.8569137323145823e-36,-3.714895257844476e-18,
   -2.6983620047884542e-17,5.922545169177514e-18,
   0.39824375457455524,0.3982438519330578,-0.8263193972913151,
   -6.727145914928453e-18,7.966342917649522e-18,1.350754723917703e-19;
 Eigen::MatrixXd domain(8,3);
 domain<<-.25,-.25,-.25, -.25,-.25,.25, .25,-.25,.25,
   .25,-.25,-.25, -.25,.25,-.25, -.25,.25,.25,
   .25,.25,.25, .25,.25,-.25;
 spin::QuadricTriangulationSettings settings;
 settings.icosphereSubdivisions=5;
 settings.maximumIcosphereSubdivisions=-1;
 settings.maximumRadiusDoublings=8;
 const auto quadric=spin::QuadricCoefficients10::fromUnnormalized(coefficients);
 const auto motion=spin::evaluateQuadricSourceMotion(quadric,domain,settings);
 require(motion.defined,"captured full-step truncation source motion unavailable: "+motion.diagnostic);
 require(motion.triangulation.closedBoundary.vertices.rows()>0,
   "captured full-step source geometry is empty");
 require(motion.velocity.size()==static_cast<std::size_t>(
   motion.triangulation.closedBoundary.vertices.rows()),
   "captured full-step derivative coverage is incomplete");
 for(const auto& jac:motion.velocity)
   require(jac.allFinite(),"captured full-step derivative is nonfinite");
 // Check the emitted crossings independently of the root solver's status.
 // Cap centers are not on the quadric; only provenance-tagged crossings are.
 const auto& triangulation=motion.triangulation;
 require(triangulation.sourceNormalEdgeByVertex.size()==
   static_cast<std::size_t>(triangulation.closedBoundary.vertices.rows()),
   "captured full-step truncation provenance is incomplete");
 int crossings=0;
 for(int i=0;i<triangulation.interfaceVertexCount;++i) {
   const auto edge=triangulation.sourceNormalEdgeByVertex[static_cast<std::size_t>(i)];
   if(edge[0]<0 && edge[1]<0)continue;
   require(edge[0]>=0 && edge[1]>edge[0],"invalid captured crossing edge identity");
   const Eigen::Vector3d point=triangulation.closedBoundary.vertices.row(i).transpose();
   require(point.allFinite(),"captured crossing coordinates are nonfinite");
   require(std::abs(point.norm()-triangulation.truncationRadius)<=
     1e-10*std::max(1.0,triangulation.truncationRadius),
     "captured crossing is not on its truncation sphere");
   require(std::abs(quadric.evaluate(point))/std::max(1.0,point.squaredNorm())<=1e-12,
     "captured crossing is not on its quadric");
   ++crossings;
 }
 require(crossings>0,"captured trial did not exercise truncation crossings");
 std::cout<<"[PASS] captured full Newton step retains certified truncation source motion\n";
}
spin::TriangleMesh cube(double minimum,double maximum) {
 spin::TriangleMesh mesh;
 mesh.vertices.resize(8,3);
 mesh.vertices<<minimum,minimum,minimum, maximum,minimum,minimum,
     maximum,maximum,minimum, minimum,maximum,minimum,
     minimum,minimum,maximum, maximum,minimum,maximum,
     maximum,maximum,maximum, minimum,maximum,maximum;
 mesh.faces.resize(12,3);
 mesh.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,
     3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
 return mesh;
}
void check(const spin::QuadricCoefficients10::Vector& v,bool hyperboloid=true) {
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(v);
 Eigen::MatrixXd domain(2,3);domain<<-1,-1,-1,1,1,1;
 spin::QuadricTriangulationSettings settings;settings.icosphereSubdivisions=2;settings.truncationRadius=8;
 const auto base=spin::evaluateQuadricSourceMotion(q,domain,settings);
 require(base.defined,"source motion unavailable: "+base.diagnostic);
 if(hyperboloid && v(0)<0.0 && v(6)<0.0) {
  // Root precision refinement changes represented coordinates and therefore
  // the scientific identity (the pre-refinement hash remains in history).
  // The multi-step differences and independent re-evaluation below must also
  // pass; the hash alone cannot certify this numerical change.
  require(base.scientificHash=="9f3025e36befbcf65e6b9f3a724649c52d9fde06752a18df3ef1a77c30872ba6",
          "current inverse-Gauss/source-motion scientific hash changed: " +
              base.scientificHash);
 }
 require(base.triangulation.type==(hyperboloid?spin::QuadricSurfaceType::Hyperboloid:spin::QuadricSurfaceType::Ellipsoid),"unexpected source type");
 if(hyperboloid) {
  require(base.triangulation.sourceNormalEdgeByVertex.size()==
              static_cast<std::size_t>(base.triangulation.closedBoundary.vertices.rows()),
          "hyperboloid triangulation must retain the normal-sphere edge that created each truncation vertex");
  int edgeVertices=0;
  for(int vertex=0;vertex<base.triangulation.interfaceVertexCount;++vertex) {
   const auto edge=base.triangulation.sourceNormalEdgeByVertex[static_cast<std::size_t>(vertex)];
   if(edge[0]>=0 || edge[1]>=0) {
    require(edge[0]>=0 && edge[1]>edge[0],"truncation edge provenance must be canonical and complete");
    ++edgeVertices;
   }
  }
  require(edgeVertices>0,"hyperboloid fixture must expose truncation-edge provenance");
 }
 if(hyperboloid)require(base.triangulation.capFaceCount>0 && base.triangulation.boundaryLoopCount==2,"cap coverage missing");
 int edges=0,centers=0;
 for(const auto& s:base.provenance) {edges+=s.rfind("edge:",0)==0;centers+=s.rfind("cap:",0)==0;}
 require(hyperboloid?(edges>0 && centers==2):(edges==0 && centers==0),"source identity coverage wrong");
 double previous=1;
 for(double h:{3e-4,1e-4,3e-5,1e-5,3e-6}) {
  double error2=0,norm2=0;
  for(int j=0;j<10;++j) {
   auto p=q.values(),m=q.values();p(j)+=h;m(j)-=h;
   const auto plus=spin::triangulateQuadric(spin::QuadricCoefficients10::fromUnnormalized(p),domain,settings);
   const auto minus=spin::triangulateQuadric(spin::QuadricCoefficients10::fromUnnormalized(m),domain,settings);
   require(plus.closedBoundary.faces==base.triangulation.closedBoundary.faces && minus.closedBoundary.faces==base.triangulation.closedBoundary.faces,"FD topology changed");
   require(plus.truncationRadius==base.triangulation.truncationRadius && minus.truncationRadius==base.triangulation.truncationRadius,"FD radius branch changed");
   const Eigen::MatrixXd fd=(plus.closedBoundary.vertices-minus.closedBoundary.vertices)/(2*h);
   for(int i=0;i<fd.rows();++i) {error2+=(fd.row(i).transpose()-base.velocity[i].col(j)).squaredNorm();norm2+=base.velocity[i].col(j).squaredNorm();}
  }
  const double relative=std::sqrt(error2/norm2);
  std::cout<<"source h="<<h<<" relative="<<relative<<" edges="<<edges<<" caps="<<centers<<"\n";
  if(h<4e-6)require(relative<1e-6,"source finite difference error");
  if(h<2e-5 && previous>1e-8)require(relative<.2*previous,"FD refinement did not converge");
  previous=relative;
 }
 for(const auto& jac:base.velocity)require((jac*q.values()).norm()<1e-9,"normalization radial-null failed");
 const auto again=spin::evaluateQuadricSourceMotion(q,domain,settings);
 require(again.defined && again.scientificHash==base.scientificHash,"nondeterministic source binding");
 std::cout<<"[PASS] native source all coefficients, hyperboloid="<<hyperboloid<<"\n";
}

void checkNearConeSourceMotionUsesCertifiedExtendedSolve() {
 spin::QuadricCoefficients10::Vector coefficients;
 coefficients<<-1.1089457960913767e-10,3.1374659509406226e-08,
   -1.2306930695891306e-07,-1.4403414898221336e-07,
   .40315348093351583,.40315365134925046,-.8215439149561166,
   -3.6688686058282268e-06,1.8991362047213163e-07,-4.039908646624158e-09;
 Eigen::MatrixXd normalizedDesign(2,3);
 normalizedDesign<<-.5,-.5,-.5,.5,.5,.5;
 spin::QuadricTriangulationSettings settings;
 settings.icosphereSubdivisions=5;
 settings.maximumRadiusDoublings=4;
 const auto quadric=spin::QuadricCoefficients10::fromUnnormalized(coefficients);
 const auto motion=spin::evaluateQuadricSourceMotion(quadric,normalizedDesign,settings);
 require(motion.defined,
         "captured near-cone source motion must pass the original 1e-10 residual gate using extended precision: "+motion.diagnostic);
 require(motion.velocity.size()==
         static_cast<std::size_t>(motion.triangulation.closedBoundary.vertices.rows()),
         "near-cone source motion must provide a Jacobian for every emitted vertex");
 for(const auto& jacobian:motion.velocity)
   require(jacobian.allFinite(),"near-cone source motion contains nonfinite Jacobian entries");
 std::cout<<"[PASS] captured near-cone source motion with extended-precision residual certification\n";
}

void checkCherchiNewtonNearConeRetainsOriginalNormalIdentity() {
 spin::QuadricCoefficients10::Vector coefficients;
 coefficients<<-2.8364972803646595e-21,-8.11768576886572e-15,
   3.615646720485601e-16,-4.875767860908935e-16,
   .4031536902205084,.4031535225168989,-.821543875483107,
   -2.3102175102938004e-13,4.8400302740371675e-15,
   -2.4278465740028797e-17;
 Eigen::MatrixXd domain(8,3);
 domain<<-.25,-.25,-.25, .25,-.25,-.25, .25,.25,-.25, -.25,.25,-.25,
         -.25,-.25,.25, .25,-.25,.25, .25,.25,.25, -.25,.25,.25;
 spin::QuadricTriangulationSettings settings;
 settings.icosphereSubdivisions=4;
 settings.maximumRadiusDoublings=4;
 settings.retainFixedNormalProvenance=true;
 const auto geometry=spin::triangulateQuadric(
   spin::QuadricCoefficients10::fromUnnormalized(coefficients),domain,settings);
 spin::QuadricCoefficients10::Vector sphereCoefficients;
 sphereCoefficients<<1,0,0,0,-1,-1,-1,0,0,0;
 const auto referenceSphere=spin::triangulateQuadric(
   spin::QuadricCoefficients10::fromUnnormalized(sphereCoefficients),
   domain,settings);
 require(geometry.type==spin::QuadricSurfaceType::Hyperboloid &&
         geometry.boundaryLoopCount==2 &&
         geometry.sourceNormalVertexByVertex.size()==
           static_cast<std::size_t>(geometry.closedBoundary.vertices.rows()),
         "captured near-cone must preserve a source-normal vertex identity "
         "for every output vertex");
 require(geometry.sourceNormalEdgeByVertex.size()==
           geometry.sourceNormalVertexByVertex.size(),
         "normal-vertex and normal-edge provenance domains must agree");
 int originalCount=0,edgeCount=0;
 for(int vertex=0;vertex<geometry.interfaceVertexCount;++vertex) {
   const int normalVertex=geometry.sourceNormalVertexByVertex[
     static_cast<std::size_t>(vertex)];
   const auto edge=geometry.sourceNormalEdgeByVertex[
     static_cast<std::size_t>(vertex)];
   if(normalVertex>=0) {
     require(normalVertex<referenceSphere.fixedNormals.rows() &&
             edge[0]==-1 && edge[1]==-1,
             "original interface vertex has conflicting source identities");
     ++originalCount;
   } else {
     require(edge[0]>=0 && edge[1]>edge[0],
             "non-original interface vertex lacks a canonical source edge");
     ++edgeCount;
   }
 }
 require(originalCount>0 && edgeCount>0,
         "near-cone source mesh must contain both original and crossing vertices");
 for(std::size_t vertex=static_cast<std::size_t>(geometry.interfaceVertexCount);
     vertex<geometry.sourceNormalVertexByVertex.size();++vertex) {
   require(geometry.sourceNormalVertexByVertex[vertex]==-1 &&
           geometry.sourceNormalEdgeByVertex[vertex][0]==-1 &&
           geometry.sourceNormalEdgeByVertex[vertex][1]==-1,
           "cap center must not claim a normal-sphere source identity");
 }
 std::cout<<"[PASS] captured Cherchi Newton near-cone source-normal vertex identity\n";
}

void unavailable() {
 Eigen::MatrixXd domain(2,3);domain<<-1,-1,-1,1,1,1;
 spin::QuadricTriangulationSettings settings;settings.icosphereSubdivisions=2;settings.truncationRadius=1;
 // Keep the first inverse-Gauss sample strictly away from the radius-2
 // membership transition while still forcing the requested radius 1 to double.
 spin::QuadricCoefficients10::Vector p;p<<4.3,0,0,0,1,1,-1,0,0,0;
 auto result=spin::evaluateQuadricSourceMotion(spin::QuadricCoefficients10::fromUnnormalized(p),domain,settings);
 require(result.defined,"adaptive radius doubling source motion unavailable: "+result.diagnostic);
 require(result.triangulation.truncationRadius>settings.truncationRadius,
         "fixture must exercise adaptive radius doubling");
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
 double previous=1;
 for(double h:{3e-4,1e-4,3e-5,1e-5,3e-6}) {
  double error2=0,norm2=0;
  for(int j=0;j<10;++j) {
   auto plus=q.values(),minus=q.values();plus(j)+=h;minus(j)-=h;
   const auto a=spin::triangulateQuadric(spin::QuadricCoefficients10::fromUnnormalized(plus),domain,settings);
   const auto b=spin::triangulateQuadric(spin::QuadricCoefficients10::fromUnnormalized(minus),domain,settings);
   require(a.closedBoundary.faces==result.triangulation.closedBoundary.faces &&
           b.closedBoundary.faces==result.triangulation.closedBoundary.faces &&
           a.truncationRadius==result.triangulation.truncationRadius &&
           b.truncationRadius==result.triangulation.truncationRadius,
           "adaptive radius branch/topology changed during source finite difference");
   const Eigen::MatrixXd fd=(a.closedBoundary.vertices-b.closedBoundary.vertices)/(2*h);
   for(int i=0;i<fd.rows();++i) {
    error2+=(fd.row(i).transpose()-result.velocity[i].col(j)).squaredNorm();
    norm2+=result.velocity[i].col(j).squaredNorm();
   }
  }
  const double relative=std::sqrt(error2/norm2);
  std::cout<<"adaptive source h="<<h<<" relative="<<relative
           <<" effective_radius="<<result.triangulation.truncationRadius<<"\n";
  if(h<4e-6)require(relative<1e-6,"adaptive-radius source derivative error");
  if(h<2e-5 && previous>1e-8)require(relative<.2*previous,"adaptive-radius FD refinement did not converge");
  previous=relative;
 }
 std::cout<<"[PASS] adaptive truncation radius source-motion branch and all coefficients\n";
 Eigen::MatrixXd smallDomain(8,3);
 smallDomain<<-.5,-.5,-.5, .5,-.5,-.5, .5,.5,-.5, -.5,.5,-.5,
              -.5,-.5,.5, .5,-.5,.5, .5,.5,.5, -.5,.5,.5;
 spin::QuadricCoefficients10::Vector crossingValues;
 crossingValues<<.15,0,0,0,1,1,-1,0,0,0;
 spin::QuadricTriangulationSettings crossingSettings;
 crossingSettings.icosphereSubdivisions=2;
 crossingSettings.truncationRadius=.3;
 const auto crossing=spin::evaluateQuadricSourceMotion(
   spin::QuadricCoefficients10::fromUnnormalized(crossingValues),
   smallDomain,crossingSettings);
 require(crossing.defined,"source-motion failed for cap-separated design intersection: "+crossing.diagnostic);
 require(crossing.triangulation.truncationRadius>crossingSettings.truncationRadius,
         "design-intersecting fixture must take an adaptive radius step");
 std::cout<<"[PASS] nonempty-design hyperboloid uses effective radius "
          <<crossing.triangulation.truncationRadius<<"\n";
 p<<0,0,0,1,0,0,0,0,0,0;
 result=spin::evaluateQuadricSourceMotion(spin::QuadricCoefficients10::fromUnnormalized(p),domain,settings);
 require(!result.defined && result.velocity.empty() &&
         result.diagnostic.rfind("RANK_DEFICIENT_STRATUM_DERIVATIVE_UNDEFINED:",0)==0,
         "plane source was not explicitly classified");
 settings.icosphereSubdivisions=8;
 result=spin::evaluateQuadricSourceMotion(spin::QuadricCoefficients10::fromUnnormalized(p),domain,settings);
 require(!result.defined && result.velocity.empty(),"source budget exceeded silently");
 std::cout<<"[PASS] doubled radius, rank-deficient plane source and budget rejected\n";
}

void checkRankDeficientRejected(
    const spin::QuadricCoefficients10::Vector& values,
    const std::string& expectedType)
{
 Eigen::MatrixXd domain(2,3);domain<<-1,-1,-1,1,1,1;
 spin::QuadricTriangulationSettings settings;
 settings.truncationRadius=8;
 const auto result=spin::evaluateQuadricSourceMotion(
     spin::QuadricCoefficients10::fromUnnormalized(values),domain,settings);
 require(!result.defined && result.velocity.empty(),
         "rank-deficient source motion was silently accepted for "+expectedType);
 require(result.diagnostic.rfind(
     "RANK_DEFICIENT_STRATUM_DERIVATIVE_UNDEFINED:",0)==0,
     "rank-deficient source failure was not explicitly classified for "+expectedType+
     ": "+result.diagnostic);
 std::cout<<"[PASS] rank-deficient source explicitly rejected type="<<expectedType<<"\n";
}

void checkStratumBasis(
    const spin::QuadricCoefficients10::Vector& values,
    spin::QuadricSurfaceSubtype expectedSubtype,
    int expectedDimension)
{
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
 const auto tangent=spin::evaluateQuadricStratumTangent(q);
 require(tangent.defined && tangent.subtype==expectedSubtype,
         "fixed-stratum tangent classification mismatch");
 require(tangent.tangentDimension==expectedDimension,
         "fixed-stratum tangent dimension mismatch");
 require(tangent.basis.rows()==10 && tangent.basis.cols()==expectedDimension,
         "fixed-stratum basis has wrong shape");
 require((tangent.constraintRows*tangent.basis).norm()<1e-10,
         "fixed-stratum basis leaves coefficient tangent space");
 require((tangent.basis.transpose()*tangent.basis-
          Eigen::MatrixXd::Identity(expectedDimension,expectedDimension)).norm()<1e-10,
         "fixed-stratum basis is not orthonormal");
 std::cout<<"[PASS] explicit fixed-stratum tangent basis type="
          <<spin::quadricSurfaceSubtypeToken(expectedSubtype)
          <<" dimension="<<expectedDimension<<"\n";
}

void checkPlaneStratumSourceMotion() {
 Eigen::MatrixXd domain(8,3);
 domain<<-1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
         -1,-1,1, 1,-1,1, 1,1,1, -1,1,1;
 spin::QuadricCoefficients10::Vector values;
 values<<.17,.31,-.47,.826,0,0,0,0,0,0;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
 spin::QuadricTriangulationSettings settings;
 settings.truncationRadius=8.0;
 const auto base=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
 require(base.defined && base.coefficientBasis.rows()==10 &&
         base.coefficientBasis.cols()==3 && base.vertexVelocity.size()==8 &&
         base.provenance.size()==8,
         "plane fixed-stratum source motion unavailable: "+base.diagnostic);
 require(base.triangulation.type==spin::QuadricSurfaceType::Plane &&
         base.triangulation.interfaceVertexCount==4 &&
         base.triangulation.closedBoundary.vertices.rows()==8,
         "plane source-motion fixture did not use the expected native patch");
 double previous=1.0;
 for(double h:{1e-4,3e-5,1e-5,3e-6}) {
  double error2=0.0,norm2=0.0;
  for(int column=0;column<3;++column) {
   const Eigen::Matrix<double,10,1> direction=base.coefficientBasis.col(column);
   const auto plus=spin::QuadricCoefficients10::fromUnnormalized(q.values()+h*direction);
   const auto minus=spin::QuadricCoefficients10::fromUnnormalized(q.values()-h*direction);
   const auto a=spin::triangulateQuadric(plus,domain,settings);
   const auto b=spin::triangulateQuadric(minus,domain,settings);
   require(a.closedBoundary.faces==base.triangulation.closedBoundary.faces &&
           b.closedBoundary.faces==base.triangulation.closedBoundary.faces &&
           a.truncationRadius==base.triangulation.truncationRadius &&
           b.truncationRadius==base.triangulation.truncationRadius,
           "plane source finite difference changed topology/radius branch");
   const Eigen::MatrixXd fd=(a.closedBoundary.vertices-b.closedBoundary.vertices)/(2*h);
   for(Eigen::Index vertex=0;vertex<fd.rows();++vertex) {
    const Eigen::Vector3d expected=base.vertexVelocity[static_cast<std::size_t>(vertex)].col(column);
    error2+=(fd.row(vertex).transpose()-expected).squaredNorm();
    norm2+=expected.squaredNorm();
   }
  }
  const double relative=std::sqrt(error2/norm2);
  std::cout<<"plane stratum source h="<<h<<" relative="<<relative<<'\n';
  if(h==3e-6)require(relative<1e-6,
      "plane reduced source derivative failed independent finite differences");
  if(h<3e-5 && previous>1e-8)require(relative<.3*previous,
      "plane source finite-difference error did not refine");
  previous=relative;
 }
 const auto again=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
 require(again.defined && again.scientificHash==base.scientificHash,
         "plane stratum source motion is nondeterministic");
 std::cout<<"[PASS] plane source motion in reduced stratum coordinates\n";
}

void checkPlaneAxialChartSourceMotion() {
 Eigen::MatrixXd domain(2,3);domain<<-1,-1,-1,1,1,1;
 spin::QuadricCoefficients10::Vector values;
 values<<0.0,0.0,0.0,1.0,0,0,0,0,0,0;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
 spin::QuadricTriangulationSettings settings;settings.truncationRadius=8.0;
 const auto base=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
 require(base.defined && base.coefficientBasis.cols()==3 &&
         base.vertexVelocity.size()==8,
         "axis-aligned plane chart source motion unavailable: "+base.diagnostic);
 double previous=1.0;
 for(double h:{1e-4,3e-5,1e-5}) {
  double error2=0.0,norm2=0.0;
  for(int column=0;column<3;++column) {
   const Eigen::Matrix<double,10,1> direction=base.coefficientBasis.col(column);
   const auto plus=spin::QuadricCoefficients10::fromUnnormalized(q.values()+h*direction);
   const auto minus=spin::QuadricCoefficients10::fromUnnormalized(q.values()-h*direction);
   const auto a=spin::triangulateQuadric(plus,domain,settings);
   const auto b=spin::triangulateQuadric(minus,domain,settings);
   require(a.closedBoundary.faces==base.triangulation.closedBoundary.faces &&
           b.closedBoundary.faces==base.triangulation.closedBoundary.faces,
           "axis-aligned plane chart changed topology during finite differences");
   const Eigen::MatrixXd fd=(a.closedBoundary.vertices-b.closedBoundary.vertices)/(2*h);
   for(Eigen::Index vertex=0;vertex<fd.rows();++vertex) {
    const Eigen::Vector3d expected=base.vertexVelocity[static_cast<std::size_t>(vertex)].col(column);
    error2+=(fd.row(vertex).transpose()-expected).squaredNorm();
    norm2+=expected.squaredNorm();
   }
  }
  const double relative=std::sqrt(error2/norm2);
  std::cout<<"axis plane source h="<<h<<" relative="<<relative<<'\n';
  if(h==1e-5)require(relative<1e-6,
      "axis-aligned plane source differential failed independent finite differences");
  if(h<1e-4 && previous>1e-8)require(relative<.3*previous,
      "axis-aligned plane source finite-difference error did not refine");
  previous=relative;
 }
}

void checkPlaneChartSeamRejected() {
 Eigen::MatrixXd domain(2,3);domain<<-1,-1,-1,1,1,1;
 spin::QuadricCoefficients10::Vector values;
 values<<0.0,0.9,std::sqrt(0.19),0.0,0,0,0,0,0,0;
 const auto result=spin::evaluateQuadricStratumSourceMotion(
     spin::QuadricCoefficients10::fromUnnormalized(values),domain,{});
 require(!result.defined && result.vertexVelocity.empty() &&
         result.diagnostic.find("chart boundary")!=std::string::npos,
         "plane chart seam was assigned a local derivative");
}

Eigen::Matrix3d quadraticMatrix(const spin::QuadricCoefficients10::Vector& p) {
 Eigen::Matrix3d A=Eigen::Matrix3d::Zero();
 A.diagonal()<<p(4),p(5),p(6);
 A(1,2)=A(2,1)=p(7);A(0,2)=A(2,0)=p(8);A(0,1)=A(1,0)=p(9);
 return A;
}

spin::QuadricCoefficients10 retractCylinder(
    const spin::QuadricCoefficients10& q,
    const spin::QuadricCoefficients10::Vector& direction,double step) {
 spin::QuadricCoefficients10::Vector p=(q.values()+step*direction).eval();
 Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(quadraticMatrix(p));
 require(eigen.info()==Eigen::Success,"cylinder retraction eigensolve failed");
 int nullIndex=0;eigen.eigenvalues().cwiseAbs().minCoeff(&nullIndex);
 auto eigenvalues=eigen.eigenvalues();eigenvalues(nullIndex)=0;
 const Eigen::Matrix3d A=eigen.eigenvectors()*eigenvalues.asDiagonal()*eigen.eigenvectors().transpose();
 const Eigen::Vector3d n=eigen.eigenvectors().col(nullIndex);
 Eigen::Vector3d b=p.segment<3>(1);b-=n*n.dot(b);
 p(4)=A(0,0);p(5)=A(1,1);p(6)=A(2,2);
 p(7)=A(1,2);p(8)=A(0,2);p(9)=A(0,1);p.segment<3>(1)=b;
 const auto result=spin::QuadricCoefficients10::fromUnnormalized(p);
 const auto subtype=spin::classifyQuadricSurfaceSubtype(result);
 require(subtype==spin::QuadricSurfaceSubtype::EllipticCylinder,
         std::string("cylinder retraction left the elliptic-cylinder stratum: ")+
         spin::quadricSurfaceSubtypeToken(subtype));
 return result;
}

spin::QuadricCoefficients10 retractHyperbolicCylinder(
    const spin::QuadricCoefficients10& q,
    const spin::QuadricCoefficients10::Vector& direction,double step) {
 auto p=(q.values()+step*direction).eval();
 Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(quadraticMatrix(p));
 require(eigen.info()==Eigen::Success,"hyperbolic-cylinder retraction eigensolve failed");
 int nullIndex=0;eigen.eigenvalues().cwiseAbs().minCoeff(&nullIndex);
 auto eigenvalues=eigen.eigenvalues();eigenvalues(nullIndex)=0;
 const Eigen::Matrix3d A=eigen.eigenvectors()*eigenvalues.asDiagonal()*eigen.eigenvectors().transpose();
 const Eigen::Vector3d n=eigen.eigenvectors().col(nullIndex);
 Eigen::Vector3d b=p.segment<3>(1);b-=n*n.dot(b);
 p(4)=A(0,0);p(5)=A(1,1);p(6)=A(2,2);
 p(7)=A(1,2);p(8)=A(0,2);p(9)=A(0,1);p.segment<3>(1)=b;
 const auto result=spin::QuadricCoefficients10::fromUnnormalized(p);
 require(spin::classifyQuadricSurfaceSubtype(result)==spin::QuadricSurfaceSubtype::HyperbolicCylinder,
         "hyperbolic-cylinder retraction left its fixed stratum");
 return result;
}

spin::QuadricCoefficients10 retractParabolicCylinder(
    const spin::QuadricCoefficients10& q,
    const spin::QuadricCoefficients10::Vector& direction,double step) {
 auto p=(q.values()+step*direction).eval();
 Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(quadraticMatrix(p));
 require(eigen.info()==Eigen::Success,"parabolic-cylinder retraction eigensolve failed");
 int curvedIndex=0;eigen.eigenvalues().cwiseAbs().maxCoeff(&curvedIndex);
 Eigen::Vector3d eigenvalues=Eigen::Vector3d::Zero();
 eigenvalues(curvedIndex)=eigen.eigenvalues()(curvedIndex);
 const Eigen::Matrix3d A=eigen.eigenvectors()*eigenvalues.asDiagonal()*
     eigen.eigenvectors().transpose();
 p(4)=A(0,0);p(5)=A(1,1);p(6)=A(2,2);
 p(7)=A(1,2);p(8)=A(0,2);p(9)=A(0,1);
 const auto result=spin::QuadricCoefficients10::fromUnnormalized(p);
 require(spin::classifyQuadricSurfaceSubtype(result)==
             spin::QuadricSurfaceSubtype::ParabolicCylinder,
         "parabolic-cylinder retraction left the rank-one stratum");
 return result;
}

spin::QuadricCoefficients10 retractParaboloid(
    const spin::QuadricCoefficients10& q,
    const spin::QuadricCoefficients10::Vector& direction,double step) {
 auto p=(q.values()+step*direction).eval();
 Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(quadraticMatrix(p));
 require(eigen.info()==Eigen::Success,"paraboloid retraction eigensolve failed");
 int nullIndex=0;eigen.eigenvalues().cwiseAbs().minCoeff(&nullIndex);
 auto eigenvalues=eigen.eigenvalues();eigenvalues(nullIndex)=0;
 const Eigen::Matrix3d A=eigen.eigenvectors()*eigenvalues.asDiagonal()*eigen.eigenvectors().transpose();
 p(4)=A(0,0);p(5)=A(1,1);p(6)=A(2,2);
 p(7)=A(1,2);p(8)=A(0,2);p(9)=A(0,1);
 const auto result=spin::QuadricCoefficients10::fromUnnormalized(p);
 require(spin::classifyQuadricSurfaceSubtype(result)==spin::QuadricSurfaceSubtype::Paraboloid,
         "paraboloid retraction left the rank-two paraboloid stratum");
 return result;
}

std::vector<std::string> booleanTopologySignature(
    const spin::TriangleMesh& material,const Eigen::VectorXi& births,
    const std::vector<std::string>& provenance) {
 require(births.size()==material.faces.rows() &&
         provenance.size()==static_cast<std::size_t>(material.vertices.rows()),
         "Boolean topology signature received inconsistent provenance");
 std::vector<std::string> faces;
 faces.reserve(static_cast<std::size_t>(material.faces.rows()));
 for(int f=0;f<material.faces.rows();++f) {
  std::array<std::string,3> keys{{
      provenance[static_cast<std::size_t>(material.faces(f,0))],
      provenance[static_cast<std::size_t>(material.faces(f,1))],
      provenance[static_cast<std::size_t>(material.faces(f,2))]}};
  std::sort(keys.begin(),keys.end());
  faces.push_back(keys[0]+";"+keys[1]+";"+keys[2]+";"+std::to_string(births(f)));
 }
 std::sort(faces.begin(),faces.end());
 return faces;
}

void checkPlaneStratumBooleanDifferential() {
 spin::SpinDesignRequest request;
 request.outerSurface=cube(-1.5,1.5);
 request.designDomain=cube(-1.0,1.0);
 request.materialDensity=1.0;
 spin::QuadricCoefficients10::Vector values;
 values<<.17,.31,-.47,.826,0,0,0,0,0,0;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
 spin::QuadricMomentEvaluatorSettings settings;
 settings.derivativeMethod=spin::QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
 settings.triangulation.truncationRadius=8.0;
 const auto base=spin::evaluateQuadricMaterialMoments(request,q,settings);
 const auto differential=spin::differentiateDiscreteQuadricStratumMoments(
     request.designDomain,base.quadricTriangulation,q,base.materialBoundaryTarget,
     base.materialBirthFaces,request.materialDensity,settings.triangulation);
 require(differential.defined && differential.coefficientBasis.rows()==10 &&
         differential.coefficientBasis.cols()==3 && differential.jacobian.rows()==10 &&
         differential.jacobian.cols()==3,
         "plane fixed-stratum Boolean differential unavailable: "+differential.diagnostic);
 require(differential.constructedIntersectionVertexCount>0,
         "plane Boolean fixture did not exercise moving intersection vertices");
 const auto baseTopology=booleanTopologySignature(base.materialBoundaryTarget,
     base.materialBirthFaces,differential.vertexProvenance);
 double previous=1.0;
 for(double h:{1e-4,3e-5,1e-5}) {
  Eigen::Matrix<double,10,3> fd=Eigen::Matrix<double,10,3>::Zero();
  for(int column=0;column<3;++column) {
   const Eigen::Matrix<double,10,1> direction=differential.coefficientBasis.col(column);
   const auto plus=spin::QuadricCoefficients10::fromUnnormalized(q.values()+h*direction);
   const auto minus=spin::QuadricCoefficients10::fromUnnormalized(q.values()-h*direction);
   const auto a=spin::evaluateQuadricMaterialMoments(request,plus,settings);
   const auto b=spin::evaluateQuadricMaterialMoments(request,minus,settings);
   require(a.materialBoundaryTarget.faces.rows()==base.materialBoundaryTarget.faces.rows() &&
           b.materialBoundaryTarget.faces.rows()==base.materialBoundaryTarget.faces.rows(),
           "plane Boolean topology changed during reduced finite difference");
   const auto da=spin::differentiateDiscreteQuadricStratumMoments(
       request.designDomain,a.quadricTriangulation,plus,a.materialBoundaryTarget,
       a.materialBirthFaces,request.materialDensity,settings.triangulation);
   const auto db=spin::differentiateDiscreteQuadricStratumMoments(
       request.designDomain,b.quadricTriangulation,minus,b.materialBoundaryTarget,
       b.materialBirthFaces,request.materialDensity,settings.triangulation);
   require(da.defined && db.defined &&
           booleanTopologySignature(a.materialBoundaryTarget,a.materialBirthFaces,
                                    da.vertexProvenance)==baseTopology &&
           booleanTopologySignature(b.materialBoundaryTarget,b.materialBirthFaces,
                                    db.vertexProvenance)==baseTopology,
           "plane Boolean provenance topology changed during reduced finite difference");
   fd.col(column)=(a.variableMomentsTarget.values()-b.variableMomentsTarget.values())/(2*h);
  }
  const double relative=(fd-differential.jacobian).norm()/std::max(1e-30,fd.norm());
  std::cout<<"plane reduced Boolean moments h="<<h<<" relative="<<relative<<'\n';
  if(h==1e-5)require(relative<1e-6,
      "plane reduced Boolean-moment derivative failed independent finite differences");
  if(h<1e-4 && previous>1e-8)require(relative<.3*previous,
      "plane Boolean-moment finite-difference error did not refine");
  previous=relative;
 }
}

void checkParaboloidStratumSourceMotion() {
 spin::QuadricCoefficients10::Vector values;values<<0,0,0,-1,1,2,0,0,0,0;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
 Eigen::MatrixXd domain(8,3);
 domain<<-2,-2,-2, 2,-2,-2, 2,2,-2, -2,2,-2,
         -2,-2,2, 2,-2,2, 2,2,2, -2,2,2;
 spin::QuadricTriangulationSettings settings;
 settings.paraboloidGridResolution=16;settings.truncationRadius=2;
 const auto motion=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
 require(motion.defined,"fixed-paraboloid source motion unavailable: "+motion.diagnostic);
 require(motion.coefficientBasis.rows()==10 && motion.coefficientBasis.cols()==8,
         "paraboloid source motion must retain its reduced 8D chart");
 require(motion.triangulation.type==spin::QuadricSurfaceType::Paraboloid &&
         motion.triangulation.interfaceVertexCount==17*17 &&
         motion.triangulation.closedBoundary.vertices.rows()==2*17*17,
         "paraboloid native grid/extruded-cap vertex layout changed");
 require(motion.vertexVelocity.size()==static_cast<std::size_t>(2*17*17),
         "paraboloid source vertex differential count mismatch");
 double previous=1;
 for(double h:{1e-4,3e-5,1e-5,3e-6}) {
  double error2=0,norm2=0;
  for(int j=0;j<motion.coefficientBasis.cols();++j) {
   const auto direction=motion.coefficientBasis.col(j).eval();
   const auto plus=spin::triangulateQuadric(retractParaboloid(q,direction,h),domain,settings);
   const auto minus=spin::triangulateQuadric(retractParaboloid(q,direction,-h),domain,settings);
   require(plus.closedBoundary.faces==motion.triangulation.closedBoundary.faces &&
           minus.closedBoundary.faces==motion.triangulation.closedBoundary.faces,
           "paraboloid finite-difference topology changed");
   const Eigen::MatrixXd fd=(plus.closedBoundary.vertices-minus.closedBoundary.vertices)/(2*h);
   for(int i=0;i<fd.rows();++i) {
    const Eigen::Vector3d predicted=motion.vertexVelocity[static_cast<std::size_t>(i)].col(j);
    error2+=(fd.row(i).transpose()-predicted).squaredNorm();
    norm2+=predicted.squaredNorm();
   }
  }
  const double relative=std::sqrt(error2/std::max(1e-30,norm2));
  std::cout<<"paraboloid stratum source h="<<h<<" relative="<<relative<<'\n';
  if(h<4e-6)require(relative<1e-6,
      "paraboloid reduced source differential failed independent finite differences");
  if(h<2e-5 && previous>1e-8)require(relative<.3*previous,
      "paraboloid source finite-difference error did not refine");
  previous=relative;
 }
 // Exercise the same reduced directions through the actual CGAL Boolean
 // material map, rather than validating only source-vertex motion.
 spin::SpinDesignRequest clippedRequest;
 clippedRequest.outerSurface=cube(-2.0,2.0);
 clippedRequest.designDomain=cube(-.8,.8);
 const Eigen::Vector3d domainOffset(.113,-.071,.057);
 clippedRequest.designDomain.vertices.rowwise()+=domainOffset.transpose();
 clippedRequest.materialDensity=1.0;
 spin::QuadricCoefficients10::Vector clippedValues;
 clippedValues<<-.2,0,0,-1,1,2,0,0,0,0;
 const auto clippedQ=spin::QuadricCoefficients10::fromUnnormalized(clippedValues);
 spin::QuadricMomentEvaluatorSettings evaluatorSettings;
 evaluatorSettings.derivativeMethod=spin::QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
 evaluatorSettings.triangulation=settings;
 const auto clippedBase=spin::evaluateQuadricMaterialMoments(
     clippedRequest,clippedQ,evaluatorSettings);
 const auto clippedDifferential=spin::differentiateDiscreteQuadricStratumMoments(
     clippedRequest.designDomain,clippedBase.quadricTriangulation,clippedQ,
     clippedBase.materialBoundaryTarget,clippedBase.materialBirthFaces,
     clippedRequest.materialDensity,evaluatorSettings.triangulation);
 require(clippedDifferential.defined && clippedDifferential.coefficientBasis.rows()==10 &&
         clippedDifferential.coefficientBasis.cols()==8 &&
         clippedDifferential.jacobian.rows()==10 && clippedDifferential.jacobian.cols()==8,
         "clipped paraboloid Boolean-moment differential unavailable: "+clippedDifferential.diagnostic);
 require(clippedDifferential.constructedIntersectionVertexCount>0,
         "paraboloid Boolean fixture did not exercise moving intersection vertices");
 const auto baseTopology=booleanTopologySignature(clippedBase.materialBoundaryTarget,
     clippedBase.materialBirthFaces,clippedDifferential.vertexProvenance);
 previous=1;
 for(double h:{1e-4,3e-5,1e-5}) {
  Eigen::Matrix<double,10,8> finiteDifference=Eigen::Matrix<double,10,8>::Zero();
  for(int j=0;j<8;++j) {
   const auto direction=clippedDifferential.coefficientBasis.col(j).eval();
   const auto plus=spin::evaluateQuadricMaterialMoments(
       clippedRequest,retractParaboloid(clippedQ,direction,h),evaluatorSettings);
   const auto minus=spin::evaluateQuadricMaterialMoments(
       clippedRequest,retractParaboloid(clippedQ,direction,-h),evaluatorSettings);
   const auto plusDifferential=spin::differentiateDiscreteQuadricStratumMoments(
       clippedRequest.designDomain,plus.quadricTriangulation,
       retractParaboloid(clippedQ,direction,h),plus.materialBoundaryTarget,
       plus.materialBirthFaces,clippedRequest.materialDensity,evaluatorSettings.triangulation);
   const auto minusDifferential=spin::differentiateDiscreteQuadricStratumMoments(
       clippedRequest.designDomain,minus.quadricTriangulation,
       retractParaboloid(clippedQ,direction,-h),minus.materialBoundaryTarget,
       minus.materialBirthFaces,clippedRequest.materialDensity,evaluatorSettings.triangulation);
   require(plusDifferential.defined && minusDifferential.defined &&
           booleanTopologySignature(plus.materialBoundaryTarget,plus.materialBirthFaces,
                                   plusDifferential.vertexProvenance)==baseTopology &&
           booleanTopologySignature(minus.materialBoundaryTarget,minus.materialBirthFaces,
                                   minusDifferential.vertexProvenance)==baseTopology,
           "paraboloid Boolean topology changed during reduced finite difference");
   finiteDifference.col(j)=(plus.variableMomentsTarget.values()-
       minus.variableMomentsTarget.values())/(2*h);
  }
  const double relative=(finiteDifference-clippedDifferential.jacobian).norm()/
      std::max(1e-30,finiteDifference.norm());
  std::cout<<"clipped paraboloid Boolean moment h="<<h<<" relative="<<relative
           <<" constructed="<<clippedDifferential.constructedIntersectionVertexCount<<'\n';
  if(h==1e-5)require(relative<1e-6,
      "paraboloid reduced Boolean-moment derivative failed independent finite differences");
  if(h<1e-4 && previous>1e-8)require(relative<.3*previous,
      "paraboloid Boolean-moment finite-difference error did not refine");
  previous=relative;
 }
 std::cout<<"[PASS] clipped paraboloid Boolean moments in reduced stratum coordinates\n";
 const auto again=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
 require(again.defined && again.scientificHash==motion.scientificHash,
         "paraboloid stratum source motion is nondeterministic");
 std::cout<<"[PASS] paraboloid source motion in reduced stratum coordinates\n";
}

void checkParaboloidEigenvectorPivotTieRejected() {
 // The null eigenvector has two equally large components with opposite
 // signs.  The native sign convention selects the largest-magnitude pivot;
 // an arbitrarily small perturbation can switch that pivot and flip the
 // canonical eigenvector.  A local differential must not cross that chart
 // boundary silently.
 const double inverseSqrtTwo=1.0/std::sqrt(2.0);
 const Eigen::Vector3d u(inverseSqrtTwo,inverseSqrtTwo,0.0);
 const Eigen::Vector3d n(inverseSqrtTwo,-inverseSqrtTwo,0.0);
 const Eigen::Vector3d v(0.0,0.0,1.0);
 const Eigen::Matrix3d A=u*u.transpose()+2.0*v*v.transpose();
 spin::QuadricCoefficients10::Vector values;
 values<<0.0,n.x(),n.y(),n.z(),A(0,0),A(1,1),A(2,2),A(1,2),A(0,2),A(0,1);
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
 Eigen::MatrixXd domain(8,3);
 domain<<-2,-2,-2, 2,-2,-2, 2,2,-2, -2,2,-2,
         -2,-2,2, 2,-2,2, 2,2,2, -2,2,2;
 spin::QuadricTriangulationSettings settings;
 settings.paraboloidGridResolution=8;
 settings.truncationRadius=2.0;
 const auto motion=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
 require(!motion.defined && motion.vertexVelocity.empty() &&
         motion.diagnostic.find("eigenvector sign-pivot boundary")!=std::string::npos,
         "paraboloid derivative silently accepted an eigenvector sign-pivot tie: "+
         motion.diagnostic);
 std::cout<<"[PASS] paraboloid eigenvector sign-pivot boundary rejected\n";

 // Equal-magnitude components with the same sign do not switch the canonical
 // sign, so that harmless chart tie must remain admissible.
 const Eigen::Vector3d sameSignNull(inverseSqrtTwo,inverseSqrtTwo,0.0);
 const Eigen::Vector3d sameSignU(2.0,-2.0,3.0);
 const Eigen::Vector3d sameSignV(3.0,-3.0,-4.0);
 const Eigen::Vector3d normalizedU=sameSignU.normalized();
 const Eigen::Vector3d normalizedV=sameSignV.normalized();
 require(std::abs(normalizedU.dot(normalizedV))<1e-14 &&
         std::abs(normalizedU.dot(sameSignNull))<1e-14 &&
         std::abs(normalizedV.dot(sameSignNull))<1e-14,
         "same-sign pivot control basis is not orthonormal");
 const Eigen::Matrix3d sameSignA=normalizedU*normalizedU.transpose()+
     2.0*normalizedV*normalizedV.transpose();
 spin::QuadricCoefficients10::Vector sameSignValues;
 sameSignValues<<0.0,sameSignNull.x(),sameSignNull.y(),sameSignNull.z(),
     sameSignA(0,0),sameSignA(1,1),sameSignA(2,2),sameSignA(1,2),
     sameSignA(0,2),sameSignA(0,1);
 const auto sameSignMotion=spin::evaluateQuadricStratumSourceMotion(
     spin::QuadricCoefficients10::fromUnnormalized(sameSignValues),domain,settings);
 require(sameSignMotion.defined && sameSignMotion.coefficientBasis.cols()==8,
         "same-sign eigenvector magnitude tie was conservatively rejected: "+
         sameSignMotion.diagnostic);
 std::cout<<"[PASS] same-sign eigenvector magnitude tie remains admissible\n";

 // The same canonical-sign branch is used by an elliptic cylinder; exercise
 // its null-axis eigenvector tie independently as well.
 spin::QuadricCoefficients10::Vector cylinderValues;
 const Eigen::Vector3d cylinderB=u;
 cylinderValues<<-1.0,cylinderB.x(),cylinderB.y(),cylinderB.z(),
     A(0,0),A(1,1),A(2,2),A(1,2),A(0,2),A(0,1);
 const auto cylinder=spin::QuadricCoefficients10::fromUnnormalized(cylinderValues);
 const auto cylinderMotion=spin::evaluateQuadricStratumSourceMotion(
     cylinder,domain,settings);
 require(!cylinderMotion.defined && cylinderMotion.vertexVelocity.empty() &&
         cylinderMotion.diagnostic.find("eigenvector sign-pivot boundary")!=std::string::npos,
         "cylinder derivative silently accepted an eigenvector sign-pivot tie: "+
         cylinderMotion.diagnostic);
 std::cout<<"[PASS] cylinder eigenvector sign-pivot boundary rejected\n";
}

void checkCylinderStratumSourceMotion() {
 spin::QuadricCoefficients10::Vector values;values<<-1,0,0,0,1,2,0,0,0,0;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
 Eigen::MatrixXd domain(8,3);
 domain<<-2,-2,-2, 2,-2,-2, 2,2,-2, -2,2,-2,
         -2,-2,2, 2,-2,2, 2,2,2, -2,2,2;
 spin::QuadricTriangulationSettings settings;
 settings.cylinderAngularResolution=24;settings.truncationRadius=2;
 const auto motion=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
 require(motion.defined,"fixed-cylinder source motion unavailable: "+motion.diagnostic);
 require(motion.coefficientBasis.rows()==10 && motion.coefficientBasis.cols()==7,
         "cylinder source motion did not expose its reduced 7D chart");
 require(motion.vertexVelocity.size()==static_cast<std::size_t>(motion.triangulation.closedBoundary.vertices.rows()),
         "cylinder source vertex motion count mismatch");
 require(motion.triangulation.type==spin::QuadricSurfaceType::EllipticCylinder,
         "cylinder source motion produced another geometry type");
 spin::TriangleMesh design;
 design.vertices.resize(8,3);design.vertices<<-5,-5,-5,5,-5,-5,5,5,-5,-5,5,-5,
     -5,-5,5,5,-5,5,5,5,5,-5,5,5;
 design.faces.resize(12,3);design.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,
     3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
 Eigen::VectorXi birth(motion.triangulation.closedBoundary.faces.rows());
 for(int f=0;f<birth.size();++f)birth(f)=design.faces.rows()+f;
 const auto sourceAudit=manufacturing::auditTriangleBoundary(motion.triangulation.closedBoundary.vertices,
     motion.triangulation.closedBoundary.faces,manufacturing::BoundaryOrientation::AnyNonzero);
 require(sourceAudit.strictBoundaryValid(),
         "native cylinder source must be a finite positive closed oriented boundary");
 const auto restricted=spin::differentiateDiscreteQuadricStratumMoments(
     design,motion.triangulation,q,motion.triangulation.closedBoundary,
     birth,1.0,settings);
 require(restricted.defined && restricted.coefficientBasis.rows()==10 &&
         restricted.coefficientBasis.cols()==7 && restricted.jacobian.rows()==10 &&
         restricted.jacobian.cols()==7,
         "restricted cylinder Boolean-moment differential unavailable: "+restricted.diagnostic);
 double previous=1;
 for(double h:{1e-4,3e-5,1e-5,3e-6}) {
  double error2=0,norm2=0;
  for(int j=0;j<motion.coefficientBasis.cols();++j) {
   const auto direction=motion.coefficientBasis.col(j).eval();
   const auto plus=spin::triangulateQuadric(retractCylinder(q,direction,h),domain,settings);
   const auto minus=spin::triangulateQuadric(retractCylinder(q,direction,-h),domain,settings);
   require(plus.closedBoundary.faces==motion.triangulation.closedBoundary.faces &&
           minus.closedBoundary.faces==motion.triangulation.closedBoundary.faces,
           "cylinder finite-difference topology changed");
   const Eigen::MatrixXd fd=(plus.closedBoundary.vertices-minus.closedBoundary.vertices)/(2*h);
   for(int i=0;i<fd.rows();++i) {
    const Eigen::Vector3d predicted=motion.vertexVelocity[static_cast<std::size_t>(i)].col(j);
    error2+=(fd.row(i).transpose()-predicted).squaredNorm();
    norm2+=predicted.squaredNorm();
   }
  }
  const double relative=std::sqrt(error2/std::max(1e-30,norm2));
  std::cout<<"cylinder stratum source h="<<h<<" relative="<<relative<<"\n";
  if(h<4e-6)require(relative<1e-6,"cylinder stratum source derivative failed finite difference");
  if(h<2e-5 && previous>1e-8)require(relative<.2*previous,"cylinder stratum derivative refinement failed");
  previous=relative;
 }
 // Independently finite-difference all ten signed-tetrahedron moments of the
 // native closed source mesh along the same exact-stratum retractions.
 for(double h:{1e-4,3e-5,1e-5}) {
  Eigen::Matrix<double,10,7> fd=Eigen::Matrix<double,10,7>::Zero();
  for(int j=0;j<7;++j) {
   const auto direction=motion.coefficientBasis.col(j).eval();
   const auto plus=spin::triangulateQuadric(retractCylinder(q,direction,h),domain,settings);
   const auto minus=spin::triangulateQuadric(retractCylinder(q,direction,-h),domain,settings);
   Eigen::Matrix<double,10,1> mp=Eigen::Matrix<double,10,1>::Zero();
   Eigen::Matrix<double,10,1> mm=Eigen::Matrix<double,10,1>::Zero();
   for(const auto& item:std::array<std::pair<const spin::TriangleMesh*,Eigen::Matrix<double,10,1>*>,2>{{
       {&plus.closedBoundary,&mp},{&minus.closedBoundary,&mm}}}) {
    for(int f=0;f<item.first->faces.rows();++f) {
     const Eigen::Vector3d a=item.first->vertices.row(item.first->faces(f,0));
     const Eigen::Vector3d b=item.first->vertices.row(item.first->faces(f,1));
     const Eigen::Vector3d c=item.first->vertices.row(item.first->faces(f,2));
     const double volume=a.dot(b.cross(c))/6;const Eigen::Vector3d sum=a+b+c;
     const Eigen::Matrix3d second=volume*(sum*sum.transpose()+a*a.transpose()+
         b*b.transpose()+c*c.transpose())/20;
     (*item.second)(0)+=volume;item.second->segment<3>(1)+=volume*sum/4;
     item.second->segment<3>(4)+=second.diagonal();
     (*item.second)(7)+=second(1,2);(*item.second)(8)+=second(0,2);(*item.second)(9)+=second(0,1);
    }
   }
   fd.col(j)=(mp-mm)/(2*h);
  }
  const double relative=(fd-restricted.jacobian).norm()/std::max(1e-30,fd.norm());
  std::cout<<"cylinder stratum moment h="<<h<<" relative="<<relative<<"\n";
  if(h==1e-5)require(relative<1e-6,"restricted cylinder discrete moments failed finite difference");
 }
 // Exact contained-region Boolean case: all design vertices lie strictly in
 // beta>0 and both artificial cylinder caps lie outside the design cube.
 // The exact intersection is therefore the unchanged design boundary.
 spin::TriangleMesh containedDesign;
 containedDesign.vertices.resize(8,3);containedDesign.vertices<<-1.5,-1.5,-1.5,1.5,-1.5,-1.5,
     1.5,1.5,-1.5,-1.5,1.5,-1.5,-1.5,-1.5,1.5,1.5,-1.5,1.5,1.5,1.5,1.5,-1.5,1.5,1.5;
 containedDesign.faces.resize(12,3);containedDesign.faces<<0,2,1,0,3,2,4,5,6,4,6,7,
     0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
 spin::QuadricCoefficients10::Vector broadValues;
 broadValues<<-1,0,0,0,.01,.02,0,0,0,0;
 const auto broad=spin::QuadricCoefficients10::fromUnnormalized(broadValues);
 auto broadSettings=settings;
 broadSettings.truncationRadius=6.0;
 Eigen::VectorXi containedBirths(containedDesign.faces.rows());
 for(int f=0;f<containedBirths.size();++f)containedBirths(f)=f;
 const auto broadTri=spin::triangulateQuadric(broad,containedDesign.vertices,broadSettings);
 const auto contained=spin::differentiateDiscreteQuadricStratumMoments(
     containedDesign,broadTri,broad,containedDesign,containedBirths,1.0,broadSettings);
 require(contained.defined && contained.jacobian.rows()==10 && contained.jacobian.cols()==7,
         "strict cylinder differential rejected an exact contained Boolean result: "+contained.diagnostic);
 require(contained.constructedIntersectionVertexCount==0 && contained.jacobian.norm()<1e-12,
         "unchanged contained Boolean material acquired fictitious quadric sensitivity");

 // Nontrivial, clipped material fixture: the positive side of this elliptic
 // cylinder intersects the four corners of the design cube.  Unlike the
 // source-only checks above, this exercises the actual CGAL Boolean output,
 // its birth-face provenance, and the reduced 10x7 material-moment map.
 auto clippedRequest=spin::SpinDesignRequest{};
 clippedRequest.outerSurface=cube(-2.0,2.0);
 clippedRequest.designDomain=cube(-1.0,1.0);
 clippedRequest.materialDensity=1.0;
 const auto& clippedDesign=clippedRequest.designDomain;
 auto clippedSettings=spin::QuadricMomentEvaluatorSettings{};
 clippedSettings.triangulation=settings;
 clippedSettings.triangulation.truncationRadius=8.0;
 clippedSettings.derivativeMethod=spin::QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
 const auto clippedBase=spin::evaluateQuadricMaterialMoments(
     clippedRequest,q,clippedSettings);
 const auto clippedDifferential=spin::differentiateDiscreteQuadricStratumMoments(
     clippedDesign,clippedBase.quadricTriangulation,q,
     clippedBase.materialBoundaryTarget,clippedBase.materialBirthFaces,
     clippedRequest.materialDensity,clippedSettings.triangulation);
 require(clippedDifferential.defined,
         "clipped cylinder discrete Boolean differential unavailable: "+clippedDifferential.diagnostic);
 require(clippedDifferential.constructedIntersectionVertexCount>0,
         "clipped cylinder fixture did not create Boolean intersection vertices");
 double previousMomentError=1.0;
 for(double h:{1e-4,3e-5,1e-5}) {
  Eigen::Matrix<double,10,7> finiteDifference=Eigen::Matrix<double,10,7>::Zero();
  for(int j=0;j<7;++j) {
   const auto direction=clippedDifferential.coefficientBasis.col(j).eval();
   const auto plus=spin::evaluateQuadricMaterialMoments(
       clippedRequest,retractCylinder(q,direction,h),clippedSettings);
   const auto minus=spin::evaluateQuadricMaterialMoments(
       clippedRequest,retractCylinder(q,direction,-h),clippedSettings);
   require(plus.materialBoundaryTarget.faces.rows()==clippedBase.materialBoundaryTarget.faces.rows() &&
           minus.materialBoundaryTarget.faces.rows()==clippedBase.materialBoundaryTarget.faces.rows(),
           "clipped-cylinder Boolean topology changed during reduced finite difference");
   finiteDifference.col(j)=(plus.variableMomentsTarget.values()-
       minus.variableMomentsTarget.values())/(2*h);
  }
  const double relative=(finiteDifference-clippedDifferential.jacobian).norm()/
      std::max(1e-30,finiteDifference.norm());
  std::cout<<"clipped cylinder Boolean moment h="<<h<<" relative="<<relative
           <<" constructed="<<clippedDifferential.constructedIntersectionVertexCount<<"\n";
  if(h==1e-5) require(relative<1e-6,
      "clipped cylinder Boolean-moment Jacobian failed independent finite differences");
  if(h<1e-4 && previousMomentError>1e-8)
      require(relative<.3*previousMomentError,
          "clipped cylinder Boolean-moment finite-difference error did not refine");
  previousMomentError=relative;
 }
 const auto again=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
 require(again.defined && again.scientificHash==motion.scientificHash,
         "cylinder stratum source motion is nondeterministic");
 spin::QuadricCoefficients10::Vector circularValues;
 circularValues<<-1,0,0,0,1,1,0,0,0,0;
 const auto circular=spin::evaluateQuadricStratumSourceMotion(
     spin::QuadricCoefficients10::fromUnnormalized(circularValues),domain,settings);
 require(!circular.defined && circular.vertexVelocity.empty() &&
         circular.diagnostic.find("simple, separated eigenspectrum")!=std::string::npos,
         "repeated-eigenvalue cylinder was granted a chart-dependent derivative");
 spin::QuadricCoefficients10::Vector paraboloidValues;
 paraboloidValues<<0,0,0,-1,1,2,0,0,0,0;
 const auto paraboloid=spin::evaluateQuadricStratumSourceMotion(
     spin::QuadricCoefficients10::fromUnnormalized(paraboloidValues),domain,settings);
 require(paraboloid.defined && paraboloid.coefficientBasis.cols()==8 &&
         paraboloid.diagnostic.find("reduced paraboloid stratum source motion")==0,
         "paraboloid source motion did not use its independent reduced chart");

 std::cout<<"[PASS] elliptic-cylinder source motion in reduced stratum coordinates\n";
}

void checkHyperbolicCylinderTriangulation() {
 spin::QuadricCoefficients10::Vector values;
 values<<-1,0,0,0,1,-1,0,0,0,0; // x^2-y^2-1=0
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
 Eigen::MatrixXd domain(8,3);
 domain<<-1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
         -1,-1,1, 1,-1,1, 1,1,1, -1,1,1;
 spin::QuadricTriangulationSettings settings;
 settings.cylinderAngularResolution=12;
 settings.hyperbolicCylinderGridResolution=12;
 settings.truncationRadius=3.0;
 const auto triangulation=spin::triangulateQuadric(q,domain,settings);
 require(triangulation.type==spin::QuadricSurfaceType::HyperbolicCylinder,
         "hyperbolic-cylinder triangulation returned the wrong type");
 require(triangulation.interfaceVertexCount>0 &&
         triangulation.interfaceFaceCount>0 && triangulation.capFaceCount>0 &&
         triangulation.closedBoundary.faces.rows()==triangulation.faceIsCap.size(),
         "hyperbolic-cylinder triangulation omitted interface or cap records");
 const auto audit=manufacturing::auditTriangleBoundary(
     triangulation.closedBoundary.vertices, triangulation.closedBoundary.faces,
     manufacturing::BoundaryOrientation::AnyNonzero);
 require(audit.strictBoundaryValid(),
         "hyperbolic-cylinder finite boundary failed strict topology audit");
 std::cout<<"hyperbolic cylinder triangulation vertices="
          <<triangulation.closedBoundary.vertices.rows()<<" faces="
          <<triangulation.closedBoundary.faces.rows()<<" loops="
          <<triangulation.boundaryLoopCount<<"\n";
 const auto motion=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
 require(motion.defined,"hyperbolic-cylinder source motion unavailable: "+motion.diagnostic);
 require(motion.coefficientBasis.rows()==10 && motion.coefficientBasis.cols()==7,
         "hyperbolic-cylinder motion did not expose its reduced 7D chart");
 require(motion.triangulation.closedBoundary.faces==triangulation.closedBoundary.faces,
         "hyperbolic-cylinder motion changed native source topology");
 double previous=1.0;
 for(double h:{1e-4,3e-5,1e-5}) {
  double error2=0.0,norm2=0.0;
  for(int column=0;column<motion.coefficientBasis.cols();++column) {
   const auto direction=motion.coefficientBasis.col(column).eval();
   const auto plus=spin::triangulateQuadric(retractHyperbolicCylinder(q,direction,h),domain,settings);
   const auto minus=spin::triangulateQuadric(retractHyperbolicCylinder(q,direction,-h),domain,settings);
   if(plus.closedBoundary.faces!=motion.triangulation.closedBoundary.faces ||
      minus.closedBoundary.faces!=motion.triangulation.closedBoundary.faces) {
    std::cerr<<"hyperbolic topology mismatch column="<<column
             <<" plus="<<(plus.closedBoundary.faces==motion.triangulation.closedBoundary.faces)
             <<" minus="<<(minus.closedBoundary.faces==motion.triangulation.closedBoundary.faces)
             <<" plus-type="<<spin::quadricSurfaceTypeToken(plus.type)
             <<" minus-type="<<spin::quadricSurfaceTypeToken(minus.type)<<"\n";
    for(int face=0;face<plus.closedBoundary.faces.rows();++face) {
     if(plus.closedBoundary.faces.row(face)!=motion.triangulation.closedBoundary.faces.row(face)) {
      std::cerr<<"first plus winding mismatch face="<<face<<" plus="
               <<plus.closedBoundary.faces.row(face)<<" base="
               <<motion.triangulation.closedBoundary.faces.row(face)<<"\n";
      break;
     }
    }
    require(false,"hyperbolic-cylinder finite-difference topology changed");
   }
   const Eigen::MatrixXd finiteDifference=(plus.closedBoundary.vertices-
       minus.closedBoundary.vertices)/(2.0*h);
   for(int vertex=0;vertex<finiteDifference.rows();++vertex) {
    const Eigen::Vector3d predicted=motion.vertexVelocity[static_cast<std::size_t>(vertex)].col(column);
    error2+=(finiteDifference.row(vertex).transpose()-predicted).squaredNorm();
    norm2+=predicted.squaredNorm();
   }
  }
  const double relative=std::sqrt(error2/std::max(1e-30,norm2));
  std::cout<<"hyperbolic-cylinder source h="<<h<<" relative="<<relative<<"\n";
  if(h==1e-5)require(relative<1e-6,"hyperbolic-cylinder source derivative failed finite difference");
  if(h<5e-5)require(relative<.2*previous,"hyperbolic-cylinder derivative refinement failed");
  previous=relative;
 }
 spin::SpinDesignRequest clippedRequest;
 clippedRequest.outerSurface=cube(-3.0,3.0);
 clippedRequest.designDomain=cube(-2.0,2.0);
 const Eigen::Vector3d designOffset(.113,-.071,.057);
 clippedRequest.designDomain.vertices.rowwise()+=designOffset.transpose();
 clippedRequest.materialDensity=1.0;
 spin::QuadricMomentEvaluatorSettings evaluatorSettings;
 evaluatorSettings.triangulation=settings;
 evaluatorSettings.triangulation.truncationRadius=8.0;
 evaluatorSettings.derivativeMethod=spin::QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
 const auto clippedBase=spin::evaluateQuadricMaterialMoments(
     clippedRequest,q,evaluatorSettings);
 const auto clippedDifferential=spin::differentiateDiscreteQuadricStratumMoments(
     clippedRequest.designDomain,clippedBase.quadricTriangulation,q,
     clippedBase.materialBoundaryTarget,clippedBase.materialBirthFaces,
     clippedRequest.materialDensity,evaluatorSettings.triangulation);
 require(clippedDifferential.defined,
         "hyperbolic-cylinder Boolean moment differential unavailable: "+
         clippedDifferential.diagnostic);
 require(clippedDifferential.coefficientBasis.cols()==7 &&
         clippedDifferential.constructedIntersectionVertexCount>0,
         "hyperbolic-cylinder Boolean fixture has no constructed intersection vertices");
 double previousBooleanError=1.0;
 for(double h:{1e-4,3e-5,1e-5}) {
  Eigen::Matrix<double,10,7> finiteDifference=Eigen::Matrix<double,10,7>::Zero();
  for(int column=0;column<7;++column) {
   const auto direction=clippedDifferential.coefficientBasis.col(column).eval();
   const auto plus=spin::evaluateQuadricMaterialMoments(
       clippedRequest,retractHyperbolicCylinder(q,direction,h),evaluatorSettings);
   const auto minus=spin::evaluateQuadricMaterialMoments(
       clippedRequest,retractHyperbolicCylinder(q,direction,-h),evaluatorSettings);
   require(plus.materialBoundaryTarget.faces.rows()==clippedBase.materialBoundaryTarget.faces.rows() &&
           minus.materialBoundaryTarget.faces.rows()==clippedBase.materialBoundaryTarget.faces.rows(),
           "hyperbolic-cylinder Boolean topology changed during reduced finite difference");
   finiteDifference.col(column)=(plus.variableMomentsTarget.values()-
       minus.variableMomentsTarget.values())/(2.0*h);
  }
  const double relative=(finiteDifference-clippedDifferential.jacobian).norm()/
      std::max(1e-30,finiteDifference.norm());
  std::cout<<"clipped hyperbolic-cylinder Boolean moment h="<<h
           <<" relative="<<relative
           <<" constructed="<<clippedDifferential.constructedIntersectionVertexCount<<"\n";
  if(h==1e-5)require(relative<1e-6,
      "hyperbolic-cylinder Boolean-moment Jacobian failed independent finite differences");
  if(h<1e-4 && previousBooleanError>1e-8)
      require(relative<.3*previousBooleanError,
          "hyperbolic-cylinder Boolean-moment finite-difference error did not refine");
  previousBooleanError=relative;
 }
 std::cout<<"[PASS] hyperbolic-cylinder source motion and clipped Boolean moment Jacobian\n";
}

void checkParabolicCylinderPolynomialContract() {
 // Independent of the grid formula: all interface vertices must satisfy the
 // public p=a+b.x+x^T A x convention, not merely match its own derivative.
 const Eigen::Matrix3d rotation=Eigen::AngleAxisd(.31,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
 for(int example=0;example<3;++example) {
  const Eigen::Vector3d u=example==0?Eigen::Vector3d::UnitX().eval():rotation.col(0).eval();
  const Eigen::Vector3d v=example==0?Eigen::Vector3d::UnitY().eval():rotation.col(1).eval();
  const double a=example==0?0.:.37,bu=example==0?0.:.41;
  const double beta=example==0?1.:1.3,lambda=example==0?1.:1.7;
  const Eigen::Matrix3d A=lambda*u*u.transpose();
  const Eigen::Vector3d b=bu*u+beta*v;
  spin::QuadricCoefficients10::Vector values;
  values<<a,b.x(),b.y(),b.z(),A(0,0),A(1,1),A(2,2),A(1,2),A(0,2),A(0,1);
  if(example==2)values=-values;
  const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
  Eigen::MatrixXd domain(2,3);domain<<-1,-1,-1,1,1,1;
  spin::QuadricTriangulationSettings settings;
  settings.truncationRadius=2.;settings.parabolicCylinderGridResolution=6;
  const auto mesh=spin::triangulateQuadric(q,domain,settings);
  double maximumResidual=0.;
  for(int i=0;i<mesh.interfaceVertexCount;++i) {
   const Eigen::Vector3d x=mesh.closedBoundary.vertices.row(i);
   // Literal polynomial in the original rotated coordinates, including sign.
   const double polynomial=a+bu*u.dot(x)+beta*v.dot(x)+lambda*std::pow(u.dot(x),2);
   maximumResidual=std::max(maximumResidual,std::abs(polynomial)/(1.+x.squaredNorm()));
  }
  std::cout<<"parabolic-cylinder polynomial example="<<example<<" residual="<<maximumResidual<<"\n";
  require(maximumResidual<1e-12,"parabolic-cylinder interface does not satisfy requested polynomial");
  const auto motion=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
  require(motion.defined && motion.coefficientBasis.cols()==6,
          "polynomial-contract source differential unavailable: "+motion.diagnostic);
  double derivativeResidual=0.;
  for(int i=0;i<mesh.interfaceVertexCount;++i) {
   const Eigen::Vector3d x=mesh.closedBoundary.vertices.row(i);
   for(int j=0;j<6;++j) {
    const auto dq=motion.coefficientBasis.col(j).eval();
    const double dp=dq(0)+dq(1)*x.x()+dq(2)*x.y()+dq(3)*x.z()+
        dq(4)*x.x()*x.x()+dq(5)*x.y()*x.y()+dq(6)*x.z()*x.z()+
        2.*(dq(7)*x.y()*x.z()+dq(8)*x.x()*x.z()+dq(9)*x.x()*x.y());
    const double transported=q.gradient(x).dot(motion.vertexVelocity[static_cast<std::size_t>(i)].col(j));
    derivativeResidual=std::max(derivativeResidual,std::abs(dp+transported)/(1.+std::abs(dp)+std::abs(transported)));
   }
  }
  std::cout<<"parabolic-cylinder implicit derivative example="<<example<<" residual="<<derivativeResidual<<"\n";
  require(derivativeResidual<1e-10,"parabolic-cylinder derivative violates implicit surface identity");
 }
}

void checkParabolicCylinderClassification() {
 spin::QuadricCoefficients10::Vector values;
 values<<0,0,1,0,1,0,0,0,0,0; // x^2 + y = 0, extruded along z
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
 const auto subtype=spin::classifyQuadricSurfaceSubtype(q);
 require(std::string(spin::quadricSurfaceSubtypeToken(subtype))=="PARABOLIC_CYLINDER",
         "rank-one real parabola extrusion was not classified as a parabolic cylinder");
 checkStratumBasis(values,spin::QuadricSurfaceSubtype::ParabolicCylinder,6);
 spin::QuadricCoefficients10::Vector planePairValues;
 planePairValues<<-1,0,0,0,1,0,0,0,0,0;
 bool planePairRejected=false;
 try {
  (void)spin::classifyQuadricSurfaceSubtype(
      spin::QuadricCoefficients10::fromUnnormalized(planePairValues));
 } catch(const std::invalid_argument&) { planePairRejected=true; }
 require(planePairRejected,
         "rank-one parallel-plane degeneration was misclassified as a parabolic cylinder");
 Eigen::MatrixXd domain(8,3);
 domain<<-2,-2,-2,2,-2,-2,2,2,-2,-2,2,-2,
         -2,-2,2,2,-2,2,2,2,2,-2,2,2;
 spin::QuadricTriangulationSettings settings;
 settings.parabolicCylinderGridResolution=12;
 settings.truncationRadius=2.0;
 const auto triangulation=spin::triangulateQuadric(q,domain,settings);
 require(triangulation.type==spin::QuadricSurfaceType::ParabolicCylinder,
         "parabolic-cylinder triangulation returned the wrong type");
 const auto audit=manufacturing::auditTriangleBoundary(
     triangulation.closedBoundary.vertices,triangulation.closedBoundary.faces,
     manufacturing::BoundaryOrientation::AnyNonzero);
 require(audit.strictBoundaryValid() && audit.signedVolume>0.0,
         "parabolic-cylinder finite construction failed closed-boundary audit: closed="+
         std::to_string(audit.closed)+" edge="+std::to_string(audit.edgeManifold)+
         " vertex="+std::to_string(audit.vertexManifold)+" orientation="+
         std::to_string(audit.consistentlyOriented)+" components="+
         std::to_string(audit.componentOrientationsValid)+" degenerate="+
         std::to_string(audit.degenerateFaceCount)+" unreferenced="+
         std::to_string(audit.unreferencedVertexCount)+" intersections="+
         std::to_string(audit.selfIntersections.completed)+"/"+
         std::to_string(audit.selfIntersections.intersectionFree));
 const auto motion=spin::evaluateQuadricStratumSourceMotion(q,domain,settings);
 require(motion.defined && motion.coefficientBasis.cols()==6,
         "parabolic-cylinder reduced source motion unavailable: "+motion.diagnostic);
 double previous=1.0;
 for(double h:{1e-4,3e-5,1e-5}) {
  double error2=0.0,norm2=0.0;
  for(int column=0;column<6;++column) {
   const auto direction=motion.coefficientBasis.col(column).eval();
   const auto plus=spin::triangulateQuadric(
       retractParabolicCylinder(q,direction,h),domain,settings);
   const auto minus=spin::triangulateQuadric(
       retractParabolicCylinder(q,direction,-h),domain,settings);
   require(plus.closedBoundary.faces==motion.triangulation.closedBoundary.faces &&
           minus.closedBoundary.faces==motion.triangulation.closedBoundary.faces,
           "parabolic-cylinder triangulation topology changed under fixed-stratum perturbation");
   const Eigen::MatrixXd finiteDifference=(plus.closedBoundary.vertices-
       minus.closedBoundary.vertices)/(2.0*h);
   for(int vertex=0;vertex<finiteDifference.rows();++vertex) {
    const Eigen::Vector3d predicted=
        motion.vertexVelocity[static_cast<std::size_t>(vertex)].col(column);
    error2+=(finiteDifference.row(vertex).transpose()-predicted).squaredNorm();
    norm2+=predicted.squaredNorm();
   }
  }
  const double relative=std::sqrt(error2/std::max(1e-30,norm2));
  std::cout<<"parabolic-cylinder source h="<<h<<" relative="<<relative<<"\n";
  if(h==1e-5)require(relative<1e-6,
      "parabolic-cylinder source derivative failed independent finite differences");
  if(h<5e-5)require(relative<.2*previous,
      "parabolic-cylinder source derivative refinement failed");
  previous=relative;
 }
 spin::SpinDesignRequest clippedRequest;
 clippedRequest.outerSurface=cube(-3.0,3.0);
 clippedRequest.designDomain=cube(-2.0,2.0);
 const Eigen::Vector3d designOffset(.113,-.071,.057);
 clippedRequest.designDomain.vertices.rowwise()+=designOffset.transpose();
 clippedRequest.materialDensity=1.0;
 spin::QuadricMomentEvaluatorSettings evaluatorSettings;
 evaluatorSettings.triangulation=settings;
 evaluatorSettings.triangulation.truncationRadius=8.0;
 evaluatorSettings.derivativeMethod=spin::QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
 const auto clippedBase=spin::evaluateQuadricMaterialMoments(
     clippedRequest,q,evaluatorSettings);
 const auto clippedDifferential=spin::differentiateDiscreteQuadricStratumMoments(
     clippedRequest.designDomain,clippedBase.quadricTriangulation,q,
     clippedBase.materialBoundaryTarget,clippedBase.materialBirthFaces,
     clippedRequest.materialDensity,evaluatorSettings.triangulation);
 require(clippedDifferential.defined &&
         clippedDifferential.coefficientBasis.cols()==6 &&
         clippedDifferential.constructedIntersectionVertexCount>0,
         "parabolic-cylinder clipped Boolean differential unavailable: "+
         clippedDifferential.diagnostic);
 auto reducedSettings=evaluatorSettings;
 reducedSettings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteReducedStratum;
 const auto unified=spin::evaluateQuadricMaterialMoments(
     clippedRequest,q,reducedSettings);
 require(!unified.derivativeDefined && unified.stratumDifferential.defined &&
         unified.stratumDifferential.jacobian.isApprox(clippedDifferential.jacobian,0.0) &&
         unified.stratumDifferential.coefficientBasis.isApprox(
             clippedDifferential.coefficientBasis,0.0),
         "parabolic-cylinder reduced differential was not carried by the unified evaluator");
 previous=1.0;
 for(double h:{1e-4,3e-5,1e-5}) {
  Eigen::Matrix<double,10,6> finiteDifference=Eigen::Matrix<double,10,6>::Zero();
  for(int column=0;column<6;++column) {
   const auto direction=clippedDifferential.coefficientBasis.col(column).eval();
   const auto plus=spin::evaluateQuadricMaterialMoments(
       clippedRequest,retractParabolicCylinder(q,direction,h),evaluatorSettings);
   const auto minus=spin::evaluateQuadricMaterialMoments(
       clippedRequest,retractParabolicCylinder(q,direction,-h),evaluatorSettings);
   require(plus.materialBoundaryTarget.faces.rows()==clippedBase.materialBoundaryTarget.faces.rows() &&
           minus.materialBoundaryTarget.faces.rows()==clippedBase.materialBoundaryTarget.faces.rows(),
           "parabolic-cylinder Boolean topology changed under fixed-stratum perturbation");
   finiteDifference.col(column)=(plus.variableMomentsTarget.values()-
       minus.variableMomentsTarget.values())/(2.0*h);
  }
  const double relative=(finiteDifference-clippedDifferential.jacobian).norm()/
      std::max(1e-30,finiteDifference.norm());
  std::cout<<"clipped parabolic-cylinder Boolean moment h="<<h
           <<" relative="<<relative
           <<" constructed="<<clippedDifferential.constructedIntersectionVertexCount<<"\n";
  if(h==1e-5)require(relative<1e-6,
      "parabolic-cylinder Boolean-moment Jacobian failed independent finite differences");
  if(h<1e-4)require(relative<.3*previous,
      "parabolic-cylinder Boolean-moment finite-difference error did not refine");
  previous=relative;
 }
 std::cout<<"[PASS] parabolic-cylinder source motion and clipped Boolean moment Jacobian\n";
}
}
int main(int argc,char** argv){try{std::cout<<std::setprecision(17);
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
 if(argc==2 && std::string(argv[1])=="--reuse-rational-source") {checkReuseRationalSource();return 0;}
#endif
 checkParabolicCylinderPolynomialContract();
 if(argc==2 && std::string(argv[1])=="--parabolic-cylinder") {
  checkParabolicCylinderClassification();return 0;
 }
 checkLargeRadiusTruncationAccuracy();
 checkRepresentedCoordinateSideWitness();
    checkCapturedFullStepTruncationSourceMotion();
    checkNearConeSourceMotionUsesCertifiedExtendedSolve();
    checkCherchiNewtonNearConeRetainsOriginalNormalIdentity();
 using Q=CGAL::Gmpq;
 const spin::detail::ExactRationalPolynomial tarskiPolynomial={Q(3)/Q(16),Q(-1),Q(1)};
 const std::vector<std::pair<spin::detail::ExactRationalPolynomial,int>> tarskiCases={
   {{Q(0),Q(0),Q(1)},2},                         // x^2: positive at both roots
   {{Q(-1)/Q(4),Q(1)},1},                       // x-1/4: zero at one root, positive at the other
   {{Q(3)/Q(16),Q(-1),Q(1)},0},                 // P itself: zero at both roots
   {{Q(0),Q(0),Q(0),Q(0),Q(0),Q(1)},2},         // x^5: high-degree weighted numerator
   {{Q(0),Q(0),Q(0),Q(0),Q(0),Q(-1)},-2}        // -x^5: negative at both roots
 };
 for(std::size_t index=0;index<tarskiCases.size();++index) {
   const auto& [value,expected]=tarskiCases[index];
   const int actual=spin::detail::sturmTarskiSignSum(
     tarskiPolynomial,value,Q(0),Q(1));
   require(actual==expected,
     "exact Sturm-Tarski sign sum case "+std::to_string(index)+
     " expected "+std::to_string(expected)+" but got "+std::to_string(actual));
 }
 spin::QuadricCoefficients10::Vector p;p<<-.31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03;check(p);
 p(0)=.31;check(p);p(0)=-.31;p(6)=.72;check(p,false);unavailable();
 p<<-1.6324804504042277e-09,2.0881696625961295e-06,
    3.1141661275861445e-06,9.5198658917346903e-06,
    .40440677544965664,.40441390012607642,-.82030759916088636,
    2.6262036827739407e-06,-1.5955878322912755e-06,1.8830914019005905e-08;
 Eigen::MatrixXd capturedDomain(2,3);capturedDomain<<-1,-1,-1,1,1,1;
 spin::QuadricTriangulationSettings capturedSettings;
 capturedSettings.icosphereSubdivisions=2;capturedSettings.truncationRadius=8;
 const auto captured=spin::evaluateQuadricSourceMotion(
   spin::QuadricCoefficients10::fromUnnormalized(p),capturedDomain,capturedSettings);
 require(captured.defined,"captured unknown-seed source motion failed: "+captured.diagnostic);
 int capturedCapVertices=0;
 for(int vertex=0;vertex<captured.triangulation.interfaceVertexCount;++vertex) {
  const double radiusResidual=std::abs(
    captured.triangulation.closedBoundary.vertices.row(vertex).norm()-
    captured.triangulation.truncationRadius);
  if(radiusResidual<1e-4*captured.triangulation.truncationRadius) {
   ++capturedCapVertices;
   require(radiusResidual<=1e-10*captured.triangulation.truncationRadius,
           "captured cap crossing violates strict radius contract");
  }
 }
 require(capturedCapVertices>0,"captured case did not exercise cap crossings");
 spin::QuadricCoefficients10::Vector cherchiNewtonLevel5;
 cherchiNewtonLevel5 <<
   -.19223703824980956,.029893982546063312,-.040178002262895494,
   .0008982683924855251,.5777227428369576,.7914604728639738,
   -.0011177055319905422,.0005885695113244015,-.0018721422682393507,
   .018916762064737025;
 spin::QuadricTriangulationSettings level5Settings;
 level5Settings.icosphereSubdivisions=5;
 level5Settings.maximumRadiusDoublings=4;
 const auto level5Motion=spin::evaluateQuadricSourceMotion(
   spin::QuadricCoefficients10::fromUnnormalized(cherchiNewtonLevel5),
   capturedDomain,level5Settings);
 require(level5Motion.defined,
         "level-5 full-rank source-motion path must support the captured resolvable hyperboloid: "+
         level5Motion.diagnostic);
 require(level5Motion.triangulation.boundaryLoopCount==2 &&
         level5Motion.velocity.size()==
           static_cast<std::size_t>(level5Motion.triangulation.closedBoundary.vertices.rows()) &&
         level5Motion.provenance.size()==level5Motion.velocity.size(),
         "level-5 source-motion must bind a complete two-loop mesh and vertex Jacobians");
 std::cout<<"[PASS] level-5 captured hyperboloid source motion\n";
 auto adaptiveCapSettings=level5Settings;
 adaptiveCapSettings.icosphereSubdivisions=3;
 adaptiveCapSettings.maximumIcosphereSubdivisions=5;
 adaptiveCapSettings.maximumRadiusDoublings=16;
 const auto adaptiveCapMotion=spin::evaluateQuadricSourceMotion(
   spin::QuadricCoefficients10::fromUnnormalized(cherchiNewtonLevel5),
   capturedDomain,adaptiveCapSettings);
 require(adaptiveCapMotion.defined,
         "bounded normal-sphere refinement must rebuild source motion on the accepted graph: "+
         adaptiveCapMotion.diagnostic);
 require(adaptiveCapMotion.triangulation.initialIcosphereSubdivisions==3 &&
         adaptiveCapMotion.triangulation.effectiveIcosphereSubdivisions==5 &&
         adaptiveCapMotion.triangulation.boundaryLoopCount==2,
         "adaptive cap topology must report requested and effective normal resolutions");
 require(adaptiveCapMotion.velocity.size()==level5Motion.velocity.size(),
         "adaptive and explicit level-5 source motion must have the same vertex domain");
 for(std::size_t vertex=0;vertex<level5Motion.velocity.size();++vertex)
  require((adaptiveCapMotion.velocity[vertex]-level5Motion.velocity[vertex]).norm()==0.0,
          "adaptive source-motion Jacobian was not rebuilt on the accepted level-5 graph");
 std::cout<<"[PASS] bounded adaptive cap source motion level 3->5\n";
 spin::QuadricCoefficients10::Vector cherchiNewtonLevel7;
 cherchiNewtonLevel7 <<
   -.1998242186450357,.02382394696275607,.005243276046278874,
   .0015024879022398209,.5734432875092278,.7938771401867838,
   -.0000336431262224418,.0007933293534953297,-.0019952147391396123,
   .01975405514590195;
 spin::QuadricTriangulationSettings level7Settings;
 level7Settings.icosphereSubdivisions=7;
 level7Settings.maximumRadiusDoublings=4;
 const auto level7Motion=spin::evaluateQuadricSourceMotion(
   spin::QuadricCoefficients10::fromUnnormalized(cherchiNewtonLevel7),
   capturedDomain,level7Settings);
 require(level7Motion.defined,
         "level-7 source motion must cover the captured Newton type-boundary trial: "+
         level7Motion.diagnostic);
 require(level7Motion.triangulation.boundaryLoopCount==2 &&
         level7Motion.velocity.size()==
           static_cast<std::size_t>(level7Motion.triangulation.closedBoundary.vertices.rows()) &&
         level7Motion.provenance.size()==level7Motion.velocity.size(),
         "level-7 source motion must bind a complete two-loop mesh and vertex Jacobians");
 std::cout<<"[PASS] level-7 captured Newton-trial source motion\n";
 p<<9.7964069485708039e-10,-3.9936663645232352e-07,
    1.9384480260996023e-07,-6.4931335978636097e-07,
    .41207792841821683,.41206631482390654,-.8126457611390473,
    2.4466098949066089e-07,2.3708460028151824e-07,6.3081416358941302e-09;
 const auto nearCone=spin::evaluateQuadricSourceMotion(
   spin::QuadricCoefficients10::fromUnnormalized(p),capturedDomain,capturedSettings);
 require(nearCone.defined,
         "a zero-Sturm-root subinterval at depth 65 must be pruned, not reported as root-isolation exhaustion: "+
         nearCone.diagnostic);
 p<<3.5482915291333073e-10,-1.4465158869595825e-07,
   7.0210965917303937e-08,-2.3518551731615803e-07,
   .41207792683184907,.41206629816542273,-.81264577039083508,
   8.8619038160607571e-08,8.5873899236902657e-08,2.2848286155595496e-09;
 const auto algebraicBranch=spin::evaluateQuadricSourceMotion(
   spin::QuadricCoefficients10::fromUnnormalized(p),capturedDomain,capturedSettings);
 require(algebraicBranch.diagnostic.find("root isolation depth budget exhausted")==std::string::npos,
         "the sign of the branch product at a unique Sturm root must be decided algebraically, not by exhausting subdivision depth: "+
         algebraicBranch.diagnostic);
 require(algebraicBranch.defined,
         "the algebraically admissible unique branch root should complete source motion: "+
         algebraicBranch.diagnostic);
 p<<4.2090703565348103e-10,-1.7158928022591721e-07,
   8.3285976014279177e-08,-2.7898267080833971e-07,
   .41207792698626994,.41206629988422155,-.8126457694409599,
   1.0512192503473947e-07,1.0186565018548999e-07,2.7103196209354624e-09;
 const auto clusteredBranches=spin::evaluateQuadricSourceMotion(
   spin::QuadricCoefficients10::fromUnnormalized(p),capturedDomain,capturedSettings);
 require(clusteredBranches.defined,
         "two clustered Sturm roots with one admissible branch must be counted exactly without exceeding the isolation depth budget: "+
         clusteredBranches.diagnostic);
 spin::QuadricCoefficients10::Vector twoPositiveBranches;
 twoPositiveBranches<<2,-2,0,0,1,1,-1,0,0,0;
 const auto twoPositiveQuadric=spin::QuadricCoefficients10::fromUnnormalized(twoPositiveBranches);
 std::string twoBranchDiagnostic;
 require(spin::certifySingleTruncationRoot(
           twoPositiveQuadric,
           Eigen::Vector3d::UnitX(),Eigen::Vector3d::UnitZ(),1.7,
           &twoBranchDiagnostic),
         "the rational two-root fixture must accept its single admissible branch: "+
         twoBranchDiagnostic);
 spin::TruncationRootCertifier sharedRootCertifier(twoPositiveQuadric,1.7);
 std::string sharedBranchDiagnostic;
 require(sharedRootCertifier.certify(Eigen::Vector3d::UnitX(),
           Eigen::Vector3d::UnitZ(),&sharedBranchDiagnostic) &&
         sharedBranchDiagnostic==twoBranchDiagnostic,
         "reused exact root-certifier context must preserve the single-edge proof and diagnostic");
 spin::QuadricCoefficients10::Vector cylinder;
 cylinder<<-1,0,0,0,1,2,0,0,0,0;
 checkRankDeficientRejected(cylinder,"ELLIPTIC_CYLINDER");
 checkStratumBasis(cylinder,spin::QuadricSurfaceSubtype::EllipticCylinder,7);
 checkCylinderStratumSourceMotion();
 checkHyperbolicCylinderTriangulation();
 checkParabolicCylinderClassification();
 spin::QuadricCoefficients10::Vector plane;
 plane<<.17,.31,-.47,.826,0,0,0,0,0,0;
 checkStratumBasis(plane,spin::QuadricSurfaceSubtype::Plane,3);
 checkPlaneStratumSourceMotion();
 checkPlaneAxialChartSourceMotion();
 checkPlaneChartSeamRejected();
 checkPlaneStratumBooleanDifferential();
 spin::QuadricCoefficients10::Vector paraboloid;
 paraboloid<<0,0,0,-1,1,2,0,0,0,0;
 checkRankDeficientRejected(paraboloid,"PARABOLOID");
 checkStratumBasis(paraboloid,spin::QuadricSurfaceSubtype::Paraboloid,8);
 checkParaboloidStratumSourceMotion();
 checkParaboloidEigenvectorPivotTieRejected();
 return 0;
}catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<"\n";return 1;}}
