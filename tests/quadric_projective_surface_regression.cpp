#include "spin/QuadricSourceMotion.h"
#include "spin/QuadricInverseGaussMap.h"
#include "verification/ExactBoundaryAudit.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool value,const std::string& message) {
 if(!value)throw std::runtime_error(message);
}
Eigen::MatrixXd domain() {
 Eigen::MatrixXd v(8,3);
 for(int i=0;i<8;++i)for(int j=0;j<3;++j)v(i,j)=(i&(1<<j))?.2:-.2;
 return v;
}
void scenario(bool twoSheets,bool algebraic=false) {
 using Q=spin::QuadricCoefficients10;
 Q::Vector p;p<<(twoSheets?1.0:-1.0),.1,-.12,.08,1,1.3,-.9,.03,-.04,.02;
 const auto q=Q::fromUnnormalized(p);
 spin::QuadricTriangulationSettings s;
 s.icosphereSubdivisions=2;s.truncationRadius=16;s.maximumRadiusDoublings=8;
 s.retainFixedNormalProvenance=true;
 const auto legacy=spin::triangulateQuadric(q,domain(),s);
 s.paperProjectiveHyperboloid=true;
 s.algebraicFiniteMap=algebraic;
 const auto t=spin::triangulateQuadric(q,domain(),s);
 require(t.paperProjectiveHyperboloid,"paper request still executed legacy analytic clipping");
 require(t.hash!=legacy.hash,"distinct discrete maps share an identity");
 require(t.boundaryLoopCount==2 && t.capFaceCount>0,"paper closure lacks two loops");
 require(t.effectiveIcosphereSubdivisions==2,"paper normals changed during construction");
 require(t.fixedNormals.rows()>0 && t.sourceProjectiveRayByVertex.size()==
         static_cast<std::size_t>(t.closedBoundary.vertices.rows()),"missing paper provenance");
 const auto audit=verification::auditExactTriangleBoundary(t.closedBoundary.vertices,
     t.closedBoundary.faces,verification::ExactBoundaryOrientation::OutwardPositive);
 require(audit.facts().boundaryValid,"independent paper source boundary audit failed");
 int rays=0;double nonAnalytic=0;
 spin::QuadricInverseGaussMap map(q);
 for(int i=0;i<t.interfaceVertexCount;++i) {
  const auto id=t.sourceProjectiveRayByVertex[i];
  if(id[0]<0)continue;
  ++rays;
  const auto finite=map.point(t.fixedNormals.row(id[0]));
  const auto ideal=map.idealEdgeDifferential(t.fixedNormals.row(id[1]),t.fixedNormals.row(id[2]));
  require(finite.valid() && ideal.defined && id[3]>=0 && id[3]<int(ideal.points.size()),
          "ray provenance cannot be replayed");
  const Eigen::Vector3d x=t.closedBoundary.vertices.row(i);
  const Eigen::Vector3d d=ideal.points[id[3]].direction,offset=x-*finite.point;
  require(offset.dot(d)>0 && offset.cross(d).norm()<=1e-11*offset.norm()*d.norm(),
          "surface truncation is not on its projective straight ray");
  require(std::abs(x.norm()-t.truncationRadius)<=1e-12*t.truncationRadius,"ray radius mismatch");
  nonAnalytic=std::max(nonAnalytic,std::abs(q.evaluate(x)));
 }
 require(rays>0 && nonAnalytic>1e-4,"analytic quadric clipping substituted for projective edges");
 for(int f=t.interfaceFaceCount;f<t.closedBoundary.faces.rows();++f) {
  bool separated=false;
  for(int a=0;a<3;++a) {
   bool above=true,below=true;
   for(int k=0;k<3;++k) {
    const double x=t.closedBoundary.vertices(t.closedBoundary.faces(f,k),a);
    above&=x>.2;below&=x<-.2;
   }
   separated|=above||below;
  }
  require(separated,"cap enters the design AABB");
 }
 s.truncationRadius=t.truncationRadius;s.maximumRadiusDoublings=0;
 const auto m=spin::evaluateQuadricSourceMotion(q,domain(),s);
 require(m.defined,"paper motion undefined: "+m.diagnostic);
 require(m.triangulation.hash==t.hash,"motion and geometry maps disagree");
 if(algebraic) {
  auto original=s;original.algebraicFiniteMap=false;
  require(spin::triangulateQuadric(q,domain(),original).hash!=t.hash,
          "algebraic policy was ignored or omitted from geometry identity");
  for(int i=0;i<t.interfaceVertexCount;++i) {
   const int normal=t.sourceNormalVertexByVertex[i];
   if(normal<0)continue;
   const auto expected=map.algebraicCoefficientDifferential(t.fixedNormals.row(normal));
   require(expected.valid(),"algebraic source point unavailable");
   require((t.closedBoundary.vertices.row(i).transpose().array()==expected.point.point->array()).all(),
           "source constructor did not use actual algebraic finite point");
   require((m.velocity[i].array()==expected.coefficientJacobian->array()).all(),
           "source derivative did not use the same algebraic point");
  }
 }
 double maxRelative=0,maxVertexRelative=0;
 for(double h:{3e-6,1e-6,3e-7})for(int j=0;j<10;++j) {
  auto pp=q.values(),pm=q.values();pp(j)+=h;pm(j)-=h;
  const auto plus=spin::triangulateQuadric(Q::fromUnnormalized(pp),domain(),s);
  const auto minus=spin::triangulateQuadric(Q::fromUnnormalized(pm),domain(),s);
  require(plus.sourceProjectiveRayByVertex==t.sourceProjectiveRayByVertex &&
          minus.sourceProjectiveRayByVertex==t.sourceProjectiveRayByVertex &&
          plus.closedBoundary.faces==t.closedBoundary.faces &&
          minus.closedBoundary.faces==t.closedBoundary.faces,"FD topology changed");
  const Eigen::MatrixXd fd=(plus.closedBoundary.vertices-minus.closedBoundary.vertices)/(2*h);
  Eigen::MatrixXd derivative(fd.rows(),3);
  for(int i=0;i<fd.rows();++i)derivative.row(i)=m.velocity[i].col(j);
  const double error=(fd-derivative).norm(),scale=derivative.norm();
  const double relative=error/std::max(1e-12,scale);
  std::cout<<"two_sheets="<<twoSheets<<" h="<<h<<" column="<<j<<" relative="<<relative<<'\n';
  require(scale<1e-8 ? error<=1e-10 : relative<=1e-6,"actual projective source derivative mismatch");
  for(int i=0;i<fd.rows();++i) {
   const double localError=(fd.row(i)-derivative.row(i)).norm();
   const double localScale=derivative.row(i).norm();
   require(localScale<1e-8 ? localError<=1e-10 : localError/localScale<=1e-6,
           "per-vertex projective source differential mismatch at vertex "+std::to_string(i));
   if(localScale>=1e-8)maxVertexRelative=std::max(maxVertexRelative,localError/localScale);
  }
  maxRelative=std::max(maxRelative,relative);
 }
 std::cout<<"source="<<t.hash<<" motion="<<m.scientificHash<<" vertices="
          <<t.closedBoundary.vertices.rows()<<" faces="<<t.closedBoundary.faces.rows()
          <<" rays="<<rays<<" max_relative="<<maxRelative<<" max_vertex_relative="<<maxVertexRelative<<'\n';
}
void rejectionContracts() {
 using Q=spin::QuadricCoefficients10;
 // Narrow normal cap about the face-interior direction (1,1,1), no icosphere
 // corner inside. Dropping all zero-crossing faces must not yield an empty pass.
 const Eigen::Vector3d axis=Eigen::Vector3d::Ones().normalized();
 const Eigen::Matrix3d a=.001*Eigen::Matrix3d::Identity()-1.001*axis*axis.transpose();
 Q::Vector p;p<<1,0,0,0,a(0,0),a(1,1),a(2,2),a(1,2),a(0,2),a(0,1);
 spin::QuadricTriangulationSettings s;s.paperProjectiveHyperboloid=true;
 s.icosphereSubdivisions=0;s.truncationRadius=16;s.maximumRadiusDoublings=0;
 bool rejected=false;
 try{(void)spin::triangulateQuadric(Q::fromUnnormalized(p),domain(),s);}
 catch(const std::runtime_error& error){
  const std::string message=error.what();
  rejected=message.find("UNSUPPORTED_PAPER_IDEAL_FACE: hidden conic")!=std::string::npos;
  std::cout<<"hidden_conic_diagnostic="<<message<<'\n';
 }
 require(rejected,"unseen ideal conic region not explicitly rejected");
 s.maximumIcosphereSubdivisions=1;rejected=false;
 try{(void)spin::triangulateQuadric(Q::fromUnnormalized(p),domain(),s);}
 catch(const std::invalid_argument&){rejected=true;}
 require(rejected,"paper mapping allowed hidden normal-sphere refinement");
}
void multipleCrossingCaps() {
 // Analytic family from the rank-transition diagnostic: two valid Gauss
 // caps can cut the same normal edge twice. Neither cap may be discarded.
 using Q=spin::QuadricCoefficients10;
 Q::Vector p;p<<.125,0,0,1,1,2,-1e-4,0,0,0;
 const auto q=Q::fromUnnormalized(p);
 spin::QuadricTriangulationSettings s;
 s.icosphereSubdivisions=2;s.maximumIcosphereSubdivisions=2;
 s.truncationRadius=16;s.maximumRadiusDoublings=16;
 s.paperProjectiveHyperboloid=true;s.algebraicFiniteMap=true;
 const auto m=spin::evaluateQuadricSourceMotion(q,domain(),s);
 require(m.defined,"multiple-crossing source derivative unavailable: "+m.diagnostic);
 const auto& t=m.triangulation;
 require(t.boundaryLoopCount==2 && t.effectiveIcosphereSubdivisions==2,
         "two-cap source closure or fixed normal resolution lost");
 int finite=0,upper=0,lower=0;
 for(int i=0;i<t.interfaceVertexCount;++i) {
  const int normal=t.sourceNormalVertexByVertex[i];
  if(normal>=0) {
   ++finite;upper+=t.fixedNormals(normal,2)>0;lower+=t.fixedNormals(normal,2)<0;
  }
 }
 require(finite==146 && upper==73 && lower==73,"a valid normal-cap component was omitted");
 const auto audit=verification::auditExactTriangleBoundary(t.closedBoundary.vertices,
     t.closedBoundary.faces,verification::ExactBoundaryOrientation::OutwardPositive);
 require(audit.facts().boundaryValid,"multiple-crossing material boundary invalid");
 s.truncationRadius=t.truncationRadius;s.maximumRadiusDoublings=0;
 std::cout<<"multiple_crossing geometry_valid=true finite="<<finite<<" radius="<<s.truncationRadius
          <<" vertices="<<t.closedBoundary.vertices.rows()<<" faces="<<t.closedBoundary.faces.rows()<<'\n';
 double maximum=0.;int failures=0;
 // Fixed before the first run: the smallest eigenvalue is O(1e-4),
 // so coefficient steps are O(1e-8), not a rank-changing perturbation.
 for(double h:{1e-8,3e-9,1e-9})for(int j=0;j<10;++j) {
  auto pp=q.values(),pm=q.values();pp(j)+=h;pm(j)-=h;
  const auto plus=spin::triangulateQuadric(Q::fromUnnormalized(pp),domain(),s);
  const auto minus=spin::triangulateQuadric(Q::fromUnnormalized(pm),domain(),s);
  require(plus.closedBoundary.faces==t.closedBoundary.faces &&
          minus.closedBoundary.faces==t.closedBoundary.faces &&
          plus.sourceProjectiveRayByVertex==t.sourceProjectiveRayByVertex &&
          minus.sourceProjectiveRayByVertex==t.sourceProjectiveRayByVertex,
          "multiple-crossing difference changed its discrete branch");
  const Eigen::MatrixXd fd=(plus.closedBoundary.vertices-minus.closedBoundary.vertices)/(2*h);
  Eigen::MatrixXd derivative(fd.rows(),3);
  for(int i=0;i<fd.rows();++i)derivative.row(i)=m.velocity[i].col(j);
  const double error=(fd-derivative).norm(),scale=derivative.norm();
  const double relative=error/std::max(1e-12,scale);
  std::cout<<"multiple_crossing h="<<h<<" column="<<j<<" relative="<<relative<<'\n';
  if(!(scale<1e-8?error<=1e-10:relative<=1e-6))++failures;
  maximum=std::max(maximum,relative);
 }
 require(failures==0,"multiple-crossing source differential mismatches="+std::to_string(failures));
 std::cout<<"[PASS] two normal caps, retained146 finite points, exact boundary, radius="
          <<t.truncationRadius<<" maximum_relative="<<maximum<<'\n';
}
void regularAlgebraicScenario(bool hyperboloid) {
 using Q=spin::QuadricCoefficients10;
 Q::Vector p;p<<-1,.1,-.12,.08,1,1.3,(hyperboloid?-.9:1.7),.03,-.04,.02;
 const auto q=Q::fromUnnormalized(p);
 spin::QuadricTriangulationSettings s;
 s.icosphereSubdivisions=2;s.truncationRadius=16;s.maximumRadiusDoublings=0;
 s.retainFixedNormalProvenance=true;s.algebraicFiniteMap=true;
 const auto m=spin::evaluateQuadricSourceMotion(q,domain(),s);
 require(m.defined,"regular algebraic source motion: "+m.diagnostic);
 const auto& t=m.triangulation;
 require(t.type==(hyperboloid?spin::QuadricSurfaceType::Hyperboloid:
                            spin::QuadricSurfaceType::Ellipsoid),"wrong regular algebraic type");
 require(!t.paperProjectiveHyperboloid,"legacy analytic clipping replaced by projective map");
 const auto audit=verification::auditExactTriangleBoundary(t.closedBoundary.vertices,
     t.closedBoundary.faces,verification::ExactBoundaryOrientation::OutwardPositive);
 require(audit.facts().boundaryValid,"regular algebraic source boundary invalid");
 auto old=s;old.algebraicFiniteMap=false;
 require(spin::triangulateQuadric(q,domain(),old).hash!=t.hash,
         "regular algebraic map missing from source identity");
 spin::QuadricInverseGaussMap map(q);
 if(!hyperboloid)for(int i=0;i<t.closedBoundary.vertices.rows();++i) {
  const auto expected=map.algebraicCoefficientDifferential(t.fixedNormals.row(i));
  require(expected.valid(),"ellipsoid algebraic point invalid");
  require((t.closedBoundary.vertices.row(i).transpose().array()==expected.point.point->array()).all(),
          "ellipsoid construction uses a different point map");
  require((m.velocity[i].array()==expected.coefficientJacobian->array()).all(),
          "ellipsoid motion uses a different point map");
 }
 double maximum=0;
 for(double h:{3e-6,1e-6,3e-7})for(int j=0;j<10;++j) {
  auto pp=q.values(),pm=q.values();pp(j)+=h;pm(j)-=h;
  const auto plus=spin::triangulateQuadric(Q::fromUnnormalized(pp),domain(),s);
  const auto minus=spin::triangulateQuadric(Q::fromUnnormalized(pm),domain(),s);
  require(plus.closedBoundary.faces==t.closedBoundary.faces &&
          minus.closedBoundary.faces==t.closedBoundary.faces &&
          plus.truncationRadius==t.truncationRadius && minus.truncationRadius==t.truncationRadius,
          "regular algebraic finite difference changed discrete branch");
  const Eigen::MatrixXd fd=(plus.closedBoundary.vertices-minus.closedBoundary.vertices)/(2*h);
  Eigen::MatrixXd derivative(fd.rows(),3);
  for(int i=0;i<fd.rows();++i)derivative.row(i)=m.velocity[i].col(j);
  const double error=(fd-derivative).norm(),scale=derivative.norm();
  const double relative=error/std::max(1e-12,scale);
  require(scale<1e-8?error<=1e-10:relative<=1e-6,"regular algebraic source differential mismatch");
  for(int i=0;i<fd.rows();++i) {
   const double e=(fd.row(i)-derivative.row(i)).norm(),d=derivative.row(i).norm();
   require(d<1e-8?e<=1e-10:e/d<=1e-6,"regular algebraic per-vertex differential mismatch");
  }
  maximum=std::max(maximum,relative);
 }
 std::cout<<"regular_algebraic_hyperboloid="<<hyperboloid<<" source="<<t.hash
          <<" vertices="<<t.closedBoundary.vertices.rows()<<" max_relative="<<maximum<<'\n';
}
void algebraicContracts() {
 scenario(false,true);scenario(true,true);
 regularAlgebraicScenario(false);regularAlgebraicScenario(true);
 spin::QuadricCoefficients10::Vector p;p<<0,0,0,1,1,2,0,0,0,0;
 spin::QuadricTriangulationSettings s;s.algebraicFiniteMap=true;
 bool rejected=false;
 try{(void)spin::triangulateQuadric(spin::QuadricCoefficients10::fromUnnormalized(p),domain(),s);}
 catch(const std::invalid_argument& e){rejected=std::string(e.what()).find("algebraic finite-map source requires full rank")!=std::string::npos;}
 require(rejected,"unsupported paraboloid algebraic source silently substituted legacy grid");
 std::cout<<"[PASS] explicit algebraic source policy and its actual differential\n";
}
int capturedState() {
 // Saved restart6/evaluation12; no optimization, radius or coefficient search.
 spin::QuadricCoefficients10::Vector p;
 p<<.4091598840032713,-.36688397145905616,-.019533580749283127,
 -.3602578587235418,.14803212299816212,-.34789822850316837,
 -.25957403943858653,-.3788370952579597,.3838860742170953,.2580806535816809;
 Eigen::MatrixXd design(8,3);
 for(int i=0;i<8;++i)for(int j=0;j<3;++j)design(i,j)=(i&(1<<j))?.25:-.25;
 design.col(2).array()+=.5;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
 int failed=0;
 for(int level:{1,2}) {
  spin::QuadricTriangulationSettings s;s.paperProjectiveHyperboloid=true;
  s.icosphereSubdivisions=level;s.maximumRadiusDoublings=0;
  s.truncationRadius=221.70250336881628;
  try {
   const auto motion=spin::evaluateQuadricSourceMotion(q,design,s);
   require(motion.defined,motion.diagnostic);
   const auto& t=motion.triangulation;
   const auto audit=verification::auditExactTriangleBoundary(t.closedBoundary.vertices,
     t.closedBoundary.faces,verification::ExactBoundaryOrientation::OutwardPositive);
   require(audit.facts().boundaryValid,"captured source boundary rejected");
   std::cout<<"captured_level="<<level<<" source="<<t.hash<<" motion="<<motion.scientificHash
            <<" vertices="<<t.closedBoundary.vertices.rows()<<" faces="<<t.closedBoundary.faces.rows()
            <<" radius="<<t.truncationRadius<<" independent_boundary_valid=true\n";
  }catch(const std::exception& e){++failed;std::cout<<"captured_level="<<level<<" REJECTED="<<e.what()<<'\n';}
 }
 return failed?2:0;
}
}
int main(int argc,char** argv){try{std::cout<<std::setprecision(17);
 if(argc==2 && std::string(argv[1])=="--multiple-crossing") {
  multipleCrossingCaps();return 0;
 }
 if(argc==2 && std::string(argv[1])=="--captured")return capturedState();
 if(argc==2 && std::string(argv[1])=="--algebraic") {
  algebraicContracts();return 0;
 }
 scenario(false);scenario(true);rejectionContracts();algebraicContracts();multipleCrossingCaps();
 std::cout<<"[PASS] paper projective surface and actual source-motion map\n";return 0;
}catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<'\n';return 1;}}
