#include "spin/PolyhedralMomentDifferential.h"
#include "spin/QuadricMomentEvaluator.h"
#include "manufacturing/ManufacturingGeometry.h"
#include "spin/CanonicalFullRankHollowing.h"
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactLibiglBoolean.h"
#include "verification/ExactMaterialRegion.h"
#include "verification/ExactRationalCoordinates.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <igl/copyleft/cgal/assign_scalar.h>
#include <sstream>
#endif
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/MeshBooleanType.h>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace {
using V=Eigen::Matrix<double,10,1>;
using J=Eigen::Matrix<double,10,10>;
void require(bool ok,const std::string& s) { if(!ok) throw std::runtime_error(s); }
spin::TriangleMesh cube() {
 spin::TriangleMesh m; m.vertices.resize(8,3); m.faces.resize(12,3);
 const double a=-.65,b=.8;
 m.vertices<<a,a,a,b,a,a,b,b,a,a,b,a,a,a,b,b,a,b,b,b,b,a,b,b;
 m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
 return m;
}
spin::TriangleMesh tetra() {
 spin::TriangleMesh m; m.vertices.resize(4,3); m.faces.resize(4,3);
 m.vertices<<-.2,-.3,-.4,1.8,-.3,-.4,-.2,1.7,-.4,-.2,-.3,1.6;
 m.faces<<0,2,1,0,1,3,0,3,2,1,2,3;
 const Eigen::Matrix3d rotation=Eigen::AngleAxisd(.23,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
 for(int i=0;i<4;++i) m.vertices.row(i)=rotation*m.vertices.row(i).transpose();
 return m;
}
V integrate(const spin::TriangleMesh& m) {
 V out=V::Zero();
 for(int f=0;f<m.faces.rows();++f) {
  const Eigen::Vector3d a=m.vertices.row(m.faces(f,0)),b=m.vertices.row(m.faces(f,1)),c=m.vertices.row(m.faces(f,2));
  const double v=1.7*a.dot(b.cross(c))/6;
  const Eigen::Vector3d s=a+b+c;
  const Eigen::Matrix3d second=v*(s*s.transpose()+a*a.transpose()+b*b.transpose()+c*c.transpose())/20;
  out(0)+=v;out.segment<3>(1)+=v*s/4;out.segment<3>(4)+=second.diagonal();
  out(7)+=second(1,2);out(8)+=second(0,2);out(9)+=second(0,1);
 }
 return out;
}
struct Cut {spin::TriangleMesh m;Eigen::VectorXi birth;};
Cut cut(const spin::TriangleMesh& d,const spin::TriangleMesh& s,bool minus) {
 Cut o;
 require(igl::copyleft::cgal::mesh_boolean(d.vertices,d.faces,s.vertices,s.faces,
  minus?igl::MESH_BOOLEAN_TYPE_MINUS:igl::MESH_BOOLEAN_TYPE_INTERSECT,o.m.vertices,o.m.faces,o.birth),"Boolean failure");
 require(o.m.faces.rows()>0,"empty fixture");
 return o;
}
void check(bool minus) {
 const auto d=cube(),s=tetra();const auto c=cut(d,s,minus);
 std::vector<spin::SourceVertexJacobian> motion(4,spin::SourceVertexJacobian::Zero());
 // Ten independent actual vertex-coordinate degrees of freedom.
 for(int k=0;k<10;++k) motion[k/3](k%3,k)=1;
 const auto a=spin::differentiatePolyhedralMaterialMoments(d,s,motion,c.m,c.birth,1.7);
 require(a.defined && a.diagnostic.rfind("VALID_FIXED_BOOLEAN_TOPOLOGY: ",0)==0,"transport unavailable: "+a.diagnostic);
 require(a.constructedIntersectionVertexCount>0,"fixture lacks constructed intersections");
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
 const auto exactCut=verification::computeExactLibiglBooleanConstruction(
   d.vertices,d.faces,s.vertices,s.faces,minus?verification::ExactMaterialBooleanOperation::Difference:
   verification::ExactMaterialBooleanOperation::Intersection);
 require(exactCut.succeeded,"exact-output Boolean failed");
 spin::TriangleMesh exactMesh;exactMesh.faces=exactCut.faces;
 exactMesh.vertices.resize(exactCut.vertices.size(),3);
 for(std::size_t v=0;v<exactCut.vertices.size();++v)for(int axis=0;axis<3;++axis)
   igl::copyleft::cgal::assign_scalar(CGAL::Epeck::FT(
     verification::parseExactRationalCoordinate(exactCut.vertices[v][axis])),exactMesh.vertices(v,axis));
 // A separate construction may relabel vertices/faces. Compare the actual
 // oriented triangles, not coincidental array indices; the exact tape below
 // belongs to exactMesh and its OWN face provenance, never to the legacy cut.
 require(spin::sameOrientedTriangleBoundary(exactMesh,c.m),
   "exact-output Boolean changed represented oriented triangles");
 const auto precise=spin::differentiatePolyhedralMaterialMoments(
   d,s,motion,exactMesh,exactCut.birthFaces,1.7,true,true,true,exactCut.vertices);
 require(precise.defined,"actual rational birth coordinates refused: "+precise.diagnostic);
 // Incorrect exact attachments must not borrow the floating boundary's proof.
 auto wrong=exactCut.vertices;wrong.pop_back();
 require(!spin::differentiatePolyhedralMaterialMoments(d,s,motion,exactMesh,exactCut.birthFaces,
   1.7,true,true,true,wrong).defined,"short exact coordinate tape accepted");
 wrong=exactCut.vertices;wrong[0][0]="not-rational";
 require(!spin::differentiatePolyhedralMaterialMoments(d,s,motion,exactMesh,exactCut.birthFaces,
   1.7,true,true,true,wrong).defined,"malformed exact coordinate accepted");
 // This change is much smaller than one ULP: rounded coordinates can agree
 // while the exact point no longer belongs to its actual birth facets.
 wrong=exactCut.vertices;
 const auto x=verification::parseExactRationalCoordinate(wrong[0][0]);
 std::ostringstream altered;altered<<x+CGAL::Gmpq("1/10000000000000000000000000000000000000000");
 wrong[0][0]=altered.str();
 require(!spin::differentiatePolyhedralMaterialMoments(d,s,motion,exactMesh,exactCut.birthFaces,
   1.7,true,true,true,wrong).defined,"same-rounded but off-birth exact point accepted");
#endif
 for(double h:{1e-3,3e-4,1e-4}) {
  J fd;
  for(int k=0;k<10;++k) {
   auto plus=s,minusMesh=s;
   for(int v=0;v<4;++v) {plus.vertices.row(v)+=h*motion[v].col(k).transpose();minusMesh.vertices.row(v)-=h*motion[v].col(k).transpose();}
   fd.col(k)=(integrate(cut(d,plus,minus).m)-integrate(cut(d,minusMesh,minus).m))/(2*h);
  }
  const double error=(fd-a.jacobian).norm()/a.jacobian.norm();
  std::cout<<"FD minus="<<minus<<" h="<<h<<" relative="<<error<<"\n";
  require(error<1e-6,"actual polyhedral finite difference error");
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
  require((fd-precise.jacobian).norm()/precise.jacobian.norm()<1e-6,
    "exact birth motion does not differentiate the actual Boolean map");
#endif
 }
 auto permuted=c;
 permuted.m.faces.col(0).swap(permuted.m.faces.col(1));
 permuted.m.faces.col(1).swap(permuted.m.faces.col(2)); // cyclic orientation-preserving
 const auto reordered=spin::differentiatePolyhedralMaterialMoments(d,s,motion,permuted.m,permuted.birth,1.7);
 require(reordered.defined && (reordered.jacobian-a.jacobian).norm()<1e-10,"cyclic face ordering dependence");
 const auto zero=spin::differentiatePolyhedralMaterialMoments(d,s,
   std::vector<spin::SourceVertexJacobian>(4,spin::SourceVertexJacobian::Zero()),c.m,c.birth,1.7);
 require(zero.defined && zero.jacobian.isZero(0),"static geometry has nonzero derivative");
 for(double rho:{-1.,0.,std::numeric_limits<double>::quiet_NaN()}) {
  const auto invalidDensity=spin::differentiatePolyhedralMaterialMoments(d,s,motion,c.m,c.birth,rho);
  require(!invalidDensity.defined && invalidDensity.diagnostic.rfind("INVALID_INPUT_OR_NONFINITE: ",0)==0,"invalid density classification");
 }
 auto invalid=c;invalid.m.faces(0,0)=999999;
 require(!spin::differentiatePolyhedralMaterialMoments(d,s,motion,invalid.m,c.birth,1.7).defined,"invalid index accepted");
 auto badBirth=c.birth;badBirth(0)=999999;
 require(!spin::differentiatePolyhedralMaterialMoments(d,s,motion,c.m,badBirth,1.7).defined,"invalid birth accepted");
 auto badMotion=motion;badMotion.pop_back();
 require(!spin::differentiatePolyhedralMaterialMoments(d,s,badMotion,c.m,c.birth,1.7).defined,"short motion accepted");
 badMotion=motion;badMotion[0](0,0)=std::numeric_limits<double>::infinity();
 require(!spin::differentiatePolyhedralMaterialMoments(d,s,badMotion,c.m,c.birth,1.7).defined,"nonfinite motion accepted");
 auto badSource=s;badSource.vertices(0,0)=std::numeric_limits<double>::quiet_NaN();
 require(!spin::differentiatePolyhedralMaterialMoments(d,badSource,motion,c.m,c.birth,1.7).defined,"nonfinite source accepted");
 badSource=s;badSource.faces(0,0)=-1;
 require(!spin::differentiatePolyhedralMaterialMoments(d,badSource,motion,c.m,c.birth,1.7).defined,"negative source index accepted");
 auto badDesign=d;badDesign.vertices.conservativeResize(8,2);
 require(!spin::differentiatePolyhedralMaterialMoments(badDesign,s,motion,c.m,c.birth,1.7).defined,"non-3D design accepted");
 auto shortBirth=c.birth;shortBirth.conservativeResize(shortBirth.size()-1);
 require(!spin::differentiatePolyhedralMaterialMoments(d,s,motion,c.m,shortBirth,1.7).defined,"short provenance accepted");
 std::cout<<"[PASS] real Boolean "<<(minus?"difference":"intersection")<<" all100 derivatives and fail-closed inputs\n";
}
}
void edgeContact() {
 auto source=tetra();
 source.vertices<<2,0,0,0,2,0,0,0,.5,0,0,-.5;
 source.faces.col(1).swap(source.faces.col(2));
 auto design=cube();
 for(int i=0;i<design.vertices.rows();++i)for(int j=0;j<3;++j)
  design.vertices(i,j)=design.vertices(i,j)<0?-1:1;
 const auto clipped=cut(design,source,false);
 const auto result=spin::differentiatePolyhedralMaterialMoments(design,source,
   std::vector<spin::SourceVertexJacobian>(4,spin::SourceVertexJacobian::Zero()),clipped.m,clipped.birth,1.7);
 require(!result.defined && result.diagnostic.rfind("TOPOLOGY_TRANSITION: ",0)==0,"edge-edge contact classification");
 auto motion=std::vector<spin::SourceVertexJacobian>(4,spin::SourceVertexJacobian::Zero());
 for(auto& v:motion)v(0,0)=1;
 require(!spin::differentiatePolyhedralMaterialMoments(design,source,motion,clipped.m,clipped.birth,1.7).defined,
   "edge-edge contact granted derivative for translating motion");
 std::cout<<"[PASS] edge-edge contact rejected independently of source velocity\n";
}
void exactSourcePolicy() {
#if !MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
 const auto d=cube(),s=tetra();const auto c=cut(d,s,false);
 const auto result=spin::differentiatePolyhedralMaterialMoments(d,s,
   std::vector<spin::SourceVertexJacobian>(4,spin::SourceVertexJacobian::Zero()),
   c.m,c.birth,1.7,false,true);
 require(!result.defined,"source policy granted without exact checker support");
 return;
#else
 auto s=tetra();
 s.vertices<<-2,-2,-2,4,-2,-2,-2,4,-2,-2,-2,4;
 // Subdivide a source face OUTSIDE the clipped material, using a dyadic point
 // with strictly nonzero exact face areas. The represented source is unchanged.
 const auto face=s.faces.row(0).eval();
 s.vertices.conservativeResize(5,3);
 const double t=std::ldexp(1.,-40);
 s.vertices.row(4)=s.vertices.row(face(0))+t*(s.vertices.row(face(1))+
   s.vertices.row(face(2))-2*s.vertices.row(face(0)));
 s.faces.conservativeResize(6,3);
 s.faces.row(0)<<face(0),face(1),4;
 s.faces.row(4)<<face(1),face(2),4;
 s.faces.row(5)<<face(2),face(0),4;
 const auto d=cube();const auto c=cut(d,s,false);
 std::vector<spin::SourceVertexJacobian> motion(5,spin::SourceVertexJacobian::Zero());
 for(auto& velocity:motion)velocity(0,0)=1;
 const auto legacy=spin::differentiatePolyhedralMaterialMoments(d,s,motion,c.m,c.birth,1.7,true);
 require(!legacy.defined,"legacy source area gate silently changed");
 const auto exact=spin::differentiatePolyhedralMaterialMoments(d,s,motion,c.m,c.birth,1.7,true,true);
 require(exact.defined,"complete exact source policy failed: "+exact.diagnostic);
 for(double h:{1e-3,3e-4,1e-4}) {
  auto plus=s,minus=s;plus.vertices.col(0).array()+=h;minus.vertices.col(0).array()-=h;
  // Own the value; an unevaluated Eigen expression would dangle after its
  // temporary integration vectors are destroyed at the semicolon.
  const V fd=(integrate(cut(d,plus,false).m)-integrate(cut(d,minus,false).m))/(2*h);
  const double error=(fd-exact.jacobian.col(0)).norm()/std::max(1.,exact.jacobian.col(0).norm());
  std::cout<<"source-only FD h="<<h<<" relative="<<error<<'\n';
  require(error<1e-6,
    "exact-source policy detached Jacobian from actual Boolean map");
 }
 auto invalid=s;invalid.vertices.row(4)=invalid.vertices.row(face(0));
 require(!spin::differentiatePolyhedralMaterialMoments(d,invalid,motion,c.m,c.birth,1.7,true,true).defined,
   "exact source accepted true zero area");
 invalid=s;invalid.faces.col(1).swap(invalid.faces.col(2));
 require(!spin::differentiatePolyhedralMaterialMoments(d,invalid,motion,c.m,c.birth,1.7,true,true).defined,
   "exact source accepted inward winding");
 invalid=s;invalid.faces.row(5)=invalid.faces.row(4);
 require(!spin::differentiatePolyhedralMaterialMoments(d,invalid,motion,c.m,c.birth,1.7,true,true).defined,
   "exact source accepted nonmanifold boundary");
 auto intersecting=s;intersecting.vertices.conservativeResize(10,3);intersecting.faces.conservativeResize(12,3);
 for(int v=0;v<5;++v)intersecting.vertices.row(v+5)=s.vertices.row(v)+Eigen::RowVector3d(.1,.2,.3);
 for(int f=0;f<6;++f)intersecting.faces.row(f+6)=s.faces.row(f).array()+5;
 auto doubledMotion=motion;doubledMotion.insert(doubledMotion.end(),motion.begin(),motion.end());
 require(!spin::differentiatePolyhedralMaterialMoments(d,intersecting,doubledMotion,c.m,c.birth,1.7,true,true).defined,
   "exact source accepted intersecting closed components");
 const auto invalidDesign=spin::differentiatePolyhedralMaterialMoments(s,s,motion,c.m,c.birth,1.7,true,true);
 require(!invalidDesign.defined && invalidDesign.diagnostic.find("audited boundary")!=std::string::npos,
   "source policy leaked into design-boundary area audit");
 // SOURCE opt-in must NEVER qualify a material whose tiny faces fail the
 // existing final-material gate, even if that same mesh is valid as a source.
 Eigen::VectorXi births(6);births<<12,13,14,15,16,17;
 auto enclosing=cube();enclosing.vertices*=20;
 require(!spin::differentiatePolyhedralMaterialMoments(enclosing,s,motion,s,births,1.7,true,true).defined,
   "source policy leaked into final-material area audit");
 std::cout<<"[PASS] source-only exact boundary policy and three-step actual-map differential\n";
#endif
}
void exactMaterialPolicy() {
#if !MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
 bool refused=false;
 try {(void)spin::auditQuadricMaterialBoundary(cube(),false,true);}
 catch(const std::runtime_error&) {refused=true;}
 require(refused,"exact material proof silently unavailable");
#else
 auto m=cube();const auto face=m.faces.row(0).eval();
 m.vertices.conservativeResize(9,3);m.faces.conservativeResize(14,3);
 const double t=std::ldexp(1.,-40);
 m.vertices.row(8)=m.vertices.row(face(0))+t*(m.vertices.row(face(1))+
   m.vertices.row(face(2))-2*m.vertices.row(face(0)));
 m.faces.row(0)<<face(0),face(1),8;
 m.faces.row(12)<<face(1),face(2),8;m.faces.row(13)<<face(2),face(0),8;
 require(!spin::auditQuadricMaterialBoundary(m,true).strictBoundaryValid(),
   "legacy material heuristic silently changed");
 require(spin::auditQuadricMaterialBoundary(m,false,true).strictBoundaryValid(),
   "complete exact material proof did not qualify nonzero faces");
 require((spin::integrateCanonicalFullRankBoundaryMoments(m,1.7,false,true).values()-integrate(m)).norm()<1e-10,
   "independent final integration did not use the same material policy");
 auto bad=m;bad.vertices.row(8)=bad.vertices.row(face(0));
 require(!spin::auditQuadricMaterialBoundary(bad,false,true).strictBoundaryValid(),"true zero material area accepted");
 bad=m;bad.faces.col(1).swap(bad.faces.col(2));
 require(!spin::auditQuadricMaterialBoundary(bad,false,true).strictBoundaryValid(),"inward whole material accepted");
 bad=m;bad.faces.row(13)=bad.faces.row(12);
 require(!spin::auditQuadricMaterialBoundary(bad,false,true).strictBoundaryValid(),"nonmanifold material accepted");
 bad=m;bad.faces.conservativeResize(13,3);
 require(!spin::auditQuadricMaterialBoundary(bad,false,true).strictBoundaryValid(),"open material accepted");
 auto overlap=cube();overlap.vertices.conservativeResize(16,3);overlap.faces.conservativeResize(24,3);
 const auto original=cube();
 for(int v=0;v<8;++v)overlap.vertices.row(v+8)=original.vertices.row(v)+Eigen::RowVector3d(.1,.2,.3);
 for(int f=0;f<12;++f)overlap.faces.row(f+12)=original.faces.row(f).array()+8;
 require(!spin::auditQuadricMaterialBoundary(overlap,false,true).strictBoundaryValid(),"crossing material shells accepted");
 // A proper hollow solid needs an inward INNER component. A blanket
 // OutwardPositive-per-component audit incorrectly rejects this legal case.
 auto cavity=cube();cavity.vertices.conservativeResize(16,3);cavity.faces.conservativeResize(24,3);
 for(int v=0;v<8;++v)cavity.vertices.row(v+8)=original.vertices.row(v)*.5;
 for(int f=0;f<12;++f)cavity.faces.row(f+12)<<original.faces(f,2)+8,original.faces(f,1)+8,original.faces(f,0)+8;
 require(spin::auditQuadricMaterialBoundary(cavity,false,true).strictBoundaryValid(),"proper inward cavity boundary rejected");
 cavity.faces.bottomRows(12).col(1).swap(cavity.faces.bottomRows(12).col(2));
 require(!spin::auditQuadricMaterialBoundary(cavity,false,true).strictBoundaryValid(),"incorrect positive nested shell accepted");
 std::cout<<"[PASS] exact material geometric proof, hollow-shell semantics and unchanged legacy gate\n";
#endif
}
void rationalMaterialAuthorityReuse() {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
 // A translated tetrahedron has hand-derived moment derivatives. Its 1/3
 // coordinates are actual rational authority, not exact binary64 values.
 spin::ExactBooleanVertexCoordinates exact={{{"0/1","0/1","0/1"}},
  {{"1/3","0/1","0/1"}},{{"0/1","1/3","0/1"}},{{"0/1","0/1","1/3"}}};
 spin::TriangleMesh source;source.vertices.resize(4,3);source.faces.resize(4,3);
 source.vertices<<0,0,0,1./3,0,0,0,1./3,0,0,0,1./3;
 source.faces<<0,2,1,0,1,3,0,3,2,1,2,3;
 const auto design=cube();spin::ExactBooleanVertexCoordinates designExact;
 for(int v=0;v<design.vertices.rows();++v) {
  std::array<std::string,3> row;
  for(int axis=0;axis<3;++axis) {std::ostringstream text;text<<CGAL::Gmpq(design.vertices(v,axis));row[axis]=text.str();}
  designExact.push_back(row);
 }
 const auto proof=verification::makeExactMaterialRegionFromRationals(exact,source.faces);
 require(proof.facts().valid && proof.matchesRationalBoundary(exact,source.faces),"actual material proof does not bind its own input");
 auto alias=exact;alias[1][0]="2/6";
 require(proof.matchesRationalBoundary(alias,source.faces),"equal rational values do not bind");
 auto wrong=exact;wrong[1][0]="10000000000000000000000000000000000000001/30000000000000000000000000000000000000000";
 require(!proof.matchesRationalBoundary(wrong,source.faces),"same-rounded unequal rational point borrowed authority");
 auto wrongFaces=source.faces;wrongFaces(0,0)=1;
 require(!proof.matchesRationalBoundary(exact,wrongFaces),"changed face index borrowed authority");
 wrong=exact;wrong[0][0]="not-a-rational";
 require(!proof.matchesRationalBoundary(wrong,source.faces),"malformed coordinate borrowed authority");
 require(!verification::emptyExactMaterialRegion().matchesRationalBoundary(exact,source.faces),"empty proof borrowed by nonempty material");
 const auto wrongProof=verification::makeExactMaterialRegionFromRationals(designExact,design.faces);
 std::vector<spin::SourceVertexJacobian> velocity(4,spin::SourceVertexJacobian::Zero());
 for(auto& v:velocity)v(0,0)=1.;
 Eigen::VectorXi births(4);births<<12,13,14,15;
 int oldMaterialAudits=0,reusedMaterialAudits=0;
 auto oldObserver=[&](const std::string& stage){if(stage=="transport-material-authority-audit-begin")++oldMaterialAudits;};
 auto newObserver=[&](const std::string& stage){if(stage=="transport-material-authority-audit-begin")++reusedMaterialAudits;};
 const auto baseline=spin::differentiatePolyhedralMaterialMoments(design,source,velocity,source,births,
  1.7,true,true,true,exact,exact,designExact,oldObserver);
 const auto reused=spin::differentiatePolyhedralMaterialMoments(design,source,velocity,source,births,
  1.7,true,true,true,exact,exact,designExact,newObserver,&proof);
 require(baseline.defined && reused.defined,"real rational transport unavailable: "+baseline.diagnostic+";"+reused.diagnostic);
 require((baseline.jacobian.array()==reused.jacobian.array()).all() &&
  baseline.vertexProvenance==reused.vertexProvenance,"same authority reuse changes actual derivative/provenance");
 const double mass=1.7/162.;V expected=V::Zero();expected(1)=mass;expected(4)=mass/6.;expected(8)=mass/12.;expected(9)=mass/12.;
 require((reused.jacobian.col(0)-expected).norm()<1e-12,"translated rational tetrahedron does not match hand-derived moments");
 require(oldMaterialAudits==1 && reusedMaterialAudits==0,"matching proof did not avoid the duplicate material construction");
 require(!spin::differentiatePolyhedralMaterialMoments(design,source,velocity,source,births,
  1.7,true,true,true,exact,exact,designExact,{},&wrongProof).defined,"wrong region proof accepted");
 require(!spin::differentiatePolyhedralMaterialMoments(design,source,velocity,source,births,
  1.7,true,true,true,wrong,exact,designExact,{},&proof).defined,"malformed material tape accepted with valid proof");
 require(!spin::differentiatePolyhedralMaterialMoments(design,source,velocity,source,births,
  -1.,true,true,true,exact,exact,designExact,{},&proof).defined,"reused proof concealed invalid density");
 auto moved=source;moved.vertices(1,0)=std::nextafter(moved.vertices(1,0),1.);
 require(!spin::differentiatePolyhedralMaterialMoments(design,source,velocity,moved,births,
  1.7,true,true,true,exact,exact,designExact,{},&proof).defined,"display assignment mismatch accepted with exact proof");
 std::cout<<"[PASS] matching immutable material authority reuse; original rational derivative and fail-closed bindings\n";
#else
 throw std::runtime_error("rational material authority reuse requires exact checks");
#endif
}
int main(int argc,char** argv) {try {std::cout<<std::setprecision(17);
 if(argc==2 && std::string(argv[1])=="--rational-authority-reuse") {rationalMaterialAuthorityReuse();return 0;}
 check(false);check(true);edgeContact();exactSourcePolicy();exactMaterialPolicy();return 0;}
 catch(const std::exception& e) {std::cerr<<"[FAIL] "<<e.what()<<"\n";return 1;}}
