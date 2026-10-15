#include "spin/PolyhedralMomentDifferential.h"
#include "verification/ExactLibiglBoolean.h"
#include "verification/ExactMaterialRegion.h"
#include "verification/ExactRationalCoordinates.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <igl/copyleft/cgal/assign_scalar.h>
#include <gmpxx.h>
#include <iostream>
#include <stdexcept>

using Tape=spin::ExactBooleanVertexCoordinates;
using Q=mpq_class;
static std::string text(Q q) {q.canonicalize();return q.get_num().get_str()+"/"+q.get_den().get_str();}
static void require(bool condition,const char* why) {if(!condition)throw std::runtime_error(why);}
static spin::TriangleMesh cube(const Q& xmin,const Q& xmax,Tape* tape) {
 spin::TriangleMesh m;m.vertices.resize(8,3);m.faces.resize(12,3);
 const std::array<std::array<Q,3>,8> points{{{{xmin,-1,-1}},{{xmax,-1,-1}},{{xmax,1,-1}},{{xmin,1,-1}},
     {{xmin,-1,1}},{{xmax,-1,1}},{{xmax,1,1}},{{xmin,1,1}}}};
 tape->resize(8);
 for(int v=0;v<8;++v)for(int a=0;a<3;++a){(*tape)[v][a]=text(points[v][a]);m.vertices(v,a)=points[v][a].get_d();}
 m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
 return m;
}
static void oneCase(bool clipped) {
 Tape designTape,sourceTape;
 auto source=cube(-1,1,&sourceTape);
 auto design=cube(-2,clipped?Q(1,3):Q(2),&designTape);
 // Enlarge only the design y,z dimensions to avoid coplanar contacts.
 for(int v=0;v<8;++v)for(int a=1;a<3;++a){
   const Q position=Q(designTape[v][a])*2+Q(a==1?1:2,8);
   design.vertices(v,a)=position.get_d();designTape[v][a]=text(position);
 }
 auto boolean=verification::computeExactLibiglBooleanConstructionFromRationals(
     designTape,design.faces,sourceTape,source.faces,verification::ExactMaterialBooleanOperation::Intersection);
 require(boolean.succeeded,"actual rational Boolean failed");
 spin::TriangleMesh material;material.faces=boolean.faces;material.vertices.resize(boolean.vertices.size(),3);
 // Assign only the compatibility view, using the pinned Boolean convention.
 for(std::size_t v=0;v<boolean.vertices.size();++v)for(int a=0;a<3;++a)
   igl::copyleft::cgal::assign_scalar(CGAL::Epeck::FT(
       verification::parseExactRationalCoordinate(boolean.vertices[v][a])),material.vertices(v,a));
 const Q speed=1+Q(1,mpz_class("10000000000000000000000000"));
 std::vector<spin::RationalSourceVertexJacobian> velocities(8);
 for(auto& velocity:velocities){for(auto& row:velocity)row.fill("0/1");velocity[0][1]=text(speed);}
 const auto evaluate=[&](const auto& motion,double density,const Tape& sourceAuthority) {
  return spin::differentiateRationalPolyhedralMaterialMoments(design,source,motion,material,
      boolean.birthFaces,density,boolean.vertices,sourceAuthority,designTape);
 };
 auto result=evaluate(velocities,1,sourceTape);
 require(result.defined,result.diagnostic.c_str());
 const std::array<Q,10> expected=clipped
    ?std::array<Q,10>{{-4,4,0,0,-4,Q(-4,3),Q(-4,3),0,0,0}}
    :std::array<Q,10>{{0,8,0,0,0,0,0,0,0,0}};
 for(int m=0;m<10;++m)for(int j=0;j<10;++j)
   require(Q(result.coefficientJacobian[m][j])==(j==1?Q(expected[m]*speed):Q(0)),
       "rational Boolean transport lost sub-binary64 velocity or moment flux");
 require(!clipped || result.summary.constructedIntersectionVertexCount>0,
     "clipped case did not exercise actual intersection transport");
 const auto proof=verification::makeExactMaterialRegionFromRationals(boolean.vertices,material.faces);
 const auto reused=spin::differentiateRationalPolyhedralMaterialMoments(design,source,velocities,material,
     boolean.birthFaces,1,boolean.vertices,sourceTape,designTape,{},&proof);
 require(reused.defined && reused.coefficientJacobian==result.coefficientJacobian,
     "matching immutable material proof changed authority derivative");
 const auto wrongProof=verification::makeExactMaterialRegionFromRationals(designTape,design.faces);
 require(!spin::differentiateRationalPolyhedralMaterialMoments(design,source,velocities,material,
     boolean.birthFaces,1,boolean.vertices,sourceTape,designTape,{},&wrongProof).defined,
     "unrelated region proof concealed authority mismatch");
 auto dense=evaluate(velocities,2,sourceTape);require(dense.defined,dense.diagnostic.c_str());
 for(int m=0;m<10;++m)require(Q(dense.coefficientJacobian[m][1])==2*Q(result.coefficientJacobian[m][1]),
     "density missing or multiplied twice in rational derivative");
 auto malformed=velocities;malformed[0][0][1]="nan";
 require(!evaluate(malformed,1,sourceTape).defined,"non-rational source velocity accepted");
 auto stale=sourceTape;stale[0][0]="-999/1000";
 require(!evaluate(velocities,1,stale).defined,"mismatched source authority accepted");
 require(!evaluate(velocities,0,sourceTape).defined,"invalid density accepted");
 auto shortMotion=velocities;shortMotion.pop_back();
 require(!evaluate(shortMotion,1,sourceTape).defined,"missing source velocity accepted");
 auto offBirth=boolean.vertices;offBirth[0][0]=text(Q(offBirth[0][0])+Q(1,mpz_class("10000000000000000000000000")));
 require(!spin::differentiateRationalPolyhedralMaterialMoments(design,source,velocities,material,
     boolean.birthFaces,1,offBirth,sourceTape,designTape).defined,
     "same-rounded material coordinate off its actual birth facet accepted");
 std::cout<<(clipped?"intersection":"source-only")<<" exact dLx="<<result.coefficientJacobian[1][1]<<'\n';
}
int main(){try{oneCase(false);oneCase(true);std::cout<<"PASS: actual rational Boolean transport and exact moment derivatives\n";}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
