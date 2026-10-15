#include "spin/FullRankMomentEvaluator.h"
#include "spin/QuadricSourceMotion.h"
#include <Eigen/Geometry>
#include <iostream>
#include <iomanip>
#include <stdexcept>

namespace {
using V=spin::QuadricCoefficients10::Vector;
using J=Eigen::Matrix<double,10,10>;
void require(bool b,const std::string& s) {if(!b)throw std::runtime_error(s);}
spin::TriangleMesh box(double a,double b) {
 spin::TriangleMesh m;m.vertices.resize(8,3);m.faces.resize(12,3);
 m.vertices<<a,a,a,b,a,a,b,b,a,a,b,a,a,a,b,b,a,b,b,b,b,a,b,b;
 m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;return m;
}
void check(const V& p,bool transformed) {
 spin::SpinDesignRequest r;r.designDomain=box(-1,1);r.outerSurface=box(-2,2);
 r.minimumWallThickness=.2;r.materialDensity=1.7;
 if(transformed) {
  r.targetFrame.origin<<.3,-.2,.1;
  r.targetFrame.axes=Eigen::AngleAxisd(.43,Eigen::Vector3d(1,2,-1).normalized()).toRotationMatrix();
  for(auto* m:{&r.designDomain,&r.outerSurface})for(int i=0;i<m->vertices.rows();++i)
   m->vertices.row(i)=r.targetFrame.axes*m->vertices.row(i).transpose()+r.targetFrame.origin;
 }
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
 spin::QuadricMomentEvaluatorSettings settings;settings.triangulation.icosphereSubdivisions=2;settings.triangulation.truncationRadius=8;
 settings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteEllipsoid;
 const auto base=spin::evaluateFullRankMaterialMoments(r,q,settings);
 require(base.coefficientJacobian.has_value(),"full-rank J unavailable: "+base.diagnostic);
 const auto old=spin::evaluateQuadricMaterialMoments(r,q,settings);
 require(base.geometryHash==old.geometryHash && base.variableMomentsTarget.values()==old.variableMomentsTarget.values(),"changed native geometry or moments");
 double previous=1;
 for(double h:{1e-4,3e-5,1e-5}) {
  J fd;
  for(int j=0;j<10;++j) {
   auto plus=q.values(),minus=q.values();plus(j)+=h;minus(j)-=h;
   const auto a=spin::evaluateQuadricMaterialMoments(r,spin::QuadricCoefficients10::fromUnnormalized(plus),settings);
   const auto b=spin::evaluateQuadricMaterialMoments(r,spin::QuadricCoefficients10::fromUnnormalized(minus),settings);
   fd.col(j)=(a.variableMomentsTarget.values()-b.variableMomentsTarget.values())/(2*h);
  }
  const double error=(fd-*base.coefficientJacobian).norm()/base.coefficientJacobian->norm();
  std::cout<<"full-map frame="<<transformed<<" h="<<h<<" relative="<<error<<"\n";
  if(h<2e-5)require(error<1e-6,"native full-map derivative error");
  if(h<1e-4 && previous>1e-8)require(error<.3*previous,"full-map FD did not converge");
  previous=error;
 }
 require((*base.coefficientJacobian*q.values()).norm()<1e-9,"full-map radial null");
 const auto again=spin::evaluateFullRankMaterialMoments(r,q,settings);
 require(again.coefficientJacobian.has_value() && again.scientificHash==base.scientificHash,"scientific nondeterminism");
 auto withShell=r;
 withShell.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(
   r.outerSurface.vertices,r.outerSurface.faces,r.materialDensity))-
   physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.designDomain.vertices,r.designDomain.faces,r.materialDensity));
 const auto shell=spin::evaluateFullRankMaterialMoments(withShell,q,settings);
 const auto oldShell=spin::evaluateQuadricMaterialMoments(withShell,q,settings);
 require(shell.coefficientJacobian && *shell.coefficientJacobian==*base.coefficientJacobian &&
  shell.totalMomentsTarget.values()==oldShell.totalMomentsTarget.values() &&
  shell.totalMomentsTarget.mass()>base.totalMomentsTarget.mass(),"fixed-shell total or derivative changed");
 auto denser=r;denser.materialDensity*=2;
 const auto doubled=spin::evaluateFullRankMaterialMoments(denser,q,settings);
 require(doubled.coefficientJacobian && (*doubled.coefficientJacobian-2*(*base.coefficientJacobian)).norm()<1e-10 &&
  doubled.scientificHash!=base.scientificHash,"density scaling/binding failed");
 std::cout<<"[PASS] full-rank native material all100 derivatives and unchanged values\n";
}
void unavailablePlane() {
 spin::SpinDesignRequest r;r.designDomain=box(-1,1);r.outerSurface=box(-2,2);r.minimumWallThickness=.2;
 V p;p<<.2,0,0,-1,0,0,0,0,0,0;
 const auto result=spin::evaluateFullRankMaterialMoments(r,spin::QuadricCoefficients10::fromUnnormalized(p));
 require(!result.coefficientJacobian && result.variableMomentsTarget.mass()>0,"unsupported source did not preserve values with absent J");
 std::cout<<"[PASS] rank-deficient native values retained without fabricated Jacobian\n";
}
void strictInactiveSelector() {
 spin::SpinDesignRequest r;r.designDomain=box(-1,1);r.outerSurface=box(-2,2);r.minimumWallThickness=.2;
 Eigen::Matrix3d a;a<<.83,.03,-.04,.03,1.17,.06,-.04,.06,-.72;
 const Eigen::Vector3d center(.1,.2,.37717397599018909);
 V p;p(0)=center.dot(a*center)+.31;p.segment<3>(1)=-2*a*center;
 p.segment<3>(4)=a.diagonal();p(7)=a(1,2);p(8)=a(0,2);p(9)=a(0,1);
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
 spin::QuadricMomentEvaluatorSettings settings;settings.triangulation.icosphereSubdivisions=2;settings.triangulation.truncationRadius=8;
 const auto source=spin::evaluateQuadricSourceMotion(q,r.designDomain.vertices,settings.triangulation);
 require(source.defined,"selector fixture lacks a regular source-motion branch");
 Eigen::Vector3d centroid;
 const auto mass=computeSurfaceMassMoments(source.triangulation.closedBoundary.vertices,source.triangulation.closedBoundary.faces,1);
 require(mass.tryCenterOfMass(centroid,1e-15) && std::abs(q.evaluate(centroid))<1e-10,"selector fixture not inside inactive-centroid band");
 const auto result=spin::evaluateFullRankMaterialMoments(r,q,settings);
 require(result.coefficientJacobian.has_value(),
   "strict inactive selector with regular source incorrectly rejected derivative");
 std::cout<<"[PASS] regular source and strict inactive selector retain local derivative\n";
}
}
int main(){try{std::cout<<std::setprecision(17);V p;p<<-.31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03;
 strictInactiveSelector();check(p,false);p(0)=.31;check(p,false);check(-p,true);unavailablePlane();return 0;
}catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<"\n";return 1;}}
