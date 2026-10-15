#include "spin/QuadricInverseGaussMap.h"
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
using Q=spin::QuadricCoefficients10;
void require(bool b,const char* message) { if(!b)throw std::runtime_error(message); }
void checkIdealRoots() {
 Q::Vector p;p<<-1,0,0,0,1,1,-1,0,0,0;
 const auto q=Q::fromUnnormalized(p);
 const spin::QuadricInverseGaussMap map(q);
 const auto one=map.idealEdgeDifferential(Eigen::Vector3d::UnitX(),Eigen::Vector3d::UnitZ());
 require(one.defined && one.points.size()==1,"paper ideal edge missing its one analytic root");
 require(std::abs(one.points[0].parameter-.5)<1e-14,"wrong ideal parameter");
 require((one.points[0].direction-Eigen::Vector3d(1,0,-1)).norm()<1e-13,"wrong ideal homogeneous direction");
 const Eigen::Vector3d a(-1,0,.5),b(1,0,.5);
 const auto two=map.idealEdgeDifferential(a,b),reverse=map.idealEdgeDifferential(b,a);
 require(two.defined && reverse.defined && two.points.size()==2 && reverse.points.size()==2,
         "two ideal roots were reduced to endpoint XOR");
 require(std::abs(two.points[0].parameter-.25)<1e-14 &&
         std::abs(two.points[1].parameter-.75)<1e-14,"wrong two-root order");
 for(int i=0;i<2;++i) {
  require((two.points[i].direction-reverse.points[1-i].direction).norm()<1e-13,
          "ideal identity changes with edge reversal");
  require((two.points[i].coefficientJacobian-reverse.points[1-i].coefficientJacobian).norm()<1e-12,
          "ideal derivative changes with edge reversal");
 }
 const auto empty=map.idealEdgeDifferential(Eigen::Vector3d::UnitX(),Eigen::Vector3d::UnitY());
 require(empty.defined && empty.points.empty(),"no-root edge is not certified empty");
 const auto endpoint=map.idealEdgeDifferential(Eigen::Vector3d(1,0,1),Eigen::Vector3d::UnitX());
 require(!endpoint.defined && endpoint.points.empty(),"ideal endpoint event not rejected");
 // n=(1,t,1) touches the ideal conic at t=0.
 const auto tangent=map.idealEdgeDifferential(Eigen::Vector3d(1,-1,1),Eigen::Vector3d(1,1,1));
 require(!tangent.defined && tangent.points.empty(),"tangent ideal event not rejected");
}
void checkRay() {
 const Eigen::Vector3d P(1,0,0),D(1,0,-1);
 const auto r=spin::truncateQuadricProjectiveRay(P,D,2);
 require(r.defined,"paper projective truncation missing");
 const double t=(std::sqrt(7.)-1)/2;
 require((r.point-Eigen::Vector3d(1+t,0,-t)).norm()<1e-13,"wrong independent ray/sphere intersection");
 // The polygon edge is not the curved analytic-quadric edge. Its truncation
 // vertex generally does NOT satisfy q(x)=0; forcing that would change maps.
 const double q=-1+r.point.x()*r.point.x()-r.point.z()*r.point.z();
 require(std::abs(q-2*t)<1e-13 && q>1,"projective edge was replaced by an analytic quadric arc");
 const auto gauge=spin::truncateQuadricProjectiveRay(P,8*D,2);
 require(gauge.defined && (gauge.point-r.point).norm()<1e-13,"positive ideal gauge changes the edge");
 for(double h:{3e-6,1e-6,3e-7})for(int j=0;j<7;++j) {
  auto pp=P,pm=P,dp=D,dm=D;double rp=2,rm=2;
  if(j<3){pp(j)+=h;pm(j)-=h;}else if(j<6){dp(j-3)+=h;dm(j-3)-=h;}else{rp+=h;rm-=h;}
  const auto plus=spin::truncateQuadricProjectiveRay(pp,dp,rp);
  const auto minus=spin::truncateQuadricProjectiveRay(pm,dm,rm);
  require(plus.defined && minus.defined,"ray FD became undefined");
  require(((plus.point-minus.point)/(2*h)-r.jacobian.col(j)).norm()<=1e-6,
          "projective ray Jacobian disagrees with independent FD");
 }
 require(!spin::truncateQuadricProjectiveRay(2*P,D,2).defined,"radius endpoint event accepted");
 require(!spin::truncateQuadricProjectiveRay(3*P,D,2).defined,"outside ray origin accepted");
 require(!spin::truncateQuadricProjectiveRay(P,Eigen::Vector3d::Zero(),2).defined,"zero ideal direction accepted");
 require(!spin::truncateQuadricProjectiveRay(P,D,std::numeric_limits<double>::infinity()).defined,"infinite radius accepted");
 // Exact binary64 input lies outside: |P|^2-1 = +1.1507582378374553e-21.
 // An extended-precision sum may round this negative, so it is not a sign proof.
 const Eigen::Vector3d outside(0x1.6a4934bdf1a0cp-2,
                            0x1.dee233e05d125p-1,0x1.4f30ad0b5cdb5p-9);
 const auto rejected=spin::truncateQuadricProjectiveRay(outside,-outside,1);
 require(!rejected.defined && rejected.point.isZero() && rejected.jacobian.isZero(),
         "exactly outside origin accepted through rounded interior sign");
}
void checkComposedDerivative() {
 Q::Vector p;p<<-1,.1,-.12,.08,1,1.3,-.9,.03,-.04,.02;
 const auto q=Q::fromUnnormalized(p);
 const Eigen::Vector3d n0(1,.2,.1),n1(.1,.1,1);
 const spin::QuadricInverseGaussMap map(q);
 const auto finite=map.coefficientDifferential(n0);
 const auto ideal=map.idealEdgeDifferential(n0,n1);
 require(finite.valid() && ideal.defined && ideal.points.size()==1,"non-axis-aligned projective fixture invalid");
 const auto ray=spin::truncateQuadricProjectiveRay(*finite.point.point,ideal.points[0].direction,4);
 require(ray.defined,"composed projective ray invalid");
 const Eigen::Matrix<double,3,10> combined=ray.jacobian.leftCols<3>()*(*finite.coefficientJacobian)+
   ray.jacobian.middleCols<3>(3)*ideal.points[0].coefficientJacobian;
 for(double h:{3e-6,1e-6,3e-7})for(int j=0;j<10;++j) {
  auto pp=q.values(),pm=q.values();pp(j)+=h;pm(j)-=h;
  const spin::QuadricInverseGaussMap mp(Q::fromUnnormalized(pp)),mm(Q::fromUnnormalized(pm));
  const auto ip=mp.idealEdgeDifferential(n0,n1),im=mm.idealEdgeDifferential(n0,n1);
  require(ip.defined && im.defined && ip.points.size()==1 && im.points.size()==1,"ideal FD changed branch");
  const auto fp=mp.point(n0),fm=mm.point(n0);require(fp.valid() && fm.valid(),"finite FD invalid");
  const auto xp=spin::truncateQuadricProjectiveRay(*fp.point,ip.points[0].direction,4);
  const auto xm=spin::truncateQuadricProjectiveRay(*fm.point,im.points[0].direction,4);
  require(xp.defined && xm.defined,"composed FD invalid");
  const double idealError=((ip.points[0].direction-im.points[0].direction)/(2*h)-
                          ideal.points[0].coefficientJacobian.col(j)).norm();
  const double error=((xp.point-xm.point)/(2*h)-combined.col(j)).norm();
  const double idealScale=ideal.points[0].coefficientJacobian.col(j).norm();
  const double scale=combined.col(j).norm();
  const double idealRelative=idealError/std::max(1e-12,idealScale);
  const double relative=error/std::max(1e-12,scale);
  std::cout<<"h="<<h<<" column="<<j<<" ideal_abs="<<idealError<<" composed_abs="<<error
           <<" ideal_rel="<<idealRelative<<" composed_rel="<<relative<<'\n';
  require((idealScale<1e-8 ? idealError<=1e-10 : idealRelative<=1e-6) &&
          (scale<1e-8 ? error<=1e-10 : relative<=1e-6),
          "paper composed actual-map derivative mismatch");
 }
 require((combined*q.values()).norm()<1e-11,"coefficient radial null direction changed geometry");
}
void checkOtherRayBranches() {
 const Eigen::Vector3d origin(1,0,0);
 for(const Eigen::Vector3d& direction :
     {Eigen::Vector3d(-1,0,1),Eigen::Vector3d(0,1,-1)}) {
  const auto base=spin::truncateQuadricProjectiveRay(origin,direction,2);
  const double a=direction.squaredNorm(),b=origin.dot(direction);
  const double t=(-b+std::sqrt(b*b+3*a))/a;
  require(base.defined && (base.point-origin-t*direction).norm()<1e-13,
          "negative/zero b branch disagrees with independent quadratic exit");
  for(double h:{3e-6,1e-6,3e-7})for(int j=0;j<7;++j) {
   auto pp=origin,pm=origin,dp=direction,dm=direction;double rp=2,rm=2;
   if(j<3){pp(j)+=h;pm(j)-=h;}else if(j<6){dp(j-3)+=h;dm(j-3)-=h;}
   else{rp+=h;rm-=h;}
   const auto plus=spin::truncateQuadricProjectiveRay(pp,dp,rp);
   const auto minus=spin::truncateQuadricProjectiveRay(pm,dm,rm);
   require(plus.defined && minus.defined,"other ray FD became undefined");
   const double error=((plus.point-minus.point)/(2*h)-base.jacobian.col(j)).norm();
   const double scale=base.jacobian.col(j).norm();
   require(scale<1e-8 ? error<=1e-10 : error/scale<=1e-6,
           "negative/zero b ray differential mismatch");
  }
 }
}
void checkRankTwoIdealEvent() {
 // A finite inverse-Gauss point can have a full coefficient derivative while
 // the ideal-edge/ray chart does not. Guard against promoting the former as
 // an everywhere-smooth source derivative or dropping one merging root.
 const Eigen::Vector3d n(-0.6937804775604491,0.7020464447761631,0.16062203564002314);
 const Eigen::Vector3d opposite(n.x(),n.y(),-n.z());
 const double g=n.x()*n.x()+n.y()*n.y()/2;
 const double x=n.x()/(2*n.z()),y=n.y()/(4*n.z());
 const Eigen::Vector3d origin(x,y,-.125-x*x-2*y*y);
 const Eigen::Vector3d limit(x,y,-std::sqrt(256-x*x-y*y));
 const Eigen::Vector3d transverse(n.x()/std::sqrt(g),n.y()/(2*std::sqrt(g)),0);
 Eigen::Vector3d leading=transverse;
 leading.z()=-(limit.x()*transverse.x()+limit.y()*transverse.y())/limit.z();
 leading*=origin.z()-limit.z();
 require(leading.norm()>1,"rank-transition witness has zero leading term");
 double previousError=std::numeric_limits<double>::infinity();
 double previousQuotient=0.;
 for(double delta:{1e-2,1e-4,1e-6,1e-8,1e-10}) {
  Q::Vector p;p<<.125,0,0,1,1,2,-delta,0,0,0;
  const spin::QuadricInverseGaussMap map(Q::fromUnnormalized(p));
  const auto finite=map.algebraicPoint(n);
  const auto ideal=map.idealEdgeDifferential(n,opposite);
  require(finite.valid() && ideal.defined && ideal.points.size()==2,
          "rank-transition two ideal roots or finite anchor lost");
  const double offset=std::sqrt(delta*g)/(2*n.z());
  require(std::abs(ideal.points[0].parameter-(.5-offset))<1e-14 &&
          std::abs(ideal.points[1].parameter-(.5+offset))<1e-14,
          "rank-transition roots disagree with independent quadratic");
  const auto ray=spin::truncateQuadricProjectiveRay(*finite.point,ideal.points[0].direction,16.);
  require(ray.defined,"rank-transition near-sheet ray unavailable");
  const Eigen::Vector3d quotient=(ray.point-limit)/std::sqrt(delta);
  const double error=(quotient-leading).norm();
  const double derivativeQuotient=(ray.point-limit).norm()/delta;
  require(error<previousError && derivativeQuotient>previousQuotient,
          "rank-transition square-root asymptotic not recovered");
  previousError=error;previousQuotient=derivativeQuotient;
  std::cout<<"rank_event delta="<<delta<<" sqrt_error="<<error
           <<" difference_quotient="<<derivativeQuotient<<'\n';
  p(6)=delta;
  const spin::QuadricInverseGaussMap positive(Q::fromUnnormalized(p));
  const auto noRoots=positive.idealEdgeDifferential(n,opposite);
  require(noRoots.defined && noRoots.points.empty(),"positive-side elliptic edge has spurious ideal roots");
 }
 require(previousError<4e-4 && previousQuotient>7e5,
         "rank-transition nonzero square-root witness lost");
 Q::Vector p;p<<.125,0,0,1,1,2,0,0,0,0;
 const spin::QuadricInverseGaussMap atEvent(Q::fromUnnormalized(p));
 const auto finite=atEvent.algebraicCoefficientDifferential(n);
 const auto ideal=atEvent.idealEdgeDifferential(n,opposite);
 require(finite.valid(),"regular rank-two finite derivative lost");
 require(!ideal.defined && ideal.points.empty(),
         "non-C1 ideal event falsely given a regular differential");
 std::cout<<"[PASS] rank-two finite derivative distinguished from non-C1 ideal/ray event\n";
}
}
int main(){try{checkIdealRoots();checkComposedDerivative();checkRay();checkOtherRayBranches();checkRankTwoIdealEvent();
 std::cout<<"[PASS] paper ideal-edge/ray construction and full coefficient derivative\n";return 0;
}catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<'\n';return 1;}}
