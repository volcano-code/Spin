#include "spin/TruncationRootCertificate.h"
#include "spin/ExactRadicalSign.h"
#include "spin/FirstMonotoneTrue.h"
#include "spin/ExactSturmTarski.h"
#include "verification/ExactRationalCoordinates.h"
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
void checkPositiveBranchSigns() {
 using R=CGAL::Gmpq;
 const auto check=[](const R& u,const R& v,const R& s,int expected) {
  if(spin::detail::positiveRadicalSumSign(u,v,s)!=expected)
   throw std::runtime_error("exact positive-root sum lost its physical sign");
 };
 check(1,1,1,1); // A zero squared polynomial is NOT a physical root.
 check(-1,-1,1,-1);
 check(1,-1,1,0);check(-1,1,1,0);
 check(1,-1,4,-1);check(-1,1,4,1);
 check(1,-1,R(1,4),1);check(-1,1,R(1,4),-1);
 check(0,1,2,1);check(0,-1,2,-1);check(0,0,2,0);
 check(1,0,2,1);check(-1,0,2,-1);
 // Tiny exact separation: no binary64 rounding may decide the side.
 CGAL::Gmpz denominator(1);
 mpz_mul_2exp(denominator.mpz(),denominator.mpz(),1100);
 const R delta(CGAL::Gmpz(1),denominator);
 check(1,-1,1+delta,-1);check(1,-1,1-delta,1);
 for(const R& s:std::array<R,2>{R(0),R(-1)}) {
  bool rejected=false;
  try{(void)spin::detail::positiveRadicalSumSign(1,-1,s);}
  catch(const std::invalid_argument&){rejected=true;}
  if(!rejected)throw std::runtime_error("nonfinite/ideal radical branch accepted");
 }
}
void checkFirstMonotoneBracket() {
 // Selecting the wrong bracket would change the authority coordinates,
 // even though a later bracket would also meet the precision threshold.
 for(std::size_t first=0;first<=16;++first) {
  int calls=0;
  const auto index=spin::detail::firstMonotoneTrue(16,[&](std::size_t i) {
   ++calls;return i>=first;
  });
  if(index!=first || calls>5)
   throw std::runtime_error("monotone enclosure search skipped the first passing bracket");
 }
 if(spin::detail::firstMonotoneTrue(0,[](std::size_t){return true;})!=0)
  throw std::runtime_error("empty bracket search produced an index");
}
void checkPreparedBranchQueries() {
 using R=CGAL::Gmpq;
 using P=spin::detail::ExactRationalPolynomial;
 // Roots 1/4 and 3/4; expectations are their literal signs, not an
 // invocation of the uncached implementation under test.
 spin::detail::PreparedSturmTarskiCounter query(P{R(3,16),R(-1),R(1)});
 for(int repeat=0;repeat<2;++repeat) {
  if(query.signSum(P{R(-1,2),R(1)},R(0),R(1))!=0 ||
     query.signSum(P{R(-1,2),R(1)},R(0),R(1,2))!=-1 ||
     query.signSum(P{R(-1,2),R(1)},R(1,2),R(1))!=1 ||
     query.signSum(P{R(3,16),R(-1),R(1)},R(0),R(1))!=0 ||
     query.signSum(P{R(0),R(0),R(0),R(0),R(0),R(1)},R(0),R(1))!=2)
   throw std::runtime_error("prepared query reused stale interval/value signs");
 }
 spin::detail::PreparedSturmTarskiCounter other(P{R(-2),R(0),R(1)});
 if(other.signSum(P{R(0),R(1)},R(-2),R(0))!=-1 ||
    other.signSum(P{R(0),R(1)},R(0),R(2))!=1)
  throw std::runtime_error("prepared query reused another polynomial's roots");
 bool rejected=false;
 try{(void)query.signSum(P{R(1)},R(1,4),R(1));}
 catch(const std::invalid_argument&){rejected=true;}
 if(!rejected)throw std::runtime_error("prepared query accepted a root endpoint");
}
void checkRegularParaboloidPartition() {
 using Q=spin::QuadricCoefficients10;
 using R=CGAL::Gmpq;
 const auto profile=spin::TruncationRootProfile::RegularParaboloid;
 Q::Vector p;p<<0,0,0,-1,1,1,0,0,0,0;
 const auto q=Q::fromUnnormalized(p);
 const Eigen::Vector3d a(1,0,-1),b(1,0,1);
 const auto check=[&](const Q& quadric,double radius,const Eigen::Vector3d& n0,
                      const Eigen::Vector3d& n1,const std::vector<bool>& expected) {
  const spin::TruncationRootCertifier certifier(quadric,radius,profile);
  const auto r=certifier.partition(n0,n1),rev=certifier.partition(n1,n0);
  if(!r.certified || r.inside!=expected || r.roots.size()+1!=expected.size())
   throw std::runtime_error("rank-two partition unavailable/incorrect: "+r.diagnostic);
  if(!rev.certified || rev.inside!=std::vector<bool>(expected.rbegin(),expected.rend()) ||
     rev.roots.size()!=r.roots.size())throw std::runtime_error("rank-two reversal mismatch");
  for(size_t i=0;i<r.roots.size();++i) {
   const auto& x=r.roots[i];const auto& y=rev.roots[r.roots.size()-1-i];
   if(!(x.lower<1-y.lower && 1-y.upper<x.upper))
    throw std::runtime_error("rank-two reversed algebraic root identity mismatch");
  }
  return r;
 };
 const auto crossing=check(q,15./16,a,b,{true,false});
 // x=3/4,z=9/16 has exact radius 15/16 and t=1/6.
 if(!(crossing.roots[0].lower<R(1,6) && R(1,6)<crossing.roots[0].upper))
  throw std::runtime_error("rank-two interval misses analytic crossing");
 check(q,15./16,a,Eigen::Vector3d(1,0,0),{true,false}); // ideal endpoint
 check(q,2,Eigen::Vector3d(1,0,1),Eigen::Vector3d(2,0,1),{false});
 // Oblique exact null axis: q=(x-y)^2+2z^2-x-y. At t=2/3,
 // n(t)=(-1+2t,-1,0) maps to (1,0,0), on the unit sphere.
 p<<0,-1,-1,0,1,1,2,0,0,-1;
 const auto oblique=check(Q::fromUnnormalized(p),1,
   Eigen::Vector3d(-1,-1,0),Eigen::Vector3d(1,-1,0),{true,false});
 if(!(oblique.roots[0].lower<R(2,3) && R(2,3)<oblique.roots[0].upper))
  throw std::runtime_error("oblique null-axis crossing missed");
 // Saddle cancellation must remove the false repeated radius root at ideal t=1/2.
 p<<0,0,0,-1,1,-1,0,0,0,0;
 check(Q::fromUnnormalized(p),1,Eigen::Vector3d(1,1,-1),Eigen::Vector3d(1,1,1),{true,false});
 // One finite normal edge enters and exits the sphere before its ideal end.
 p<<-2,0,0,-1,1,1,0,0,0,0;
 check(Q::fromUnnormalized(p),1.5,Eigen::Vector3d(0,0,-1),Eigen::Vector3d(1,0,0),{false,true,false});
 check(Q::fromUnnormalized(p),1.5,Eigen::Vector3d(-2.5,0,-1),Eigen::Vector3d(4.5,0,1),
       {true,false,true,false}); // retain exit, reentry, final exit
 const auto reject=[](const spin::TruncationEdgePartition& r) {
  if(r.certified || !r.roots.empty() || !r.inside.empty() || r.diagnostic.empty())
   throw std::runtime_error("rank-two unknown leaked success/partial partition");
 };
 // On x=3/4, radius 15/16 is tangent when y=0; not a regular crossing.
 reject(spin::TruncationRootCertifier(q,15./16,profile).partition(
   Eigen::Vector3d(1.5,-1,-1),Eigen::Vector3d(1.5,1,-1)));
 reject(spin::TruncationRootCertifier(q,15./16,profile).partition(
   Eigen::Vector3d(1.5,0,-1),Eigen::Vector3d(2,0,-1))); // radius endpoint
 reject(spin::TruncationRootCertifier(q,1,profile).partition(
   Eigen::Vector3d(1,0,0),Eigen::Vector3d(0,1,0))); // wholly ideal
 reject(spin::TruncationRootCertifier(q,1).partition(a,b)); // legacy opt-in fence
 p<<-1,0,0,0,1,1,1,0,0,0;
 reject(spin::TruncationRootCertifier(Q::fromUnnormalized(p),2,profile).partition(a,b));
 p<<-1,0,0,0,1,1,0,0,0,0;
 reject(spin::TruncationRootCertifier(Q::fromUnnormalized(p),2,profile).partition(a,b));
 std::cout<<"[PASS] explicit rank-two exact edge partition; double/simple poles, ideal endpoints, multiple roots and guards\n";
}
void checkEdgePartition() {
 const auto check=[](const spin::QuadricCoefficients10& q,double radius,
                    const Eigen::Vector3d& a,const Eigen::Vector3d& b,
                    const std::vector<bool>& expected) {
  const spin::TruncationRootCertifier certifier(q,radius);
  const auto actual=certifier.partition(a,b);
  if(!actual.certified || actual.inside!=expected ||
     actual.roots.size()+1!=expected.size())
   throw std::runtime_error("complete edge membership partition is wrong: "+actual.diagnostic);
  const auto reverse=certifier.partition(b,a);
  if(!reverse.certified || reverse.inside!=std::vector<bool>(expected.rbegin(),expected.rend()) ||
     reverse.roots.size()!=actual.roots.size())
   throw std::runtime_error("edge reversal changed membership or root count");
  for(std::size_t i=0;i<actual.roots.size();++i) {
   const auto& f=actual.roots[i];
   const auto& r=reverse.roots[reverse.roots.size()-1-i];
   // Isolation need not choose identical brackets after reversal, but the
   // corresponding algebraic root's brackets must overlap strictly.
   if(!(f.lower<1-r.lower && 1-r.upper<f.upper))
    throw std::runtime_error("reversed root ordinal does not identify the same root");
  }
 };
 spin::QuadricCoefficients10::Vector p;
 p<<1,0,0,0,1,1,-1,0,0,0;
 auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
 const Eigen::Vector3d a(-1,0,1),b(1,0,1);
 check(q,2,a,b,{false,true,false}); // both endpoints ideal/outside
 check(q,.5,a,b,{false});
 p<<0,-2,0,0,1,1,1,0,0,0; // sphere centred at (1,0,0)
 q=spin::QuadricCoefficients10::fromUnnormalized(p);
 check(q,1.8,Eigen::Vector3d(1,2,0),Eigen::Vector3d(1,-2,0),{true,false,true});
 check(q,3,Eigen::Vector3d(1,2,0),Eigen::Vector3d(1,-2,0),{true});
 p<<-1,0,0,0,1,1,-1,0,0,0;
 q=spin::QuadricCoefficients10::fromUnnormalized(p);
 check(q,2,Eigen::Vector3d::UnitX(),Eigen::Vector3d::UnitZ(),{true,false});
 // A mixed-endpoint edge can cross three times, not just once.
 Eigen::Matrix3d A;A<<1.6530585186285478,-1.4703811461027434,.023949489618097675,
 -1.4703811461027434,.44900212325215494,.6030355405690727,
 .023949489618097675,.6030355405690727,2.1869169579519263;
 const Eigen::Vector3d c(3.049131259992132,3.41662257009421,2.3750832053948923);
 p(0)=c.dot(A*c)+1.8865497742300783;p.segment<3>(1)=-2*A*c;
 p.segment<3>(4)=A.diagonal();p(7)=A(1,2);p(8)=A(0,2);p(9)=A(0,1);
 const double phi=(1+std::sqrt(5.))/2;
 const auto n0=Eigen::Vector3d(-1,phi,0).normalized().eval();
 const auto n1=Eigen::Vector3d(0,1,-phi).normalized().eval();
 q=spin::QuadricCoefficients10::fromUnnormalized(p);
 check(q,3.006066543711902,n0,n1,{false,true,false,true});
 const auto reject=[](const spin::TruncationEdgePartition& result) {
  if(result.certified || !result.roots.empty() || !result.inside.empty() || result.diagnostic.empty())
   throw std::runtime_error("unknown partition leaked partial success");
 };
 reject(spin::TruncationRootCertifier(q,3.006066543711902).partition(n0,n1,0));
 p<<1,0,0,0,1,1,-1,0,0,0;q=spin::QuadricCoefficients10::fromUnnormalized(p);
 reject(spin::TruncationRootCertifier(q,1).partition(a,b)); // tangent
 p<<-1,0,0,0,1,1,-1,0,0,0;q=spin::QuadricCoefficients10::fromUnnormalized(p);
 reject(spin::TruncationRootCertifier(q,1).partition(Eigen::Vector3d::UnitX(),Eigen::Vector3d::UnitZ()));
 std::cout<<"[PASS] complete exact membership for zero/one/two/three roots and reversal\n";
}
void checkMultipleCrossingIsolation() {
 spin::QuadricCoefficients10::Vector p;p<<1,0,0,0,1,1,-1,0,0,0;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
 const Eigen::Vector3d a(-1,0,1),b(1,0,1);
 const auto roots=spin::TruncationRootCertifier(q,2).isolate(a,b);
 if(!roots.certified || roots.intervals.size()!=2)
  throw std::runtime_error("two interior crossings between inadmissible endpoints were lost: "+roots.diagnostic);
 const double expected[2]={(1-std::sqrt(.6))/2,(1+std::sqrt(.6))/2};
 for(int i=0;i<2;++i) {
  const auto& interval=roots.intervals[i];
  if(!(interval.lower<CGAL::Gmpq(expected[i]) && CGAL::Gmpq(expected[i])<interval.upper))
   throw std::runtime_error("isolating interval does not enclose the analytic radius crossing");
 }
 if(roots.intervals[0].upper>roots.intervals[1].lower)
  throw std::runtime_error("multiple root intervals overlap");
 const auto exhausted=spin::TruncationRootCertifier(q,2).isolate(a,b,0);
 if(exhausted.certified || !exhausted.intervals.empty())
  throw std::runtime_error("exhausted isolation budget returned partial roots as complete");
 std::string why;
 if(spin::TruncationRootCertifier(q,2).certify(a,b,&why))
  throw std::runtime_error("legacy single-root certificate accepted two roots");
 const auto empty=spin::TruncationRootCertifier(q,.5).isolate(a,b);
 if(!empty.certified || !empty.intervals.empty())
  throw std::runtime_error("certified empty edge confused with unresolved isolation");
 const auto tangent=spin::TruncationRootCertifier(q,1).isolate(a,b);
 if(tangent.certified || !tangent.intervals.empty())
  throw std::runtime_error("tangent root granted a regular crossing certificate");
 const auto reversed=spin::TruncationRootCertifier(q,2).isolate(b,a);
 if(!reversed.certified || reversed.intervals.size()!=2)
  throw std::runtime_error("edge reversal lost regular crossing roots");
 // Captured normal edge 13--14 from the failed fixed-state level-1 map.
 p<<.4091598840032713,-.36688397145905616,-.019533580749283127,
 -.3602578587235418,.14803212299816212,-.34789822850316837,
 -.25957403943858653,-.3788370952579597,.3838860742170953,.2580806535816809;
 const auto captured=spin::TruncationRootCertifier(
  spin::QuadricCoefficients10::fromUnnormalized(p),443.40500673763256).isolate(
   Eigen::Vector3d(-.5000000000000001,.30901699437494745,.8090169943749475),
   Eigen::Vector3d(-.30901699437494745,.8090169943749475,.5000000000000001));
 if(!captured.certified || captured.intervals.size()!=2)
  throw std::runtime_error("captured double crossing was not isolated: "+captured.diagnostic);
 std::cout<<"[PASS] complete multiple-root isolation, empty edge, tangent rejection and captured crossing\n";
}
}
int main(int argc,char** argv){try{
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
 checkPositiveBranchSigns();
 checkFirstMonotoneBracket();
 checkPreparedBranchQueries();
 if(argc==2 && std::string(argv[1])=="--near-cone-rational-point") {
  // Actual event 17, normal edge 7--67, from the frozen level-2 Cuboid.
  // Squaring the radius equation creates a nearby inadmissible branch root;
  // rejecting that root must not lose the uniquely certified physical root.
  spin::QuadricCoefficients10::Vector p;
  p<<-1.0815175509191138e-12,3.074304257986158e-8,
     5.4940100928584925e-9,3.89182874524446e-8,
     .12168625877604855,.263115530927828,-.9276754486009523,
     1.5085428667989776e-8,-.23533154032699224,3.9530246828365786e-11;
  const spin::TruncationRootCertifier certifier(
   spin::QuadricCoefficients10::fromNormalized(p),5.116961229773625);
  const Eigen::Vector3d a(0,.52573111211913359,-.85065080835203999);
  const Eigen::Vector3d b(-.25989191300775444,.43388856455269481,-.86266848041618616);
  std::string why;
  if(!certifier.certify(a,b,&why))
   throw std::runtime_error("captured root lacks a unique admissible certificate: "+why);
  const auto isolated=certifier.isolate(a,b);
  std::cout<<"[DIAGNOSTIC] unique admissible root certified; legacy isolation: "
           <<isolated.certified<<" "<<isolated.diagnostic<<std::endl;
  const auto point=certifier.rationalPoint(a,b);
  const auto compact=certifier.rationalPoint(a,b,true);
  if(!compact.defined || compact.lower!=point.lower || compact.upper!=point.upper)
   throw std::runtime_error("compact captured root changed the certified analytic interval");
  CGAL::Gmpz compactDenominator(1);
  mpz_mul_2exp(compactDenominator.mpz(),compactDenominator.mpz(),1024);
  for(int axis=0;axis<3;++axis) {
   const auto x=verification::parseExactRationalCoordinate(compact.coordinates[axis]);
   if(x<verification::parseExactRationalCoordinate(point.lower[axis]) ||
      x>verification::parseExactRationalCoordinate(point.upper[axis]) ||
      compactDenominator%x.denominator()!=0)
    throw std::runtime_error("compact truncation coordinate escapes original enclosure/grid");
  }
  if(!point.defined)
   throw std::runtime_error("captured unique physical crossing lost: "+point.diagnostic);
  using R=CGAL::Gmpq;
  CGAL::Gmpz denominator(1);mpz_mul_2exp(denominator.mpz(),denominator.mpz(),384);
  const R tolerance(CGAL::Gmpz(1),denominator);
  const auto reverse=certifier.rationalPoint(b,a);
  if(!reverse.defined)throw std::runtime_error("reversed captured crossing unavailable: "+reverse.diagnostic);
  R radiusLower=0,radiusUpper=0;
  for(int axis=0;axis<3;++axis) {
   const auto lo=verification::parseExactRationalCoordinate(point.lower[axis]);
   const auto hi=verification::parseExactRationalCoordinate(point.upper[axis]);
   const auto value=verification::parseExactRationalCoordinate(point.coordinates[axis]);
   const auto rlo=verification::parseExactRationalCoordinate(reverse.lower[axis]);
   const auto rhi=verification::parseExactRationalCoordinate(reverse.upper[axis]);
   if(lo>value || value>hi || hi-lo>tolerance || rhi-rlo>tolerance || hi<rlo || rhi<lo)
    throw std::runtime_error("captured crossing enclosure/reversal invalid");
   radiusLower+=(lo<=0 && hi>=0)?R(0):std::min(lo*lo,hi*hi);
   radiusUpper+=std::max(lo*lo,hi*hi);
  }
  const R radiusSquared=R(5.116961229773625)*R(5.116961229773625);
  if(radiusLower>radiusSquared || radiusUpper<radiusSquared)
   throw std::runtime_error("captured enclosure does not contain the actual truncation sphere");
  std::cout<<"[PASS] captured near-cone physical crossing, unchanged radius and <=2^-384 enclosures\n";
  return 0;
 }
 if(argc==2 && std::string(argv[1])=="--rational-point") {
  using R=CGAL::Gmpq;
  spin::QuadricCoefficients10::Vector p;p<<1,0,0,0,1,1,-1,0,0,0;
  const spin::TruncationRootCertifier certifier(spin::QuadricCoefficients10::fromUnnormalized(p),1.5);
  const auto point=certifier.rationalPoint(Eigen::Vector3d(0,0,-1),Eigen::Vector3d(1,0,-1));
  if(!point.defined)throw std::runtime_error("rational truncation point unavailable: "+point.diagnostic);
  CGAL::Gmpz denominator(1);mpz_mul_2exp(denominator.mpz(),denominator.mpz(),384);
  const R tolerance(CGAL::Gmpz(1),denominator);
  for(int axis=0;axis<3;++axis) {
   const auto lo=verification::parseExactRationalCoordinate(point.lower[axis]);
   const auto hi=verification::parseExactRationalCoordinate(point.upper[axis]);
   const auto value=verification::parseExactRationalCoordinate(point.coordinates[axis]);
   if(lo>value || value>hi || hi-lo>tolerance)
    throw std::runtime_error("invalid or too wide rational truncation coordinate enclosure");
   const R squared=axis==0?R(5,8):(axis==2?R(13,8):R(0));
   if(lo<0 || lo*lo>squared || hi*hi<squared)
    throw std::runtime_error("rational truncation enclosure misses analytic crossing");
  }
  const auto repeated=certifier.rationalPoint(Eigen::Vector3d(0,0,-1),Eigen::Vector3d(1,0,-1));
  if(point.coordinates!=repeated.coordinates || point.lower!=repeated.lower || point.upper!=repeated.upper)
   throw std::runtime_error("rational root enclosure is nondeterministic");
  if(certifier.rationalPoint(Eigen::Vector3d(0,0,-1),Eigen::Vector3d(.01,0,-1)).defined)
   throw std::runtime_error("missing regular unique crossing granted a rational point");
  if(spin::TruncationRootCertifier(spin::QuadricCoefficients10::fromUnnormalized(p),2)
      .rationalPoint(Eigen::Vector3d(-1,0,1),Eigen::Vector3d(1,0,1)).defined)
   throw std::runtime_error("two physical roots granted a unique rational point");
  if(spin::TruncationRootCertifier(spin::QuadricCoefficients10::fromUnnormalized(p),1)
      .rationalPoint(Eigen::Vector3d(-1,0,1),Eigen::Vector3d(1,0,1)).defined)
   throw std::runtime_error("tangent/repeated root granted a regular rational point");
  p<<-1,0,0,0,1,1,-1,0,0,0;
  if(spin::TruncationRootCertifier(spin::QuadricCoefficients10::fromUnnormalized(p),1)
      .rationalPoint(Eigen::Vector3d::UnitX(),Eigen::Vector3d::UnitZ()).defined)
   throw std::runtime_error("endpoint root granted an interior rational point");
  std::cout<<"[PASS] analytic hyperboloid truncation point rational enclosure <=2^-384\n";
  return 0;
 }
#endif
 checkRegularParaboloidPartition();
 checkEdgePartition();
 checkMultipleCrossingIsolation();
 spin::QuadricCoefficients10::Vector p;p<<-1,0,0,0,1,1,-1,0,0,0;
 std::string why;
 if(!spin::certifySingleTruncationRoot(spin::QuadricCoefficients10::fromUnnormalized(p),
  Eigen::Vector3d::UnitX(),Eigen::Vector3d::UnitZ(),2,&why))throw std::runtime_error("unique analytic root not certified: "+why);
 Eigen::Matrix3d a;a<<1.6530585186285478,-1.4703811461027434,.023949489618097675,
 -1.4703811461027434,.44900212325215494,.6030355405690727,
 .023949489618097675,.6030355405690727,2.1869169579519263;
 Eigen::Vector3d c(3.049131259992132,3.41662257009421,2.3750832053948923);
 p(0)=c.dot(a*c)+1.8865497742300783;p.segment<3>(1)=-2*a*c;
 p.segment<3>(4)=a.diagonal();p(7)=a(1,2);p(8)=a(0,2);p(9)=a(0,1);
 const double phi=(1+std::sqrt(5.))/2;
 const Eigen::Vector3d n0=Eigen::Vector3d(-1,phi,0).normalized(),n1=Eigen::Vector3d(0,1,-phi).normalized();
 if(spin::certifySingleTruncationRoot(spin::QuadricCoefficients10::fromUnnormalized(p),n0,n1,3.006066543711902,&why))
  throw std::runtime_error("three admissible roots falsely certified unique");
 if(why!="multiple admissible positive-branch roots")
  throw std::runtime_error("multiple-root fixture rejected for wrong reason: "+why);
 std::cout<<"[PASS] exact single-root witness and multiple-root rejection: "<<why<<"\n";
 p<<-1,0,0,0,1,1,-1,0,0,0;
 if(spin::certifySingleTruncationRoot(spin::QuadricCoefficients10::fromUnnormalized(p),
  Eigen::Vector3d::UnitX(),Eigen::Vector3d::UnitZ(),1,&why))
  throw std::runtime_error("endpoint root certified interior");
 // Sphere centered at (1,0,0), tangent to radius 2 at the edge midpoint.
 p<<0,-2,0,0,1,1,1,0,0,0;
 if(spin::certifySingleTruncationRoot(spin::QuadricCoefficients10::fromUnnormalized(p),
  Eigen::Vector3d(1,-1,0),Eigen::Vector3d(1,1,0),2,&why))
  throw std::runtime_error("tangent/repeated root certified regular");
 std::cout<<"[PASS] endpoint and tangent/repeated roots rejected\n";return 0;
}catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<"\n";return 1;}}
