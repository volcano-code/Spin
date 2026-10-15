#include "spin/TruncationRootCertificate.h"
#include "spin/RationalSourceRepresentative.h"
#include "spin/ExactSturmTarski.h"
#include "spin/ExactRadicalSign.h"
#include "spin/FirstMonotoneTrue.h"
#include <CGAL/Gmpq.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <locale>
#include <stdexcept>
#include <vector>

namespace {
using Q=CGAL::Gmpq;
using P=std::vector<Q>;
void trim(P& p){while(!p.empty() && p.back()==0)p.pop_back();}
P add(P a,const P& b){a.resize(std::max(a.size(),b.size()),Q(0));for(size_t i=0;i<b.size();++i)a[i]+=b[i];trim(a);return a;}
P scale(P a,const Q& s){for(auto& v:a)v*=s;trim(a);return a;}
P mul(const P& a,const P& b){
 if(a.empty()||b.empty())return {};
 P p(a.size()+b.size()-1,Q(0));
 for(size_t i=0;i<a.size();++i)for(size_t j=0;j<b.size();++j)p[i+j]+=a[i]*b[j];
 trim(p);return p;
}
P derivative(const P& p){P d;for(size_t i=1;i<p.size();++i)d.push_back(p[i]*Q(static_cast<int>(i)));trim(d);return d;}
Q eval(const P& p,const Q& x){Q y(0);for(auto i=p.rbegin();i!=p.rend();++i)y=y*x+*i;return y;}
P remainder(P a,const P& b){
 while(!a.empty() && a.size()>=b.size()){
  const size_t offset=a.size()-b.size();const Q factor=a.back()/b.back();
  for(size_t j=0;j<b.size();++j)a[offset+j]-=factor*b[j];trim(a);
 }return a;
}
P exactQuotient(P a,const P& b){
 if(b.empty())throw std::invalid_argument("zero polynomial divisor");
 P quotient(a.size()>=b.size()?a.size()-b.size()+1:0,Q(0));
 while(!a.empty() && a.size()>=b.size()){
  const size_t offset=a.size()-b.size();const Q factor=a.back()/b.back();
  quotient[offset]+=factor;
  for(size_t j=0;j<b.size();++j)a[offset+j]-=factor*b[j];trim(a);
 }
 if(!a.empty())throw std::logic_error("nonexact rank-two pole cancellation");
 trim(quotient);return quotient;
}
int sign(const Q& x){return x>0?1:(x<0?-1:0);}
int variations(const std::vector<P>& chain,const Q& x){
 int previous=0,count=0;
 for(const auto& p:chain){const int s=sign(eval(p,x));if(s){if(previous && previous!=s)++count;previous=s;}}
 return count;
}
std::vector<P> tarskiChain(const P& polynomial,const P& value){
 if(value.empty())return {};
 P weighted=mul(derivative(polynomial),value);
 if(weighted.empty())return {};
 std::vector<P> chain={polynomial,weighted};
 while(!chain.back().empty()){
  P rest=scale(remainder(chain[chain.size()-2],chain.back()),Q(-1));
  if(rest.empty())break;
  // Only positive rescaling is permitted: this is a signed remainder chain.
  Q magnitude=rest.back();if(magnitude<0)magnitude=-magnitude;
  chain.push_back(scale(rest,1/magnitude));
 }
 return chain;
}
int tarskiQueryUnchecked(const P& polynomial,const P& value,const Q& lo,const Q& hi){
 const auto chain=tarskiChain(polynomial,value);
 return variations(chain,lo)-variations(chain,hi);
}
int branchSignSum(const P& polynomial,const P& value,const Q& lo,const Q& hi,
                   spin::detail::PreparedSturmTarskiCounter* prepared) {
 return prepared?prepared->signSum(value,lo,hi):
     spin::detail::sturmTarskiSignSum(polynomial,value,lo,hi);
}
struct Interval {Q lo,hi;};
Interval sum(const Interval& a,const Interval& b){return {a.lo+b.lo,a.hi+b.hi};}
Interval subtract(const Interval& a,const Interval& b){return {a.lo-b.hi,a.hi-b.lo};}
Interval product(const Interval& a,const Interval& b){
 const std::array<Q,4> v={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
 return {*std::min_element(v.begin(),v.end()),*std::max_element(v.begin(),v.end())};
}
Interval quotient(const Interval& a,const Interval& b){
 if(b.lo<=0 && b.hi>=0)throw std::runtime_error("root derivative determinant enclosure contains zero");
 return product(a,{1/b.hi,1/b.lo});
}
using IntervalMatrix=std::array<std::array<Interval,3>,3>;
Interval determinant(const IntervalMatrix& a){
 return sum(subtract(
  product(a[0][0],subtract(product(a[1][1],a[2][2]),product(a[1][2],a[2][1]))),
  product(a[0][1],subtract(product(a[1][0],a[2][2]),product(a[1][2],a[2][0])))),
  product(a[0][2],subtract(product(a[1][0],a[2][1]),product(a[1][1],a[2][0]))));
}
bool containsZero(const Interval& a){return a.lo<=0 && a.hi>=0;}
Interval bound(const P& p,const Q& lo,const Q& hi){
 Interval result{Q(0),Q(0)};
 for(auto i=p.rbegin();i!=p.rend();++i){result=product(result,{lo,hi});result.lo+=*i;result.hi+=*i;}
 return result;
}
int countPositiveAtRoots(const P& polynomial,const P& value,
                         const Q& lo,const Q& hi,int roots,
                         spin::detail::PreparedSturmTarskiCounter* prepared){
 const int first=branchSignSum(polynomial,value,lo,hi,prepared);
 const int second=branchSignSum(polynomial,mul(value,value),lo,hi,prepared);
 const int numerator=first+second;
 if(numerator%2!=0)return -1;
 const int count=numerator/2;
 return count<0||count>roots?-1:count;
}
int countNegativeAtRoots(const P& polynomial,const P& value,
                         const Q& lo,const Q& hi,int roots,
                         spin::detail::PreparedSturmTarskiCounter* prepared){
 const int first=branchSignSum(polynomial,value,lo,hi,prepared);
 const int second=branchSignSum(polynomial,mul(value,value),lo,hi,prepared);
 const int numerator=second-first;
 if(numerator%2!=0)return -1;
 const int count=numerator/2;
 return count<0||count>roots?-1:count;
}
int countPositiveNegativeAtRoots(const P& polynomial,const P& positive,
                                 const P& negative,const Q& lo,const Q& hi,
                                 int roots,spin::detail::PreparedSturmTarskiCounter* prepared){
 const P positiveSquared=mul(positive,positive);
 const P negativeSquared=mul(negative,negative);
 const int q22=branchSignSum(polynomial,mul(positiveSquared,negativeSquared),lo,hi,prepared);
 const int q21=branchSignSum(polynomial,mul(positiveSquared,negative),lo,hi,prepared);
 const int q12=branchSignSum(polynomial,mul(positive,negativeSquared),lo,hi,prepared);
 const int q11=branchSignSum(polynomial,mul(positive,negative),lo,hi,prepared);
 const int numerator=q22-q21+q12-q11;
 if(numerator%4!=0)return -1;
 const int count=numerator/4;
 return count<0||count>roots?-1:count;
}
}
namespace spin::detail {
PreparedSturmTarskiCounter::PreparedSturmTarskiCounter(ExactRationalPolynomial polynomial)
 : polynomial_(std::move(polynomial)) {
 if(polynomial_.size()<2)
  throw std::invalid_argument("prepared Sturm-Tarski query requires a nonconstant polynomial");
}
int PreparedSturmTarskiCounter::signSum(const ExactRationalPolynomial& value,
                                       const Q& lo,const Q& hi) {
 if(lo>=hi || eval(polynomial_,lo)==0 || eval(polynomial_,hi)==0)
  throw std::invalid_argument("Sturm-Tarski query requires nonroot interval endpoints");
 auto found=chains_.find(value);
 if(found==chains_.end())found=chains_.emplace(value,tarskiChain(polynomial_,value)).first;
 return variations(found->second,lo)-variations(found->second,hi);
}
int sturmTarskiSignSum(const ExactRationalPolynomial& polynomial,
                       const ExactRationalPolynomial& value,
                       const CGAL::Gmpq& lo,const CGAL::Gmpq& hi) {
 if(polynomial.size()<2 || lo>=hi ||
    eval(polynomial,lo)==0 || eval(polynomial,hi)==0)
  throw std::invalid_argument("Sturm-Tarski query requires a nonconstant polynomial and nonroot interval endpoints");
 return tarskiQueryUnchecked(polynomial,value,lo,hi);
}
}
namespace spin {
namespace {
std::array<Q,10> representedCoefficients(const QuadricCoefficients10& quadric) {
 std::array<Q,10> out;
 for(int i=0;i<10;++i)out[i]=Q(quadric.values()(i));
 return out;
}
std::array<Q,10> authoritativeCoefficients(const QuadricCoefficientAuthority10& authority) {
 std::array<Q,10> out;
 for(int i=0;i<10;++i)out[i]=Q(authority.geometry()[i]);
 return out;
}
}
struct TruncationRootCertifier::Impl {
 const std::array<Q,10> coefficients;
 std::array<std::array<Q,3>,3> inverse{};
 std::array<std::array<Q,3>,3> quadratic{};
 std::array<Q,3> linear{},axis{};
 Q constant=0,axisLinear=0,radiusSquared=0;
 bool rankTwo=false;
 std::array<Q,3> center{};
 Q denominator=0;
 Q centerMinusRadius=0;
 bool valid=false;
 std::string diagnostic;

 Impl(const std::array<Q,10>& p,double radius,TruncationRootProfile profile) : coefficients(p) {
  if(!std::isfinite(radius)||radius<=0) {
   diagnostic="invalid normal edge or radius";
   return;
  }
  radiusSquared=Q(radius)*Q(radius);
  Q a[3][3]={{p[4],p[9],p[8]},{p[9],p[5],p[7]},{p[8],p[7],p[6]}};
  const Q determinant=a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
   -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
  if(profile==TruncationRootProfile::RegularParaboloid) {
   if(determinant!=0){diagnostic="regular paraboloid certificate requires exact rank two";return;}
   bool found=false;
   for(int column=0;column<3 && !found;++column) {
    for(int row=0;row<3;++row)axis[row]=
     a[(column+1)%3][(row+1)%3]*a[(column+2)%3][(row+2)%3]-
     a[(column+1)%3][(row+2)%3]*a[(column+2)%3][(row+1)%3];
    found=std::any_of(axis.begin(),axis.end(),[](const Q& x){return x!=0;});
   }
   if(!found){diagnostic="regular paraboloid certificate rejects rank below two";return;}
   Q pivot=0;for(const auto& x:axis)if(x!=0){pivot=x;break;}
   for(auto& x:axis)x/=pivot;
   for(int i=0;i<3;++i) {
    linear[i]=p[1+i];axisLinear+=axis[i]*linear[i];
    for(int j=0;j<3;++j)quadratic[i][j]=a[i][j];
   }
   if(axisLinear==0){diagnostic="rank-two certificate requires a nonzero linear null-axis term";return;}
   Q bordered[3][3];for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    bordered[i][j]=a[i][j]+axis[i]*axis[j];
   const Q bd=bordered[0][0]*(bordered[1][1]*bordered[2][2]-bordered[1][2]*bordered[2][1])
    -bordered[0][1]*(bordered[1][0]*bordered[2][2]-bordered[1][2]*bordered[2][0])
    +bordered[0][2]*(bordered[1][0]*bordered[2][1]-bordered[1][1]*bordered[2][0]);
   if(bd==0){diagnostic="singular regular-paraboloid nullspace solve";return;}
   for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    inverse[i][j]=(bordered[(j+1)%3][(i+1)%3]*bordered[(j+2)%3][(i+2)%3]-
                   bordered[(j+1)%3][(i+2)%3]*bordered[(j+2)%3][(i+1)%3])/bd;
   constant=p[0];radiusSquared=Q(radius)*Q(radius);
   rankTwo=true;valid=true;return;
  }
  if(determinant==0) {
   diagnostic="singular rational quadratic matrix";
   return;
  }
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)
   inverse[i][j]=(a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]-
                  a[(j+1)%3][(i+2)%3]*a[(j+2)%3][(i+1)%3])/determinant;
  const std::array<Q,3> linear={p[1],p[2],p[3]};
  denominator=-4*p[0];
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
   center[i]-=inverse[i][j]*linear[j]/2;
   denominator+=linear[i]*inverse[i][j]*linear[j];
  }
  if(denominator==0) {
   diagnostic="zero rational inverse-Gauss scale";
   return;
  }
  centerMinusRadius=-Q(radius)*Q(radius);
  for(const auto& coordinate:center)centerMinusRadius+=coordinate*coordinate;
  valid=true;
 }
 bool prepareEdge(const Eigen::Vector3d& n0,const Eigen::Vector3d& n1,
                  P& polynomial,P& s,P& u,P& v,std::vector<P>& chain,
                  std::string& why) const;
};

TruncationRootCertifier::TruncationRootCertifier(
 const QuadricCoefficients10& quadric,double radius)
 : TruncationRootCertifier(quadric,radius,TruncationRootProfile::FullRank) {}

TruncationRootCertifier::TruncationRootCertifier(
 const QuadricCoefficients10& quadric,double radius,TruncationRootProfile profile)
 : impl_(std::make_shared<Impl>(representedCoefficients(quadric),radius,profile)) {}

TruncationRootCertifier::TruncationRootCertifier(
 const QuadricCoefficientAuthority10& authority,double radius)
 : impl_(std::make_shared<Impl>(authoritativeCoefficients(authority),radius,
                             TruncationRootProfile::FullRank)) {}

bool TruncationRootCertifier::Impl::prepareEdge(
 const Eigen::Vector3d& n0,const Eigen::Vector3d& n1,
 P& polynomial,P& s,P& u,P& v,std::vector<P>& chain,std::string& why) const {
 const auto fail=[&](const std::string& message){why=message;return false;};
 if(!valid)return fail(diagnostic);
 if(!n0.allFinite()||!n1.allFinite()||
    n0.norm()<1e-12||n1.norm()<1e-12||n0.cross(n1).norm()<1e-12)
  return fail("invalid normal edge or radius");
 std::array<P,3> n,bn;
 for(int i=0;i<3;++i)n[i]={Q(n0(i)),Q(n1(i))-Q(n0(i))};
 if(rankTwo) {
  P w;for(int i=0;i<3;++i)w=add(w,scale(n[i],axis[i]));
  if(w.empty())return fail("normal edge is entirely ideal for the regular paraboloid");
  // y=Y/w, Y=B*(ub*n-b*w)/2; x=N/w² after imposing q(x)=0.
  std::array<P,3> y,numerator;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)
   y[i]=add(y[i],scale(add(scale(n[j],axisLinear),scale(w,-linear[j])),inverse[i][j]/2));
  P potential=scale(mul(w,w),constant);
  for(int i=0;i<3;++i) {
   potential=add(potential,scale(mul(y[i],w),linear[i]));
   for(int j=0;j<3;++j)potential=add(potential,scale(mul(y[i],y[j]),quadratic[i][j]));
  }
  for(int i=0;i<3;++i)numerator[i]=add(mul(y[i],w),scale(potential,-axis[i]/axisLinear));
  P coordinateDenominator=mul(w,w);
  // Saddle asymptotic normals can have a simple pole. Remove only exact
  // factors common to ALL coordinates before forming the radius polynomial;
  // otherwise w=0 becomes a spurious repeated sphere root.
  for(int cancel=0;w.size()==2 && cancel<2;++cancel) {
   bool common=true;for(const auto& p:numerator)common&=remainder(p,w).empty();
   if(!common)break;
   for(auto& p:numerator)p=exactQuotient(p,w);
   coordinateDenominator=exactQuotient(coordinateDenominator,w);
  }
  if(coordinateDenominator.size()==1 && w.size()==2)
   return fail("regular paraboloid normal edge unexpectedly has no ideal pole");
  s=scale(w,axisLinear); // positive-gradient admissibility, not a squared norm
  for(const auto& p:numerator)polynomial=add(polynomial,mul(p,p));
  polynomial=add(polynomial,scale(mul(coordinateDenominator,coordinateDenominator),-radiusSquared));
  u=polynomial; // denominator²>0 on every finite admissible region
 } else {
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)
   bn[i]=add(bn[i],scale(n[j],inverse[i][j]));
  for(int i=0;i<3;++i){s=add(s,scale(mul(n[i],bn[i]),4/denominator));u=add(u,mul(bn[i],bn[i]));v=add(v,scale(bn[i],2*center[i]));}
  u=add(u,scale(s,centerMinusRadius));
  const bool centered=v.empty();
  polynomial=centered?u:add(mul(u,u),scale(mul(mul(v,v),s),Q(-1)));
 }
 trim(polynomial);
 if(polynomial.size()<2)return fail("degenerate radius polynomial");
 if(eval(polynomial,Q(0))==0||eval(polynomial,Q(1))==0)return fail("radius polynomial has endpoint root");
 chain={polynomial,derivative(polynomial)};
 while(!chain.back().empty()) {
  P rest=scale(remainder(chain[chain.size()-2],chain.back()),Q(-1));
  if(rest.empty()) {
   if(chain.back().size()>1)return fail("radius polynomial has repeated root");
   break;
  }
  Q magnitude=rest.back();if(magnitude<0)magnitude=-magnitude;
  chain.push_back(scale(rest,1/magnitude));
 }
 return true;
}

namespace {
int admissibleRoots(const P& polynomial,const P& s,const P& u,const P& v,
                    const std::vector<P>& chain,const Q& lo,const Q& hi,
                    int& roots,spin::detail::PreparedSturmTarskiCounter* prepared=nullptr) {
 roots=variations(chain,lo)-variations(chain,hi);
 if(roots<=0)return roots;
 const auto scaleBounds=bound(s,lo,hi);
 if(scaleBounds.hi<=0)return 0;
 if(v.empty())return scaleBounds.lo>0?roots:countPositiveAtRoots(polynomial,s,lo,hi,roots,prepared);
 if(scaleBounds.lo>0)return countNegativeAtRoots(polynomial,mul(u,v),lo,hi,roots,prepared);
 return countPositiveNegativeAtRoots(polynomial,s,mul(u,v),lo,hi,roots,prepared);
}
}

bool TruncationRootCertifier::certify(
 const Eigen::Vector3d& n0,const Eigen::Vector3d& n1,std::string* why) const {
 const auto fail=[&](const std::string& message){if(why)*why=message;return false;};
 if(!impl_)return fail("missing exact certificate context");
 P polynomial,s,u,v;std::vector<P> chain;std::string diagnostic;
 if(!impl_->prepareEdge(n0,n1,polynomial,s,u,v,chain,diagnostic))return fail(diagnostic);
 const Q lo(0),hi(1);
 const int roots=variations(chain,lo)-variations(chain,hi);
 if(roots<0)return fail("invalid Sturm root count");
 if(roots==0)return fail("no uniquely admissible positive-branch root");
 int countedRoots=0;
 const int admissible=admissibleRoots(polynomial,s,u,v,chain,lo,hi,countedRoots);
 if(admissible<0)return fail("inconsistent exact Sturm-Tarski branch-sign count");
 if(admissible>1)return fail("multiple admissible positive-branch roots");
 if(admissible!=1)return fail("no uniquely admissible positive-branch root");
 if(why)*why="exact rational Sturm certificate: one simple admissible open-edge root";
 return true;
}

TruncationRootIsolation TruncationRootCertifier::isolate(
 const Eigen::Vector3d& n0,const Eigen::Vector3d& n1,int maximumDepth) const {
 const auto fail=[](const std::string& why){TruncationRootIsolation out;out.diagnostic=why;return out;};
 if(!impl_)return fail("missing exact certificate context");
 if(maximumDepth<0 || maximumDepth>64)return fail("invalid root isolation depth budget");
 P polynomial,s,u,v;std::vector<P> chain;std::string diagnostic;
 if(!impl_->prepareEdge(n0,n1,polynomial,s,u,v,chain,diagnostic))return fail(diagnostic);
 struct Node {Q lo,hi;int depth;};
 detail::PreparedSturmTarskiCounter prepared(polynomial);
 std::vector<Node> pending{{Q(0),Q(1),0}};
 TruncationRootIsolation result;
 while(!pending.empty()) {
  const Node node=pending.back();pending.pop_back();
  int roots=0;
  const int admissible=admissibleRoots(polynomial,s,u,v,chain,node.lo,node.hi,roots,&prepared);
  if(roots<0 || admissible<0 || admissible>roots)
   return fail("inconsistent exact Sturm-Tarski branch-sign count");
  if(admissible==0)continue;
  if(roots==1) {result.intervals.push_back({node.lo,node.hi});continue;}
  if(node.depth>=maximumDepth)return fail("root isolation depth budget exhausted");
  // A radius polynomial has degree at most four. Five distinct candidates
  // guarantee a nonroot split without perturbing the represented polynomial.
  Q split;bool found=false;
  for(const Q& fraction:std::array<Q,5>{Q(1,2),Q(1,4),Q(3,4),Q(1,8),Q(3,8)}) {
   split=node.lo+(node.hi-node.lo)*fraction;
   if(eval(polynomial,split)!=0){found=true;break;}
  }
  if(!found)return fail("no nonroot isolation split");
  pending.push_back({split,node.hi,node.depth+1});
  pending.push_back({node.lo,split,node.depth+1});
 }
 result.certified=true;
 result.diagnostic="exact rational Sturm-Tarski isolation: all simple admissible open-edge roots";
 return result;
}

QuadricRationalPointResult TruncationRootCertifier::rationalPoint(
 const Eigen::Vector3d& n0,const Eigen::Vector3d& n1,bool compactRationalCoordinates) const {
 return rationalPointImpl(n0,n1,compactRationalCoordinates,384);
}

QuadricRationalPointResult TruncationRootCertifier::rationalPointImpl(
 const Eigen::Vector3d& n0,const Eigen::Vector3d& n1,bool compactRationalCoordinates,
 unsigned enclosureBits) const {
 QuadricRationalPointResult out;
 const auto fail=[&](const std::string& why){out.diagnostic=why;return out;};
 if(!impl_ || impl_->rankTwo)return fail("rational truncation point requires regular full rank");
 P polynomial,s,u,v;std::vector<P> chain;std::string why;
 if(!impl_->prepareEdge(n0,n1,polynomial,s,u,v,chain,why))return fail(why);
 detail::PreparedSturmTarskiCounter prepared(polynomial);
 Q lo(0),hi(1);int polynomialRoots=0;
 const int admissible=admissibleRoots(polynomial,s,u,v,chain,lo,hi,polynomialRoots,&prepared);
 if(polynomialRoots<0 || admissible!=1)
  return fail("rational truncation point requires exactly one certified regular crossing");
 // Reuse the isolator's exact branch-sign oracle, but refine the one
 // PHYSICAL root within the existing coordinate-refinement budget. The
 // squared radius equation can have an arbitrarily nearby inadmissible root;
 // the legacy 64-depth all-polynomial-root isolator is not a prerequisite
 // for enclosing this uniquely certified analytic coordinate. Its public
 // contract and budget remain unchanged.
 std::array<P,3> inverseNormal;
 for(int axis=0;axis<3;++axis)for(int k=0;k<3;++k)
  inverseNormal[axis]=add(inverseNormal[axis],scale(
   P{Q(n0(k)),Q(n1(k))-Q(n0(k))},impl_->inverse[axis][k]));
 CGAL::Gmpz toleranceDenominator(1),sqrtDenominator(1);
 mpz_mul_2exp(toleranceDenominator.mpz(),toleranceDenominator.mpz(),enclosureBits);
 mpz_mul_2exp(sqrtDenominator.mpz(),sqrtDenominator.mpz(),1024);
 const Q tolerance(CGAL::Gmpz(1),toleranceDenominator);
 const auto sqrtBounds=[&](const Q& x) {
  CGAL::Gmpz scaled,root;
  mpz_mul_2exp(scaled.mpz(),x.numerator().mpz(),2048);
  mpz_fdiv_q(scaled.mpz(),scaled.mpz(),x.denominator().mpz());
  mpz_sqrt(root.mpz(),scaled.mpz());
  const Q low(root,sqrtDenominator);
  return Interval{low,low*low==x?low:Q(root+1,sqrtDenominator)};
 };
 const auto text=[](const Q& value) {
  std::ostringstream stream;stream.imbue(std::locale::classic());stream<<value;return stream.str();
 };
 const auto enclosedCoordinates=[&](const Q& first,const Q& last,
                                    std::array<Interval,3>* coordinates) {
  const auto scale=bound(s,first,last);
  if(scale.lo>Q(1e-14)) {
   const auto lowerSqrt=sqrtBounds(scale.lo),upperSqrt=sqrtBounds(scale.hi);
   if(lowerSqrt.lo>0) {
    const Interval reciprocal{1/upperSqrt.hi,1/lowerSqrt.lo};
    bool narrow=true;
    for(int axis=0;axis<3;++axis) {
     (*coordinates)[axis]=product(bound(inverseNormal[axis],first,last),reciprocal);
     (*coordinates)[axis].lo+=impl_->center[axis];(*coordinates)[axis].hi+=impl_->center[axis];
     narrow&=(*coordinates)[axis].hi-(*coordinates)[axis].lo<=tolerance;
    }
    return narrow;
   }
  }
  return false;
 };
 // Horner interval arithmetic, exact sqrt floor/ceil, positive reciprocal
 // and interval products all preserve interval inclusion. On nested root
 // brackets, the availability of these enclosures and their narrowness
 // predicate are therefore monotone. Check a bounded block's last bracket
 // first; once it is narrow, find the FIRST passing bracket in that block.
 // This preserves the original exact endpoints/representative rather than
 // accepting a later, different approximation. The original 2048 budget,
 // sqrt precision and all root refinements are unchanged.
 std::vector<std::pair<Q,Q>> pending;
 pending.reserve(16);
 for(int refinement=0;refinement<=2048;++refinement) {
  pending.emplace_back(lo,hi);
  if(refinement==0 || pending.size()==16 || refinement==2048) {
   std::array<Interval,3> coordinates;
   if(enclosedCoordinates(lo,hi,&coordinates)) {
    const auto first=detail::firstMonotoneTrue(pending.size(),[&](std::size_t index) {
     return enclosedCoordinates(pending[index].first,pending[index].second,&coordinates);
    });
    if(first==pending.size() || !enclosedCoordinates(
        pending[first].first,pending[first].second,&coordinates))
     return fail("nested coordinate enclosure lost its first passing bracket");
     for(int axis=0;axis<3;++axis) {
      out.lower[axis]=text(coordinates[axis].lo);out.upper[axis]=text(coordinates[axis].hi);
      out.coordinates[axis]=text(rationalSourceRepresentative(
          coordinates[axis].lo,coordinates[axis].hi,compactRationalCoordinates));
     }
     out.defined=true;
     out.diagnostic=enclosureBits==384
      ? "same isolated full-rank positive-root crossing; rational coordinate enclosure <=2^-384; not an exact analytic root"
      : "same isolated full-rank positive-root crossing; independent differential proof enclosure <=2^-512; not a new source point";
     return out;
   }
   pending.clear();
  }
  const auto scale=bound(s,lo,hi);
  Q mid;bool found=false;
  for(const Q& fraction:std::array<Q,5>{Q(1,2),Q(1,4),Q(3,4),Q(1,8),Q(3,8)}) {
   mid=lo+(hi-lo)*fraction;
   if(eval(polynomial,mid)!=0){found=true;break;}
  }
  if(!found)return fail("no nonroot rational coordinate refinement split");
  int leftRoots=0,left=-1;
  // The initial exact Sturm-Tarski proof establishes ONE simple physical
  // root. On an entirely finite bracket F=u+v*sqrt(s) is continuous, and
  // opposite nonzero endpoint signs enclose that same unique root. Exact
  // sign comparison then selects precisely the half that the branch-root
  // count would select, without re-counting nearby inadmissible roots.
  // Nonroot splits above ensure the midpoint is not a zero of F. Ideal or
  // uncertified-side brackets retain the original counting path.
  if(scale.lo>0) {
   const int before=detail::positiveRadicalSumSign(eval(u,lo),eval(v,lo),eval(s,lo));
   const int after=detail::positiveRadicalSumSign(eval(u,hi),eval(v,hi),eval(s,hi));
   if(before*after<0) {
    const int middle=detail::positiveRadicalSumSign(eval(u,mid),eval(v,mid),eval(s,mid));
    if(middle!=0)left=middle==before?0:1;
   }
  }
  if(left<0) {
   const int totalRoots=variations(chain,lo)-variations(chain,hi);
   left=totalRoots==1 ? variations(chain,lo)-variations(chain,mid)
    : admissibleRoots(polynomial,s,u,v,chain,lo,mid,leftRoots,&prepared);
  }
  if(left==1)hi=mid;
  else if(left==0)lo=mid;
  else return fail("isolated truncation root lost its unique branch during refinement");
 }
 return fail("rational truncation point resource limit before coordinate enclosure");
}

TruncationRationalDifferentialResult TruncationRootCertifier::rationalDifferential(
 const Eigen::Vector3d& n0,const Eigen::Vector3d& n1,bool compactRationalCoordinates) const {
 TruncationRationalDifferentialResult out;
 const auto fail=[&](const std::string& why){out.defined=false;out.diagnostic=why;return out;};
 // Preserve the actual frozen coordinate construction. The second enclosure
 // is proof-only, refers to the same unique certified analytic root, and never
 // substitutes its narrower representative into the Boolean operand.
 out.point=rationalPoint(n0,n1,compactRationalCoordinates);
 if(!out.point.defined)return fail(out.point.diagnostic);
 const auto proof=rationalPointImpl(n0,n1,compactRationalCoordinates,512);
 if(!proof.defined)return fail(proof.diagnostic);
 const auto& q=impl_->coefficients;
 const std::array<std::array<int,3>,3> indices{{{{4,9,8}},{{9,5,7}},{{8,7,6}}}};
 std::array<Interval,3> x;
 std::array<Q,3> edge;
 const auto constant=[](const Q& value){return Interval{value,value};};
 const auto scaled=[&](const Interval& value,const Q& factor){return product(value,constant(factor));};
 for(int i=0;i<3;++i) {
  x[i]={Q(proof.lower[i]),Q(proof.upper[i])};
  if(x[i].lo<Q(out.point.lower[i]) || x[i].hi>Q(out.point.upper[i]))
   return fail("root derivative proof is not nested in the frozen point enclosure");
  out.derivativePointLower[i]=proof.lower[i];out.derivativePointUpper[i]=proof.upper[i];
  edge[i]=Q(n0((i+1)%3))*Q(n1((i+2)%3))-Q(n0((i+2)%3))*Q(n1((i+1)%3));
 }
 IntervalMatrix system{};
 Interval surface=constant(q[0]),radius=constant(-impl_->radiusSquared),edgeGradient=constant(0);
 for(int i=0;i<3;++i) {
  system[0][i]=constant(q[i+1]);system[1][i]=scaled(x[i],Q(2));system[2][i]=constant(0);
  surface=sum(surface,scaled(x[i],q[i+1]));radius=sum(radius,product(x[i],x[i]));
  for(int j=0;j<3;++j) {
   system[0][i]=sum(system[0][i],scaled(x[j],2*q[indices[i][j]]));
   system[2][i]=sum(system[2][i],constant(2*edge[j]*q[indices[j][i]]));
   surface=sum(surface,scaled(product(x[i],x[j]),q[indices[i][j]]));
  }
  edgeGradient=sum(edgeGradient,scaled(system[0][i],edge[i]));
 }
 if(!containsZero(surface) || !containsZero(radius) || !containsZero(edgeGradient))
  return fail("root differential proof enclosure fails the actual implicit equations");
 const auto denominator=determinant(system);
 if(containsZero(denominator))return fail("root derivative determinant enclosure contains zero");
 CGAL::Gmpz scaleDenominator(1);mpz_mul_2exp(scaleDenominator.mpz(),scaleDenominator.mpz(),384);
 const Q tolerance(CGAL::Gmpz(1),scaleDenominator);
 const auto text=[](const Q& value){std::ostringstream stream;stream.imbue(std::locale::classic());stream<<value;return stream.str();};
 try {
  for(int column=0;column<10;++column) {
   std::array<Interval,3> rhs{{constant(column==0?Q(-1):Q(0)),constant(0),constant(0)}};
   if(column>=1 && column<=3) {
    rhs[0]=scaled(x[column-1],Q(-1));rhs[2]=constant(-edge[column-1]);
   } else if(column>=4 && column<=6) {
    const int axis=column-4;
    rhs[0]=scaled(product(x[axis],x[axis]),Q(-1));
    rhs[2]=scaled(x[axis],-2*edge[axis]);
   } else if(column>=7) {
    const int first=column==7?1:0,last=column==9?1:2;
    rhs[0]=scaled(product(x[first],x[last]),Q(-2));
    rhs[2]=scaled(sum(scaled(x[last],edge[first]),scaled(x[first],edge[last])),Q(-2));
   }
   for(int axis=0;axis<3;++axis) {
    auto replaced=system;for(int row=0;row<3;++row)replaced[row][axis]=rhs[row];
    const auto value=quotient(determinant(replaced),denominator);
    if(value.hi-value.lo>tolerance)
     return fail("root differential enclosure exceeds 2^-384 at the fixed proof precision");
    out.lower[axis][column]=text(value.lo);out.upper[axis][column]=text(value.hi);
    out.coefficientJacobian[axis][column]=text(rationalSourceRepresentative(
     value.lo,value.hi,compactRationalCoordinates));
   }
  }
 } catch(const std::runtime_error& error){return fail(error.what());}
 out.defined=true;
 out.diagnostic="unique positive-branch radius root; rational interval implicit differential <=2^-384; original source coordinates unchanged";
 return out;
}

TruncationEdgePartition TruncationRootCertifier::partition(
 const Eigen::Vector3d& n0,const Eigen::Vector3d& n1,int maximumDepth) const {
 const auto fail=[](const std::string& why) {
  TruncationEdgePartition out;out.diagnostic=why;return out;
 };
 const auto isolated=isolate(n0,n1,maximumDepth);
 if(!isolated.certified)return fail(isolated.diagnostic);
 P polynomial,s,u,v;std::vector<P> chain;std::string why;
 if(!impl_->prepareEdge(n0,n1,polynomial,s,u,v,chain,why))return fail(why);
 // Full rank: x=c+B*n/sqrt(s); for s>0, sign(|x|^2-R^2) is
 // sign(u+v*sqrt(s)). Compare exactly without taking a floating square root.
 // At s<=0 there is no finite positive branch. Across an ideal point the
 // finite branch diverges (B nonsingular, n nonzero), so finite-radius
 // membership can change only at one of the certified admissible roots.
 // Rank two: s=ub*w, v=0, u is the reduced radius polynomial. Its finite
 // denominator is squared and positive. At an ideal boundary the reduced
 // numerator is nonzero, so the radius polynomial is strictly positive.
 const auto insideAt=[&](const Q& t) {
  const Q st=eval(s,t);
  if(st<=0)return 0;
  const Q ut=eval(u,t),vt=eval(v,t);
  const int radiusSign=detail::positiveRadicalSumSign(ut,vt,st);
  // A bracket endpoint on the physical radius is not a valid side witness.
  return radiusSign==0 ? -1 : (radiusSign<0 ? 1 : 0);
 };
 TruncationEdgePartition result;
 int previous=insideAt(Q(0));
 if(previous<0)return fail("edge endpoint lies on the truncation sphere");
 result.inside.push_back(previous!=0);
 for(const auto& root:isolated.intervals) {
  const int before=insideAt(root.lower),after=insideAt(root.upper);
  if(before<0 || after<0 || before!=previous || before==after)
   return fail("isolated radius root lacks consistent exact crossing sides");
  result.roots.push_back(root);
  result.inside.push_back(after!=0);
  previous=after;
 }
 if(insideAt(Q(1))!=previous)
  return fail("final edge membership disagrees with the complete root sequence");
 result.certified=true;
 result.diagnostic="exact finite-radius normal-edge partition; face connectivity not certified";
 return result;
}

bool certifySingleTruncationRoot(const QuadricCoefficients10& quadric,
 const Eigen::Vector3d& n0,const Eigen::Vector3d& n1,double radius,std::string* why) {
 return TruncationRootCertifier(quadric,radius).certify(n0,n1,why);
}
}
