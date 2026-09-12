#include "spin2014/SourceConformingBoundary.h"
#include "spin2014/SurfaceRawMomentDifferential.h"
#include "source_surface_test_utils.h"
#include <iostream>
#include <iomanip>
#include <map>
template<class T> void seq(const T& xs){for(auto x:xs)std::cout<<' '<<x;}
template<class T> void row(const T& xs){for(Eigen::Index i=0;i<xs.size();++i)std::cout<<' '<<xs(i);}
using Boundary=spin2014::SourceConformingBoundary;
physics::RawMassMoments10 moments(const Eigen::MatrixXd& points,
 const std::vector<Boundary::Triangle>& triangles,const std::vector<Boundary::OrientedTriangle>& refs){
 std::map<std::size_t,int> index;
 for(const auto& ref:refs)for(auto v:triangles[ref.triangle])if(!index.count(v))index.emplace(v,int(index.size()));
 Eigen::MatrixXd v(index.size(),3);for(auto entry:index)v.row(entry.second)=points.row(entry.first);
 Eigen::MatrixXi f(refs.size(),3);
 for(std::size_t i=0;i<refs.size();++i){auto t=triangles[refs[i].triangle];if(refs[i].orientation<0)std::swap(t[1],t[2]);for(int a=0;a<3;++a)f(i,a)=index.at(t[a]);}
 return spin2014::integrateSurfaceRawMoments(v,f);
}
int main(){
 using namespace spin2014;
 std::cout<<std::setprecision(17);
 for(int fixture=0;fixture<3;++fixture){
  const double scale=fixture==1?1.7:1.;
  const Eigen::Vector3d offset=fixture==1?Eigen::Vector3d(-.2,.15,.6):Eigen::Vector3d::Zero();
  auto outer=fixture==1?surface_test::box({.31,.32,.3},{.69,.65,.7}):surface_test::cube(.3,fixture==2?.4:.7);
  auto design=fixture==1?surface_test::box({.35,.36,.38},{.61,.6,.65}):surface_test::cube(.34,fixture==2?.35:.66);
  outer.v=(scale*outer.v).rowwise()+offset.transpose();design.v=(scale*design.v).rowwise()+offset.transpose();
  SourceSurfaceOctreeConfig sc;sc.maximumDepth=fixture==2?9:5;sc.rootOrigin=offset;sc.rootSide=scale;
  auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
  SourceAdaptiveDomain domain(source,2,sc.maximumDepth);
  SourceSurfaceCageConfig cc;cc.cageLevel=2;cc.rootOrigin=offset;cc.rootSide=scale;
  SourceCageGridConfig gc;gc.fineLevel=3;
  SourceCageFemConfig fc;SourceCageBbwConfig bc;
  const SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),fc),bc);
  if(!bbw.success())throw std::runtime_error(bbw.terminalReason());
  SourceCageBindingConfig bindingConfig;
  const SourceCageBinding binding(domain,bbw,bindingConfig);
  SourceConformingBoundaryConfig config;
  const Boundary b(binding,config);
  const auto& grid=bbw.fem().grid();const auto& cage=grid.cage();
  Eigen::MatrixXd moved=b.points();
  for(Eigen::Index start=0;start<b.points().rows();start+=7){
   const Eigen::Index count=std::min<Eigen::Index>(7,b.points().rows()-start);
   const auto w=grid.interpolate(bbw.finalWeights(),b.points().middleRows(start,count));
   for(Eigen::Index i=0;i<count;++i){
    Eigen::Vector3d result=Eigen::Vector3d::Zero();
    for(Eigen::Index h=0;h<w.cols();++h){
     const Eigen::Vector3d center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]).transpose();
     const double hs=1+.01*(int(h%5)-2);
     const Eigen::Vector3d t(.006*(int(h%3)-1),.004*(int(h%5)-2),.005*(int(h%7)-3));
     result+=w(i,h)*(center+hs*(b.points().row(start+i).transpose()-center)+scale*t);
    }
    moved.row(start+i)=result.transpose();
   }
  }
  std::cout<<"CASE "<<fixture<<' '<<sc.maximumDepth<<' '<<b.latticeLevel()<<' '<<binding.primitives().size()<<' '
   <<b.points().rows()<<' '<<b.patches().size()<<' '<<b.triangles().size()<<' '<<b.outerPointCount()<<' '<<binding.handleCount()<<'\n';
  std::cout<<"IDENT "<<binding.outputSha256()<<' '<<b.configSha256()<<' '<<b.outputSha256()<<'\n';
  std::cout<<"ROOT";row(offset);std::cout<<' '<<scale<<'\n';
  const std::array<int,3> lo=fixture==0?std::array<int,3>{11,11,11}:(fixture==1?std::array<int,3>{12,12,13}:std::array<int,3>{175,175,175});
  const std::array<int,3> hi=fixture==0?std::array<int,3>{21,21,21}:(fixture==1?std::array<int,3>{19,19,20}:std::array<int,3>{179,179,179});
  std::cout<<"EXPECTED_LO";seq(lo);std::cout<<'\n';std::cout<<"EXPECTED_HI";seq(hi);std::cout<<'\n';
  for(Eigen::Index i=0;i<grid.nodePositions().rows();++i){std::cout<<"NODE "<<i;row(grid.nodePositions().row(i));row(bbw.finalWeights().row(i));std::cout<<'\n';}
  for(std::size_t i=0;i<grid.cells().size();++i){std::cout<<"GRIDCELL "<<i;seq(grid.cells()[i]);seq(grid.cellNodes()[i]);std::cout<<'\n';}
  for(std::size_t h=0;h<binding.handleCount();++h){std::cout<<"HANDLE "<<h;row(cage.boundaryVertices().row(cage.handleVertexIndices()[h]));std::cout<<' '<<1+.01*(int(h%5)-2)<<' '<<scale*.006*(int(h%3)-1)<<' '<<scale*.004*(int(h%5)-2)<<' '<<scale*.005*(int(h%7)-3)<<'\n';}
  for(Eigen::Index i=0;i<outer.v.rows();++i){std::cout<<"OUTER "<<i;row(outer.v.row(i));std::cout<<'\n';}
  for(Eigen::Index i=0;i<b.points().rows();++i){std::cout<<"POINT "<<i;row(b.points().row(i));std::cout<<'\n';std::cout<<"MOVED "<<i;row(moved.row(i));std::cout<<'\n';}
  for(std::size_t i=0;i<b.materialPointKeys().size();++i){std::cout<<"KEY "<<i;seq(b.materialPointKeys()[i]);std::cout<<'\n';}
  for(std::size_t i=0;i<binding.primitives().size();++i){const auto& p=binding.primitives()[i];std::cout<<"PRIM "<<i<<' '<<p.betaIndex<<' '<<p.sourceKey.level;seq(p.sourceKey.coordinates);row(p.bounds.lower);row(p.bounds.upper);std::cout<<'\n';}
  for(std::size_t i=0;i<b.patches().size();++i){const auto& p=b.patches()[i];std::cout<<"PATCH "<<i<<' '<<p.axis<<' '<<p.plane<<' '<<p.level;seq(p.coordinates);std::cout<<' '<<(p.negativeOwner==Boundary::noOwner?-1:static_cast<long long>(p.negativeOwner))<<' '<<(p.positiveOwner==Boundary::noOwner?-1:static_cast<long long>(p.positiveOwner))<<' '<<p.triangleBegin<<' '<<p.triangleCount<<'\n';}
  for(std::size_t i=0;i<b.triangles().size();++i){std::cout<<"TRI "<<i;seq(b.triangles()[i]);std::cout<<'\n';}
  for(std::size_t i=0;i<b.primitiveBoundaries().size();++i){const auto& refs=b.primitiveBoundaries()[i];for(const auto& ref:refs)std::cout<<"REF "<<i<<' '<<ref.triangle<<' '<<ref.orientation<<'\n';std::cout<<"MOMREST "<<i;row(moments(b.points(),b.triangles(),refs).values());std::cout<<'\n';std::cout<<"MOMMOVED "<<i;row(moments(moved,b.triangles(),refs).values());std::cout<<'\n';}
  std::vector<Boundary::OrientedTriangle> exterior;
  for(const auto& p:b.patches())if(p.negativeOwner==Boundary::noOwner||p.positiveOwner==Boundary::noOwner)
   for(std::size_t j=0;j<p.triangleCount;++j)exterior.push_back({p.triangleBegin+j,p.positiveOwner==Boundary::noOwner?1:-1});
  std::cout<<"UNIONREST";row(moments(b.points(),b.triangles(),exterior).values());std::cout<<'\n';
  std::cout<<"UNIONMOVED";row(moments(moved,b.triangles(),exterior).values());std::cout<<"\nEND\n";
 }
}
