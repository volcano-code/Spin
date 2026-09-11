#include "spin2014/SourceCageBinding.h"
#include "source_surface_test_utils.h"
#include <iomanip>
#include <iostream>
template<class T> void values(const T& a){for(const auto& x:a)std::cout<<' '<<x;}
template<class T> void row(const T& a){for(Eigen::Index i=0;i<a.size();++i)std::cout<<' '<<a(i);}
int main(){
 using namespace spin2014;
 std::cout<<std::setprecision(17);
 for(int fixture=0;fixture<2;++fixture)for(int protocol=0;protocol<2;++protocol){
  const double scale=fixture?1.7:1.;
  const Eigen::Vector3d offset=fixture?Eigen::Vector3d(-.2,.15,.6):Eigen::Vector3d::Zero();
  auto outer=fixture?surface_test::box({.31,.32,.3},{.69,.65,.7}):surface_test::cube(.3,.7);
  auto design=fixture?surface_test::box({.35,.36,.38},{.61,.6,.65}):surface_test::cube(.34,.66);
  outer.v=(scale*outer.v).rowwise()+offset.transpose();design.v=(scale*design.v).rowwise()+offset.transpose();
  SourceSurfaceOctreeConfig sc;sc.maximumDepth=5;sc.rootOrigin=offset;sc.rootSide=scale;
  auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
  SourceAdaptiveDomain domain(source,2,5);
  SourceSurfaceCageConfig cc;cc.cageLevel=2;cc.rootOrigin=offset;cc.rootSide=scale;
  SourceCageGridConfig gc;gc.fineLevel=3;
  SourceCageFemConfig fc;
  SourceCageBbwConfig bc;bc.protocol=protocol?SourceCageBbwProtocol::INDEPENDENT_NORMALIZED:SourceCageBbwProtocol::COUPLED_PARTITION;
  const SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),fc),bc);
  if(!bbw.success())throw std::runtime_error(bbw.terminalReason());
  const SourceCageBinding binding(domain,bbw,{});
  const auto& grid=bbw.fem().grid();
  std::cout<<"CASE "<<fixture<<' '<<protocol<<' '<<sc.maximumDepth<<' '<<domain.leaves().size()<<' '
    <<binding.primitives().size()<<' '<<binding.points().rows()<<' '<<binding.handleCount()<<' '
    <<grid.nodePositions().rows()<<' '<<grid.cells().size()<<' '<<binding.outerPointCount()<<'\n';
  std::cout<<"IDENT "<<binding.outputSha256()<<' '<<binding.configSha256()<<' '<<source->sourceSha256()<<' '<<bbw.outputSha256()<<'\n';
  std::cout<<"ROOT";row(offset);std::cout<<' '<<scale<<'\n';
  std::cout<<"EXPECTED_LO";values(fixture?std::array<int,3>{12,12,13}:std::array<int,3>{11,11,11});std::cout<<'\n';
  std::cout<<"EXPECTED_HI";values(fixture?std::array<int,3>{19,19,20}:std::array<int,3>{21,21,21});std::cout<<'\n';
  for(Eigen::Index i=0;i<outer.v.rows();++i){std::cout<<"OUTER "<<i;row(outer.v.row(i));std::cout<<'\n';}
  for(std::size_t i=0;i<domain.leaves().size();++i){const auto& l=domain.leaves()[i];std::cout<<"LEAF "<<i<<' '<<l.level;values(l.coordinates);row(l.moments.values());std::cout<<'\n';}
  for(std::size_t i=0;i<binding.primitives().size();++i){const auto& p=binding.primitives()[i];std::cout<<"PRIM "<<i<<' '<<p.betaIndex<<' '<<p.sourceKey.level;values(p.sourceKey.coordinates);row(p.bounds.lower);row(p.bounds.upper);values(p.pointIndices);std::cout<<'\n';}
  for(Eigen::Index i=0;i<binding.points().rows();++i){std::cout<<"POINT "<<i;row(binding.points().row(i));std::cout<<'\n';const auto& e=binding.embeddings()[i];std::cout<<"EMBED "<<i<<' '<<e.cellIndex;values(e.cellKey);values(e.nodeIndices);values(e.coefficients);std::cout<<'\n';}
  for(std::size_t i=0;i<binding.materialCornerKeys().size();++i){std::cout<<"KEY "<<i;values(binding.materialCornerKeys()[i]);std::cout<<'\n';}
  for(Eigen::Index i=0;i<grid.nodePositions().rows();++i){std::cout<<"NODE "<<i;row(grid.nodePositions().row(i));row(bbw.finalWeights().row(i));std::cout<<'\n';}
  for(std::size_t i=0;i<grid.cells().size();++i){std::cout<<"GRIDCELL "<<i;values(grid.cells()[i]);values(grid.cellNodes()[i]);std::cout<<'\n';}
  for(std::size_t begin=0;begin<std::size_t(binding.points().rows());begin+=7){const auto count=std::min<std::size_t>(7,binding.points().rows()-begin);const auto w=binding.pointWeights(begin,count);for(Eigen::Index i=0;i<w.rows();++i){std::cout<<"WEIGHT "<<begin+i;row(w.row(i));std::cout<<'\n';}}
  std::cout<<"END\n";
 }
}
