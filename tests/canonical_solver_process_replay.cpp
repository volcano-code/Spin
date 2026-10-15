#define main frozen_regression_main
#include M1_FROZEN_REGRESSION
#undef main
#include <vector>
int main(int argc,char** argv){try {
 const int noise=argc>1?std::stoi(argv[1]):0;
 require(noise>=0 && noise<=12,"noise must be0..12");
 auto r=request();auto es=settings();r.designDomain=box(-.8,.8);r.outerSurface=box(-1,1);
 r.designDomain.vertices.col(2).array()+=2;r.outerSurface.vertices.col(2).array()+=2;
 r.minimumWallThickness=.15;r.materialDensity=1;r.mode=spin::SpinMode::SpinningTop;r.centerOfMassWeight=1;
 r.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.outerSurface.vertices,r.outerSurface.faces,1))-
  physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(r.designDomain.vertices,r.designDomain.faces,1));
 r.settings.maximumRestarts=1;r.settings.maximumFunctionEvaluationsPerRestart=20;es.triangulation.icosphereSubdivisions=1;
 const auto seed=fixture();const Eigen::Vector3d center(0,0,2);V p=seed.values();
 p(0)=seed.constant()-seed.linearVector().dot(center)+center.dot(seed.quadraticMatrix()*center);
 p.segment<3>(1)=seed.linearVector()-2*seed.quadraticMatrix()*center;
 const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
 std::vector<Eigen::MatrixXd> retained;
 for(int i=0;i<noise;++i) {
  retained.emplace_back(Eigen::MatrixXd::Constant(17+i*19,23+i,static_cast<double>(i)));
  require(spin::evaluateCanonicalFullRankMaterialMoments(r,q,es).coefficientJacobian.has_value(),"warmup unavailable");
 }
 auto valueSettings=es;valueSettings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteEllipsoid;
 const auto native=spin::evaluateQuadricMaterialMoments(r,q,valueSettings);
 const auto canonical=spin::evaluateCanonicalFullRankMaterialMoments(r,q,es);
 require(native.geometryHash==canonical.geometryHash && native.totalMomentsTarget.values()==canonical.totalMomentsTarget.values(),"canonical derivative changed native material values/geometry");
 require(!canonical.scientificHash.empty(),"canonical evaluation policy lacks identity");
 const auto result=spin::solveCanonicalFullRankQuadricHollowing(r,es,{q});
 require(result.finalGeometryVerified,"process replay final audit failed");
 std::cout<<result.inputHash<<' '<<result.outputHash<<' '<<result.numerical.eventHash<<' '<<result.numerical.traceHash<<' '<<result.finalEvaluation->geometryHash<<'\n';
 return 0;
}catch(const std::exception& e){std::cerr<<"[FAIL] "<<e.what()<<'\n';return 1;}}
