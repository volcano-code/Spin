#include "spin2014/SourceCageBinding.h"
#include "source_surface_test_utils.h"
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <type_traits>
using namespace spin2014;
using namespace surface_test;
template<class E,class F> void throws(F f) {
    bool caught=false; try {f();} catch(const E& e) {
        caught=true;
        if constexpr(std::is_same<E,std::length_error>::value)
            require(std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos,"resource diagnostic");
    }
    require(caught,"expected exception type");
}
SourceCageBbw solve(const Mesh& outer,Eigen::Vector3d origin=Eigen::Vector3d::Zero(),double side=1,
                    SourceCageBbwConfig bc={}) {
    SourceSurfaceCageConfig cc; cc.cageLevel=2; cc.rootOrigin=origin; cc.rootSide=side;
    SourceCageGridConfig gc; gc.fineLevel=3;
    return SourceCageBbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),bc);
}
std::shared_ptr<const SourceSurfaceOctree> oracle(const Mesh& outer,const Mesh& design,
    Eigen::Vector3d origin=Eigen::Vector3d::Zero(),double side=1) {
    SourceSurfaceOctreeConfig c; c.maximumDepth=5;c.rootOrigin=origin;c.rootSide=side;
    return std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,c);
}
void checkGeometry(const SourceCageBinding& b) {
    const auto& source=b.source(); const auto& leaves=b.domain().leaves();
    require(b.outerPointCount()==std::size_t(source.outerVertices().rows()),"original outer count");
    require(b.points().topRows(b.outerPointCount())==source.outerVertices(),"original outer order");
    const auto& keys=b.materialCornerKeys();
    require(std::is_sorted(keys.begin(),keys.end()) && std::set<SourceCageBinding::Key>(keys.begin(),keys.end()).size()==keys.size(),"canonical unique material order");
    for(std::size_t i=0;i<keys.size();++i) for(int a=0;a<3;++a)
        require(b.points()(b.outerPointCount()+i,a)==std::fma(std::ldexp(double(keys[i][a]),-source.config().maximumDepth),source.config().rootSide,source.config().rootOrigin[a]),"canonical point coordinates");
    std::vector<Ten> sums(leaves.size(),Ten::Zero());
    std::vector<std::size_t> owned(leaves.size(),0);
    for(const auto& p:b.primitives()) {
        require(source.cell(p.sourceKey.level,p.sourceKey.coordinates).kind==SourceSurfaceCellKind::Full,"actual Full source primitive");
        std::size_t owners=0;
        for(std::size_t i=0;i<leaves.size();++i) {
            bool match=leaves[i].level<=p.sourceKey.level;
            if(match) for(int a=0;a<3;++a) match=match && (p.sourceKey.coordinates[a]>>(p.sourceKey.level-leaves[i].level))==leaves[i].coordinates[a];
            if(match) {++owners;require(p.betaIndex==i,"owner preserves leaf order");}
        }
        require(owners==1,"unique ancestor owner"); ++owned.at(p.betaIndex);
        sums.at(p.betaIndex)+=boxMoments(p.bounds.lower,p.bounds.upper);
        const auto canonical=source.cellBounds(p.sourceKey.level,p.sourceKey.coordinates);
        require(p.bounds.lower==canonical.lower && p.bounds.upper==canonical.upper,"canonical primitive bounds");
        for(int k=0;k<8;++k) for(int a=0;a<3;++a) {
            require(p.pointIndices[k]>=b.outerPointCount(),"outer never deduplicated into corners");
            require(b.points()(p.pointIndices[k],a)==((k>>a)&1?p.bounds.upper[a]:p.bounds.lower[a]),"x-bit-first connectivity");
        }
    }
    Ten weighted=Ten::Zero(); Eigen::VectorXd beta(leaves.size());
    for(std::size_t i=0;i<leaves.size();++i) {
        require(owned[i]>0,"every adaptive leaf bound");
        close(sums[i],leaves[i].moments.values(),"all ten independently integrated owner moments");
        beta[i]=double(i+1)/(leaves.size()+1);weighted+=beta[i]*sums[i];
    }
    close(weighted,b.domain().removedMoments(beta).values(),"nonuniform beta moments");
    const auto& grid=b.bbw().fem().grid();
    for(std::size_t begin=0;begin<std::size_t(b.points().rows());begin+=7) {
        const auto count=std::min(std::size_t(7),std::size_t(b.points().rows())-begin);
        const auto w=b.pointWeights(begin,count); require(w.cols()==26 && w.allFinite(),"all handles finite");
        for(std::size_t i=0;i<count;++i) {
            const auto& e=b.embeddings().at(begin+i);
            const auto actual=grid.embed(b.points().row(begin+i).transpose());
            require(e.nodeIndices==actual.nodeIndices && e.coefficients==actual.coefficients && e.cellKey==actual.cellKey && e.cellIndex==actual.cellIndex,"actual cached embed");
            Eigen::Vector3d p=Eigen::Vector3d::Zero();Eigen::RowVectorXd expected=Eigen::RowVectorXd::Zero(26);
            for(int k=0;k<8;++k) {p+=e.coefficients[k]*grid.nodePositions().row(e.nodeIndices[k]).transpose();expected+=e.coefficients[k]*b.bbw().finalWeights().row(e.nodeIndices[k]);}
            require((p-b.points().row(begin+i).transpose()).norm()<1e-14,"embedding reconstructs unsnapped point");
            require((w.row(i)-expected).norm()<1e-15,"manual eight-node all-handle weights");
        }
    }
}
int main() {
    try {
        std::cout<<std::setprecision(17);
        const auto outer=cube(.3,.7),design=cube(.34,.66);
        const auto src=oracle(outer,design); SourceAdaptiveDomain domain(src,2,5);
        const auto bbw=solve(outer);require(bbw.success(),"actual coupled BBW success");
        const auto sourceHash=src->partitionSha256(),bbwHash=bbw.outputSha256(),partitionHash=domain.partitionHash();
        SourceCageBindingConfig cfg;cfg.maximumWeightValues=7*26;
        const SourceCageBinding b(domain,bbw,cfg),again(domain,bbw,cfg);
        require(domain.leaves().size()==8 && b.primitives().size()==496 && b.handleCount()==26,"unaltered full fixture");
        checkGeometry(b);
        require(std::abs(src->carveableMoments().mass()-.030517578125)<1e-16,"actual volume excludes bounding-cell excess");
        Ten total=Ten::Zero();for(const auto& p:b.primitives())total+=boxMoments(p.bounds.lower,p.bounds.upper);
        close(total,boxMoments(Eigen::Vector3d::Constant(11./32),Eigen::Vector3d::Constant(21./32)),"literal source-grid material extent and all nonzero cross moments");
        require(b.points()==again.points() && b.configSha256()==again.configSha256() && b.outputSha256()==again.outputSha256(),"deterministic complete binding");
        require(b.source().sourceSha256()!=b.bbw().fem().grid().cage().sourceSha256(),"different source digests accepted for exact common outer");
        require(b.pointWeights(b.points().rows(),0).rows()==0 && b.pointWeights(0,0).cols()==26,"empty chunks");
        throws<std::length_error>([&]{b.pointWeights(0,b.points().rows());});
        throws<std::out_of_range>([&]{b.pointWeights(b.points().rows()+1,0);});
        throws<std::out_of_range>([&]{b.pointWeights(1,std::numeric_limits<std::size_t>::max());});
        for(int which=0;which<3;++which) for(auto value:{std::size_t(0),std::numeric_limits<std::size_t>::max()}) {
            auto bad=cfg;if(which==0)bad.maximumPrimitives=value;if(which==1)bad.maximumPoints=value;if(which==2)bad.maximumWeightValues=value;
            throws<std::length_error>([&]{SourceCageBinding rejected(domain,bbw,bad);});
        }
        auto bad=cfg;bad.maximumPrimitives=495;throws<std::length_error>([&]{SourceCageBinding rejected(domain,bbw,bad);});
        bad=cfg;bad.maximumPoints=b.points().rows()-1;throws<std::length_error>([&]{SourceCageBinding rejected(domain,bbw,bad);});
        bad=cfg;bad.maximumPoints=7;throws<std::length_error>([&]{SourceCageBinding rejected(domain,bbw,bad);});
        for(int which=0;which<3;++which) {
            auto changed=cfg;if(which==0)++changed.maximumPrimitives;if(which==1)++changed.maximumPoints;if(which==2)++changed.maximumWeightValues;
            SourceCageBinding other(domain,bbw,changed);require(other.configSha256()!=b.configSha256() && other.outputSha256()!=b.outputSha256(),"every budget binds identity");
        }
        Eigen::VectorXd beta=Eigen::VectorXd::Constant(8,.5);domain.splitFractionalLeaves(beta,1e-6);
        SourceCageBinding split(domain,bbw,cfg);checkGeometry(split);
        require(split.outputSha256()!=b.outputSha256() && split.domain().leaves().size()>8,"split remaps new owner snapshot");
        beta.setOnes();domain.mergeBinarySiblings(beta,1e-6);
        SourceCageBinding merged(domain,bbw,cfg);checkGeometry(merged);
        require(b.domain().partitionHash()==partitionHash && b.domain().leaves().size()==8 && b.outputSha256()==again.outputSha256(),"old binding survives split/merge");
        require(src->partitionSha256()==sourceHash && bbw.outputSha256()==bbwHash && src->outerVertices()==outer.v && src->outerFaces()==outer.f,"immutable upstream arrays and identity");
        SourceAdaptiveDomain otherDesign(oracle(outer,cube(.36,.64)),2,5);SourceCageBinding accepted(otherDesign,bbw,cfg);
        require(accepted.outputSha256()!=b.outputSha256(),"different design changes binding");
        auto reordered=outer;reordered.v.row(0).swap(reordered.v.row(1));for(Eigen::Index i=0;i<reordered.f.size();++i){int& x=reordered.f.data()[i];if(x==0)x=1;else if(x==1)x=0;}
        SourceAdaptiveDomain wrongOuter(oracle(reordered,design),2,5);
        throws<std::invalid_argument>([&]{SourceCageBinding rejected(wrongOuter,bbw,cfg);});
        auto reorderedFaces=outer;reorderedFaces.f.row(0).swap(reorderedFaces.f.row(1));
        SourceAdaptiveDomain wrongFaces(oracle(reorderedFaces,design),2,5);
        throws<std::invalid_argument>([&]{SourceCageBinding rejected(wrongFaces,bbw,cfg);});
        SourceAdaptiveDomain wrongBytes(oracle(cube(std::nextafter(.3,0.),.7),design),2,5);
        throws<std::invalid_argument>([&]{SourceCageBinding rejected(wrongBytes,bbw,cfg);});
        SourceAdaptiveDomain wrongRoot(oracle(outer,design,Eigen::Vector3d(-0.,0,0)),2,5);
        throws<std::invalid_argument>([&]{SourceCageBinding rejected(wrongRoot,bbw,cfg);});
        SourceAdaptiveDomain wrongSide(oracle(outer,design,Eigen::Vector3d::Zero(),2),2,5);
        throws<std::invalid_argument>([&]{SourceCageBinding rejected(wrongSide,bbw,cfg);});
        SourceCageBbwConfig failConfig;failConfig.maximumVariables=1;const SourceCageBbw failed(bbw.fem(),failConfig);
        require(!failed.success(),"actual resource-failed BBW");throws<std::invalid_argument>([&]{SourceCageBinding rejected(b.domain(),failed,cfg);});
        failConfig=SourceCageBbwConfig{};failConfig.primalTolerance=0;const SourceCageBbw invalid(bbw.fem(),failConfig);
        throws<std::invalid_argument>([&]{SourceCageBinding rejected(b.domain(),invalid,cfg);});
        SourceCageBbwConfig independentConfig;independentConfig.protocol=SourceCageBbwProtocol::INDEPENDENT_NORMALIZED;
        const SourceCageBbw independent(bbw.fem(),independentConfig);SourceCageBinding independentBinding(b.domain(),independent,cfg);
        require(independent.success() && independentBinding.bbw().finalAssessment().primalPassed && !independentBinding.bbw().finalAssessment().optimalityPassed,"normalized acceptance retains noncoupled optimality");
        const Eigen::Vector3d shift(1.25,-.75,.125);const double scale=2;
        auto transformedOuter=outer,transformedDesign=design;
        transformedOuter.v=(outer.v*scale).rowwise()+shift.transpose();transformedDesign.v=(design.v*scale).rowwise()+shift.transpose();
        SourceAdaptiveDomain transformedDomain(oracle(transformedOuter,transformedDesign,shift,scale),2,5);
        const auto transformedBbw=solve(transformedOuter,shift,scale);SourceCageBinding transformed(transformedDomain,transformedBbw,cfg);checkGeometry(transformed);
        // Numerically equal outer vertices with distinct signed-zero bytes must
        // be rejected. Both actual audited geometries are otherwise identical.
        const auto zeroOuter=cube(0.,.4),zeroDesign=cube(.04,.36);
        const Eigen::Vector3d zeroRoot=Eigen::Vector3d::Constant(-.3);
        const auto zeroBbw=solve(zeroOuter,zeroRoot);
        auto negativeZeroOuter=zeroOuter;
        for(Eigen::Index i=0;i<negativeZeroOuter.v.size();++i)if(negativeZeroOuter.v.data()[i]==0)negativeZeroOuter.v.data()[i]=-0.;
        SourceAdaptiveDomain signedZeroDomain(oracle(negativeZeroOuter,zeroDesign,zeroRoot),2,5);
        require(negativeZeroOuter.v==zeroOuter.v,"signed zero is numerically equal");
        throws<std::invalid_argument>([&]{SourceCageBinding rejected(signedZeroDomain,zeroBbw,cfg);});
        SourceAdaptiveDomain matchingSignedDomain(oracle(negativeZeroOuter,zeroDesign,zeroRoot),2,5);
        const auto signedBbw=solve(negativeZeroOuter,zeroRoot);SourceCageBinding signedBinding(matchingSignedDomain,signedBbw,cfg);
        for(Eigen::Index i=0;i<8;++i)for(int a=0;a<3;++a)
            require(std::memcmp(&signedBinding.points()(i,a),&negativeZeroOuter.v(i,a),sizeof(double))==0,"original signed-zero point bytes preserved");
        std::cout<<"FACT source="<<sourceHash<<" bbw="<<bbwHash<<" binding_config="<<b.configSha256()<<" binding_output="<<b.outputSha256()<<" primitives="<<b.primitives().size()<<" points="<<b.points().rows()<<" volume="<<src->carveableMoments().mass()<<'\n';
        std::cout<<"PASS SYNTHETIC_REGRESSION actual immutable source material/full BBW binding\n";
    } catch(const std::exception& e) {std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
    return 0;
}
