#include "spin2014/SourceCageBinding.h"
#include "spin2014/SourceBinaryMoments.h"
#include "verification/Sha256.h"
#include <cstring>
#include <limits>
#include <map>
#include <new>
#include <unordered_map>

namespace spin2014 {
namespace {
using Key=SourceCageBinding::Key;
using Primitive=SourceCageBinding::Primitive;
struct Corner { Eigen::Vector3d position; std::size_t index=0; };
struct PartitionHash {
    std::size_t operator()(const SourcePartitionKey& k) const {
        std::size_t h=k.level;
        for(auto c:k.coordinates) h=h*std::size_t(16777619)^c;
        return h;
    }
};
using Owners=std::unordered_map<SourcePartitionKey,std::size_t,PartitionHash>;
class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {
        unsigned char bytes[8];for(int i=0;i<8;++i)bytes[7-i]=static_cast<unsigned char>(x>>(8*i));
        hash_.update(bytes,8);
    }
    void string(const std::string& x) {
        integer(x.size());hash_.update(reinterpret_cast<const unsigned char*>(x.data()),x.size());
    }
    void real(double x) {
        static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,"IEEE754 binary64 required");
        std::uint64_t bits;std::memcpy(&bits,&x,8);integer(bits);
    }
    template<class T,std::size_t N> void integers(const std::array<T,N>& x) {for(auto v:x)integer(v);}
    template<class T> void matrix(const T& x) {
        integer(x.rows());integer(x.cols());
        for(Eigen::Index r=0;r<x.rows();++r)for(Eigen::Index c=0;c<x.cols();++c) {
            if constexpr(std::is_floating_point<typename T::Scalar>::value)real(x(r,c));else integer(x(r,c));
        }
    }
    std::string finish() {return hash_.finish();}
};
bool sameReal(double a,double b) {return std::memcmp(&a,&b,sizeof(double))==0;}
template<class T> bool sameMatrix(const T& a,const T& b) {
    if(a.rows()!=b.rows() || a.cols()!=b.cols())return false;
    for(Eigen::Index r=0;r<a.rows();++r)for(Eigen::Index c=0;c<a.cols();++c)
        if(std::memcmp(&a(r,c),&b(r,c),sizeof(typename T::Scalar))!=0)return false;
    return true;
}
void validate(const SourceAdaptiveDomain& domain,const SourceCageBbw& bbw,const SourceCageBindingConfig& c) {
    const auto indexMax=std::size_t(std::numeric_limits<Eigen::Index>::max());
    const auto bytesMax=std::numeric_limits<std::size_t>::max();
    if(!c.maximumPrimitives || c.maximumPrimitives>std::vector<Primitive>().max_size() ||
       c.maximumPrimitives>std::vector<SourceRemovedCuboid>().max_size() ||
       c.maximumPrimitives>Owners().max_size() || c.maximumPrimitives>indexMax ||
       c.maximumPrimitives>bytesMax/sizeof(std::size_t) ||
       !c.maximumPoints || c.maximumPoints>std::vector<Key>().max_size() ||
       c.maximumPoints>std::map<Key,Corner>().max_size() ||
       c.maximumPoints>std::vector<SourceCageGrid::Embedding>().max_size() ||
       c.maximumPoints>indexMax/3 || c.maximumPoints>bytesMax/(3*sizeof(double)) ||
       !c.maximumWeightValues || c.maximumWeightValues>indexMax || c.maximumWeightValues>bytesMax/sizeof(double))
        throw std::length_error("RESOURCE_LIMIT: invalid source cage binding primitive/point/weight budget");
    if(!bbw.success() || !bbw.finalAssessment().available || !bbw.finalAssessment().primalPassed)
        throw std::invalid_argument("source cage binding requires actual successful BBW and final primal validity");
    const auto& source=domain.source();const auto& cage=bbw.fem().grid().cage();
    if(!sameMatrix(source.outerVertices(),cage.sourceVertices()) || !sameMatrix(source.outerFaces(),cage.sourceFaces()) ||
       !sameMatrix(source.config().rootOrigin,cage.config().rootOrigin) || !sameReal(source.config().rootSide,cage.config().rootSide))
        throw std::invalid_argument("source cage binding original outer arrays/root differ in dimensions or scalar bytes");
    if(domain.leaves().size()>c.maximumPrimitives || std::size_t(source.outerVertices().rows())>c.maximumPoints)
        throw std::length_error("RESOURCE_LIMIT: source cage binding leaf/original point budget");
    const auto& grid=bbw.fem().grid();const auto& weights=bbw.finalWeights();
    if(weights.rows()!=grid.nodePositions().rows() || weights.cols()<=0 ||
       std::size_t(weights.cols())!=grid.handleNodeIndices().size() || !weights.allFinite())
        throw std::invalid_argument("source cage binding invalid actual final BBW field");
}
Key cornerKey(const SourcePartitionKey& p,int corner,int depth) {
    Key key;for(int a=0;a<3;++a)key[a]=(p.coordinates[a]+((corner>>a)&1))<<(depth-p.level);
    return key;
}
}

struct SourceCageBinding::Data {
    SourceAdaptiveDomain domain;
    SourceCageBbw bbw;
    SourceCageBindingConfig config;
    Eigen::MatrixXd points;
    std::vector<Key> keys;
    std::vector<Primitive> primitives;
    std::vector<SourceCageGrid::Embedding> embeddings;
    std::string configHash,outputHash;
    Data(const SourceAdaptiveDomain& d,const SourceCageBbw& b,const SourceCageBindingConfig& c):domain(d),bbw(b),config(c) {
        domain.validateCoverage();
        const auto& source=domain.source();const auto& leaves=domain.leaves();
        const auto outerCount=std::size_t(source.outerVertices().rows());
        const auto actual=sourceBinaryCuboids(domain,Eigen::VectorXd::Ones(static_cast<Eigen::Index>(leaves.size())),config.maximumPrimitives);
        Owners owners;owners.reserve(leaves.size());
        for(std::size_t i=0;i<leaves.size();++i)
            if(!owners.emplace(leaves[i],i).second)throw std::invalid_argument("source cage binding duplicate adaptive owner");
        std::vector<std::size_t> owned(leaves.size(),0);
        std::map<Key,Corner> corners;
        primitives.reserve(actual.cuboids.size());
        for(const auto& cuboid:actual.cuboids) {
            Primitive p;p.sourceKey=cuboid.cell;p.bounds={cuboid.lower,cuboid.upper};
            // Expected O(depth) hash lookups; never scan all adaptive leaves.
            auto ancestor=p.sourceKey;std::size_t count=0;
            for(;;) {
                const auto found=owners.find(ancestor);
                if(found!=owners.end()){p.betaIndex=found->second;++count;}
                if(ancestor.level==0)break;
                --ancestor.level;for(auto& coordinate:ancestor.coordinates)coordinate>>=1;
            }
            if(count!=1)throw std::invalid_argument("source cage binding requires exactly one primitive owner");
            ++owned[p.betaIndex];
            for(int k=0;k<8;++k) {
                const auto key=cornerKey(p.sourceKey,k,source.config().maximumDepth);
                Eigen::Vector3d position;for(int a=0;a<3;++a)position[a]=(k>>a)&1?p.bounds.upper[a]:p.bounds.lower[a];
                const auto found=corners.find(key);
                if(found!=corners.end()) {
                    if(!sameMatrix(found->second.position,position))throw std::invalid_argument("source cage binding inconsistent canonical endpoint bytes");
                    continue;
                }
                if(corners.size()>=config.maximumPoints-outerCount)
                    throw std::length_error("RESOURCE_LIMIT: source cage binding material point insertion budget");
                corners.emplace(key,Corner{position,0});
            }
            primitives.push_back(p);
        }
        for(auto count:owned)if(!count)throw std::invalid_argument("source cage binding adaptive leaf has no Full primitives");
        points.resize(static_cast<Eigen::Index>(outerCount+corners.size()),3);
        points.topRows(source.outerVertices().rows())=source.outerVertices();keys.reserve(corners.size());
        for(auto& entry:corners) {
            entry.second.index=outerCount+keys.size();keys.push_back(entry.first);
            points.row(static_cast<Eigen::Index>(entry.second.index))=entry.second.position.transpose();
        }
        for(auto& p:primitives)for(int k=0;k<8;++k)
            p.pointIndices[k]=corners.at(cornerKey(p.sourceKey,k,source.config().maximumDepth)).index;
        embeddings.reserve(static_cast<std::size_t>(points.rows()));
        for(Eigen::Index i=0;i<points.rows();++i) {
            try {embeddings.push_back(bbw.fem().grid().embed(points.row(i).transpose()));}
            catch(const std::out_of_range&) {throw std::invalid_argument("source cage binding point outside actual BBW grid");}
        }
        Digest cfg;cfg.string("source-cage-binding-config-v1:actual-Full-rest-cuboids:outer-first:canonical-maxDepth-corners:all-handles:cached-eight-node");
        cfg.string(source.configSha256());cfg.string(bbw.configSha256());
        cfg.integer(config.maximumPrimitives);cfg.integer(config.maximumPoints);cfg.integer(config.maximumWeightValues);
        configHash=cfg.finish();
        Digest out;out.string("source-cage-binding-output-v1");out.string(configHash);
        out.string(source.sourceSha256());out.string(source.partitionSha256());out.string(domain.partitionHash());
        out.string(bbw.fem().grid().cage().sourceSha256());out.string(bbw.outputSha256());
        out.matrix(source.outerVertices());out.matrix(source.outerFaces());out.matrix(source.config().rootOrigin);out.real(source.config().rootSide);
        out.matrix(bbw.finalWeights());out.integer(bbw.fem().grid().handleNodeIndices().size());
        out.integer(leaves.size());for(const auto& leaf:leaves){out.integer(leaf.level);out.integers(leaf.coordinates);out.matrix(leaf.moments.values());}
        out.integer(outerCount);out.matrix(points);out.integer(keys.size());for(const auto& key:keys)out.integers(key);
        out.integer(primitives.size());for(const auto& p:primitives) {
            out.integer(p.sourceKey.level);out.integers(p.sourceKey.coordinates);out.integer(p.betaIndex);
            out.matrix(p.bounds.lower);out.matrix(p.bounds.upper);out.integers(p.pointIndices);
        }
        out.integer(embeddings.size());for(const auto& e:embeddings) {
            out.integer(e.cellIndex);out.integers(e.cellKey);out.integers(e.nodeIndices);for(auto w:e.coefficients)out.real(w);
        }
        outputHash=out.finish();
    }
};

SourceCageBinding::SourceCageBinding(const SourceAdaptiveDomain& d,const SourceCageBbw& b,const SourceCageBindingConfig& c) {
    try {validate(d,b,c);data_=std::make_shared<const Data>(d,b,c);}
    catch(const std::bad_alloc&) {throw std::length_error("RESOURCE_LIMIT: source cage binding allocation failed");}
}
const SourceAdaptiveDomain& SourceCageBinding::domain() const{return data_->domain;}
const SourceSurfaceOctree& SourceCageBinding::source() const{return data_->domain.source();}
const SourceCageBbw& SourceCageBinding::bbw() const{return data_->bbw;}
const SourceCageBindingConfig& SourceCageBinding::config() const{return data_->config;}
const Eigen::MatrixXd& SourceCageBinding::points() const{return data_->points;}
const std::vector<SourceCageBinding::Key>& SourceCageBinding::materialCornerKeys() const{return data_->keys;}
const std::vector<SourceCageBinding::Primitive>& SourceCageBinding::primitives() const{return data_->primitives;}
const std::vector<SourceCageGrid::Embedding>& SourceCageBinding::embeddings() const{return data_->embeddings;}
std::size_t SourceCageBinding::outerPointCount() const{return static_cast<std::size_t>(source().outerVertices().rows());}
std::size_t SourceCageBinding::handleCount() const{return bbw().fem().grid().handleNodeIndices().size();}
const std::string& SourceCageBinding::configSha256() const{return data_->configHash;}
const std::string& SourceCageBinding::outputSha256() const{return data_->outputHash;}
Eigen::MatrixXd SourceCageBinding::pointWeights(std::size_t begin,std::size_t count) const {
    const auto total=static_cast<std::size_t>(points().rows()),handles=handleCount();
    if(begin>total || count>total-begin)throw std::out_of_range("source cage binding point chunk range");
    if(count>config().maximumWeightValues/handles)
        throw std::length_error("RESOURCE_LIMIT: source cage binding weight chunk budget");
    // Validated maximumWeightValues also bounds Eigen index and allocation bytes.
    try {
        Eigen::MatrixXd output=Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(count),static_cast<Eigen::Index>(handles));
        for(std::size_t i=0;i<count;++i) {
            const auto& e=embeddings()[begin+i];
            for(int k=0;k<8;++k)output.row(static_cast<Eigen::Index>(i))+=e.coefficients[k]*bbw().finalWeights().row(e.nodeIndices[k]);
        }
        if(!output.allFinite())throw std::overflow_error("source cage binding nonfinite interpolated weights");
        return output;
    } catch(const std::bad_alloc&) {throw std::length_error("RESOURCE_LIMIT: source cage binding weight chunk allocation failed");}
}
}
