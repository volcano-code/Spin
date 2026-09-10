#include "spin2014/SourceCageGrid.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <set>

namespace spin2014 {
namespace {
using Key=SourceCageGrid::Key;
using CellNodes=SourceCageGrid::CellNodes;
using Edge=SourceCageGrid::Edge;

class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {
        unsigned char bytes[8];
        for(int i=0;i<8;++i) bytes[7-i]=static_cast<unsigned char>(x>>(8*i));
        hash_.update(bytes,8);
    }
    void string(const std::string& x) {
        integer(x.size()); hash_.update(reinterpret_cast<const unsigned char*>(x.data()),x.size());
    }
    void real(double x) {
        static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,"IEEE754 binary64 required");
        std::uint64_t bits; std::memcpy(&bits,&x,8); integer(bits);
    }
    template<class T,std::size_t N> void array(const std::array<T,N>& x) {for(auto v:x) integer(v);}
    std::string finish() {return hash_.finish();}
};

double plane(const SourceSurfaceCageConfig& root,int level,int axis,std::uint32_t index) {
    return std::fma(std::ldexp(double(index),-level),root.rootSide,root.rootOrigin(axis));
}
Key cornerKey(Key key,int b) {for(int a=0;a<3;++a) key[a]+=(b>>a)&1; return key;}

void validate(const SourceSurfaceCageConfig& root,const SourceCageGridConfig& c) {
    if(c.fineLevel<=root.cageLevel || c.fineLevel>20)
        throw std::invalid_argument("source cage grid requires cageLevel < fineLevel <= 20");
    const auto indexMax=static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max());
    const auto bytesMax=std::numeric_limits<std::size_t>::max();
    if(!c.maximumCells || c.maximumCells>std::vector<Key>().max_size() ||
       c.maximumCells>std::vector<CellNodes>().max_size() ||
       !c.maximumNodes || c.maximumNodes>std::size_t(std::numeric_limits<int>::max()) ||
       c.maximumNodes>std::vector<Key>().max_size() || c.maximumNodes>indexMax/3 ||
       c.maximumNodes>bytesMax/(3*sizeof(double)) ||
       !c.maximumEdges || c.maximumEdges>std::vector<Edge>().max_size() ||
       c.maximumEdges>std::set<Edge>().max_size() ||
       !c.maximumInterpolationValues || c.maximumInterpolationValues>indexMax ||
       c.maximumInterpolationValues>bytesMax/sizeof(double))
        throw std::length_error("RESOURCE_LIMIT: invalid source cage grid cell/node/edge/interpolation budget");
    // Binary searches below rely on strictly increasing represented planes.
    // Check finest spacing against the largest root-endpoint ULP, as the
    // upstream oracle does, but at the requested finer level.
    const double step=std::ldexp(root.rootSide,-c.fineLevel);
    for(int a=0;a<3;++a) {
        const double magnitude=std::max(std::abs(plane(root,0,a,0)),std::abs(plane(root,0,a,1)));
        const double spacing=std::nextafter(magnitude,std::numeric_limits<double>::infinity())-magnitude;
        if(!std::isfinite(step) || step<=0 || step<spacing)
            throw std::invalid_argument("source cage grid finer coordinates numerically collapse");
    }
}
} // namespace

struct SourceCageGrid::Data {
    SourceSurfaceCage cage;
    SourceCageGridConfig config;
    std::vector<Key> cells,keys;
    Eigen::MatrixXd positions;
    std::vector<CellNodes> connectivity;
    std::vector<Edge> edges;
    std::vector<int> handles;
    std::string configHash,outputHash;

    Data(const SourceSurfaceCage& c,const SourceCageGridConfig& f):cage(c),config(f) {
        const auto& root=cage.config(); validate(root,config);
        const int delta=config.fineLevel-root.cageLevel;
        const std::uint64_t perCoarse=std::uint64_t(1)<<(3*delta); // <= 2^60
        if(cage.selectedCells().size()>config.maximumCells/perCoarse)
            throw std::length_error("RESOURCE_LIMIT: fine cell budget exceeded before expansion");
        const auto count=cage.selectedCells().size()*static_cast<std::size_t>(perCoarse);
        cells.reserve(count);
        const auto factor=std::uint32_t(1)<<delta;
        for(const auto& coarse:cage.selectedCells()) {
            Key base; for(int a=0;a<3;++a) base[a]=coarse[a]*factor;
            for(std::uint32_t x=0;x<factor;++x) for(std::uint32_t y=0;y<factor;++y) for(std::uint32_t z=0;z<factor;++z)
                cells.push_back({base[0]+x,base[1]+y,base[2]+z});
        }
        std::sort(cells.begin(),cells.end());
        // Unique containers keep temporary storage proportional to actual
        // nodes/edges. No unchecked 8*Ncell or 12*Ncell intermediate vectors.
        std::map<Key,int> nodeIndex;
        for(const auto& cell:cells) for(int b=0;b<8;++b) {
            const auto key=cornerKey(cell,b);
            if(nodeIndex.find(key)!=nodeIndex.end()) continue;
            if(nodeIndex.size()>=config.maximumNodes)
                throw std::length_error("RESOURCE_LIMIT: fine node budget exceeded");
            nodeIndex.emplace(key,0);
        }
        keys.reserve(nodeIndex.size()); positions.resize(static_cast<Eigen::Index>(nodeIndex.size()),3);
        for(auto& entry:nodeIndex) {
            entry.second=static_cast<int>(keys.size()); keys.push_back(entry.first);
            for(int a=0;a<3;++a) positions(entry.second,a)=plane(root,config.fineLevel,a,entry.first[a]);
        }
        connectivity.reserve(cells.size());
        std::set<Edge> uniqueEdges;
        for(const auto& cell:cells) {
            CellNodes row; for(int b=0;b<8;++b) row[b]=nodeIndex.at(cornerKey(cell,b));
            connectivity.push_back(row);
            for(int b=0;b<8;++b) for(int a=0;a<3;++a) if(!(b&(1<<a))) {
                const Edge edge{row[b],row[b|(1<<a)]};
                if(uniqueEdges.find(edge)!=uniqueEdges.end()) continue;
                if(uniqueEdges.size()>=config.maximumEdges)
                    throw std::length_error("RESOURCE_LIMIT: fine edge budget exceeded");
                uniqueEdges.insert(edge);
            }
        }
        edges.assign(uniqueEdges.begin(),uniqueEdges.end());
        handles.reserve(cage.handleVertexIndices().size());
        for(int handle:cage.handleVertexIndices()) {
            Key key=cage.boundaryVertexKeys()[handle]; for(auto& x:key) x*=factor;
            handles.push_back(nodeIndex.at(key));
        }
        Digest configuration;
        configuration.string("source-cage-grid-selected-coarse-union-refinement-IMPLEMENTATION_DEFAULT-closed-represented-lex-ties-v1");
        configuration.string(cage.configSha256());
        configuration.integer(config.fineLevel); configuration.integer(config.maximumCells);
        configuration.integer(config.maximumNodes); configuration.integer(config.maximumEdges);
        configuration.integer(config.maximumInterpolationValues); configHash=configuration.finish();
        Digest output;
        output.string("source-cage-grid-lex-keys-FMA-positions-x-bit-first-axis-edges-complete-handles-v1");
        output.string(cage.sourceSha256()); output.string(cage.configSha256());
        output.string(cage.outputSha256()); output.string(configHash);
        output.integer(cells.size()); for(const auto& cell:cells) output.array(cell);
        output.integer(keys.size()); for(const auto& key:keys) output.array(key);
        output.integer(positions.rows());
        for(Eigen::Index i=0;i<positions.rows();++i) for(int a=0;a<3;++a) output.real(positions(i,a));
        output.integer(connectivity.size()); for(const auto& row:connectivity) output.array(row);
        output.integer(edges.size()); for(const auto& edge:edges) output.array(edge);
        output.integer(handles.size()); for(auto handle:handles) output.integer(handle);
        outputHash=output.finish();
    }
};

SourceCageGrid::SourceCageGrid(const SourceSurfaceCage& c,const SourceCageGridConfig& f)
    :data_(std::make_shared<const Data>(c,f)) {}
const SourceSurfaceCage& SourceCageGrid::cage() const {return data_->cage;}
const SourceCageGridConfig& SourceCageGrid::config() const {return data_->config;}
const std::vector<SourceCageGrid::Key>& SourceCageGrid::cells() const {return data_->cells;}
const std::vector<SourceCageGrid::Key>& SourceCageGrid::nodeKeys() const {return data_->keys;}
const Eigen::MatrixXd& SourceCageGrid::nodePositions() const {return data_->positions;}
const std::vector<SourceCageGrid::CellNodes>& SourceCageGrid::cellNodes() const {return data_->connectivity;}
const std::vector<SourceCageGrid::Edge>& SourceCageGrid::edges() const {return data_->edges;}
const std::vector<int>& SourceCageGrid::handleNodeIndices() const {return data_->handles;}
const std::string& SourceCageGrid::sourceSha256() const {return cage().sourceSha256();}
const std::string& SourceCageGrid::cageConfigSha256() const {return cage().configSha256();}
const std::string& SourceCageGrid::cageOutputSha256() const {return cage().outputSha256();}
const std::string& SourceCageGrid::configSha256() const {return data_->configHash;}
const std::string& SourceCageGrid::outputSha256() const {return data_->outputHash;}

SourceCageGrid::Embedding SourceCageGrid::embed(const Eigen::Vector3d& point) const {
    if(!point.allFinite()) throw std::invalid_argument("source cage grid embedding requires a finite point");
    const auto& root=cage().config(); const int level=config().fineLevel;
    const auto side=std::uint32_t(1)<<level;
    std::array<std::array<std::uint32_t,2>,3> candidates{};
    std::array<int,3> counts{};
    for(int a=0;a<3;++a) {
        if(point(a)<plane(root,level,a,0) || point(a)>plane(root,level,a,side))
            throw std::out_of_range("point outside source cage grid domain");
        // First represented plane >= point; inverse floating lattice indices
        // can round across non-dyadic planes and cannot decide closed ties.
        std::uint32_t lo=0,hi=side;
        while(lo<hi) {
            const auto mid=lo+(hi-lo)/2;
            if(plane(root,level,a,mid)<point(a)) lo=mid+1; else hi=mid;
        }
        if(lo>0) candidates[a][counts[a]++]=lo-1;
        if(lo<side && plane(root,level,a,lo)==point(a)) candidates[a][counts[a]++]=lo;
    }
    for(int x=0;x<counts[0];++x) for(int y=0;y<counts[1];++y) for(int z=0;z<counts[2];++z) {
        Key key{candidates[0][x],candidates[1][y],candidates[2][z]};
        const auto it=std::lower_bound(cells().begin(),cells().end(),key);
        if(it==cells().end() || *it!=key) continue;
        Embedding e; e.cellIndex=static_cast<std::size_t>(it-cells().begin());
        e.cellKey=key; e.nodeIndices=cellNodes()[e.cellIndex];
        std::array<double,3> t;
        for(int a=0;a<3;++a) {
            const double lower=plane(root,level,a,key[a]),upper=plane(root,level,a,key[a]+1);
            t[a]=static_cast<double>((static_cast<long double>(point(a))-lower)/(static_cast<long double>(upper)-lower));
        }
        for(int b=0;b<8;++b) {
            double w=1;
            for(int a=0;a<3;++a) w*=((b>>a)&1)?t[a]:1-t[a];
            e.coefficients[b]=w;
        }
        return e;
    }
    throw std::out_of_range("point outside source cage grid selected-cell union (including holes)");
}

Eigen::MatrixXd SourceCageGrid::interpolate(const Eigen::MatrixXd& nodalField,const Eigen::MatrixXd& queries) const {
    if(nodalField.rows()!=nodePositions().rows() || nodalField.cols()<=0 || queries.cols()!=3)
        throw std::invalid_argument("source cage grid interpolation requires Nnodes x K (K>0) and Nqueries x 3");
    const auto columns=static_cast<std::size_t>(nodalField.cols());
    const auto rows=static_cast<std::size_t>(queries.rows());
    if(rows>config().maximumInterpolationValues/columns)
        throw std::length_error("RESOURCE_LIMIT: interpolation output value budget exceeded before allocation");
    if(!nodalField.allFinite() || !queries.allFinite())
        throw std::invalid_argument("source cage grid interpolation requires finite field and query values");
    Eigen::MatrixXd output(queries.rows(),nodalField.cols());
    for(Eigen::Index i=0;i<queries.rows();++i) {
        const auto e=embed(queries.row(i).transpose());
        for(Eigen::Index k=0;k<nodalField.cols();++k) {
            long double value=0;
            for(int b=0;b<8;++b) value+=static_cast<long double>(e.coefficients[b])*nodalField(e.nodeIndices[b],k);
            output(i,k)=static_cast<double>(value);
            if(!std::isfinite(output(i,k))) throw std::overflow_error("source cage grid interpolated field overflowed");
        }
    }
    return output;
}
} // namespace spin2014
