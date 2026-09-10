#include "spin2014/SourceSurfaceCage.h"
#include "spin2014/SourceSurfaceOctree.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
namespace spin2014 {
namespace {
using Key=SourceSurfaceCage::Key;
// Outward cyclic corners: -x,+x,-y,+y,-z,+z, x bit first.
// Every quad uses the diagonal corner[0]--corner[2].
constexpr int faceCorners[6][4]={{0,4,6,2},{1,3,7,5},{0,1,5,4},
                               {2,6,7,3},{0,2,3,1},{4,5,7,6}};
class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {
        unsigned char bytes[8]; for(int i=0;i<8;++i) bytes[7-i]=static_cast<unsigned char>(x>>(8*i));
        hash_.update(bytes,8);
    }
    void string(const std::string& x) {
        integer(x.size()); hash_.update(reinterpret_cast<const unsigned char*>(x.data()),x.size());
    }
    void real(double x) {
        static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,"IEEE754 binary64 required");
        std::uint64_t bits; std::memcpy(&bits,&x,8); integer(bits);
    }
    void key(const Key& k) {for(auto x:k) integer(x);}
    std::string finish() {return hash_.finish();}
};
void validateBudgets(const SourceSurfaceCageConfig& c) {
    const auto matrixLimit=static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()/3);
    if(c.maximumSelectedCells==0 || c.maximumSelectedCells>std::vector<Key>().max_size() ||
       c.maximumBoundaryVertices==0 || c.maximumBoundaryVertices>std::size_t(std::numeric_limits<int>::max()) ||
       c.maximumBoundaryVertices>matrixLimit || c.maximumBoundaryVertices>std::vector<Key>().max_size() ||
       c.maximumBoundaryVertices>std::numeric_limits<std::size_t>::max()/(3*sizeof(double)) ||
       c.maximumBoundaryTriangles==0 || c.maximumBoundaryTriangles>matrixLimit ||
       c.maximumBoundaryTriangles>std::size_t(std::numeric_limits<int>::max()/2) ||
       c.maximumBoundaryTriangles>std::numeric_limits<std::size_t>::max()/(3*sizeof(int)))
        throw std::length_error("RESOURCE_LIMIT: invalid source cage cell/boundary budget");
}
Key cornerKey(Key cell,int corner) {
    for(int a=0;a<3;++a) cell[a]+=(corner>>a)&1;
    return cell;
}
// Work scales with selected cells, never the root lattice's empty volume.
template<class Emit> void boundaryQuads(const std::vector<Key>& cells,std::uint32_t side,Emit emit) {
    for(const auto& cell:cells) for(int face=0;face<6;++face) {
        const int axis=face/2; const bool positive=face%2;
        Key neighbor=cell;
        const bool outside=positive ? cell[axis]+1==side : cell[axis]==0;
        if(!outside) {
            if(positive) ++neighbor[axis]; else --neighbor[axis];
            if(std::binary_search(cells.begin(),cells.end(),neighbor)) continue;
        }
        std::array<Key,4> quad;
        for(int j=0;j<4;++j) quad[j]=cornerKey(cell,faceCorners[face][j]);
        emit(quad);
    }
}
} // namespace
struct SourceSurfaceCage::Data {
    Eigen::MatrixXd sourceV,vertices;
    Eigen::MatrixXi sourceF,faces;
    SourceSurfaceCageConfig config;
    std::vector<Key> cells,vertexKeys;
    std::vector<int> handles;
    verification::ExactBoundaryAudit audit;
    std::size_t nodes=0;
    std::string sourceHash,configHash,outputHash;
    Data(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f,const SourceSurfaceCageConfig& c)
        :sourceV(v),sourceF(f),config(c) {
        validateBudgets(config);
        SourceSurfaceOctreeConfig oracleConfig;
        oracleConfig.rootOrigin=config.rootOrigin; oracleConfig.rootSide=config.rootSide;
        oracleConfig.maximumDepth=config.cageLevel; oracleConfig.maximumNodes=config.maximumSourceNodes;
        // The original signed source is BOTH outer and design. No carved design,
        // prior lattice or occupancy may substitute for original object volume.
        const SourceSurfaceOctree oracle(sourceV,sourceF,sourceV,sourceF,oracleConfig);
        nodes=oracle.cells().size(); sourceHash=oracle.sourceSha256();
        std::size_t selectedCount=0;
        for(const auto& node:oracle.cells()) {
            std::uint64_t count=0;
            if(node.kind==SourceSurfaceCellKind::Boundary) count=1;
            if(node.kind==SourceSurfaceCellKind::Full)
                count=std::uint64_t(1)<<(3*(config.cageLevel-node.level)); // depth<=20: at most 2^60
            if(count>config.maximumSelectedCells-selectedCount)
                throw std::length_error("RESOURCE_LIMIT: source cage selected-cell budget exceeded before Full expansion");
            selectedCount+=static_cast<std::size_t>(count);
        }
        cells.reserve(selectedCount);
        for(const auto& node:oracle.cells()) {
            if(node.kind==SourceSurfaceCellKind::Boundary) cells.push_back(node.coordinates);
            if(node.kind!=SourceSurfaceCellKind::Full) continue;
            const auto side=std::uint32_t(1)<<(config.cageLevel-node.level);
            Key base; for(int a=0;a<3;++a) base[a]=node.coordinates[a]*side;
            for(std::uint32_t x=0;x<side;++x) for(std::uint32_t y=0;y<side;++y) for(std::uint32_t z=0;z<side;++z)
                cells.push_back({base[0]+x,base[1]+y,base[2]+z});
        }
        std::sort(cells.begin(),cells.end());
        const auto side=std::uint32_t(1)<<config.cageLevel;
        std::map<Key,int> vertexIndices;
        std::size_t triangleCount=0;
        boundaryQuads(cells,side,[&](const std::array<Key,4>& quad) {
            if(config.maximumBoundaryTriangles-triangleCount<2)
                throw std::length_error("RESOURCE_LIMIT: source cage boundary triangle budget exceeded");
            triangleCount+=2;
            for(const auto& key:quad) if(vertexIndices.find(key)==vertexIndices.end()) {
                if(vertexIndices.size()>=config.maximumBoundaryVertices)
                    throw std::length_error("RESOURCE_LIMIT: source cage boundary vertex/handle budget exceeded");
                vertexIndices.emplace(key,0);
            }
        });
        vertexKeys.reserve(vertexIndices.size());
        vertices.resize(vertexIndices.size(),3); handles.resize(vertexIndices.size());
        for(auto& entry:vertexIndices) {
            entry.second=static_cast<int>(vertexKeys.size()); vertexKeys.push_back(entry.first);
            Key cell=entry.first;
            for(auto& x:cell) x=std::min(x,side-1);
            const auto bounds=oracle.cellBounds(config.cageLevel,cell);
            for(int a=0;a<3;++a)
                vertices(entry.second,a)=entry.first[a]==side?bounds.upper(a):bounds.lower(a);
        }
        std::iota(handles.begin(),handles.end(),0);
        faces.resize(triangleCount,3); Eigen::Index row=0;
        boundaryQuads(cells,side,[&](const std::array<Key,4>& quad) {
            const int a=vertexIndices.at(quad[0]),b=vertexIndices.at(quad[1]),
                      c=vertexIndices.at(quad[2]),d=vertexIndices.at(quad[3]);
            faces.row(row++)<<a,b,c; faces.row(row++)<<a,c,d;
        });
        // AnyNonzero retains negative cavity shells. Their material-relative
        // direction follows the exposed-cell face winding above, not repair.
        audit=verification::auditExactTriangleBoundary(vertices,faces,verification::ExactBoundaryOrientation::AnyNonzero);
        if(!audit.facts().boundaryValid) throw SourceSurfaceCageGeometryError(audit);
        Digest configuration;
        configuration.string("source-surface-cage-CLOSED-contact-IMPLEMENTATION_DEFAULT-all-boundary-components-v1");
        configuration.string(oracle.configSha256());
        for(int a=0;a<3;++a) configuration.real(config.rootOrigin(a));
        configuration.real(config.rootSide); configuration.integer(config.cageLevel);
        configuration.integer(config.maximumSourceNodes); configuration.integer(config.maximumSelectedCells);
        configuration.integer(config.maximumBoundaryVertices); configuration.integer(config.maximumBoundaryTriangles);
        configHash=configuration.finish();
        Digest output; output.string("source-surface-cage-lexicographic-keys-x-bit-first-outward-fixed-diagonal-v1");
        output.string(sourceHash); output.string(configHash);
        output.integer(cells.size()); for(const auto& key:cells) output.key(key);
        output.integer(vertexKeys.size()); for(const auto& key:vertexKeys) output.key(key);
        output.integer(vertices.rows());
        for(Eigen::Index i=0;i<vertices.rows();++i) for(int a=0;a<3;++a) output.real(vertices(i,a));
        output.integer(faces.rows());
        for(Eigen::Index i=0;i<faces.rows();++i) for(int a=0;a<3;++a) output.integer(faces(i,a));
        output.integer(handles.size()); for(auto index:handles) output.integer(index);
        outputHash=output.finish();
    }
};
SourceSurfaceCageGeometryError::SourceSurfaceCageGeometryError(verification::ExactBoundaryAudit audit)
    :std::runtime_error("INVALID_GEOMETRY: source cage boundary failed exact manifold/orientation/intersection audit"),audit_(std::move(audit)) {}
SourceSurfaceCage::SourceSurfaceCage(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f,const SourceSurfaceCageConfig& c)
    :data_(std::make_shared<const Data>(v,f,c)) {}
const Eigen::MatrixXd& SourceSurfaceCage::sourceVertices() const {return data_->sourceV;}
const Eigen::MatrixXi& SourceSurfaceCage::sourceFaces() const {return data_->sourceF;}
const SourceSurfaceCageConfig& SourceSurfaceCage::config() const {return data_->config;}
const std::vector<SourceSurfaceCage::Key>& SourceSurfaceCage::selectedCells() const {return data_->cells;}
const std::vector<SourceSurfaceCage::Key>& SourceSurfaceCage::boundaryVertexKeys() const {return data_->vertexKeys;}
const Eigen::MatrixXd& SourceSurfaceCage::boundaryVertices() const {return data_->vertices;}
const Eigen::MatrixXi& SourceSurfaceCage::boundaryFaces() const {return data_->faces;}
const std::vector<int>& SourceSurfaceCage::handleVertexIndices() const {return data_->handles;}
const verification::ExactBoundaryAudit& SourceSurfaceCage::boundaryAudit() const {return data_->audit;}
std::size_t SourceSurfaceCage::sourceNodeCount() const {return data_->nodes;}
const std::string& SourceSurfaceCage::sourceSha256() const {return data_->sourceHash;}
const std::string& SourceSurfaceCage::configSha256() const {return data_->configHash;}
const std::string& SourceSurfaceCage::outputSha256() const {return data_->outputHash;}
} // namespace spin2014
