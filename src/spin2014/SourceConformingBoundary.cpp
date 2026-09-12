#include "spin2014/SourceConformingBoundary.h"
#include "verification/Sha256.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <unordered_map>
#include <new>

namespace spin2014 {
namespace {
using Boundary=SourceConformingBoundary;
using Key=Boundary::Key;
using Patch=Boundary::Patch;
using Triangle=Boundary::Triangle;
using Reference=Boundary::OrientedTriangle;
constexpr auto absent=Boundary::noOwner;
[[noreturn]] void limit(const char* what) {throw std::length_error(std::string("RESOURCE_LIMIT: source conforming boundary ")+what);}
struct Node {
    std::array<std::size_t,4> children{{absent,absent,absent,absent}};
    std::array<std::size_t,2> owners{{absent,absent}};
};
using PlaneRoots=std::unordered_map<std::uint64_t,std::size_t>;
using PointMap=std::map<Key,std::size_t>;
using Lines=std::map<Key,std::vector<std::uint32_t>>;
using Edges=std::map<std::array<std::size_t,2>,std::pair<unsigned,int>>;
class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=static_cast<unsigned char>(x>>(8*i));hash_.update(b,8);}
    void string(const std::string& s) {integer(s.size());hash_.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());}
    void real(double x) {static_assert(sizeof(double)==8&&std::numeric_limits<double>::is_iec559,"IEEE754 required");std::uint64_t b;std::memcpy(&b,&x,8);integer(b);}
    template<class T,std::size_t N> void array(const std::array<T,N>& a){for(auto x:a)integer(x);}
    std::string finish(){return hash_.finish();}
};
void validate(const SourceCageBinding& b,const SourceConformingBoundaryConfig& c) {
    const auto bytes=std::numeric_limits<std::size_t>::max();
    const auto indices=std::size_t(std::numeric_limits<Eigen::Index>::max());
    // Bound aggregate derived entry counts as well as each container's bytes.
    if(!c.maximumFaceTreeNodes||c.maximumFaceTreeNodes>std::vector<Node>().max_size()||c.maximumFaceTreeNodes>PlaneRoots().max_size()||
       !c.maximumPatches||c.maximumPatches>std::vector<Patch>().max_size()||c.maximumPatches>bytes/(5*sizeof(Key))||
       !c.maximumPoints||c.maximumPoints>std::size_t(std::numeric_limits<int>::max())||c.maximumPoints>indices/3||
       c.maximumPoints>bytes/(3*sizeof(double))||c.maximumPoints>std::vector<Key>().max_size()||c.maximumPoints>PointMap().max_size()||
       c.maximumPoints>Lines().max_size()/3||c.maximumPoints>std::vector<std::uint32_t>().max_size()/3||
       !c.maximumTriangles||c.maximumTriangles>std::vector<Triangle>().max_size()||c.maximumTriangles>indices/3||
       c.maximumTriangles>std::vector<Reference>().max_size()/2||c.maximumTriangles>Edges().max_size()/6||
       b.primitives().size()>std::vector<std::vector<Reference>>().max_size()||b.primitives().size()>bytes/(6*sizeof(std::uint64_t)))
        limit("invalid/overflow budget");
    if(b.outerPointCount()>c.maximumPoints)limit("original point budget");
    const auto& root=b.source().config();
    if(root.maximumDepth<0||root.maximumDepth>20)throw std::invalid_argument("source conforming boundary requires source depth <=20");
    const double step=std::ldexp(root.rootSide,-root.maximumDepth-1);
    for(int a=0;a<3;++a) {
        const double endpoint=std::fma(1.,root.rootSide,root.rootOrigin[a]);
        const double magnitude=std::max(std::abs(root.rootOrigin[a]),std::abs(endpoint));
        const double ulp=std::nextafter(magnitude,std::numeric_limits<double>::infinity())-magnitude;
        if(!std::isfinite(endpoint)||!std::isfinite(step))throw std::overflow_error("source conforming boundary nonfinite root arithmetic");
        if(!(step>ulp))throw std::invalid_argument("source conforming boundary D+1 spacing is numerically unsupported");
    }
}
std::size_t addNode(std::vector<Node>& nodes,std::size_t cap) {
    if(nodes.size()>=cap)limit("face tree node insertion budget");
    nodes.emplace_back();return nodes.size()-1;
}
std::size_t combine(std::size_t inherited,std::size_t local) {
    if(inherited!=absent&&local!=absent&&inherited!=local)throw std::invalid_argument("source conforming boundary same-side primitive overlap");
    return local==absent?inherited:local;
}
void emitPatches(const std::vector<Node>& nodes,std::size_t node,int axis,std::uint32_t plane,int level,
                 std::array<std::uint32_t,2> uv,std::array<std::size_t,2> owners,
                 std::vector<Patch>& patches,std::size_t cap) {
    bool descend=false;
    if(node!=absent) {
        for(int s=0;s<2;++s)owners[s]=combine(owners[s],nodes[node].owners[s]);
        for(auto child:nodes[node].children)descend|=child!=absent;
    }
    if(!descend) {
        if(owners[0]==absent&&owners[1]==absent)return;
        if(owners[0]!=absent&&owners[0]==owners[1])throw std::invalid_argument("source conforming boundary same primitive on both sides");
        if(patches.size()>=cap)limit("patch insertion budget");
        patches.push_back({axis,plane,level,uv,owners[0],owners[1],0,0});return;
    }
    for(int q=0;q<4;++q)emitPatches(nodes,nodes[node].children[q],axis,plane,level+1,
        {2*uv[0]+unsigned(q&1),2*uv[1]+unsigned((q>>1)&1)},owners,patches,cap);
}
std::array<Key,5> patchPoints(const Patch& p,int latticeLevel) {
    std::array<Key,5> result{};const auto side=1u<<(latticeLevel-p.level);
    const int u=(p.axis+1)%3,v=(p.axis+2)%3;
    for(int k=0;k<5;++k) {
        auto& key=result[k];key[p.axis]=2*p.plane;
        key[u]=p.coordinates[0]*side;key[v]=p.coordinates[1]*side;
        if(k==1||k==2)key[u]+=side;if(k==2||k==3)key[v]+=side;
        if(k==4){key[u]+=side/2;key[v]+=side/2;}
    }
    return result;
}
Key lineKey(int axis,const Key& point) {return {unsigned(axis),point[(axis+1)%3],point[(axis+2)%3]};}
} // namespace

struct SourceConformingBoundary::Data {
    SourceCageBinding binding;
    SourceConformingBoundaryConfig config;
    Eigen::MatrixXd points;
    std::vector<Key> keys;
    std::vector<Patch> patches;
    std::vector<Triangle> triangles;
    std::vector<std::vector<Reference>> boundaries;
    int level;
    std::string configHash,outputHash;
    Data(const SourceCageBinding& b,const SourceConformingBoundaryConfig& c):binding(b),config(c),level(b.source().config().maximumDepth+1) {
        const int depth=level-1;
        std::vector<Node> nodes;PlaneRoots roots;
        // Expected O(D*faces) insertion. Only actual paths allocate nodes;
        // ordered plane traversal below makes hash-table iteration irrelevant.
        for(std::size_t owner=0;owner<binding.primitives().size();++owner) {
            const auto& p=binding.primitives()[owner].sourceKey;
            if(p.level<0||p.level>depth)throw std::invalid_argument("source conforming boundary invalid primitive depth");
            for(auto x:p.coordinates)if(x>=(1u<<p.level))throw std::invalid_argument("source conforming boundary invalid primitive coordinates");
            for(int axis=0;axis<3;++axis)for(int upper=0;upper<2;++upper) {
                const auto plane=(p.coordinates[axis]+upper)<<(depth-p.level);
                const auto planeId=(std::uint64_t(axis)<<32)|plane;
                auto found=roots.find(planeId);
                if(found==roots.end()){const auto root=addNode(nodes,config.maximumFaceTreeNodes);found=roots.emplace(planeId,root).first;}
                auto node=found->second;const int u=(axis+1)%3,v=(axis+2)%3;
                for(int l=1;l<=p.level;++l) {
                    const int shift=p.level-l;const auto q=((p.coordinates[u]>>shift)&1)|(((p.coordinates[v]>>shift)&1)<<1);
                    if(nodes[node].children[q]==absent){const auto child=addNode(nodes,config.maximumFaceTreeNodes);nodes[node].children[q]=child;}
                    node=nodes[node].children[q];
                }
                auto& owned=nodes[node].owners[upper?0:1];owned=combine(owned,owner);
            }
        }
        std::vector<std::uint64_t> planes;planes.reserve(roots.size());for(const auto& r:roots)planes.push_back(r.first);
        std::sort(planes.begin(),planes.end());
        for(auto p:planes)emitPatches(nodes,roots.at(p),int(p>>32),std::uint32_t(p),0,{0,0},{absent,absent},patches,config.maximumPatches);
        PointMap index;const auto outer=binding.outerPointCount();
        for(const auto& p:patches)for(const auto& key:patchPoints(p,level)) {
            if(index.find(key)!=index.end())continue;
            if(index.size()>=config.maximumPoints-outer)limit("material point insertion budget");
            index.emplace(key,0);
        }
        // Primitive corners must survive even when face overlays reindex them.
        for(auto key:binding.materialCornerKeys()) {
            for(auto& x:key)x*=2;
            if(index.find(key)==index.end())throw std::invalid_argument("source conforming boundary omitted primitive corner");
        }
        points.resize(static_cast<Eigen::Index>(outer+index.size()),3);points.topRows(outer)=binding.source().outerVertices();keys.reserve(index.size());
        const auto& root=binding.source().config();
        Lines lines;
        for(auto& entry:index) {
            entry.second=outer+keys.size();keys.push_back(entry.first);
            for(int axis=0;axis<3;++axis) {
                const double value=std::fma(std::ldexp(double(entry.first[axis]),-level),root.rootSide,root.rootOrigin[axis]);
                if(!std::isfinite(value))throw std::overflow_error("source conforming boundary nonfinite point");
                points(entry.second,axis)=value;
                // <=3*maximumPoints entries, validated before any insertion.
                lines[lineKey(axis,entry.first)].push_back(entry.first[axis]);
            }
        }
        for(auto& line:lines)std::sort(line.second.begin(),line.second.end());
        boundaries.resize(binding.primitives().size());
        std::size_t references=0;
        for(auto& patch:patches) {
            const auto corners=patchPoints(patch,level);std::vector<std::size_t> polygon;
            // Query ALL-plane points on every edge, with no point/edge scan.
            for(int e=0;e<4;++e) {
                const auto& from=corners[e];const auto& to=corners[(e+1)%4];
                int axis=0;while(from[axis]==to[axis])++axis;
                const auto& line=lines.at(lineKey(axis,from));
                const auto lo=std::lower_bound(line.begin(),line.end(),std::min(from[axis],to[axis]));
                const auto hi=std::upper_bound(line.begin(),line.end(),std::max(from[axis],to[axis]));
                auto append=[&](std::uint32_t coordinate) {
                    if(coordinate==to[axis])return;
                    if(polygon.size()>=config.maximumTriangles-triangles.size())limit("triangle perimeter budget");
                    auto key=from;key[axis]=coordinate;polygon.push_back(index.at(key));
                };
                if(from[axis]<to[axis])for(auto it=lo;it!=hi;++it)append(*it);
                else for(auto it=hi;it!=lo;){--it;append(*it);}
            }
            if(polygon.size()<4)throw std::invalid_argument("source conforming boundary invalid patch polygon");
            patch.triangleBegin=triangles.size();patch.triangleCount=polygon.size();const auto center=index.at(corners[4]);
            const auto ownerCount=std::size_t(patch.negativeOwner!=absent)+std::size_t(patch.positiveOwner!=absent);
            if(polygon.size()>(2*config.maximumTriangles-references)/ownerCount)limit("derived oriented reference budget");
            references+=ownerCount*polygon.size();
            for(std::size_t k=0;k<polygon.size();++k) {
                const Triangle t{center,polygon[k],polygon[(k+1)%polygon.size()]};
                const Eigen::Vector3d a=(points.row(t[1])-points.row(t[0])).transpose(),d=(points.row(t[2])-points.row(t[0])).transpose();
                const Eigen::Vector3d cross=a.cross(d);
                if(!cross.allFinite())throw std::overflow_error("source conforming boundary nonfinite triangle arithmetic");
                if(t[0]==t[1]||t[0]==t[2]||t[1]==t[2]||!(cross[patch.axis]>0))throw std::invalid_argument("source conforming boundary degenerate/reversed rest triangle");
                const auto id=triangles.size();triangles.push_back(t);
                if(patch.negativeOwner!=absent)boundaries.at(patch.negativeOwner).push_back({id,1});
                if(patch.positiveOwner!=absent)boundaries.at(patch.positiveOwner).push_back({id,-1});
            }
        }
        // Fail closed if any perpendicular hanging point or owner convention
        // was lost. Per-primitive temporary edges <=6*maximumTriangles.
        for(const auto& boundary:boundaries) {
            Edges edges;
            for(const auto& ref:boundary) {
                auto t=triangles[ref.triangle];if(ref.orientation<0)std::swap(t[1],t[2]);
                for(int e=0;e<3;++e) {
                    auto a=t[e],b=t[(e+1)%3];const int sign=a<b?1:-1;if(a>b)std::swap(a,b);
                    auto& count=edges[{a,b}];if(++count.first>2)throw std::invalid_argument("source conforming boundary nonmanifold primitive edge");count.second+=sign;
                }
            }
            if(edges.empty())throw std::invalid_argument("source conforming boundary empty primitive boundary");
            for(const auto& e:edges)if(e.second.first!=2||e.second.second!=0)throw std::invalid_argument("source conforming boundary open/inconsistently oriented primitive");
        }
        Digest cfg;cfg.string("source-conforming-boundary-config-v1:IMPLEMENTATION_DEFAULT:sparse-dyadic-face-overlay:all-plane-edge-points:shared-center-fan:cyclic-uv:+axis:outer-first:sorted-D+1-keys:plane-sorted-Morton-patches");
        cfg.string(binding.outputSha256());cfg.integer(level);
        cfg.integer(config.maximumFaceTreeNodes);cfg.integer(config.maximumPatches);cfg.integer(config.maximumPoints);cfg.integer(config.maximumTriangles);configHash=cfg.finish();
        Digest out;out.string("source-conforming-boundary-output-v1");out.string(configHash);out.string(binding.outputSha256());out.integer(level);out.integer(outer);
        const auto& faces=binding.source().outerFaces();out.integer(faces.rows());for(Eigen::Index r=0;r<faces.rows();++r)for(int a=0;a<3;++a)out.integer(faces(r,a));
        out.integer(points.rows());for(Eigen::Index r=0;r<points.rows();++r)for(int a=0;a<3;++a)out.real(points(r,a));
        out.integer(keys.size());for(const auto& k:keys)out.array(k);
        out.integer(patches.size());for(const auto& p:patches){out.integer(p.axis);out.integer(p.plane);out.integer(p.level);out.array(p.coordinates);out.integer(p.negativeOwner);out.integer(p.positiveOwner);out.integer(p.triangleBegin);out.integer(p.triangleCount);}
        out.integer(triangles.size());for(const auto& t:triangles)out.array(t);
        out.integer(boundaries.size());for(const auto& b:boundaries){out.integer(b.size());for(const auto& r:b){out.integer(r.triangle);out.integer(r.orientation==1?1:0);}}
        outputHash=out.finish();
    }
};
SourceConformingBoundary::SourceConformingBoundary(const SourceCageBinding& b,const SourceConformingBoundaryConfig& c) {
    try{validate(b,c);data_=std::make_shared<const Data>(b,c);}catch(const std::bad_alloc&){limit("allocation failed");}
}
const SourceCageBinding& SourceConformingBoundary::binding() const{return data_->binding;}
const SourceConformingBoundaryConfig& SourceConformingBoundary::config() const{return data_->config;}
const Eigen::MatrixXd& SourceConformingBoundary::points() const{return data_->points;}
const std::vector<SourceConformingBoundary::Key>& SourceConformingBoundary::materialPointKeys() const{return data_->keys;}
const std::vector<SourceConformingBoundary::Patch>& SourceConformingBoundary::patches() const{return data_->patches;}
const std::vector<SourceConformingBoundary::Triangle>& SourceConformingBoundary::triangles() const{return data_->triangles;}
const std::vector<std::vector<SourceConformingBoundary::OrientedTriangle>>& SourceConformingBoundary::primitiveBoundaries() const{return data_->boundaries;}
std::size_t SourceConformingBoundary::outerPointCount() const{return binding().outerPointCount();}
int SourceConformingBoundary::latticeLevel() const{return data_->level;}
const std::string& SourceConformingBoundary::configSha256() const{return data_->configHash;}
const std::string& SourceConformingBoundary::outputSha256() const{return data_->outputHash;}
} // namespace spin2014
