#include "spin2014/SourceSurfaceOctree.h"
#include "spin2014/SurfaceRawMomentDifferential.h"
#include "verification/ExactMaterialRegion.h"
#include "verification/Sha256.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/intersections.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <stdexcept>

namespace spin2014 {
namespace {
using Kernel=CGAL::Simple_cartesian<CGAL::Gmpq>;
using Point=Kernel::Point_3;
using Triangle=Kernel::Triangle_3;
using Mesh=CGAL::Surface_mesh<Point>;
struct Bounds { Eigen::Vector3d lo,hi; };
struct OriginalTriangle { Triangle triangle; Bounds bounds; };

// Fixed big-endian field encoding, independent of matrix storage order/padding.
class Digest {
    verification::Sha256 hash;
public:
    void string(const std::string& value) {
        integer(value.size()); hash.update(reinterpret_cast<const unsigned char*>(value.data()),value.size());
    }
    void integer(std::uint64_t value) {
        unsigned char bytes[8]; for(int i=0;i<8;++i) bytes[7-i]=static_cast<unsigned char>(value>>(8*i));
        hash.update(bytes,8);
    }
    void real(double value) {
        static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,"IEEE754 binary64 required");
        std::uint64_t bits; std::memcpy(&bits,&value,8); integer(bits);
    }
    void surface(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f) {
        integer(v.rows()); integer(f.rows());
        for(Eigen::Index i=0;i<v.rows();++i) for(int j=0;j<3;++j) real(v(i,j));
        for(Eigen::Index i=0;i<f.rows();++i) for(int j=0;j<3;++j) integer(f(i,j));
    }
    void moments(const physics::RawMassMoments10& m) {for(int i=0;i<10;++i) real(m.values()(i));}
    std::string finish() {return hash.finish();}
};

Bounds bounds(const SourceSurfaceOctreeConfig& c,int level,const std::array<std::uint32_t,3>& xyz) {
    Bounds b;
    for(int a=0;a<3;++a) {
        // The same dyadic plane gets identical coordinates at every level.
        b.lo(a)=std::fma(std::ldexp(double(xyz[a]),-level),c.rootSide,c.rootOrigin(a));
        b.hi(a)=std::fma(std::ldexp(double(xyz[a])+1,-level),c.rootSide,c.rootOrigin(a));
    }
    return b;
}
void validateCellKey(const SourceSurfaceOctreeConfig& config,int level,
    const std::array<std::uint32_t,3>& xyz) {
    if(level<0 || level>config.maximumDepth)
        throw std::out_of_range("source octree query level outside configured depth");
    for(auto x:xyz) if(x>=(std::uint32_t(1)<<level))
        throw std::out_of_range("source octree query coordinate outside root");
}
physics::RawMassMoments10 cuboid(const Bounds& b) {
    const Eigen::Vector3d length=b.hi-b.lo;
    const double volume=length.prod();
    if(!length.allFinite() || (length.array()<=0).any() || !std::isfinite(volume) || volume<=0)
        throw std::overflow_error("source octree cell volume overflow/underflow or collapsed coordinates");
    const Eigen::Vector3d center=b.lo/2+b.hi/2;
    MassMoments m; m.mass=volume; m.firstMoment=volume*center;
    m.secondMoment=volume*(center*center.transpose());
    m.secondMoment.diagonal().array()+=volume*length.array().square()/12;
    const auto result=physics::RawMassMoments10::fromMassMoments(m);
    if(!result.isFinite()) throw std::overflow_error("source octree cuboid moments overflowed");
    return result;
}
void validateConfig(const SourceSurfaceOctreeConfig& c) {
    if(!c.rootOrigin.allFinite() || !std::isfinite(c.rootSide) || c.rootSide<=0 || c.maximumDepth<0 || c.maximumDepth>20)
        throw std::invalid_argument("source octree requires finite positive root and depth 0..20");
    if(c.maximumNodes==0 || c.maximumNodes>std::vector<SourceSurfaceCell>().max_size())
        throw std::length_error("RESOURCE_LIMIT: invalid source octree node budget");
    const auto root=bounds(c,0,{0,0,0});
    cuboid(root);
    const double step=std::ldexp(c.rootSide,-c.maximumDepth);
    for(int a=0;a<3;++a) {
        const double magnitude=std::max(std::abs(root.lo(a)),std::abs(root.hi(a)));
        const double spacing=std::nextafter(magnitude,std::numeric_limits<double>::infinity())-magnitude;
        if(step<spacing || !std::isfinite(step) || step<=0)
            throw std::invalid_argument("source octree dyadic coordinates numerically collapse");
    }
    cuboid(bounds(c,c.maximumDepth,{0,0,0}));
}
void contain(const Eigen::MatrixXd& v,const Bounds& root) {
    if(v.cols()!=3 || v.rows()==0 || !v.allFinite()) throw std::invalid_argument("invalid source vertices");
    for(Eigen::Index i=0;i<v.rows();++i) for(int a=0;a<3;++a)
        if(v(i,a)<root.lo(a) || v(i,a)>root.hi(a)) throw std::invalid_argument("original boundary outside root cube");
}
Point point(const Eigen::MatrixXd& v,int i) {return {v(i,0),v(i,1),v(i,2)};}
bool overlap(const Bounds& a,const Bounds& b) {
    return (a.lo.array()<=b.hi.array()).all() && (b.lo.array()<=a.hi.array()).all();
}
} // namespace

struct SourceSurfaceOctree::Data {
    Eigen::MatrixXd outerV,designV;
    Eigen::MatrixXi outerF,designF;
    SourceSurfaceOctreeConfig config;
    physics::RawMassMoments10 outerMoments,designMoments;
    SourceSurfaceOctreeStats stats;
    std::vector<SourceSurfaceCell> cells;
    std::string sourceHash,configHash,partitionHash;

    std::size_t allocate(std::size_t count) {
        if(count>config.maximumNodes-cells.size())
            throw std::length_error("RESOURCE_LIMIT: source octree node budget exceeded");
        const auto first=cells.size(), required=first+count;
        if(required>cells.capacity()) {
            const auto growth=cells.capacity()>config.maximumNodes/2?config.maximumNodes:2*cells.capacity();
            cells.reserve(std::min(config.maximumNodes,std::max(required,growth)));
        }
        cells.resize(required); return first;
    }

    void visit(std::size_t index,SourceSurfaceCell current,
               const std::vector<std::size_t>& candidates,
               const std::vector<OriginalTriangle>& triangles,
               const CGAL::Side_of_triangle_mesh<Mesh,Kernel>& side) {
        current.index=index;
        ++stats.levels[current.level];
        stats.maximumVisitedLevel=std::max(stats.maximumVisitedLevel,current.level);
        const auto b=bounds(config,current.level,current.coordinates);
        const Kernel::Iso_cuboid_3 exactBox(Point(b.lo.x(),b.lo.y(),b.lo.z()),Point(b.hi.x(),b.hi.y(),b.hi.z()));
        std::vector<std::size_t> touching;
        for(const auto candidate:candidates) {
            const auto& t=triangles[candidate];
            // Bounds formed by exact min/max of represented source coordinates;
            // inclusive broad phase cannot miss enclosed or touching triangles.
            // This pinned CGAL revision exposes triangle/closed-Bbox_3 (not
            // triangle/Iso_cuboid_3). Its SAT converts the binary64 endpoints
            // to Kernel::FT before projection; Gmpq makes all predicates exact.
            const CGAL::Bbox_3 closedBox(b.lo.x(),b.lo.y(),b.lo.z(),b.hi.x(),b.hi.y(),b.hi.z());
            if(overlap(t.bounds,b) && CGAL::do_intersect(t.triangle,closedBox)) touching.push_back(candidate);
        }
        if(touching.empty()) {
            const Point center((exactBox.xmin()+exactBox.xmax())/2,
                               (exactBox.ymin()+exactBox.ymax())/2,
                               (exactBox.zmin()+exactBox.zmax())/2);
            const auto location=side(center);
            if(location==CGAL::ON_BOUNDARY) throw std::logic_error("boundary-free cell center on original surface");
            if(location==CGAL::ON_BOUNDED_SIDE) {
                current.kind=SourceSurfaceCellKind::Full; current.moments=cuboid(b); ++stats.fullNodes;
            } else {current.kind=SourceSurfaceCellKind::Empty; ++stats.emptyNodes;}
        } else if(current.level==config.maximumDepth) {
            current.kind=SourceSurfaceCellKind::Boundary;
            current.unresolvedVolume=cuboid(b).mass(); ++stats.boundaryNodes;
        } else {
            current.kind=SourceSurfaceCellKind::Internal; ++stats.internalNodes;
            current.firstChild=allocate(8);
            for(int child=0;child<8;++child) {
                SourceSurfaceCell next; next.parent=index; next.level=current.level+1;
                for(int a=0;a<3;++a) next.coordinates[a]=2*current.coordinates[a]+((child>>a)&1);
                visit(current.firstChild+child,next,touching,triangles,side);
                const auto& stored=cells[current.firstChild+child];
                current.moments+=stored.moments; current.unresolvedVolume+=stored.unresolvedVolume;
            }
            if(!current.moments.isFinite() || !std::isfinite(current.unresolvedVolume))
                throw std::overflow_error("source octree parent aggregate overflowed");
        }
        cells[index]=current;
    }

    Data(const Eigen::MatrixXd& ov,const Eigen::MatrixXi& of,const Eigen::MatrixXd& dv,const Eigen::MatrixXi& df,const SourceSurfaceOctreeConfig& c)
        :outerV(ov),designV(dv),outerF(of),designF(df),config(c) {
        validateConfig(config);
        const auto root=bounds(config,0,{0,0,0}); contain(outerV,root); contain(designV,root);
        // These existing audits compute actual boundary/manifold/intersection and
        // signed nesting facts. Nef sets validate geometry only; no Nef moments.
        const auto outer=verification::makeExactMaterialRegion(outerV,outerF);
        const auto design=verification::makeExactMaterialRegion(designV,designF);
        if(!outer.facts().valid || !outer.facts().completed || !design.facts().valid || !design.facts().completed)
            throw std::invalid_argument("invalid original material region: outer="+outer.facts().error+" design="+design.facts().error);
        const auto relation=verification::relateExactMaterialRegions(design,outer);
        if(!relation.completed || relation.aSubsetB!=true)
            throw std::invalid_argument("original design is not an exact subset of outer material: "+relation.error);
        outerMoments=integrateSurfaceRawMoments(outerV,outerF);
        designMoments=integrateSurfaceRawMoments(designV,designF);
        if(outerMoments.mass()<=0 || designMoments.mass()<=0)
            throw std::overflow_error("positive original surface volume lost in binary64 integration");

        Mesh mesh;
        for(Eigen::Index i=0;i<designV.rows();++i) mesh.add_vertex(point(designV,i));
        std::vector<OriginalTriangle> triangles; triangles.reserve(designF.rows());
        for(Eigen::Index i=0;i<designF.rows();++i) {
            if(mesh.add_face(Mesh::Vertex_index(designF(i,0)),Mesh::Vertex_index(designF(i,1)),Mesh::Vertex_index(designF(i,2)))==Mesh::null_face())
                throw std::logic_error("validated source failed exact query-mesh construction");
            Bounds b{designV.row(designF(i,0)).transpose(),designV.row(designF(i,0)).transpose()};
            for(int j=1;j<3;++j) {
                b.lo=b.lo.cwiseMin(designV.row(designF(i,j)).transpose());
                b.hi=b.hi.cwiseMax(designV.row(designF(i,j)).transpose());
            }
            triangles.push_back({Triangle(point(designV,designF(i,0)),point(designV,designF(i,1)),point(designV,designF(i,2))),b});
        }
        const CGAL::Side_of_triangle_mesh<Mesh,Kernel> side(mesh);
        std::vector<std::size_t> candidates(triangles.size()); std::iota(candidates.begin(),candidates.end(),0);
        stats.levels.resize(config.maximumDepth+1);
        allocate(1); visit(0,SourceSurfaceCell{},candidates,triangles,side);

        Digest source; source.string("source-surface-original-binary64-v1"); source.surface(outerV,outerF); source.surface(designV,designF); sourceHash=source.finish();
        Digest configuration; configuration.string("source-surface-octree-closed-boundary-exclusion-IMPLEMENTATION_DEFAULT-v1");
        configuration.string("CGAL-f7c3c8212b56c0d6dae63787efc99093f4383415/Simple_cartesian<Gmpq>");
        for(int a=0;a<3;++a) configuration.real(config.rootOrigin(a));
        configuration.real(config.rootSide); configuration.integer(config.maximumDepth);
        configuration.integer(config.maximumNodes); configHash=configuration.finish();
        Digest partition; partition.string("source-surface-octree-scientific-tree-v1");
        partition.string(sourceHash); partition.string(configHash); partition.moments(outerMoments); partition.moments(designMoments);
        partition.integer(cells.size());
        for(const auto& cell:cells) {
            partition.integer(cell.index); partition.integer(cell.parent); partition.integer(cell.firstChild);
            partition.integer(cell.level); for(auto x:cell.coordinates) partition.integer(x);
            partition.integer(static_cast<std::uint8_t>(cell.kind)); partition.moments(cell.moments); partition.real(cell.unresolvedVolume);
        }
        partitionHash=partition.finish();
    }
};

SourceSurfaceOctree::SourceSurfaceOctree(const Eigen::MatrixXd& ov,const Eigen::MatrixXi& of,const Eigen::MatrixXd& dv,const Eigen::MatrixXi& df,const SourceSurfaceOctreeConfig& c)
    :data_(std::make_shared<const Data>(ov,of,dv,df,c)) {}
const Eigen::MatrixXd& SourceSurfaceOctree::outerVertices() const {return data_->outerV;}
const Eigen::MatrixXi& SourceSurfaceOctree::outerFaces() const {return data_->outerF;}
const Eigen::MatrixXd& SourceSurfaceOctree::designVertices() const {return data_->designV;}
const Eigen::MatrixXi& SourceSurfaceOctree::designFaces() const {return data_->designF;}
const SourceSurfaceOctreeConfig& SourceSurfaceOctree::config() const {return data_->config;}
const physics::RawMassMoments10& SourceSurfaceOctree::outerSurfaceMoments() const {return data_->outerMoments;}
const physics::RawMassMoments10& SourceSurfaceOctree::designSurfaceMoments() const {return data_->designMoments;}
const physics::RawMassMoments10& SourceSurfaceOctree::carveableMoments() const {return data_->cells.front().moments;}
double SourceSurfaceOctree::unresolvedBoundaryVolume() const {return data_->cells.front().unresolvedVolume;}
const std::string& SourceSurfaceOctree::sourceSha256() const {return data_->sourceHash;}
const std::string& SourceSurfaceOctree::configSha256() const {return data_->configHash;}
const std::string& SourceSurfaceOctree::partitionSha256() const {return data_->partitionHash;}
const SourceSurfaceOctreeStats& SourceSurfaceOctree::stats() const {return data_->stats;}
const std::vector<SourceSurfaceCell>& SourceSurfaceOctree::cells() const {return data_->cells;}
SourceSurfaceCell SourceSurfaceOctree::cell(int level,std::array<std::uint32_t,3> xyz) const {
    validateCellKey(data_->config,level,xyz);
    std::size_t index=0;
    while(true) {
        const auto& node=data_->cells[index];
        if(node.level==level) return node;
        if(node.kind!=SourceSurfaceCellKind::Internal) {
            SourceSurfaceCell result; result.parent=index; result.level=level; result.coordinates=xyz; result.kind=node.kind;
            if(node.kind==SourceSurfaceCellKind::Full) result.moments=cuboid(bounds(data_->config,level,xyz));
            return result;
        }
        const int shift=level-node.level-1;
        unsigned child=0; for(int a=0;a<3;++a) child|=((xyz[a]>>shift)&1)<<a;
        index=node.firstChild+child;
    }
}
SourceSurfaceCellBounds SourceSurfaceOctree::cellBounds(int level,
    std::array<std::uint32_t,3> xyz) const {
    validateCellKey(data_->config,level,xyz);
    const auto b=bounds(data_->config,level,xyz);
    return {b.lo,b.hi};
}
} // namespace spin2014
