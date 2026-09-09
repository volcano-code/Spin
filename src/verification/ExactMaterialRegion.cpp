#include "verification/ExactMaterialRegion.h"
#include "verification/ExactCuboidUnion.h"
#include "verification/Sha256.h"
#include "ExactBoundaryDistance.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <map>
#include <functional>

namespace verification {
namespace {
using Kernel = CGAL::Simple_cartesian<CGAL::Gmpq>;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;
using Nef = CGAL::Nef_polyhedron_3<Kernel>;

// Domain-separated NUL-delimited strings; digests are lowercase fixed-width hex.
std::string digest(std::initializer_list<std::string> fields) {
    Sha256 hash;
    const unsigned char zero = 0;
    for (const auto& field : fields) {
        hash.update(reinterpret_cast<const unsigned char*>(field.data()),field.size());
        hash.update(&zero,1);
    }
    return hash.finish();
}

void checked(Nef& set, const char* operation) {
    if (!set.is_valid()) throw std::runtime_error(std::string(operation)+": invalid Nef structure");
}

void finiteRegularSolid(Nef& set, bool allowEmpty, const char* operation) {
    checked(set,operation);
    // Pinned standard-kernel is_bounded() always returns true. The first SNC
    // volume is the unbounded exterior (Mark_bounded_volumes.h); its mark is
    // therefore the actual finite-set check for this rational Cartesian kernel.
    if (set.volumes_begin()==set.volumes_end() || set.volumes_begin()->mark())
        throw std::runtime_error(std::string(operation)+": marked unbounded exterior");
    if (!allowEmpty && set.is_empty())
        throw std::runtime_error(std::string(operation)+": unexpectedly empty solid");
    Nef regular = set.regularization();
    checked(regular,"regularity check");
    Nef difference = set.symmetric_difference(regular);
    checked(difference,"regularity equality");
    if (!difference.is_empty()) throw std::runtime_error(std::string(operation)+": non-regular solid");
}

class ComponentBuilder : public CGAL::Modifier_base<Polyhedron::HalfedgeDS> {
    const Eigen::MatrixXd& vertices_;
    const Eigen::MatrixXi& faces_;
    const std::vector<int>& components_;
    int component_;
    bool reverse_;
public:
    ComponentBuilder(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f,
                     const std::vector<int>& c,int component,bool reverse)
        : vertices_(v),faces_(f),components_(c),component_(component),reverse_(reverse) {}
    void operator()(Polyhedron::HalfedgeDS& hds) override {
        std::vector<int> local(vertices_.rows(),-1);
        int count=0,faceCount=0;
        for (Eigen::Index i=0;i<faces_.rows();++i) if (components_[i]==component_) {
            ++faceCount;
            for (int j=0;j<3;++j) local[faces_(i,j)]=0;
        }
        for (auto& index : local) if (index==0) index=count++;
        CGAL::Polyhedron_incremental_builder_3<Polyhedron::HalfedgeDS> builder(hds,true);
        builder.begin_surface(count,faceCount);
        for (Eigen::Index i=0;i<vertices_.rows();++i) if (local[i]>=0)
            builder.add_vertex(Kernel::Point_3(CGAL::Gmpq(vertices_(i,0)),
                                              CGAL::Gmpq(vertices_(i,1)),CGAL::Gmpq(vertices_(i,2))));
        for (Eigen::Index i=0;i<faces_.rows();++i) if (components_[i]==component_) {
            builder.begin_facet();
            for (int j=0;j<3;++j) builder.add_vertex_to_facet(local[faces_(i,reverse_?2-j:j)]);
            builder.end_facet();
        }
        builder.end_surface();
        if (builder.error()) throw std::runtime_error("component Polyhedron construction failed");
    }
};

std::string boundaryError(const ExactBoundaryFacts& b) {
    std::ostringstream out;
    out << "invalid exact boundary: closed=" << b.closed << ", edge_manifold=" << b.edgeManifold
        << ", vertex_manifold=" << b.vertexManifold << ", oriented=" << b.consistentlyOriented
        << ", nonzero_components=" << b.componentOrientationsValid
        << ", unreferenced_vertices=" << b.unreferencedVertices.size()
        << ", degenerate_faces=" << b.degenerateFaces.size()
        << ", intersection_completed=" << b.intersectionCompleted
        << ", intersecting_pairs=" << b.intersectingPairs.size();
    if (!b.intersectionError.empty()) out << ", " << b.intersectionError;
    return out.str();
}
} // namespace

struct ExactMaterialRegion::Data {
    // Original arrays retained exactly, independent of caller lifetime/mutation.
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    Nef solid;
};

ExactMaterialRegion makeExactMaterialRegion(const Eigen::MatrixXd& vertices,const Eigen::MatrixXi& faces) {
    ExactMaterialRegion result;
    result.boundary_=auditExactTriangleBoundary(vertices,faces,ExactBoundaryOrientation::AnyNonzero).facts();
    auto& r=result.facts_;
    r.inputDigest=result.boundary_.inputDigest;
    r.regionDigest=digest({r.algorithm,r.backend,r.semantics,"indexed-boundary",r.inputDigest});
    r.originalSigns=result.boundary_.componentVolumeSigns;
    auto data=std::make_shared<ExactMaterialRegion::Data>();
    data->vertices=vertices; data->faces=faces;
    result.data_=data;
    if (!result.boundary_.boundaryValid) {
        r.completed=result.boundary_.intersectionCompleted;
        r.error=boundaryError(result.boundary_);
        return result;
    }
    try {
        const int n=static_cast<int>(r.originalSigns.size());
        std::vector<Nef> solids;
        for (int i=0;i<n;++i) {
            const bool reverse=r.originalSigns[i]<0;
            if (reverse) r.checkerReversedComponents.push_back(i);
            Polyhedron surface;
            ComponentBuilder builder(vertices,faces,result.boundary_.faceComponents,i,reverse);
            surface.delegate(builder);
            if (surface.empty() || !surface.is_valid() || !surface.is_closed())
                throw std::runtime_error("component Polyhedron is not a valid closed surface");
            // Pinned constructor marks bounded volumes, not an oriented halfspace.
            // Positive orientation is provided on the checker-only surface copy.
            Nef solid(surface);
            finiteRegularSolid(solid,false,"component construction");
            solids.push_back(std::move(solid));
        }
        std::vector<std::vector<bool>> inside(n,std::vector<bool>(n,false));
        for (int i=0;i<n;++i) for (int j=i+1;j<n;++j) {
            Nef intersection=solids[i].intersection(solids[j]);
            checked(intersection,"component intersection");
            if (intersection.is_empty()) continue;
            Nef ij=solids[i].difference(solids[j]); checked(ij,"component containment i-j");
            Nef ji=solids[j].difference(solids[i]); checked(ji,"component containment j-i");
            inside[i][j]=ij.is_empty(); inside[j][i]=ji.is_empty();
            if (inside[i][j]==inside[j][i]) {
                r.error="component solids are duplicate or intersect without strict nesting";
                r.completed=true;
                return result;
            }
        }
        r.parents.assign(n,-1); r.depths.assign(n,0);
        for (int i=0;i<n;++i) {
            for (int j=0;j<n;++j) if (inside[i][j]) {
                ++r.depths[i];
                if (r.parents[i]<0 || inside[j][r.parents[i]]) r.parents[i]=j;
            }
            // All containers must form a strict chain, independent of volumes.
            for (int j=0;j<n;++j) if (inside[i][j])
                for (int k=j+1;k<n;++k) if (inside[i][k] && !inside[j][k] && !inside[k][j])
                    throw std::runtime_error("component containers do not form a nesting chain");
        }
        for (int i=0;i<n;++i) if (r.originalSigns[i] != (r.depths[i]%2==0?1:-1)) {
            r.error="component "+std::to_string(i)+" orientation violates positive roots and alternating nesting";
            r.completed=true;
            return result;
        }
        std::vector<int> order(n); std::iota(order.begin(),order.end(),0);
        std::stable_sort(order.begin(),order.end(),[&](int a,int b){return r.depths[a]<r.depths[b];});
        for (int i : order) {
            Nef next=r.originalSigns[i]>0?data->solid.join(solids[i]):data->solid.difference(solids[i]);
            checked(next,"material signed operation");
            data->solid=next.regularization();
            checked(data->solid,"material regularization");
        }
        finiteRegularSolid(data->solid,false,"assembled region");
        // SNC volumes are connected 3D cells. The first is unbounded exterior.
        // Disjoint embedded boundary components and strict alternating nesting
        // ensure distinct material cells cannot connect merely at a boundary.
        int material=0,voids=0;
        auto volume=data->solid.volumes_begin();
        ++volume;
        for (;volume!=data->solid.volumes_end();++volume) {
            if (volume->mark()) ++material; else ++voids;
        }
        // Independent topology consistency check: every even-depth connected
        // shell starts one material component; each odd-depth shell one cavity.
        const int positives=static_cast<int>(std::count(r.originalSigns.begin(),r.originalSigns.end(),1));
        if (material!=positives || voids!=n-positives)
            throw std::runtime_error("Nef volume topology disagrees with strict shell hierarchy");
        r.materialComponents=material; r.boundedVoidComponents=voids;
        r.valid=true; r.completed=true;
    } catch (const std::exception& e) { r.error=e.what(); }
    catch (...) { r.error="unknown exact material operation failure"; }
    return result;
}

ExactMaterialRegion emptyExactMaterialRegion() {
    ExactMaterialRegion result;
    auto& r=result.facts_;
    r.explicitEmpty=true;
    r.regionDigest=digest({r.algorithm,r.backend,r.semantics,"explicit-empty"});
    try {
        auto data=std::make_shared<ExactMaterialRegion::Data>();
        finiteRegularSolid(data->solid,true,"explicit empty construction");
        if (!data->solid.is_empty()) throw std::runtime_error("explicit empty constructor returned nonempty set");
        result.data_=data;
        r.materialComponents=0; r.boundedVoidComponents=0;
        r.valid=true; r.completed=true;
    } catch (const std::exception& e) { r.error=e.what(); }
    catch (...) { r.error="unknown explicit empty construction failure"; }
    return result;
}

ExactMaterialRelationFacts relateExactMaterialRegions(const ExactMaterialRegion& a,const ExactMaterialRegion& b) {
    for (const auto* operand : {&a,&b})
        if (!operand->facts_.valid || !operand->facts_.completed || !operand->data_)
            throw std::invalid_argument("invalid exact material operand: "+operand->facts_.error);
    ExactMaterialRelationFacts r;
    r.aDigest=a.facts_.regionDigest; r.bDigest=b.facts_.regionDigest;
    r.aExplicitEmpty=a.facts_.explicitEmpty; r.bExplicitEmpty=b.facts_.explicitEmpty;
    r.relationDigest=digest({r.algorithm,r.semantics,r.aDigest,r.bDigest,
                            r.aExplicitEmpty?"empty":"nonempty",r.bExplicitEmpty?"empty":"nonempty"});
    try {
        Nef intersection=a.data_->solid.intersection(b.data_->solid);
        checked(intersection,"relation intersection");
        Nef volume=intersection.regularization(); checked(volume,"relation intersection regularization");
        // Unregularized differences authenticate inclusive set containment,
        // retaining boundary-only differences. Do not reuse regularized emptiness.
        Nef ab=a.data_->solid.difference(b.data_->solid); checked(ab,"relation difference A-B");
        Nef ba=b.data_->solid.difference(a.data_->solid); checked(ba,"relation difference B-A");
        const bool contact=!intersection.is_empty(), overlap=!volume.is_empty();
        const bool aSubsetB=ab.is_empty(), bSubsetA=ba.is_empty();
        // Publish atomically only after every exact operation/predicate succeeds.
        r.intersectionNonempty=contact; r.positiveVolumeOverlap=overlap;
        r.aSubsetB=aSubsetB; r.bSubsetA=bSubsetA;
        r.equal=aSubsetB && bSubsetA;
        r.completed=true;
    } catch (const std::exception& e) { r.error=e.what(); }
    catch (...) { r.error="unknown exact relation operation failure"; }
    return r;
}

ExactMaterialPartitionFacts checkExactMaterialPartition(
    const ExactMaterialRegion& whole,const std::vector<ExactMaterialRegion>& parts) {
    if (parts.empty()) throw std::invalid_argument("partition requires a nonempty ordered list of parts");
    const auto validate=[](const ExactMaterialRegion& region) {
        if (!region.facts_.valid || !region.facts_.completed || !region.data_)
            throw std::invalid_argument("invalid exact material operand: "+region.facts_.error);
    };
    validate(whole);
    for (const auto& part : parts) validate(part);
    ExactMaterialPartitionFacts r;
    r.backend=whole.facts_.backend;
    r.wholeDigest=whole.facts_.regionDigest; r.wholeExplicitEmpty=whole.facts_.explicitEmpty;
    // Fixed-field prefix plus count binds the ordered operand chain, including
    // repetitions and explicit empties; no caller facts are trusted.
    std::string identity=digest({r.algorithm,r.backend,r.semantics,r.emptyProtocol,r.wholeDigest,
        r.wholeExplicitEmpty?"empty":"nonempty",std::to_string(parts.size())});
    for (const auto& part : parts) {
        r.partDigests.push_back(part.facts_.regionDigest);
        r.partExplicitEmpty.push_back(part.facts_.explicitEmpty);
        identity=digest({"partition-ordered-operand-v1",identity,part.facts_.regionDigest,
                         part.facts_.explicitEmpty?"empty":"nonempty"});
    }
    r.partitionDigest=identity;
    try {
        Nef united;
        std::vector<bool> subsets;
        std::vector<std::array<std::size_t,2>> overlaps;
        for (std::size_t i=0;i<parts.size();++i) {
            Nef next=united.join(parts[i].data_->solid); checked(next,"partition union");
            united=std::move(next);
            Nef excess=parts[i].data_->solid.difference(whole.data_->solid);
            checked(excess,"partition inclusive part subset");
            subsets.push_back(excess.is_empty());
            for (std::size_t j=0;j<i;++j) {
                Nef intersection=parts[j].data_->solid.intersection(parts[i].data_->solid);
                checked(intersection,"partition pair intersection");
                Nef volume=intersection.regularization(); checked(volume,"partition positive volume");
                if (!volume.is_empty()) overlaps.push_back({j,i});
            }
        }
        // Do not regularize differences: boundary-only missing/excess content
        // is significant for inclusive equality.
        Nef missing=whole.data_->solid.difference(united); checked(missing,"partition missing whole content");
        Nef excess=united.difference(whole.data_->solid); checked(excess,"partition excess union content");
        const bool hasMissing=!missing.is_empty(),hasExcess=!excess.is_empty();
        const bool equal=!hasMissing && !hasExcess,disjoint=overlaps.empty();
        std::sort(overlaps.begin(),overlaps.end());
        // Construct a completed copy before publishing optional facts: even an
        // allocation failure leaves the returned incomplete result all unknown.
        auto complete=r;
        complete.partSubsetWhole=std::move(subsets); complete.positiveVolumeOverlapPairs=std::move(overlaps);
        complete.missingWholeContent=hasMissing; complete.excessUnionContent=hasExcess;
        complete.unionEqualsWhole=equal; complete.noPositiveVolumeOverlap=disjoint;
        complete.completePartition=equal && disjoint; complete.completed=true;
        return complete;
    } catch (const std::exception& e) { r.error=e.what(); }
    catch (...) { r.error="unknown exact partition operation failure"; }
    return r;
}

ExactBoundaryClearanceFacts checkExactBoundaryClearance(
    const ExactMaterialRegion& a,const ExactMaterialRegion& b,double minimumClearance) {
    if (!std::isfinite(minimumClearance) || minimumClearance<0)
        throw std::invalid_argument("minimum clearance must be a finite nonnegative binary64 value");
    for (const auto* operand : {&a,&b})
        if (!operand->facts_.valid || !operand->facts_.completed || !operand->data_)
            throw std::invalid_argument("invalid exact material operand: "+operand->facts_.error);
    ExactBoundaryClearanceFacts r;
    r.backend=a.facts_.backend+"/Distance_3/Intersections_3";
    r.aDigest=a.facts_.regionDigest; r.bDigest=b.facts_.regionDigest;
    r.aInputDigest=a.facts_.inputDigest; r.bInputDigest=b.facts_.inputDigest;
    r.aExplicitEmpty=a.facts_.explicitEmpty; r.bExplicitEmpty=b.facts_.explicitEmpty;
    std::ostringstream threshold; threshold << CGAL::Gmpq(minimumClearance);
    r.minimumClearanceExact=threshold.str();
    r.clearanceDigest=digest({r.algorithm,r.backend,r.semantics,r.aDigest,r.bDigest,
        r.aInputDigest,r.bInputDigest,r.minimumClearanceExact,
        r.aExplicitEmpty?"empty":"nonempty",r.bExplicitEmpty?"empty":"nonempty"});
    if (r.aExplicitEmpty || r.bExplicitEmpty) {
        r.state="NOT_APPLICABLE"; r.applicable=false; r.completed=true;
        return r;
    }
    try {
        const auto measured=detail::exactBoundaryDistance(a.data_->vertices,a.data_->faces,
            b.data_->vertices,b.data_->faces,minimumClearance);
        auto complete=r;
        complete.squaredDistanceExact=measured.squaredExact; complete.witnessFaces=measured.witness;
        complete.atLeastMinimum=measured.atLeastMinimum; complete.contact=measured.contact;
        complete.applicable=true; complete.state="COMPLETE"; complete.completed=true;
        return complete;
    } catch (const std::exception& e) { r.error=e.what(); }
    catch (...) { r.error="unknown exact boundary distance failure"; }
    return r;
}

void validateExactCuboidUnionSettings(const ExactCuboidUnionSettings& s) {
    const auto limit=static_cast<std::size_t>(std::numeric_limits<int>::max());
    for(auto value:{s.maximumInputCuboids,s.maximumOutputVertices,s.maximumOutputTriangles,
            s.maximumIntermediateComplexity})
        if(value==0 || value>limit)
            throw std::invalid_argument("cuboid union budgets must be positive and fit signed native indices");
}

std::string hashExactCuboidUnionSettings(const ExactCuboidUnionSettings& s) {
    validateExactCuboidUnionSettings(s);
    return digest({"exact-cuboid-union-settings-v1","IMPLEMENTATION_DEFAULT",
        std::to_string(s.maximumInputCuboids),std::to_string(s.maximumOutputVertices),
        std::to_string(s.maximumOutputTriangles),std::to_string(s.maximumIntermediateComplexity)});
}

namespace {
ExactNefComplexity complexity(const Nef& n) {
    return {n.number_of_vertices(),n.number_of_halfedges(),n.number_of_halffacets(),n.number_of_volumes()};
}

// Checks explicit live structures and conservative known-growth preflights.
// CGAL's internal temporary allocations are not an allocator/cgroup hard cap.
class UnionBudget {
    const ExactCuboidUnionSettings& settings_;
    ExactCuboidUnionFacts& facts_;
public:
    UnionBudget(const ExactCuboidUnionSettings& s,ExactCuboidUnionFacts& f):settings_(s),facts_(f) {}
    void preflight(std::size_t requested) const {
        if(requested>settings_.maximumIntermediateComplexity)
            throw std::length_error("RESOURCE_LIMIT: exact union intermediate complexity");
    }
    void observed(std::size_t actual) {
        facts_.peakComplexity=std::max(facts_.peakComplexity,actual);
        preflight(actual);
    }
    void validate(Nef& set,std::size_t otherLive) {
        observed(otherLive+complexity(set).total());
        checked(set,"cuboid union structure");
        if(set.volumes_begin()==set.volumes_end() || set.volumes_begin()->mark())
            throw std::runtime_error("cuboid union has marked unbounded exterior");
        preflight(otherLive+3*complexity(set).total());
        Nef regular=set.regularization();
        observed(otherLive+complexity(set).total()+complexity(regular).total());
        checked(regular,"cuboid union regularization");
        Nef delta=set.symmetric_difference(regular);
        observed(otherLive+complexity(set).total()+complexity(regular).total()+complexity(delta).total());
        checked(delta,"cuboid union regularity equality");
        if(!delta.is_empty()) throw std::runtime_error("cuboid union is not regular closed");
    }
};

Nef boxSolid(const ExactCuboid& box) {
    Eigen::MatrixXd vertices(8,3);
    Eigen::MatrixXi faces(12,3);
    for(int i=0;i<8;++i) for(int axis=0;axis<3;++axis)
        vertices(i,axis)=(i&(1<<axis))?box.upper(axis):box.lower(axis);
    faces<<0,2,1,1,2,3,4,5,6,5,7,6,0,1,4,1,5,4,
        2,6,3,3,6,7,0,4,2,2,4,6,1,3,5,3,7,5;
    std::vector<int> components(12,0);
    Polyhedron surface;
    ComponentBuilder builder(vertices,faces,components,0,false);
    surface.delegate(builder);
    return Nef(surface);
}

struct Nonrepresentable : std::runtime_error {using std::runtime_error::runtime_error;};

std::size_t triangulatedBoundarySize(const Nef& solid) {
    std::size_t bothOrientations=0;
    for(auto facet=solid.halffacets_begin();facet!=solid.halffacets_end();++facet) {
        std::size_t vertices=0,cycles=0;
        for(auto cycle=facet->facet_cycles_begin();cycle!=facet->facet_cycles_end();++cycle) {
            if(!cycle.is_shalfedge()) throw std::runtime_error("simple union facet has non-polygon cycle");
            Nef::SHalfedge_around_facet_const_circulator current(cycle),end(current);
            std::size_t ringSize=0;
            do {++ringSize;++current;} while(current!=end);
            if(ringSize<3) throw std::runtime_error("simple union facet has short polygon cycle");
            vertices+=ringSize;++cycles;
        }
        if(cycles==0) throw std::runtime_error("simple union facet has no polygon cycle");
        // Planar constrained triangulation with n boundary vertices and h holes
        // has n+2h-2 triangles; simple exact constraints need no Steiner points.
        bothOrientations+=vertices+2*cycles-4;
    }
    if(bothOrientations%2) throw std::runtime_error("union facet twins disagree on triangulation size");
    return bothOrientations/2;
}

// Exact Surface_mesh conversion visits the outer shell of every bounded SNC
// volume, including unmarked voids; the pinned converter uses mark-dependent
// orientation and constrained exact triangulation for polygon facets/holes.
std::pair<Eigen::MatrixXd,Eigen::MatrixXi> exportUnion(
    const Nef& solid,const ExactCuboidUnionSettings& settings) {
    const auto counts=complexity(solid);
    // Count actual polygon cycles, including holes, before triangulation.
    // Output-cap equality must not fail due to a loose growth estimate.
    const auto triangleCount=triangulatedBoundarySize(solid);
    if(counts.vertices>settings.maximumOutputVertices ||
            triangleCount>settings.maximumOutputTriangles)
        throw std::length_error("RESOURCE_LIMIT: exact union conversion preflight");
    CGAL::Surface_mesh<Kernel::Point_3> mesh;
    CGAL::convert_nef_polyhedron_to_polygon_mesh(solid,mesh,true);
    if(mesh.number_of_vertices()>settings.maximumOutputVertices ||
            mesh.number_of_faces()>settings.maximumOutputTriangles)
        throw std::length_error("RESOURCE_LIMIT: exact union triangulated output");
    if(mesh.number_of_vertices()!=counts.vertices || mesh.number_of_faces()!=triangleCount)
        throw std::runtime_error("exact converter output disagrees with polygon-cycle counts");
    using Vertex=std::array<double,3>;
    std::map<Vertex,int> ordered;
    std::vector<Vertex> converted(mesh.number_of_vertices());
    for(auto v:mesh.vertices()) {
        const auto& point=mesh.point(v);
        Vertex value;
        for(int axis=0;axis<3;++axis) {
            double d=CGAL::to_double(point[axis]);
            if(!std::isfinite(d) || CGAL::Gmpq(d)!=point[axis])
                throw Nonrepresentable("exact union vertex does not roundtrip through binary64");
            value[axis]=(d==0)?0:d;
        }
        converted[static_cast<std::size_t>(v)]=value;
        ordered.emplace(value,0);
    }
    Eigen::MatrixXd vertices(ordered.size(),3);
    int index=0;
    for(auto& entry:ordered) {
        entry.second=index;
        for(int axis=0;axis<3;++axis) vertices(index,axis)=entry.first[axis];
        ++index;
    }
    std::vector<std::array<int,3>> triangles;
    for(auto f:mesh.faces()) {
        std::array<int,3> triangle;
        int count=0;
        for(auto v:CGAL::vertices_around_face(mesh.halfedge(f),mesh)) {
            if(count==3) throw std::runtime_error("CGAL export returned non-triangle");
            triangle[count++]=ordered.at(converted[static_cast<std::size_t>(v)]);
        }
        if(count!=3) throw std::runtime_error("CGAL export returned incomplete triangle");
        std::rotate(triangle.begin(),std::min_element(triangle.begin(),triangle.end()),triangle.end());
        triangles.push_back(triangle);
    }
    std::sort(triangles.begin(),triangles.end());
    Eigen::MatrixXi faces(triangles.size(),3);
    for(std::size_t i=0;i<triangles.size();++i)
        for(int j=0;j<3;++j) faces(i,j)=triangles[i][j];
    return {std::move(vertices),std::move(faces)};
}
}

ExactCuboidUnion makeExactCuboidUnion(const std::vector<ExactCuboid>& boxes,
    const ExactCuboidUnionSettings& settings) {
    validateExactCuboidUnionSettings(settings);
    for(const auto& b:boxes)
        if(!b.lower.allFinite() || !b.upper.allFinite() || (b.upper.array()<=b.lower.array()).any())
            throw std::invalid_argument("cuboid union requires finite strictly positive binary64 boxes");
    ExactCuboidUnion result;
    auto& facts=result.facts_;
    facts.settings=settings;facts.inputCuboids=boxes.size();
    std::ostringstream input;
    input.imbue(std::locale::classic());input<<std::hexfloat<<boxes.size()<<';';
    for(const auto& b:boxes) for(int axis=0;axis<3;++axis) input<<b.lower(axis)<<','<<b.upper(axis)<<';';
    facts.inputDigest=digest({facts.algorithm,facts.backend,hashExactCuboidUnionSettings(settings),input.str()});
    try {
        if(boxes.size()>settings.maximumInputCuboids)
            throw std::length_error("RESOURCE_LIMIT: exact union input cuboids");
        UnionBudget budget(settings,facts);
        // Depth-first balanced reduction holds only logarithmically many
        // sub-results; source input order remains part of the input identity.
        std::function<Nef(std::size_t,std::size_t,std::size_t)> reduce;
        reduce=[&](std::size_t begin,std::size_t end,std::size_t otherLive)->Nef {
            if(begin==end) {Nef empty;budget.validate(empty,otherLive);return empty;}
            if(end-begin==1) {
                budget.preflight(otherLive+64); // triangulated cube: 8V+36H+12F+2 volumes
                Nef box=boxSolid(boxes[begin]);budget.validate(box,otherLive);return box;
            }
            const auto mid=begin+(end-begin)/2;
            Nef left=reduce(begin,mid,otherLive);
            Nef right=reduce(mid,end,otherLive+complexity(left).total());
            const auto operands=complexity(left).total()+complexity(right).total();
            budget.preflight(otherLive+2*operands);
            Nef united=left.join(right);
            ++facts.booleanOperations;
            budget.validate(united,otherLive+operands);
            return united;
        };
        Nef united=reduce(0,boxes.size(),0);
        facts.finalComplexity=complexity(united);
        facts.simple=united.is_simple();
        if(!*facts.simple) {
            facts.status="NONMANIFOLD_GEOMETRY";
            facts.diagnostic="exact union is not simple; edge/point contacts are not repaired";
        } else {
            auto arrays=boxes.empty()?std::make_pair(Eigen::MatrixXd(0,3),Eigen::MatrixXi(0,3)):exportUnion(united,settings);
            auto region=boxes.empty()?emptyExactMaterialRegion():makeExactMaterialRegion(arrays.first,arrays.second);
            if(!region.facts().valid || !region.facts().completed)
                throw std::runtime_error("exported signed boundary invalid: "+region.facts().error);
            const auto live=complexity(united).total()+complexity(region.data_->solid).total();
            budget.observed(live);
            budget.preflight(2*live);
            Nef difference=united.symmetric_difference(region.data_->solid);
            budget.observed(live+complexity(difference).total());
            checked(difference,"exported union exact inclusive roundtrip");
            facts.roundtripEqual=difference.is_empty();
            if(!*facts.roundtripEqual) throw std::runtime_error("exported union differs from exact original set");
            facts.meshDigest=region.facts().inputDigest;
            facts.status="EXACT_UNION_VALIDATED";
            facts.diagnostic=boxes.empty()?"explicit empty algebra union":"exact signed boundary and inclusive roundtrip validated";
            result.region_=std::make_shared<const ExactMaterialRegion>(std::move(region));
            result.vertices_=std::move(arrays.first);result.faces_=std::move(arrays.second);
        }
    } catch(const Nonrepresentable& e) {facts.status="NONREPRESENTABLE_GEOMETRY";facts.diagnostic=e.what();}
    catch(const std::length_error& e) {facts.status="RESOURCE_LIMIT_UNKNOWN";facts.diagnostic=e.what();}
    catch(const std::bad_alloc&) {facts.status="RESOURCE_LIMIT_UNKNOWN";facts.diagnostic="exact union allocation failure";}
    catch(const std::exception& e) {facts.status="GEOMETRY_UNKNOWN";facts.diagnostic=e.what();}
    catch(...) {facts.status="GEOMETRY_UNKNOWN";facts.diagnostic="unknown exact union failure";}
    if(facts.status!="EXACT_UNION_VALIDATED") {
        result.region_.reset();
        result.vertices_.resize(0,3);result.faces_.resize(0,3);
    }
    const auto& c=facts.finalComplexity;
    facts.outputDigest=digest({facts.inputDigest,facts.status,facts.diagnostic,facts.meshDigest,
        facts.simple?(*facts.simple?"simple":"nonsimple"):"unknown",
        facts.roundtripEqual?(*facts.roundtripEqual?"equal":"unequal"):"unknown",
        std::to_string(c.vertices),std::to_string(c.halfedges),std::to_string(c.halffacets),std::to_string(c.volumes),
        std::to_string(facts.booleanOperations),std::to_string(facts.peakComplexity),facts.complexityDefinition});
    return result;
}
} // namespace verification
