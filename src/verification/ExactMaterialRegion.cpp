#include "verification/ExactMaterialRegion.h"
#include "verification/Sha256.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <stdexcept>

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
} // namespace verification
