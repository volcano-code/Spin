#include "verification/ExactMaterialRegion.h"
#include "verification/ExactRationalCoordinates.h"
#include "verification/ExactCuboidUnion.h"
#include "verification/Sha256.h"
#include "ExactBoundaryDistance.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <map>
#include <functional>
#include <set>
#include <limits>
#include <chrono>

namespace verification {
ExactNefBudget::ExactNefBudget(std::size_t maximumComplexity) {
    if(maximumComplexity==0 || maximumComplexity>static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("exact Nef budget must be positive and fit native signed indices");
    facts_.maximumComplexity=maximumComplexity;
}
void ExactNefBudget::preflight(std::size_t additional,const std::string& operation) {
    ++facts_.preflightChecks;
    const auto requested=facts_.currentComplexity+additional;
    if(requested>facts_.maximumComplexity) {
        facts_.lastOperation=operation;facts_.lastRequestedComplexity=requested;facts_.exceeded=true;
        throw std::length_error("RESOURCE_LIMIT: exact Nef preflight "+operation);
    }
}
void ExactNefBudget::acquire(std::size_t count,const std::string& operation) {
    // Throwing string/map work precedes lease admission. A rejected result is
    // observed in the peak, but never left charged after constructor unwinding.
    ++facts_.observationsByOperation[operation];
    facts_.lastOperation=operation;
    const auto requested=facts_.currentComplexity+count;
    facts_.peakComplexity=std::max(facts_.peakComplexity,requested);
    if(requested>facts_.maximumComplexity) {
        facts_.lastRequestedComplexity=requested;facts_.exceeded=true;
        throw std::length_error("RESOURCE_LIMIT: exact Nef observed "+operation);
    }
    facts_.currentComplexity=requested;
}
void ExactNefBudget::release(std::size_t count) noexcept {facts_.currentComplexity-=count;}

class ExactNefBudgetLease {
    ExactNefBudgetPtr budget_;
    std::size_t count_=0;
public:
    ExactNefBudgetLease()=default;
    ExactNefBudgetLease(const ExactNefBudgetPtr& budget,std::size_t count,const std::string& operation):budget_(budget) {
        if(budget_) budget_->acquire(count,operation);
        count_=budget_?count:0;
    }
    ExactNefBudgetLease(const ExactNefBudgetLease&)=delete;
    ExactNefBudgetLease& operator=(const ExactNefBudgetLease&)=delete;
    ExactNefBudgetLease(ExactNefBudgetLease&& other) noexcept
        :budget_(std::move(other.budget_)),count_(other.count_) {other.count_=0;}
    ExactNefBudgetLease& operator=(ExactNefBudgetLease&& other) noexcept {
        if(this!=&other) {
            if(budget_) budget_->release(count_);
            budget_=std::move(other.budget_);count_=other.count_;other.count_=0;
        }
        return *this;
    }
    ~ExactNefBudgetLease() {if(budget_) budget_->release(count_);}
    static void preflight(const ExactNefBudgetPtr& budget,std::size_t additional,const std::string& operation) {
        if(budget) budget->preflight(additional,operation);
    }
};

namespace {
using Kernel = CGAL::Simple_cartesian<CGAL::Gmpq>;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;
using NativeNef = CGAL::Nef_polyhedron_3<Kernel>;

std::size_t nativeComplexity(const NativeNef& n) {
    return n.number_of_vertices()+n.number_of_halfedges()+n.number_of_halffacets()+n.number_of_volumes();
}

// Local ownership adapter: every explicit native result acquires a lease
// immediately after its exact operation, before any successful facts. Moving
// the owner moves the native object and its charge together without copying it.
class Nef {
    std::unique_ptr<NativeNef> value_;
    ExactNefBudgetPtr budget_;
    ExactNefBudgetLease lease_;
    template<class Build> static std::unique_ptr<NativeNef> construct(
        const ExactNefBudgetPtr& budget,const char* operation,Build build) {
        // Every Nef has at least its exterior volume. More complex growth is
        // measured after the operation, never rejected using a guessed mesh bound.
        ExactNefBudgetLease::preflight(budget,1,operation);
        return std::unique_ptr<NativeNef>(new NativeNef(build()));
    }
    template<class Build> Nef(const ExactNefBudgetPtr& budget,const char* operation,Build build)
        :value_(construct(budget,operation,build)),budget_(budget),lease_(budget,nativeComplexity(*value_),operation) {}
public:
    explicit Nef(const ExactNefBudgetPtr& budget={})
        :Nef(budget,"empty construction",[]{return NativeNef();}) {}
    Nef(Polyhedron& surface,const ExactNefBudgetPtr& budget={})
        :Nef(budget,"component construction",[&]{return NativeNef(surface);}) {}
    Nef(const Nef& other):Nef(other.budget_,"solid copy",[&]{return NativeNef(other.native());}) {}
    Nef(Nef&&)=default;
    Nef& operator=(Nef&&)=default;
    Nef& operator=(const Nef& other) {Nef copy(other);*this=std::move(copy);return *this;}
    const NativeNef& native() const {return *value_;}
    const ExactNefBudgetPtr& budget() const {return budget_;}
    bool is_valid() const {return value_->is_valid();}
    bool is_empty() const {return value_->is_empty();}
    bool is_simple() const {return value_->is_simple();}
    auto volumes_begin() const {return value_->volumes_begin();}
    auto volumes_end() const {return value_->volumes_end();}
    auto halffacets_begin() const {return value_->halffacets_begin();}
    auto halffacets_end() const {return value_->halffacets_end();}
    auto number_of_vertices() const {return value_->number_of_vertices();}
    auto number_of_halfedges() const {return value_->number_of_halfedges();}
    auto number_of_halffacets() const {return value_->number_of_halffacets();}
    auto number_of_volumes() const {return value_->number_of_volumes();}
    using SHalfedge_around_facet_const_circulator=NativeNef::SHalfedge_around_facet_const_circulator;
    Nef join(const Nef& b) const {return join(b,budget_);}
    Nef join(const Nef& b,const ExactNefBudgetPtr& budget) const {
        return Nef(budget,"union",[&]{return value_->join(b.native());});
    }
    Nef intersection(const Nef& b) const {return intersection(b,budget_);}
    Nef intersection(const Nef& b,const ExactNefBudgetPtr& budget) const {
        return Nef(budget,"intersection",[&]{return value_->intersection(b.native());});
    }
    Nef difference(const Nef& b) const {return difference(b,budget_);}
    Nef difference(const Nef& b,const ExactNefBudgetPtr& budget) const {
        return Nef(budget,"difference",[&]{return value_->difference(b.native());});
    }
    Nef symmetric_difference(const Nef& b) const {return symmetric_difference(b,budget_);}
    Nef symmetric_difference(const Nef& b,const ExactNefBudgetPtr& budget) const {
        return Nef(budget,"symmetric difference",[&]{return value_->symmetric_difference(b.native());});
    }
    Nef regularization() const {
        return Nef(budget_,"regularization",[&]{return value_->regularization();});
    }
};

class InputBudgetScope {
    ExactNefBudgetPtr budget_;
    std::set<const NativeNef*> seen_;
    std::vector<ExactNefBudgetLease> leases_;
public:
    explicit InputBudgetScope(const ExactNefBudgetPtr& budget):budget_(budget) {}
    void add(const Nef& solid) {
        if(budget_ && solid.budget()!=budget_ && seen_.insert(&solid.native()).second)
            leases_.emplace_back(budget_,nativeComplexity(solid.native()),"existing operand");
    }
};

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

CGAL::Gmpq exactRegularizedVolume(const Nef& region) {
    if(region.is_empty()) return CGAL::Gmpq(0);
    Nef regular=region.regularization();
    checked(regular,"symmetric-difference regularization");
    if(regular.is_empty()) return CGAL::Gmpq(0);
    CGAL::Gmpq signedVolume(0);
    const Kernel::Point_3 origin(0,0,0);
    for(auto facet=regular.halffacets_begin();facet!=regular.halffacets_end();++facet) {
        if(!facet->incident_volume()->mark()) continue;
        // The marked-side twin carries the outward-oriented cycles of this
        // boundary. Summing signed fan tetrahedra over every cycle naturally
        // subtracts polygon holes and internal cavity shells without meshing.
        const auto orientedFacet=facet->twin();
        for(auto cycle=orientedFacet->facet_cycles_begin();cycle!=orientedFacet->facet_cycles_end();++cycle) {
            if(!cycle.is_shalfedge())
                throw std::runtime_error("finite Nef volume has a non-polygon facet cycle");
            Nef::SHalfedge_around_facet_const_circulator current(cycle),end(current);
            std::vector<Kernel::Point_3> ring;
            do {
                ring.push_back(current->source()->center_vertex()->point());
                ++current;
            } while(current!=end);
            if(ring.size()<3)
                throw std::runtime_error("finite Nef volume has a short facet boundary cycle");
            for(std::size_t index=1;index+1<ring.size();++index)
                signedVolume+=CGAL::volume(origin,ring[0],ring[index],ring[index+1]);
        }
    }
    return CGAL::sign(signedVolume)<0 ? -signedVolume : signedVolume;
}

std::string exactString(const CGAL::Gmpq& value) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out<<value;
    return out.str();
}

class ComponentBuilder : public CGAL::Modifier_base<Polyhedron::HalfedgeDS> {
    const std::vector<Kernel::Point_3>& vertices_;
    const Eigen::MatrixXi& faces_;
    const std::vector<int>& components_;
    int component_;
    bool reverse_;
public:
    ComponentBuilder(const std::vector<Kernel::Point_3>& v,const Eigen::MatrixXi& f,
                     const std::vector<int>& c,int component,bool reverse)
        : vertices_(v),faces_(f),components_(c),component_(component),reverse_(reverse) {}
    void operator()(Polyhedron::HalfedgeDS& hds) override {
        std::vector<int> local(vertices_.size(),-1);
        int count=0,faceCount=0;
        for (Eigen::Index i=0;i<faces_.rows();++i) if (components_[i]==component_) {
            ++faceCount;
            for (int j=0;j<3;++j) local[faces_(i,j)]=0;
        }
        for (auto& index : local) if (index==0) index=count++;
        CGAL::Polyhedron_incremental_builder_3<Polyhedron::HalfedgeDS> builder(hds,true);
        builder.begin_surface(count,faceCount);
        for (std::size_t i=0;i<vertices_.size();++i) if (local[i]>=0)
            builder.add_vertex(vertices_[i]);
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

ExactMaterialBoundaryFacts auditExactRationalMaterialBoundary(
    const std::vector<std::array<std::string,3>>& vertices,
    const Eigen::MatrixXi& faces,const VerifiedFixedBoundary* fixedProof) {
    using Clock=std::chrono::steady_clock;
    ExactMaterialBoundaryFacts result;
    const auto started=Clock::now();
    if(fixedProof) {
        const auto audit=auditRationalBoundaryWithFixedProof(vertices,faces,
            ExactBoundaryOrientation::AnyNonzero,*fixedProof);
        result.boundary=audit.audit.facts();
        result.reusedPairTests=audit.reusedPairTests;
        result.fixedProofIdentity=fixedProof->facts().inputDigest;
    } else result.boundary=auditExactTriangleBoundaryRational(vertices,faces,
        ExactBoundaryOrientation::AnyNonzero).facts();
    result.boundarySeconds=std::chrono::duration<double>(Clock::now()-started).count();
    if(!result.boundary.boundaryValid) {
        result.completed=result.boundary.intersectionCompleted;
        result.error=boundaryError(result.boundary);
        return result;
    }
    const auto nestingStarted=Clock::now();
    try {
        std::vector<Kernel::Point_3> points;
        points.reserve(vertices.size());
        for(const auto& row:vertices)
            points.emplace_back(parseExactRationalCoordinate(row[0]),
                parseExactRationalCoordinate(row[1]),parseExactRationalCoordinate(row[2]));
        result.originalSigns=result.boundary.componentVolumeSigns;
        const int n=static_cast<int>(result.originalSigns.size());
        result.witnessVertices.assign(n,-1);
        for(Eigen::Index face=0;face<faces.rows();++face) {
            const int component=result.boundary.faceComponents[face];
            if(result.witnessVertices[component]<0) result.witnessVertices[component]=faces(face,0);
        }
        using Side=CGAL::Side_of_triangle_mesh<Polyhedron,Kernel>;
        // Trees borrow their component surfaces. Keep all surfaces alive until
        // every exact query completes, and destroy trees before the surfaces.
        std::vector<std::unique_ptr<Polyhedron>> surfaces;
        std::vector<std::unique_ptr<Side>> sides;
        for(int i=0;i<n;++i) {
            auto surface=std::make_unique<Polyhedron>();
            ComponentBuilder builder(points,faces,result.boundary.faceComponents,i,result.originalSigns[i]<0);
            surface->delegate(builder);
            if(surface->empty() || !surface->is_valid() || !surface->is_closed())
                throw std::runtime_error("point-location component is not a closed valid surface");
            surfaces.push_back(std::move(surface));
            sides.push_back(std::make_unique<Side>(*surfaces.back()));
        }
        result.inside.assign(n,std::vector<bool>(n,false));
        // A connected embedded boundary disjoint from another boundary lies
        // wholly on one side. Exact location of its actual vertex determines
        // that side; complete pair intersection checks above are essential.
        for(int i=0;i<n;++i) for(int j=0;j<n;++j) if(i!=j) {
            const auto side=(*sides[j])(points[result.witnessVertices[i]]);
            if(side==CGAL::ON_BOUNDARY)
                throw std::runtime_error("exact component witness contacts another component");
            result.inside[i][j]=side==CGAL::ON_BOUNDED_SIDE;
        }
        result.parents.assign(n,-1);result.depths.assign(n,0);
        for(int i=0;i<n;++i) {
            for(int j=0;j<n;++j) if(result.inside[i][j]) {
                if(result.inside[j][i]) throw std::runtime_error("mutual component containment");
                ++result.depths[i];
                if(result.parents[i]<0 || result.inside[j][result.parents[i]]) result.parents[i]=j;
            }
            for(int j=0;j<n;++j) if(result.inside[i][j])
                for(int k=j+1;k<n;++k)
                    if(result.inside[i][k] && !result.inside[j][k] && !result.inside[k][j])
                        throw std::runtime_error("component containers do not form a nesting chain");
        }
        result.completed=true;
        for(int i=0;i<n;++i) {
            if(result.originalSigns[i]!=(result.depths[i]%2==0?1:-1)) {
                result.error="component "+std::to_string(i)+" orientation violates positive roots and alternating nesting";
                result.nestingSeconds=std::chrono::duration<double>(Clock::now()-nestingStarted).count();
                return result;
            }
        }
        result.materialComponents=0;result.boundedVoidComponents=0;
        for(int depth:result.depths) {
            if(depth%2==0)++*result.materialComponents;else ++*result.boundedVoidComponents;
        }
        result.classificationValid=true;
    } catch(const std::exception& error) {result.error=error.what();result.completed=false;}
      catch(...) {result.error="unknown exact point-location material classification failure";result.completed=false;}
    result.nestingSeconds=std::chrono::duration<double>(Clock::now()-nestingStarted).count();
    return result;
}

std::string hashExactNefBudgetFacts(const ExactNefBudgetFacts& f) {
    std::string operations;
    for(const auto& entry:f.observationsByOperation)
        operations+=std::to_string(entry.first.size())+":"+entry.first+":"+std::to_string(entry.second)+";";
    return digest({"exact-nef-budget-snapshot-v1",f.definition,std::to_string(f.maximumComplexity),
        std::to_string(f.currentComplexity),std::to_string(f.peakComplexity),std::to_string(f.preflightChecks),
        std::to_string(f.lastRequestedComplexity),f.exceeded?"exceeded":"within",f.lastOperation,operations});
}

struct ExactMaterialRegion::Data {
    // Original arrays retained exactly, independent of caller lifetime/mutation.
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    bool rationalCoordinates=false;
    std::vector<Kernel::Point_3> exactVertices;
    Nef solid;
    std::vector<ExactMaterialRegion> operationOperands;
    explicit Data(const ExactNefBudgetPtr& budget={}):solid(budget) {}
};

bool ExactMaterialRegion::matchesRationalBoundary(
    const std::vector<std::array<std::string,3>>& vertices,
    const Eigen::MatrixXi& faces) const {
    if(!data_ || !facts_.completed || !facts_.valid || facts_.explicitEmpty ||
       vertices.size()!=data_->exactVertices.size() ||
       faces.rows()!=data_->faces.rows() || faces.cols()!=data_->faces.cols() ||
       faces.cols()!=3 || !(faces.array()==data_->faces.array()).all()) return false;
    try {
        for(std::size_t v=0;v<vertices.size();++v) {
            const auto& actual=data_->exactVertices[v];
            if(parseExactRationalCoordinate(vertices[v][0])!=actual.x() ||
               parseExactRationalCoordinate(vertices[v][1])!=actual.y() ||
               parseExactRationalCoordinate(vertices[v][2])!=actual.z()) return false;
        }
    } catch(const std::invalid_argument&) { return false; }
    return true;
}

std::vector<std::array<std::string,3>> ExactMaterialRegion::rationalBoundaryVertices() const {
    if(!data_ || !facts_.completed || !facts_.valid || facts_.explicitEmpty ||
       data_->exactVertices.empty() || data_->faces.rows()==0)
        throw std::invalid_argument("exact boundary export requires a valid nonempty owned region");
    std::vector<std::array<std::string,3>> rows;
    rows.reserve(data_->exactVertices.size());
    const auto token=[](const CGAL::Gmpq& q) {
        std::ostringstream out;
        out.imbue(std::locale::classic());
        out<<q.numerator()<<"/"<<q.denominator();
        return out.str();
    };
    for(const auto& point:data_->exactVertices)
        rows.push_back({token(point.x()),token(point.y()),token(point.z())});
    return rows;
}

Eigen::MatrixXi ExactMaterialRegion::boundaryFaces() const {
    if(!data_ || !facts_.completed || !facts_.valid || facts_.explicitEmpty ||
       data_->exactVertices.empty() || data_->faces.rows()==0)
        throw std::invalid_argument("exact boundary export requires a valid nonempty owned region");
    return data_->faces;
}

bool ExactMaterialRegion::boundaryUsesRationals() const {
    if(!data_ || !facts_.completed || !facts_.valid || facts_.explicitEmpty)
        throw std::invalid_argument("coordinate policy requires a valid nonempty owned region");
    return data_->rationalCoordinates;
}

Eigen::MatrixXd ExactMaterialRegion::binary64BoundaryVertices() const {
    if(boundaryUsesRationals() || data_->vertices.rows()!=static_cast<Eigen::Index>(data_->exactVertices.size()))
        throw std::invalid_argument("region does not own a binary64 indexed boundary");
    return data_->vertices;
}

std::vector<ExactMaterialRegion> ExactMaterialRegion::constructionOperands() const {
    if(!data_ || !facts_.completed || !facts_.valid || facts_.explicitEmpty)
        throw std::invalid_argument("construction operands require a valid nonempty owned region");
    return data_->operationOperands;
}

ExactMaterialRegion makeExactInteriorDifferenceRegion(
    const ExactMaterialRegion& outer,const ExactMaterialRegion& carved) {
    for(const auto* operand:{&outer,&carved})
        if(!operand->data_ || !operand->facts_.completed || !operand->facts_.valid || operand->facts_.explicitEmpty)
            throw std::invalid_argument("interior difference requires nonempty owned exact regions");
    if(outer.facts_.originalSigns!=std::vector<int>{1} ||
       outer.facts_.parents!=std::vector<int>{-1} ||
       carved.facts_.originalSigns.empty() ||
       !std::all_of(carved.facts_.originalSigns.begin(),carved.facts_.originalSigns.end(),[](int s){return s==1;}) ||
       !std::all_of(carved.facts_.parents.begin(),carved.facts_.parents.end(),[](int p){return p==-1;}))
        throw std::invalid_argument("NOT_SUPPORTED: interior difference requires one outer shell and disjoint carved root shells");
    const auto contained=checkExactMaterialSubset(carved,outer);
    if(!contained.completed || !contained.aSubsetB || !*contained.aSubsetB)
        throw std::invalid_argument("interior difference carved region is not certified inside outer");
    auto rows=outer.rationalBoundaryVertices();
    const auto removedRows=carved.rationalBoundaryVertices();
    const std::size_t offset=rows.size();
    if(offset+removedRows.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
       outer.data_->faces.rows()+carved.data_->faces.rows()>std::numeric_limits<int>::max())
        throw std::length_error("interior difference boundary exceeds signed indices");
    rows.insert(rows.end(),removedRows.begin(),removedRows.end());
    Eigen::MatrixXi faces(outer.data_->faces.rows()+carved.data_->faces.rows(),3);
    faces.topRows(outer.data_->faces.rows())=outer.data_->faces;
    for(Eigen::Index i=0;i<carved.data_->faces.rows();++i)
        for(int j=0;j<3;++j)
            faces(outer.data_->faces.rows()+i,j)=carved.data_->faces(i,2-j)+static_cast<int>(offset);
    ExactMaterialRegion result;
    result.boundary_=auditExactTriangleBoundaryRational(rows,faces,ExactBoundaryOrientation::AnyNonzero).facts();
    // Containment alone permits contact. The exact embedded boundary audit
    // rejects contact/intersection between outer and carved surfaces as well
    // as any defect in the actual concatenated indexed boundary.
    if(!result.boundary_.boundaryValid)
        throw std::invalid_argument("interior difference boundary is not strictly separated: "+boundaryError(result.boundary_));
    auto& r=result.facts_;
    r.algorithm="strict-interior-nef-difference-shared-boundary-v1";
    r.inputDigest=result.boundary_.inputDigest;
    r.regionDigest=digest({r.algorithm,r.backend,r.semantics,r.inputDigest,
        outer.facts_.regionDigest,carved.facts_.regionDigest});
    r.originalSigns.assign(carved.facts_.originalSigns.size()+1,-1);
    r.originalSigns[0]=1;
    if(result.boundary_.componentVolumeSigns!=r.originalSigns)
        throw std::invalid_argument("interior difference boundary component correspondence changed");
    r.parents.assign(r.originalSigns.size(),0);r.parents[0]=-1;
    r.depths.assign(r.originalSigns.size(),1);r.depths[0]=0;
    for(std::size_t i=1;i<r.originalSigns.size();++i)r.checkerReversedComponents.push_back(static_cast<int>(i));
    auto data=std::make_shared<ExactMaterialRegion::Data>();
    data->rationalCoordinates=true;data->vertices.resize(0,3);data->faces=faces;
    data->operationOperands={outer,carved};
    data->exactVertices=outer.data_->exactVertices;
    data->exactVertices.insert(data->exactVertices.end(),carved.data_->exactVertices.begin(),carved.data_->exactVertices.end());
    result.data_=data;
    try {
        Nef raw=outer.data_->solid.difference(carved.data_->solid);
        checked(raw,"shared-operand interior difference raw set");
        // Difference of closed sets excludes the carved boundary. Material
        // regions use regular-closed semantics, as in assemble() above.
        Nef actual=raw.regularization();
        finiteRegularSolid(actual,false,"shared-operand interior difference");
        int material=0,voids=0;
        auto volume=actual.volumes_begin();++volume;
        for(;volume!=actual.volumes_end();++volume) {
            if(volume->mark())++material;else ++voids;
        }
        if(material!=1 || voids!=static_cast<int>(carved.facts_.originalSigns.size()))
            throw std::runtime_error("interior difference Nef topology disagrees with certified shell hierarchy");
        data->solid=std::move(actual);
        r.materialComponents=material;r.boundedVoidComponents=voids;
        r.valid=true;r.completed=true;
    } catch(const std::exception& error) {r.error=error.what();}
    catch(...) {r.error="unknown shared-operand interior difference failure";}
    return result;
}

ExactMaterialRegion makeExactMaterialRegion(const Eigen::MatrixXd& vertices,const Eigen::MatrixXi& faces) {
    return makeExactMaterialRegionWithBudget(vertices,faces,{});
}
ExactMaterialRegion makeExactMaterialRegionWithBudget(const Eigen::MatrixXd& vertices,const Eigen::MatrixXi& faces,const ExactNefBudgetPtr& budget) {
    return ExactMaterialRegion::assemble(vertices,nullptr,faces,budget);
}

ExactMaterialRegion makeExactMaterialRegionFromRationals(
    const std::vector<std::array<std::string,3>>& vertices,const Eigen::MatrixXi& faces) {
    return ExactMaterialRegion::assemble(Eigen::MatrixXd(),&vertices,faces,{});
}

ExactMaterialRegion ExactMaterialRegion::assemble(const Eigen::MatrixXd& vertices,
    const std::vector<std::array<std::string,3>>* rationalVertices,
    const Eigen::MatrixXi& faces,const ExactNefBudgetPtr& budget) {
    ExactMaterialRegion result;
    result.boundary_=(rationalVertices?
        auditExactTriangleBoundaryRational(*rationalVertices,faces,ExactBoundaryOrientation::AnyNonzero):
        auditExactTriangleBoundary(vertices,faces,ExactBoundaryOrientation::AnyNonzero)).facts();
    auto& r=result.facts_;
    if(rationalVertices) r.algorithm="signed-shell-nef-exact-rational-construction-v1";
    r.inputDigest=result.boundary_.inputDigest;
    r.regionDigest=digest({r.algorithm,r.backend,r.semantics,"indexed-boundary",r.inputDigest});
    r.originalSigns=result.boundary_.componentVolumeSigns;
    auto data=std::make_shared<ExactMaterialRegion::Data>(budget);
    data->rationalCoordinates=rationalVertices!=nullptr;
    if(!rationalVertices) data->vertices=vertices;
    else data->vertices.resize(0,3);
    data->faces=faces;
    data->exactVertices.reserve(rationalVertices?rationalVertices->size():static_cast<std::size_t>(vertices.rows()));
    if(rationalVertices) for(const auto& row:*rationalVertices)
        data->exactVertices.emplace_back(parseExactRationalCoordinate(row[0]),
            parseExactRationalCoordinate(row[1]),parseExactRationalCoordinate(row[2]));
    else for(Eigen::Index i=0;i<vertices.rows();++i)
        data->exactVertices.emplace_back(CGAL::Gmpq(vertices(i,0)),CGAL::Gmpq(vertices(i,1)),CGAL::Gmpq(vertices(i,2)));
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
            ComponentBuilder builder(data->exactVertices,faces,result.boundary_.faceComponents,i,reverse);
            surface.delegate(builder);
            if (surface.empty() || !surface.is_valid() || !surface.is_closed())
                throw std::runtime_error("component Polyhedron is not a valid closed surface");
            // Pinned constructor marks bounded volumes, not an oriented halfspace.
            // Positive orientation is provided on the checker-only surface copy.
            Nef solid(surface,budget);
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
        if (n==1) {
            // The first finiteRegularSolid proof certifies this exact owned
            // component; the sign/depth check above requires a positive root.
            // regularize(empty union S) == S. Move the same certified object
            // and its budget lease, without another Boolean or geometry change.
            data->solid=std::move(solids[0]);
        } else {
            std::vector<int> order(n); std::iota(order.begin(),order.end(),0);
            std::stable_sort(order.begin(),order.end(),[&](int a,int b){return r.depths[a]<r.depths[b];});
            for (int i : order) {
                Nef next=r.originalSigns[i]>0?data->solid.join(solids[i]):data->solid.difference(solids[i]);
                checked(next,"material signed operation");
                data->solid=next.regularization();
                checked(data->solid,"material regularization");
            }
            finiteRegularSolid(data->solid,false,"assembled region");
        }
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
    return emptyExactMaterialRegionWithBudget({});
}
ExactMaterialRegion emptyExactMaterialRegionWithBudget(const ExactNefBudgetPtr& budget) {
    ExactMaterialRegion result;
    auto& r=result.facts_;
    r.explicitEmpty=true;
    r.regionDigest=digest({r.algorithm,r.backend,r.semantics,"explicit-empty"});
    try {
        auto data=std::make_shared<ExactMaterialRegion::Data>(budget);
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
    return relateExactMaterialRegionsWithBudget(a,b,{});
}
ExactMaterialRelationFacts relateExactBooleanOperationToMaterial(
    const ExactMaterialRegion& design,const ExactMaterialRegion& cutter,
    const ExactMaterialRegion& material,ExactMaterialBooleanOperation operation) {
    for(const auto* operand:{&design,&cutter,&material})
        if(!operand->facts_.valid || !operand->facts_.completed || !operand->data_)
            throw std::invalid_argument("invalid exact Boolean operand: "+operand->facts_.error);
    if(operation!=ExactMaterialBooleanOperation::Intersection &&
       operation!=ExactMaterialBooleanOperation::Difference)
        throw std::invalid_argument("unsupported exact Boolean operation");
    ExactMaterialRelationFacts r;
    r.algorithm="nef-exact-input-boolean-vs-material-v1";
    const std::string token=operation==ExactMaterialBooleanOperation::Intersection
        ? "intersection":"difference";
    r.aDigest=digest({r.algorithm,token,design.facts_.regionDigest,
        cutter.facts_.regionDigest});
    r.bDigest=material.facts_.regionDigest;
    r.bExplicitEmpty=material.facts_.explicitEmpty;
    r.relationDigest=digest({r.algorithm,r.semantics,r.aDigest,r.bDigest,
        r.bExplicitEmpty?"empty":"nonempty"});
    try {
        Nef raw=operation==ExactMaterialBooleanOperation::Intersection
            ? design.data_->solid.intersection(cutter.data_->solid)
            : design.data_->solid.difference(cutter.data_->solid);
        checked(raw,"exact input Boolean");
        Nef expected=raw.regularization();
        finiteRegularSolid(expected,true,"exact input Boolean result");
        r.aExplicitEmpty=expected.is_empty();
        Nef intersection=expected.intersection(material.data_->solid);
        checked(intersection,"exact Boolean/material intersection");
        Nef volume=intersection.regularization();
        checked(volume,"exact Boolean/material intersection regularization");
        Nef ab=expected.difference(material.data_->solid);
        checked(ab,"exact Boolean/material difference A-B");
        Nef ba=material.data_->solid.difference(expected);
        checked(ba,"exact Boolean/material difference B-A");
        const bool equal=ab.is_empty() && ba.is_empty();
        CGAL::Gmpq symmetricDifferenceVolume(0);
        if(!equal) {
            Nef symmetric=expected.symmetric_difference(material.data_->solid);
            checked(symmetric,"exact Boolean/material symmetric difference");
            symmetricDifferenceVolume=exactRegularizedVolume(symmetric);
        }
        r.intersectionNonempty=!intersection.is_empty();
        r.positiveVolumeOverlap=!volume.is_empty();
        r.aSubsetB=ab.is_empty();r.bSubsetA=ba.is_empty();r.equal=equal;
        r.symmetricDifferenceVolumeExact=exactString(symmetricDifferenceVolume);
        r.completed=true;
    } catch(const std::exception& e) {r.error=e.what();}
    catch(...) {r.error="unknown exact input Boolean operation failure";}
    return r;
}
ExactMaterialSubsetFacts checkExactMaterialSubset(
    const ExactMaterialRegion& a,const ExactMaterialRegion& b) {
    return checkExactMaterialSubsetWithBudget(a,b,{});
}
ExactMaterialOverlapFacts checkExactPositiveVolumeOverlap(
    const ExactMaterialRegion& a,const ExactMaterialRegion& b) {
    return checkExactPositiveVolumeOverlapWithBudget(a,b,{});
}
ExactMaterialOverlapFacts checkExactPositiveVolumeOverlapWithBudget(
    const ExactMaterialRegion& a,const ExactMaterialRegion& b,const ExactNefBudgetPtr& budget) {
    for(const auto* operand:{&a,&b})
        if(!operand->facts_.valid||!operand->facts_.completed||!operand->data_)
            throw std::invalid_argument("invalid exact material operand: "+operand->facts_.error);
    ExactMaterialOverlapFacts r;
    r.aDigest=a.facts_.regionDigest;r.bDigest=b.facts_.regionDigest;
    r.operationDigest=digest({r.algorithm,r.semantics,r.aDigest,r.bDigest,
        a.facts_.explicitEmpty?"empty":"nonempty",b.facts_.explicitEmpty?"empty":"nonempty"});
    if(a.facts_.explicitEmpty||b.facts_.explicitEmpty) {
        r.positiveVolumeOverlap=false;r.completed=true;return r;
    }
    try {
        InputBudgetScope inputs(budget);inputs.add(a.data_->solid);inputs.add(b.data_->solid);
        Nef intersection=a.data_->solid.intersection(b.data_->solid,budget);
        checked(intersection,"positive-volume overlap intersection");
        Nef volume=intersection.regularization();
        checked(volume,"positive-volume overlap regularization");
        const bool overlaps=!volume.is_empty();
        r.positiveVolumeOverlap=overlaps;r.completed=true;
    } catch(const std::exception& e) {r.error=e.what();}
    catch(...) {r.error="unknown exact positive-volume overlap failure";}
    return r;
}
ExactMaterialRelationFacts relateExactMaterialRegionsWithBudget(const ExactMaterialRegion& a,const ExactMaterialRegion& b,const ExactNefBudgetPtr& budget) {
    for (const auto* operand : {&a,&b})
        if (!operand->facts_.valid || !operand->facts_.completed || !operand->data_)
            throw std::invalid_argument("invalid exact material operand: "+operand->facts_.error);
    ExactMaterialRelationFacts r;
    r.aDigest=a.facts_.regionDigest; r.bDigest=b.facts_.regionDigest;
    r.aExplicitEmpty=a.facts_.explicitEmpty; r.bExplicitEmpty=b.facts_.explicitEmpty;
    r.relationDigest=digest({r.algorithm,r.semantics,r.aDigest,r.bDigest,
                            r.aExplicitEmpty?"empty":"nonempty",r.bExplicitEmpty?"empty":"nonempty"});
    try {
        InputBudgetScope inputs(budget);inputs.add(a.data_->solid);inputs.add(b.data_->solid);
        Nef intersection=a.data_->solid.intersection(b.data_->solid,budget);
        checked(intersection,"relation intersection");
        Nef volume=intersection.regularization(); checked(volume,"relation intersection regularization");
        // Unregularized differences authenticate inclusive set containment,
        // retaining boundary-only differences. Do not reuse regularized emptiness.
        Nef ab=a.data_->solid.difference(b.data_->solid,budget); checked(ab,"relation difference A-B");
        Nef ba=b.data_->solid.difference(a.data_->solid,budget); checked(ba,"relation difference B-A");
        const bool contact=!intersection.is_empty(), overlap=!volume.is_empty();
        const bool aSubsetB=ab.is_empty(), bSubsetA=ba.is_empty();
        const bool equal=aSubsetB && bSubsetA;
        CGAL::Gmpq symmetricDifferenceVolume(0);
        if(!equal) {
            Nef symmetricDifference=a.data_->solid.symmetric_difference(b.data_->solid,budget);
            checked(symmetricDifference,"relation symmetric difference");
            symmetricDifferenceVolume=exactRegularizedVolume(symmetricDifference);
        }
        // Publish atomically only after every exact operation/predicate succeeds.
        r.intersectionNonempty=contact; r.positiveVolumeOverlap=overlap;
        r.aSubsetB=aSubsetB; r.bSubsetA=bSubsetA;
        r.equal=equal;
        r.symmetricDifferenceVolumeExact=exactString(symmetricDifferenceVolume);
        r.completed=true;
    } catch (const std::exception& e) { r.error=e.what(); }
    catch (...) { r.error="unknown exact relation operation failure"; }
    return r;
}

ExactMaterialSubsetFacts checkExactMaterialSubsetWithBudget(
    const ExactMaterialRegion& a,const ExactMaterialRegion& b,
    const ExactNefBudgetPtr& budget) {
    for (const auto* operand : {&a,&b})
        if (!operand->facts_.valid || !operand->facts_.completed || !operand->data_)
            throw std::invalid_argument("invalid exact material operand: "+operand->facts_.error);
    ExactMaterialSubsetFacts result;
    result.aDigest=a.facts_.regionDigest;
    result.bDigest=b.facts_.regionDigest;
    try {
        if (a.facts_.explicitEmpty) {
            result.aSubsetB=true;
            result.completed=true;
            return result;
        }
        if (b.facts_.explicitEmpty) {
            result.aSubsetB=false;
            result.completed=true;
            return result;
        }
        InputBudgetScope inputs(budget);
        inputs.add(a.data_->solid);
        inputs.add(b.data_->solid);
        // Do not regularize this difference: boundary-only excess is still
        // content outside the container for inclusive set containment.
        Nef excess=a.data_->solid.difference(b.data_->solid,budget);
        checked(excess,"inclusive material subset difference");
        result.aSubsetB=excess.is_empty();
        result.completed=true;
    } catch (const std::exception& error) {
        result.error=error.what();
    } catch (...) {
        result.error="unknown exact material subset failure";
    }
    return result;
}

ExactMaterialPartitionFacts checkExactMaterialPartition(
    const ExactMaterialRegion& whole,const std::vector<ExactMaterialRegion>& parts) {
    return checkExactMaterialPartitionWithBudget(whole,parts,{});
}
ExactMaterialPartitionFacts checkExactMaterialPartitionWithBudget(
    const ExactMaterialRegion& whole,const std::vector<ExactMaterialRegion>& parts,const ExactNefBudgetPtr& budget) {
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
        InputBudgetScope inputs(budget);inputs.add(whole.data_->solid);
        for(const auto& part:parts) inputs.add(part.data_->solid);
        Nef united(budget);
        std::vector<bool> subsets;
        std::vector<std::array<std::size_t,2>> overlaps;
        for (std::size_t i=0;i<parts.size();++i) {
            Nef next=united.join(parts[i].data_->solid); checked(next,"partition union");
            united=std::move(next);
            Nef excess=parts[i].data_->solid.difference(whole.data_->solid,budget);
            checked(excess,"partition inclusive part subset");
            subsets.push_back(excess.is_empty());
            for (std::size_t j=0;j<i;++j) {
                Nef intersection=parts[j].data_->solid.intersection(parts[i].data_->solid,budget);
                checked(intersection,"partition pair intersection");
                Nef volume=intersection.regularization(); checked(volume,"partition positive volume");
                if (!volume.is_empty()) overlaps.push_back({j,i});
            }
        }
        // Do not regularize differences: boundary-only missing/excess content
        // is significant for inclusive equality.
        Nef missing=whole.data_->solid.difference(united,budget); checked(missing,"partition missing whole content");
        Nef excess=united.difference(whole.data_->solid,budget); checked(excess,"partition excess union content");
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
    r.backend=a.facts_.backend+"/RationalBoundaryDistance_3/Intersections_3";
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
        auto rationalVertices=[](const auto& data) {
            std::vector<std::array<std::string,3>> vertices;
            vertices.reserve(data.exactVertices.size());
            for(const auto& point:data.exactVertices) {
                std::array<std::string,3> row;
                for(int axis=0;axis<3;++axis) {
                    std::ostringstream coordinate; coordinate<<point[axis];
                    row[axis]=coordinate.str();
                }
                vertices.push_back(std::move(row));
            }
            return vertices;
        };
        detail::BoundaryDistance measured;
        if(a.data_->rationalCoordinates || b.data_->rationalCoordinates) {
            // Exact vertices are retained for both construction paths, so mixed
            // rational/binary64 operands are measured in one exact coordinate field.
            const auto av=rationalVertices(*a.data_),bv=rationalVertices(*b.data_);
            measured=detail::exactBoundaryDistance(av,a.data_->faces,bv,b.data_->faces,minimumClearance);
        } else {
            measured=detail::exactBoundaryDistance(a.data_->vertices,a.data_->faces,
                b.data_->vertices,b.data_->faces,minimumClearance);
        }
        auto complete=r;
        complete.squaredDistanceExact=measured.squaredExact; complete.witnessFaces=measured.witness;
        complete.atLeastMinimum=measured.atLeastMinimum; complete.contact=measured.contact;
        complete.applicable=true; complete.state="COMPLETE"; complete.completed=true;
        return complete;
    } catch (const std::exception& e) { r.error=e.what(); }
    catch (...) { r.error="unknown exact boundary distance failure"; }
    return r;
}

ExactBoundaryClearanceFacts checkExactBoundaryClearance(
    const ExactMaterialRegion& a,const ExactMaterialRegion& b,double minimumClearance,
    const detail::ExactBoundaryDistanceIndex& aIndex,
    const detail::ExactBoundaryDistanceIndex& bIndex) {
    if(!std::isfinite(minimumClearance)||minimumClearance<0)
        throw std::invalid_argument("minimum clearance must be a finite nonnegative binary64 value");
    for(const auto* operand:{&a,&b})
        if(!operand->facts_.valid||!operand->facts_.completed||!operand->data_)
            throw std::invalid_argument("invalid exact material operand: "+operand->facts_.error);
    if(a.data_->rationalCoordinates || b.data_->rationalCoordinates)
        throw std::invalid_argument("binary64 boundary clearance does not accept rational construction operands");
    ExactBoundaryClearanceFacts r;
    r.backend=a.facts_.backend+"/Distance_3/Intersections_3";
    r.aDigest=a.facts_.regionDigest;r.bDigest=b.facts_.regionDigest;
    r.aInputDigest=a.facts_.inputDigest;r.bInputDigest=b.facts_.inputDigest;
    r.aExplicitEmpty=a.facts_.explicitEmpty;r.bExplicitEmpty=b.facts_.explicitEmpty;
    std::ostringstream threshold;threshold<<CGAL::Gmpq(minimumClearance);
    r.minimumClearanceExact=threshold.str();
    r.clearanceDigest=digest({r.algorithm,r.backend,r.semantics,r.aDigest,r.bDigest,
        r.aInputDigest,r.bInputDigest,r.minimumClearanceExact,
        r.aExplicitEmpty?"empty":"nonempty",r.bExplicitEmpty?"empty":"nonempty"});
    if(r.aExplicitEmpty||r.bExplicitEmpty) {
        r.state="NOT_APPLICABLE";r.applicable=false;r.completed=true;return r;
    }
    if(!aIndex.matches(a.data_->vertices,a.data_->faces)||
       !bIndex.matches(b.data_->vertices,b.data_->faces))
        throw std::invalid_argument("exact boundary distance index does not match its material region");
    try {
        const auto measured=detail::exactBoundaryDistance(aIndex,bIndex,minimumClearance);
        auto complete=r;
        complete.squaredDistanceExact=measured.squaredExact;complete.witnessFaces=measured.witness;
        complete.atLeastMinimum=measured.atLeastMinimum;complete.contact=measured.contact;
        complete.applicable=true;complete.state="COMPLETE";complete.completed=true;return complete;
    } catch(const std::exception& e) {r.error=e.what();}
    catch(...) {r.error="unknown exact boundary distance failure";}
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

Nef boxSolid(const ExactCuboid& box,const ExactNefBudgetPtr& budget) {
    Eigen::MatrixXd vertices(8,3);
    Eigen::MatrixXi faces(12,3);
    for(int i=0;i<8;++i) for(int axis=0;axis<3;++axis)
        vertices(i,axis)=(i&(1<<axis))?box.upper(axis):box.lower(axis);
    faces<<0,2,1,1,2,3,4,5,6,5,7,6,0,1,4,1,5,4,
        2,6,3,3,6,7,0,4,2,2,4,6,1,3,5,3,7,5;
    std::vector<int> components(12,0);
    std::vector<Kernel::Point_3> exactVertices;
    exactVertices.reserve(8);
    for(int i=0;i<8;++i)
        exactVertices.emplace_back(CGAL::Gmpq(vertices(i,0)),CGAL::Gmpq(vertices(i,1)),
                                   CGAL::Gmpq(vertices(i,2)));
    Polyhedron surface;
    ComponentBuilder builder(exactVertices,faces,components,0,false);
    surface.delegate(builder);
    return Nef(surface,budget);
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
    CGAL::convert_nef_polyhedron_to_polygon_mesh(solid.native(),mesh,true);
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
    return makeExactCuboidUnionWithBudget(boxes,settings,
        std::make_shared<ExactNefBudget>(settings.maximumIntermediateComplexity));
}
ExactCuboidUnion makeExactCuboidUnionWithBudget(const std::vector<ExactCuboid>& boxes,
    const ExactCuboidUnionSettings& settings,const ExactNefBudgetPtr& budget) {
    validateExactCuboidUnionSettings(settings);
    if(!budget || budget->facts().maximumComplexity!=settings.maximumIntermediateComplexity)
        throw std::invalid_argument("union context must match intermediate complexity setting");
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
        // Depth-first balanced reduction holds only logarithmically many
        // sub-results; source input order remains part of the input identity.
        std::function<Nef(std::size_t,std::size_t)> reduce;
        reduce=[&](std::size_t begin,std::size_t end)->Nef {
            if(begin==end) {Nef empty(budget);finiteRegularSolid(empty,true,"empty union");return empty;}
            if(end-begin==1) {
                Nef box=boxSolid(boxes[begin],budget);finiteRegularSolid(box,false,"box union");return box;
            }
            const auto mid=begin+(end-begin)/2;
            Nef left=reduce(begin,mid);
            Nef right=reduce(mid,end);
            Nef united=left.join(right);
            ++facts.booleanOperations;
            finiteRegularSolid(united,false,"balanced union");
            return united;
        };
        Nef united=reduce(0,boxes.size());
        facts.finalComplexity=complexity(united);
        facts.simple=united.is_simple();
        if(!*facts.simple) {
            facts.status="NONMANIFOLD_GEOMETRY";
            facts.diagnostic="exact union is not simple; edge/point contacts are not repaired";
        } else {
            auto arrays=boxes.empty()?std::make_pair(Eigen::MatrixXd(0,3),Eigen::MatrixXi(0,3)):exportUnion(united,settings);
            auto region=boxes.empty()?emptyExactMaterialRegionWithBudget(budget):makeExactMaterialRegionWithBudget(arrays.first,arrays.second,budget);
            if(!region.facts().valid || !region.facts().completed)
                throw std::runtime_error("exported signed boundary invalid: "+region.facts().error);
            Nef difference=united.symmetric_difference(region.data_->solid);
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
    if(budget->facts().exceeded) facts.status="RESOURCE_LIMIT_UNKNOWN";
    if(facts.status!="EXACT_UNION_VALIDATED") {
        result.region_.reset();
        result.vertices_.resize(0,3);result.faces_.resize(0,3);
    }
    facts.budget=budget->facts();facts.peakComplexity=facts.budget.peakComplexity;
    const auto& c=facts.finalComplexity;
    facts.outputDigest=digest({facts.inputDigest,facts.status,facts.diagnostic,facts.meshDigest,
        facts.simple?(*facts.simple?"simple":"nonsimple"):"unknown",
        facts.roundtripEqual?(*facts.roundtripEqual?"equal":"unequal"):"unknown",
        std::to_string(c.vertices),std::to_string(c.halfedges),std::to_string(c.halffacets),std::to_string(c.volumes),
        std::to_string(facts.booleanOperations),std::to_string(facts.peakComplexity),facts.complexityDefinition,
        hashExactNefBudgetFacts(facts.budget)});
    return result;
}
} // namespace verification
