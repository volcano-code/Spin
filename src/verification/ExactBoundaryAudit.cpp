#include "verification/ExactBoundaryAudit.h"
#include "verification/ExactRationalCoordinates.h"
#ifdef MAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS
#include "verification/ConservativeDiagonalSlabs.h"
#endif
#include "verification/ExactSharedEdgeIntersection.h"
#ifdef MAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE
#include "verification/ExactSharedVertexIntersection.h"
#endif
#include "verification/Sha256.h"
#ifdef MAKE_IT_STAND_RANKED_BOUNDARY_BOXES
#include "verification/ExactCoordinateRankBox.h"
#endif

#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#ifdef MAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL
#include <CGAL/Filtered_kernel.h>
#endif
#include <CGAL/intersections.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <chrono>
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
#include <chrono>
#include <iostream>
#endif

namespace verification {
CGAL::Gmpq parseExactRationalCoordinate(const std::string& value) {
    const auto slash=value.find('/');
    if(slash==std::string::npos || value.find('/',slash+1)!=std::string::npos)
        throw std::invalid_argument("exact coordinate must be decimal numerator/positive denominator");
    const auto numerator=value.substr(0,slash),denominator=value.substr(slash+1);
    auto decimal=[](const std::string& token,bool signedValue) {
        const auto start=signedValue && !token.empty() && token[0]=='-'?1U:0U;
        return start<token.size() && std::all_of(token.begin()+start,token.end(),
            [](unsigned char c){return c>='0' && c<='9';});
    };
    if(!decimal(numerator,true) || !decimal(denominator,false) ||
       std::all_of(denominator.begin(),denominator.end(),[](char c){return c=='0';}))
        throw std::invalid_argument("exact coordinate has invalid or zero denominator");
    return CGAL::Gmpq(CGAL::Gmpz(numerator.c_str()),CGAL::Gmpz(denominator.c_str()));
}
namespace {
// Coordinates and constructions are exact rationals in both variants. The
// opt-in predicate filter uses certified intervals and exact fallback; neither
// variant uses a geometric epsilon or corner trimming.
using BaseKernel = CGAL::Simple_cartesian<CGAL::Gmpq>;
#ifdef MAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL
using Kernel = CGAL::Filtered_kernel<BaseKernel>;
#else
using Kernel = BaseKernel;
#endif
using Point = Kernel::Point_3;
using Triangle = Kernel::Triangle_3;
using Segment = Kernel::Segment_3;
using Edge = std::array<int,2>;
struct BoundaryStageClock {
    BoundaryPerformance* output;
    std::chrono::steady_clock::time_point last;
    explicit BoundaryStageClock(BoundaryPerformance* value):output(value) {
        if(output) last=std::chrono::steady_clock::now();
    }
    void mark(const char* name) {
        if(!output) return;
        const auto now=std::chrono::steady_clock::now();
        output->exclusiveSeconds[name]+=std::chrono::duration<double>(now-last).count();
        last=std::chrono::steady_clock::now();
    }
};
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
struct AuditPhaseClock {
    std::chrono::steady_clock::time_point last=std::chrono::steady_clock::now();
    void mark(const char* name) {
        const auto now=std::chrono::steady_clock::now();
        std::cerr<<"AUDIT_PHASE "<<name<<' '
                 <<std::chrono::duration_cast<std::chrono::nanoseconds>(now-last).count()<<'\n';
        last=std::chrono::steady_clock::now();
    }
};
struct AuditPairTotals {
    std::array<std::uint64_t,4> count{},nanoseconds{};
};
thread_local AuditPairTotals pairProfile;
struct AuditPairClock {
    std::size_t category;
    std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
    explicit AuditPairClock(std::size_t value):category(value){}
    ~AuditPairClock() {
        ++pairProfile.count[category];
        pairProfile.nanoseconds[category]+=static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now()-start).count());
    }
};
#endif

std::string inputDigest(const Eigen::MatrixXd& v, const Eigen::MatrixXi& f,
                        ExactBoundaryOrientation orientation) {
    static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,
                  "canonical identity requires IEEE754 binary64");
    Sha256 hash;
    constexpr char prefix[] = "exact-boundary-input-v1";
    hash.update(reinterpret_cast<const unsigned char*>(prefix),sizeof(prefix));
    const auto mode=static_cast<unsigned char>(orientation);
    hash.update(&mode,1);
    auto integer = [&](std::uint64_t x) {
        unsigned char bytes[8];
        for(int k=0;k<8;++k) bytes[7-k]=static_cast<unsigned char>(x>>(8*k));
        hash.update(bytes,8);
    };
    integer(v.rows()); integer(f.rows());
    for(Eigen::Index i=0;i<v.rows();++i) for(int j=0;j<3;++j) {
        std::uint64_t bits;
        const double x=v(i,j);
        std::memcpy(&bits,&x,8);
        integer(bits);
    }
    for(Eigen::Index i=0;i<f.rows();++i) for(int j=0;j<3;++j) integer(f(i,j));
    return hash.finish();
}

bool disallowedIntersection(int a, int b, const Eigen::MatrixXi& faces,
                           const std::vector<Point>& points,
                           const std::vector<Triangle>& triangles) {
    std::vector<int> shared;
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        if(faces(a,i)==faces(b,j)) shared.push_back(faces(a,i));
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
    AuditPairClock pairClock(shared.size());
#endif
    if(shared.size()==3) return true;
#ifdef MAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE
    if(shared.size()==1) {
        const int common=shared[0];
        std::array<int,2> first{},second{};
        int nextA=0,nextB=0;
        for(int i=0;i<3;++i) {
            if(faces(a,i)!=common)first[nextA++]=faces(a,i);
            if(faces(b,i)!=common)second[nextB++]=faces(b,i);
        }
        return detail::sharedVertexHasExcessIntersection<Kernel>(points[common],
            points[first[0]],points[first[1]],points[second[0]],points[second[1]]);
    }
#endif
#ifndef SPIN_DISABLE_EXACT_SHARED_EDGE_SHORTCUT
    if(shared.size()==2) {
        int oppositeA=-1,oppositeB=-1;
        for(int i=0;i<3;++i) {
            if(faces(a,i)!=shared[0] && faces(a,i)!=shared[1]) oppositeA=faces(a,i);
            if(faces(b,i)!=shared[0] && faces(b,i)!=shared[1]) oppositeB=faces(b,i);
        }
        return detail::sharedEdgeHasExcessIntersection<Kernel>(
            points[shared[0]],points[shared[1]],points[oppositeA],points[oppositeB]);
    }
#endif
    // Without an indexed shared simplex, any contact is forbidden. Exact
    // decision suffices; constructing rational intersection vertices is wasted
    // work. Shared-index contacts still require the full containment proof below.
    if(shared.empty()) return CGAL::do_intersect(triangles[a],triangles[b]);
    const auto intersection=CGAL::intersection(triangles[a],triangles[b]);
    if(!intersection) return false;
    auto allowed = [&](const Point& p) {
        if(shared.size()==1) return p==points[shared[0]];
        return Segment(points[shared[0]],points[shared[1]]).has_on(p);
    };
    // A triangle intersection is convex. Checking every extreme point against
    // the shared indexed point/segment proves containment of the whole result.
    if(const auto* p=boost::get<Point>(&*intersection)) return !allowed(*p);
    if(const auto* s=boost::get<Segment>(&*intersection))
        return !allowed(s->source()) || !allowed(s->target());
    if(const auto* t=boost::get<Triangle>(&*intersection))
        return !allowed(t->vertex(0)) || !allowed(t->vertex(1)) || !allowed(t->vertex(2));
    if(const auto* polygon=boost::get<std::vector<Point>>(&*intersection)) {
        if(polygon->empty()) throw std::runtime_error("unexpected empty intersection polygon");
        return std::any_of(polygon->begin(),polygon->end(),[&](const Point& p){return !allowed(p);});
    }
    throw std::runtime_error("unrecognized exact intersection variant");
}

#ifdef MAKE_IT_STAND_RANKED_BOUNDARY_BOXES
using Box=detail::ExactRankBox3;
#else
struct Box { std::array<CGAL::Gmpq,3> low,high; };
#endif
bool overlaps(const Box& a,const Box& b) {
    for(int k=0;k<3;++k) if(a.high[k]<b.low[k] || b.high[k]<a.low[k]) return false;
    return true;
}

// Balanced deterministic median BVH. Bounds use ONLY input-coordinate min/max,
// never centroids or rounded expansion/shrink arithmetic. Closed boxes preserve
// point/edge contact. Same-node recursion partitions unordered pairs exactly once.
class PairScan {
    struct Node { Box box; int begin,end,left=-1,right=-1; };
    const Eigen::MatrixXi& faces_;
    const std::vector<Point>& points_;
    const std::vector<Triangle>& triangles_;
    ExactBoundaryFacts& facts_;
    const std::vector<unsigned char>* fixedFaces_;
    std::uint64_t* reusedPairs_;
    std::vector<Box> boxes_;
#ifdef MAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS
    using Slabs=detail::ConservativeDiagonalSlabs<Point>;
    std::vector<Slabs::Box> diagonalBoxes_;
#endif
    std::vector<int> order_;
    std::vector<Node> nodes_;
    int build(int begin,int end,int depth) {
        Box bounds=boxes_[order_[begin]];
        for(int i=begin+1;i<end;++i) for(int k=0;k<3;++k) {
            bounds.low[k]=std::min(bounds.low[k],boxes_[order_[i]].low[k]);
            bounds.high[k]=std::max(bounds.high[k],boxes_[order_[i]].high[k]);
        }
        const int id=static_cast<int>(nodes_.size());
        nodes_.push_back({bounds,begin,end});
        if(end-begin>8) {
            int axis=depth%3,mid=begin+(end-begin)/2;
            std::nth_element(order_.begin()+begin,order_.begin()+mid,order_.begin()+end,
                [&](int a,int b) {
                    if(boxes_[a].low[axis]!=boxes_[b].low[axis]) return boxes_[a].low[axis]<boxes_[b].low[axis];
                    return a<b;
                });
            const int left=build(begin,mid,depth+1),right=build(mid,end,depth+1);
            nodes_[id].left=left; nodes_[id].right=right;
        }
        return id;
    }
    void scan(int a,int b) {
        const Node& x=nodes_[a]; const Node& y=nodes_[b];
        if(!overlaps(x.box,y.box)) return;
        if(x.left<0 && y.left<0) {
            for(int i=x.begin;i<x.end;++i) for(int j=(a==b?i+1:y.begin);j<y.end;++j) {
                int f=order_[i],g=order_[j];
                if(!overlaps(boxes_[f],boxes_[g])) continue;
                if(f>g) std::swap(f,g);
                if(fixedFaces_ && (*fixedFaces_)[f] && (*fixedFaces_)[g]) {
                    ++*reusedPairs_;
                    continue;
                }
#ifdef MAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS
                if(Slabs::separated(diagonalBoxes_[f],diagonalBoxes_[g])) continue;
#endif
                ++facts_.exactPairTests;
                if(disallowedIntersection(f,g,faces_,points_,triangles_)) facts_.intersectingPairs.push_back({f,g});
            }
        } else if(a==b) {
            scan(x.left,x.left); scan(x.left,x.right); scan(x.right,x.right);
        } else if(x.left>=0 && (y.left<0 || x.end-x.begin>=y.end-y.begin)) {
            scan(x.left,b); scan(x.right,b);
        } else { scan(a,y.left); scan(a,y.right); }
    }
public:
    PairScan(const Eigen::MatrixXi& f,const std::vector<Point>& p,
             const std::vector<Triangle>& t,ExactBoundaryFacts& facts,
             const std::vector<unsigned char>* fixedFaces=nullptr,std::uint64_t* reusedPairs=nullptr)
        :faces_(f),points_(p),triangles_(t),facts_(facts),fixedFaces_(fixedFaces),reusedPairs_(reusedPairs) {
        boxes_.resize(f.rows()); order_.resize(f.rows());
#ifdef MAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS
        const Slabs slabs(p);
        diagonalBoxes_.reserve(f.rows());
        for(int i=0;i<f.rows();++i) diagonalBoxes_.push_back(slabs.box({f(i,0),f(i,1),f(i,2)}));
#endif
        std::iota(order_.begin(),order_.end(),0);
#ifdef MAKE_IT_STAND_RANKED_BOUNDARY_BOXES
        const detail::ExactCoordinateRanks<Point> ranks(p);
        for(int i=0;i<f.rows();++i)boxes_[i]=ranks.box({f(i,0),f(i,1),f(i,2)});
#else
        for(int i=0;i<f.rows();++i) for(int k=0;k<3;++k) {
            boxes_[i].low[k]=p[f(i,0)][k]; boxes_[i].high[k]=p[f(i,0)][k];
            for(int j=1;j<3;++j) {
                boxes_[i].low[k]=std::min(boxes_[i].low[k],p[f(i,j)][k]);
                boxes_[i].high[k]=std::max(boxes_[i].high[k],p[f(i,j)][k]);
            }
        }
#endif
        build(0,static_cast<int>(f.rows()),0);
    }
    void run() { scan(0,0); }
};

} // namespace

static ExactBoundaryFacts auditPoints(const std::vector<Point>& points,const Eigen::MatrixXi& f,
                                     ExactBoundaryOrientation orientation,const std::string& digest,
                                     const std::vector<unsigned char>* fixedFaces=nullptr,
                                     std::uint64_t* reusedPairs=nullptr,
                                     BoundaryPerformance* performance=nullptr) {
    BoundaryStageClock measured(performance);
    ExactBoundaryFacts r;
#ifdef MAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL
    r.backend="CGAL-f7c3c8212b56c0d6dae63787efc99093f4383415/Filtered_kernel<Simple_cartesian<Gmpq>>/certified-interval-exact-fallback-v1";
#endif
#ifdef MAKE_IT_STAND_RANKED_BOUNDARY_BOXES
    r.backend+="/exact-coordinate-rank-boxes-v1";
#endif
#ifdef MAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE
    r.backend+="/exact-shared-vertex-cone-v1";
#endif
    if(fixedFaces) r.backend+="/verified-original-triangle-pairs-v1";
#ifdef MAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS
    r.backend+="/conservative-exact-diagonal-slabs-v1";
#endif
    r.orientation=orientation==ExactBoundaryOrientation::OutwardPositive?"OutwardPositive":
                  orientation==ExactBoundaryOrientation::InwardNegative?"InwardNegative":"AnyNonzero";
    r.inputDigest=digest;
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
    AuditPhaseClock phase;
    pairProfile={};
#endif
    std::vector<Triangle> triangles;
    triangles.reserve(f.rows());
    std::vector<std::vector<int>> incident(points.size()),adjacent(f.rows());
    std::map<Edge,std::vector<std::pair<int,int>>> edges;
    for(int i=0;i<f.rows();++i) {
        triangles.emplace_back(points[f(i,0)],points[f(i,1)],points[f(i,2)]);
        if(triangles.back().is_degenerate()) r.degenerateFaces.push_back(i);
        for(int k=0;k<3;++k) {
            int a=f(i,k),b=f(i,(k+1)%3);
            incident[a].push_back(i);
            edges[{std::min(a,b),std::max(a,b)}].push_back({i,a<b?1:-1});
        }
    }
    for(const auto& entry:edges) {
        const auto& e=entry.second;
        if(e.size()==1) r.boundaryEdges.push_back(entry.first);
        if(e.size()>2) r.nonManifoldEdges.push_back(entry.first);
        if(e.size()==2 && e[0].second==e[1].second) r.inconsistentEdges.push_back(entry.first);
        // Connect a fan in linear space even for grossly nonmanifold edges.
        for(std::size_t i=1;i<e.size();++i) {
            adjacent[e[0].first].push_back(e[i].first);
            adjacent[e[i].first].push_back(e[0].first);
        }
    }
    for(int vertex=0;vertex<static_cast<int>(points.size());++vertex) {
        if(incident[vertex].empty()) {r.unreferencedVertices.push_back(vertex);continue;}
        std::map<int,std::vector<int>> link;
        for(int face:incident[vertex]) {
            int k=0; while(f(face,k)!=vertex) ++k;
            int a=f(face,(k+1)%3),b=f(face,(k+2)%3);
            link[a].push_back(b); link[b].push_back(a);
        }
        bool cycle=link.size()>=3;
        for(const auto& entry:link) cycle=cycle && entry.second.size()==2;
        std::set<int> seen;
        std::vector<int> stack{link.begin()->first};
        while(!stack.empty()) {
            int p=stack.back();stack.pop_back();
            if(!seen.insert(p).second) continue;
            for(int q:link.at(p)) if(!seen.count(q)) stack.push_back(q);
        }
        if(!cycle || seen.size()!=link.size()) r.nonManifoldVertices.push_back(vertex);
    }
    r.closed=r.boundaryEdges.empty();
    r.edgeManifold=r.nonManifoldEdges.empty();
    r.vertexManifold=r.nonManifoldVertices.empty();
    r.consistentlyOriented=r.inconsistentEdges.empty() && r.edgeManifold;
    r.faceComponents.assign(f.rows(),-1);
    measured.mark("topology_and_degeneracy");
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
    phase.mark("topology_and_triangle_degeneracy");
#endif
    std::vector<CGAL::Gmpq> volumes;
    const Point origin(0,0,0);
    for(int seed=0;seed<f.rows();++seed) {
        if(r.faceComponents[seed]>=0) continue;
        int component=static_cast<int>(volumes.size());
        CGAL::Gmpq volume(0);
        std::vector<int> stack{seed}; r.faceComponents[seed]=component;
        while(!stack.empty()) {
            int face=stack.back();stack.pop_back();
            volume+=CGAL::volume(origin,points[f(face,0)],points[f(face,1)],points[f(face,2)]);
            for(int next:adjacent[face]) if(r.faceComponents[next]<0) {
                r.faceComponents[next]=component;stack.push_back(next);
            }
        }
        volumes.push_back(volume);
    }
    r.componentOrientationsValid=true;
    for(const auto& volume:volumes) {
        int sign=static_cast<int>(CGAL::sign(volume));
        r.componentVolumeSigns.push_back(sign);
        std::ostringstream exact;
        exact.imbue(std::locale::classic()); exact << volume;
        r.componentSignedVolumesExact.push_back(exact.str());
        r.componentOrientationsValid=r.componentOrientationsValid && sign!=0 &&
            (orientation==ExactBoundaryOrientation::AnyNonzero ||
             (orientation==ExactBoundaryOrientation::OutwardPositive?sign>0:sign<0));
    }
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
    phase.mark("exact_component_volumes");
#endif
    measured.mark("component_volumes");
    if(performance) {
        performance->exclusiveSeconds["bvh_build"] = 0;
        performance->exclusiveSeconds["intersection_query"] = 0;
    }
    if(r.degenerateFaces.empty()) {
        try {
            PairScan scan(f,points,triangles,r,fixedFaces,reusedPairs);
            measured.mark("bvh_build");
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
            phase.mark("exact_bvh_build");
#endif
            scan.run();
            measured.mark("intersection_query");
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
            phase.mark("exact_bvh_and_intersection_scan");
#endif
            r.intersectionCompleted=true;
        } catch(const std::exception& error) {r.intersectionError=error.what();}
          catch(...) {r.intersectionError="unknown exact intersection scan exception";}
        std::sort(r.intersectingPairs.begin(),r.intersectingPairs.end());
    } else r.intersectionError="exact zero-area input faces; intersection predicates not run";
    r.intersectionFree=r.intersectionCompleted && r.intersectingPairs.empty();
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
    for(int category=0;category<4;++category)
        std::cerr<<"AUDIT_PAIR "<<category<<' '<<pairProfile.count[category]<<' '
                 <<pairProfile.nanoseconds[category]<<'\n';
#endif
    r.boundaryValid=r.closed && r.edgeManifold && r.vertexManifold && r.consistentlyOriented &&
        r.componentOrientationsValid && r.unreferencedVertices.empty() && r.degenerateFaces.empty() && r.intersectionFree;
    measured.mark("finalize");
    return r;
}

static void validateIndexedFaces(std::size_t vertexCount,const Eigen::MatrixXi& f,
                                 ExactBoundaryOrientation orientation) {
    if(vertexCount==0 || f.cols()!=3 || f.rows()==0)
        throw std::invalid_argument("exact boundary requires nonempty (N,3) vertices and (F,3) faces");
    if(vertexCount>static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
       f.rows()>std::numeric_limits<int>::max()/2)
        throw std::invalid_argument("exact boundary index capacity exceeded");
    if(orientation!=ExactBoundaryOrientation::OutwardPositive && orientation!=ExactBoundaryOrientation::InwardNegative &&
       orientation!=ExactBoundaryOrientation::AnyNonzero)
        throw std::invalid_argument("invalid exact boundary orientation");
    for(int i=0;i<f.rows();++i) for(int k=0;k<3;++k) {
        if(f(i,k)<0 || static_cast<std::size_t>(f(i,k))>=vertexCount)
            throw std::invalid_argument("face index out of range");
        if(f(i,k)==f(i,(k+1)%3)) throw std::invalid_argument("face indices must be distinct");
    }
}

ExactBoundaryAudit auditExactTriangleBoundary(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f,
                                              ExactBoundaryOrientation orientation) {
    if(v.cols()!=3) throw std::invalid_argument("exact boundary requires (N,3) vertices");
    validateIndexedFaces(static_cast<std::size_t>(v.rows()),f,orientation);
    if(!v.allFinite()) throw std::invalid_argument("exact boundary coordinates must be finite binary64");
    std::vector<Point> points;
    points.reserve(v.rows());
    for(int i=0;i<v.rows();++i) points.emplace_back(v(i,0),v(i,1),v(i,2));
    ExactBoundaryAudit audit;
    audit.facts_=auditPoints(points,f,orientation,inputDigest(v,f,orientation));
    return audit;
}

static std::pair<std::vector<Point>,std::string> rationalBoundaryInput(
    const std::vector<std::array<std::string,3>>& v,const Eigen::MatrixXi& f,
    ExactBoundaryOrientation orientation, BoundaryPerformance* performance=nullptr) {
    BoundaryStageClock measured(performance);
    validateIndexedFaces(v.size(),f,orientation);
    measured.mark("index_validation");
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
    AuditPhaseClock inputPhase;
#endif
    std::vector<Point> points;
    points.reserve(v.size());
    for(const auto& row:v) {
        const CGAL::Gmpq x=parseExactRationalCoordinate(row[0]);
        const CGAL::Gmpq y=parseExactRationalCoordinate(row[1]);
        const CGAL::Gmpq z=parseExactRationalCoordinate(row[2]);
        points.emplace_back(x,y,z);
    }
    measured.mark("rational_parse");
    Sha256 hash;
    constexpr char prefix[]="exact-rational-boundary-input-v1";
    hash.update(reinterpret_cast<const unsigned char*>(prefix),sizeof(prefix));
    const auto mode=static_cast<unsigned char>(orientation);
    hash.update(&mode,1);
    auto integer=[&](std::uint64_t value) {
        unsigned char bytes[8];
        for(int k=0;k<8;++k) bytes[7-k]=static_cast<unsigned char>(value>>(8*k));
        hash.update(bytes,8);
    };
    auto rational=[&](const CGAL::Gmpq& value) {
        std::ostringstream out;
        out.imbue(std::locale::classic()); out<<value;
        const auto text=out.str();
        integer(text.size());
        hash.update(reinterpret_cast<const unsigned char*>(text.data()),text.size());
    };
    integer(v.size());integer(f.rows());
    for(const auto& point:points) {rational(point.x());rational(point.y());rational(point.z());}
    for(Eigen::Index i=0;i<f.rows();++i) for(int j=0;j<3;++j) integer(f(i,j));
#ifdef MAKE_IT_STAND_PROFILE_EXACT_BOUNDARY
    inputPhase.mark("rational_decode_and_identity");
#endif
    auto digest=hash.finish();
    measured.mark("canonical_digest");
    return {std::move(points),std::move(digest)};
}

ExactBoundaryAudit auditExactTriangleBoundaryRational(
    const std::vector<std::array<std::string,3>>& v,const Eigen::MatrixXi& f,
    ExactBoundaryOrientation orientation) {
    const auto input=rationalBoundaryInput(v,f,orientation);
    ExactBoundaryAudit audit;
    audit.facts_=auditPoints(input.first,f,orientation,input.second);
    return audit;
}

using CoordinateKey=std::tuple<CGAL::Gmpq,CGAL::Gmpq,CGAL::Gmpq>;
static CoordinateKey coordinateKey(const Point& p) {return {p.x(),p.y(),p.z()};}

struct VerifiedFixedBoundary::Data {
    ExactBoundaryFacts facts;
    std::map<CoordinateKey,int> vertexByCoordinates;
    std::set<std::array<int,3>> faces;
    std::size_t vertexCount=0;
};

const ExactBoundaryFacts& VerifiedFixedBoundary::facts() const {return data_->facts;}

VerifiedFixedBoundary verifyFixedRationalBoundary(
    const std::vector<std::array<std::string,3>>& v,const Eigen::MatrixXi& f) {
    const auto input=rationalBoundaryInput(v,f,ExactBoundaryOrientation::AnyNonzero);
    const auto facts=auditPoints(input.first,f,ExactBoundaryOrientation::AnyNonzero,input.second);
    if(!facts.boundaryValid) throw std::invalid_argument("fixed boundary proof requires a valid exact boundary");
    auto data=std::make_shared<VerifiedFixedBoundary::Data>();
    data->facts=facts; data->vertexCount=v.size();
    for(std::size_t i=0;i<input.first.size();++i) {
        const auto inserted=data->vertexByCoordinates.emplace(coordinateKey(input.first[i]),static_cast<int>(i));
        if(!inserted.second) inserted.first->second=-1; // Ambiguous coordinates never establish indexed adjacency.
    }
    for(int i=0;i<f.rows();++i) {
        std::array<int,3> face{f(i,0),f(i,1),f(i,2)};
        std::sort(face.begin(),face.end()); data->faces.insert(face);
    }
    return VerifiedFixedBoundary(std::move(data));
}

FixedSubsetBoundaryAudit auditRationalBoundaryWithFixedProof(
    const std::vector<std::array<std::string,3>>& v,const Eigen::MatrixXi& f,
    ExactBoundaryOrientation orientation,const VerifiedFixedBoundary& proof,BoundaryPerformance* performance) {
    const auto start=performance?std::chrono::steady_clock::now():std::chrono::steady_clock::time_point{};
    if(performance) {
        *performance={}; performance->vertexCount=v.size(); performance->faceCount=f.rows();
    }
    if(!proof.data_ || !proof.data_->facts.boundaryValid)
        throw std::invalid_argument("invalid immutable fixed boundary proof");
    const auto input=rationalBoundaryInput(v,f,orientation,performance);
    BoundaryStageClock measured(performance);
    const auto& data=*proof.data_;
    std::vector<int> labels(v.size(),-1),seen(data.vertexCount,-1);
    std::vector<unsigned char> ambiguous(data.vertexCount,0),fixedFaces(f.rows(),0);
    for(std::size_t i=0;i<input.first.size();++i) {
        const auto found=data.vertexByCoordinates.find(coordinateKey(input.first[i]));
        if(found==data.vertexByCoordinates.end() || found->second<0) continue;
        const int label=found->second; labels[i]=label;
        if(seen[label]>=0) ambiguous[label]=1;
        else seen[label]=static_cast<int>(i);
    }
    std::map<std::array<int,3>,std::vector<int>> matched;
    for(int i=0;i<f.rows();++i) {
        std::array<int,3> key{labels[f(i,0)],labels[f(i,1)],labels[f(i,2)]};
        if(std::any_of(key.begin(),key.end(),[&](int label){return label<0 || ambiguous[label];})) continue;
        std::sort(key.begin(),key.end());
        if(data.faces.count(key)) matched[key].push_back(i);
    }
    for(const auto& entry:matched) if(entry.second.size()==1) fixedFaces[entry.second[0]]=1;
    if(performance) performance->matchedOriginalFaces=std::count(fixedFaces.begin(),fixedFaces.end(),1);
    measured.mark("original_face_matching");
    FixedSubsetBoundaryAudit result;
    result.audit.facts_=auditPoints(input.first,f,orientation,input.second,&fixedFaces,&result.reusedPairTests,performance);
    if(performance) performance->totalSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return result;
}
} // namespace verification
