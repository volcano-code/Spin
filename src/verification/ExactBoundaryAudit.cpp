#include "verification/ExactBoundaryAudit.h"
#include "verification/Sha256.h"

#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>

namespace verification {
namespace {
// Exact rational arithmetic on the actual binary64 inputs, including all
// constructions. No floating filter, geometric epsilon, or corner trimming.
using Kernel = CGAL::Simple_cartesian<CGAL::Gmpq>;
using Point = Kernel::Point_3;
using Triangle = Kernel::Triangle_3;
using Segment = Kernel::Segment_3;
using Edge = std::array<int,2>;

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
    if(shared.size()==3) return true;
    const auto intersection=CGAL::intersection(triangles[a],triangles[b]);
    if(!intersection) return false;
    if(shared.empty()) return true;
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

struct Box { std::array<double,3> low,high; };
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
    std::vector<Box> boxes_;
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
    PairScan(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f,
             const std::vector<Point>& p,const std::vector<Triangle>& t,ExactBoundaryFacts& facts)
        :faces_(f),points_(p),triangles_(t),facts_(facts) {
        boxes_.resize(f.rows()); order_.resize(f.rows());
        std::iota(order_.begin(),order_.end(),0);
        for(int i=0;i<f.rows();++i) for(int k=0;k<3;++k) {
            boxes_[i].low[k]=std::min({v(f(i,0),k),v(f(i,1),k),v(f(i,2),k)});
            boxes_[i].high[k]=std::max({v(f(i,0),k),v(f(i,1),k),v(f(i,2),k)});
        }
        build(0,static_cast<int>(f.rows()),0);
    }
    void run() { scan(0,0); }
};

} // namespace

ExactBoundaryAudit auditExactTriangleBoundary(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f,
                                              ExactBoundaryOrientation orientation) {
    if(v.cols()!=3 || f.cols()!=3 || v.rows()==0 || f.rows()==0)
        throw std::invalid_argument("exact boundary requires nonempty (N,3) vertices and (F,3) faces");
    if(v.rows()>std::numeric_limits<int>::max() || f.rows()>std::numeric_limits<int>::max()/2)
        throw std::invalid_argument("exact boundary index capacity exceeded");
    if(!v.allFinite()) throw std::invalid_argument("exact boundary coordinates must be finite binary64");
    if(orientation!=ExactBoundaryOrientation::OutwardPositive && orientation!=ExactBoundaryOrientation::InwardNegative &&
       orientation!=ExactBoundaryOrientation::AnyNonzero) throw std::invalid_argument("invalid exact boundary orientation");
    for(int i=0;i<f.rows();++i) for(int k=0;k<3;++k) {
        if(f(i,k)<0 || f(i,k)>=v.rows()) throw std::invalid_argument("face index out of range");
        if(f(i,k)==f(i,(k+1)%3)) throw std::invalid_argument("face indices must be distinct");
    }
    ExactBoundaryAudit audit;
    auto& r=audit.facts_;
    r.orientation=orientation==ExactBoundaryOrientation::OutwardPositive?"OutwardPositive":
                  orientation==ExactBoundaryOrientation::InwardNegative?"InwardNegative":"AnyNonzero";
    r.inputDigest=inputDigest(v,f,orientation);
    std::vector<Point> points;
    points.reserve(v.rows());
    for(int i=0;i<v.rows();++i) points.emplace_back(v(i,0),v(i,1),v(i,2));
    std::vector<Triangle> triangles;
    triangles.reserve(f.rows());
    std::vector<std::vector<int>> incident(v.rows()),adjacent(f.rows());
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
    for(int vertex=0;vertex<v.rows();++vertex) {
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
    if(r.degenerateFaces.empty()) {
        try {
            PairScan(v,f,points,triangles,r).run();
            r.intersectionCompleted=true;
        } catch(const std::exception& error) {r.intersectionError=error.what();}
          catch(...) {r.intersectionError="unknown exact intersection scan exception";}
        std::sort(r.intersectingPairs.begin(),r.intersectingPairs.end());
    } else r.intersectionError="exact zero-area input faces; intersection predicates not run";
    r.intersectionFree=r.intersectionCompleted && r.intersectingPairs.empty();
    r.boundaryValid=r.closed && r.edgeManifold && r.vertexManifold && r.consistentlyOriented &&
        r.componentOrientationsValid && r.unreferencedVertices.empty() && r.degenerateFaces.empty() && r.intersectionFree;
    return audit;
}
} // namespace verification
