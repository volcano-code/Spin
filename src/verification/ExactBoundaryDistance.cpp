#include "ExactBoundaryDistance.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangle_3_Triangle_3_do_intersect.h>
#include <CGAL/squared_distance_3.h>
#include <algorithm>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace verification::detail {
namespace {
using Q = CGAL::Gmpq;
using K = CGAL::Simple_cartesian<Q>;
using Triangle = K::Triangle_3;

struct Box {
    std::array<Q,3> lo, hi;
    explicit Box(const Triangle& t) {
        for (int axis=0;axis<3;++axis) {
            lo[axis]=hi[axis]=t[0][axis];
            for (int i=1;i<3;++i) {
                lo[axis]=std::min(lo[axis],t[i][axis]);
                hi[axis]=std::max(hi[axis],t[i][axis]);
            }
        }
    }
    void include(const Box& b) {
        for (int i=0;i<3;++i) {
            lo[i]=std::min(lo[i],b.lo[i]); hi[i]=std::max(hi[i],b.hi[i]);
        }
    }
};

Q lowerBound(const Box& a,const Box& b) {
    Q result(0);
    for (int axis=0;axis<3;++axis) {
        const Q gap=std::max(Q(0),std::max(a.lo[axis]-b.hi[axis],b.lo[axis]-a.hi[axis]));
        result+=gap*gap;
    }
    return result;
}

Q triangleDistance(const Triangle& a,const Triangle& b) {
    // Crossing edge/face interiors are not captured by vertex/face or
    // edge/edge candidates, so exact intersection must be tested first.
    if (CGAL::do_intersect(a,b)) return Q(0);
    // For disjoint compact triangles a closest pair exists. The minimal
    // features are vertex/face or edge/edge; parallel face/face and edge/face
    // minimizer continua extend to one of those features. Pinned CGAL's
    // rational point/triangle and segment/segment implementations cover face
    // projections, endpoints, parallel and skew edges without square roots.
    Q best=CGAL::squared_distance(a[0],b);
    for (int i=0;i<3;++i) {
        best=std::min(best,CGAL::squared_distance(a[i],b));
        best=std::min(best,CGAL::squared_distance(b[i],a));
        for (int j=0;j<3;++j)
            best=std::min(best,CGAL::squared_distance(
                K::Segment_3(a[i],a[(i+1)%3]),K::Segment_3(b[j],b[(j+1)%3])));
    }
    return best;
}

struct Node {
    Box box;
    std::size_t begin,end,left=0,right=0;
    bool leaf() const { return left==right; }
};

class Tree {
    std::size_t build(std::size_t begin,std::size_t end) {
        Box box=boxes[order[begin]];
        for (std::size_t i=begin+1;i<end;++i) box.include(boxes[order[i]]);
        const std::size_t node=nodes.size();
        nodes.push_back({box,begin,end});
        if (end-begin>4) {
            int axis=0;
            for (int i=1;i<3;++i)
                if (box.hi[i]-box.lo[i]>box.hi[axis]-box.lo[axis]) axis=i;
            const auto middle=begin+(end-begin)/2;
            // Total exact ordering makes median membership deterministic.
            std::nth_element(order.begin()+begin,order.begin()+middle,order.begin()+end,
                [&](std::size_t a,std::size_t b) {
                    const Q ac=boxes[a].lo[axis]+boxes[a].hi[axis];
                    const Q bc=boxes[b].lo[axis]+boxes[b].hi[axis];
                    return ac<bc || (ac==bc && a<b);
                });
            const auto left=build(begin,middle),right=build(middle,end);
            nodes[node].left=left; nodes[node].right=right;
        }
        return node;
    }
public:
    std::vector<Triangle> triangles;
    std::vector<Box> boxes;
    std::vector<std::size_t> order;
    std::vector<Node> nodes;
    Tree(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f) {
        if (f.rows()==0) throw std::invalid_argument("nonempty exact boundary required");
        std::vector<K::Point_3> points;
        points.reserve(v.rows());
        for (Eigen::Index i=0;i<v.rows();++i) points.emplace_back(Q(v(i,0)),Q(v(i,1)),Q(v(i,2)));
        triangles.reserve(f.rows()); boxes.reserve(f.rows());
        for (Eigen::Index i=0;i<f.rows();++i) {
            triangles.emplace_back(points[f(i,0)],points[f(i,1)],points[f(i,2)]);
            boxes.emplace_back(triangles.back());
        }
        order.resize(f.rows()); std::iota(order.begin(),order.end(),0);
        build(0,order.size());
    }
};

class Search {
    const Tree &a_,&b_;
    void visit(std::size_t ai,std::size_t bi) {
        const auto& a=a_.nodes[ai]; const auto& b=b_.nodes[bi];
        // Strict > only: equal bounds may contain a lexicographically earlier
        // witness. These boxes are rational extrema, never double Bbox_3.
        if (best && lowerBound(a.box,b.box)>*best) return;
        if (a.leaf() && b.leaf()) {
            for (auto i=a.begin;i<a.end;++i) for (auto j=b.begin;j<b.end;++j) {
                const auto x=a_.order[i],y=b_.order[j];
                if (best && lowerBound(a_.boxes[x],b_.boxes[y])>*best) continue;
                const Q distance=triangleDistance(a_.triangles[x],b_.triangles[y]);
                const std::array<std::size_t,2> pair={x,y};
                if (!best || distance<*best || (distance==*best && pair<witness)) {
                    best=distance; witness=pair;
                }
            }
        } else {
            std::array<std::size_t,2> first,second;
            if (!a.leaf() && (b.leaf() || a.end-a.begin>=b.end-b.begin)) {
                first={a.left,bi}; second={a.right,bi};
            } else {
                first={ai,b.left}; second={ai,b.right};
            }
            if (lowerBound(a_.nodes[second[0]].box,b_.nodes[second[1]].box)<
                lowerBound(a_.nodes[first[0]].box,b_.nodes[first[1]].box)) std::swap(first,second);
            visit(first[0],first[1]); visit(second[0],second[1]);
        }
    }
public:
    std::optional<Q> best;
    std::array<std::size_t,2> witness{};
    Search(const Tree& a,const Tree& b):a_(a),b_(b) { visit(0,0); }
};
} // namespace

BoundaryDistance exactBoundaryDistance(const Eigen::MatrixXd& av,const Eigen::MatrixXi& af,
    const Eigen::MatrixXd& bv,const Eigen::MatrixXi& bf,double minimumClearance) {
    const Tree a(av,af),b(bv,bf);
    const Search search(a,b);
    if (!search.best || *search.best<Q(0)) throw std::runtime_error("incomplete exact boundary distance");
    std::ostringstream out; out << *search.best;
    const Q threshold(minimumClearance);
    return {out.str(),search.witness,*search.best>=threshold*threshold,*search.best==Q(0)};
}
} // namespace verification::detail
