#include "ExactBoundaryDistance.h"
#include "verification/ExactRationalCoordinates.h"
#include "verification/Sha256.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangle_3_Triangle_3_do_intersect.h>
#include <CGAL/squared_distance_3.h>
#include <algorithm>
#include <cmath>
#include <cstring>
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
    Tree(const std::vector<std::array<std::string,3>>& v,const Eigen::MatrixXi& f) {
        if(f.rows()==0)throw std::invalid_argument("nonempty exact boundary required");
        std::vector<K::Point_3> points;points.reserve(v.size());
        for(const auto& row:v)points.emplace_back(parseExactRationalCoordinate(row[0]),
            parseExactRationalCoordinate(row[1]),parseExactRationalCoordinate(row[2]));
        triangles.reserve(f.rows());boxes.reserve(f.rows());
        for(Eigen::Index i=0;i<f.rows();++i) {
            for(int j=0;j<3;++j)if(f(i,j)<0||static_cast<std::size_t>(f(i,j))>=points.size())
                throw std::invalid_argument("exact boundary face index out of range");
            triangles.emplace_back(points[f(i,0)],points[f(i,1)],points[f(i,2)]);
            if(triangles.back().is_degenerate())
                throw std::invalid_argument("degenerate rational boundary triangle");
            boxes.emplace_back(triangles.back());
        }
        order.resize(f.rows());std::iota(order.begin(),order.end(),0);
        build(0,order.size());
    }
};

std::string meshIdentity(const Eigen::MatrixXd& vertices,const Eigen::MatrixXi& faces) {
    Sha256 hash;const std::string domain="exact-boundary-distance-index-binary64-v1";
    hash.update(reinterpret_cast<const unsigned char*>(domain.data()),domain.size());
    auto integer=[&](std::uint64_t value) {
        unsigned char bytes[8];for(int i=0;i<8;++i)bytes[7-i]=static_cast<unsigned char>(value>>(8*i));
        hash.update(bytes,sizeof(bytes));
    };
    integer(vertices.rows());integer(vertices.cols());
    for(Eigen::Index i=0;i<vertices.rows();++i)for(Eigen::Index j=0;j<vertices.cols();++j) {
        std::uint64_t bits=0;std::memcpy(&bits,&vertices(i,j),sizeof(bits));integer(bits);
    }
    integer(faces.rows());integer(faces.cols());
    for(Eigen::Index i=0;i<faces.rows();++i)for(Eigen::Index j=0;j<faces.cols();++j)
        integer(static_cast<std::uint32_t>(faces(i,j)));
    return hash.finish();
}

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

class NearSelfSearch {
    const Tree& tree_;
    const Eigen::MatrixXi& faces_;
    const Q limit_;
    void visit(std::size_t x,std::size_t nodeIndex) {
        ++result.nodesVisited;
        const auto& node=tree_.nodes[nodeIndex];
        if (lowerBound(tree_.boxes[x],node.box)>limit_) return;
        if (!node.leaf()) {
            visit(x,node.left); visit(x,node.right); return;
        }
        for (auto j=node.begin;j<node.end;++j) {
            const auto y=tree_.order[j];
            if (y<=x || lowerBound(tree_.boxes[x],tree_.boxes[y])>limit_) continue;
            ++result.boxCandidates;
            bool incident=false;
            for (int a=0;a<3;++a) for (int b=0;b<3;++b)
                if (faces_(x,a)==faces_(y,b)) incident=true;
            if (incident) { ++result.incidentExcluded; continue; }
            ++result.distancesEvaluated;
            const auto& a=tree_.triangles[x]; const auto& b=tree_.triangles[y];
            const Q distance=triangleDistance(a,b);
            if (distance>limit_) continue;
            const auto an=CGAL::cross_product(a[1]-a[0],a[2]-a[0]);
            const auto bn=CGAL::cross_product(b[1]-b[0],b[2]-b[0]);
            const bool parallel=CGAL::cross_product(an,bn)==CGAL::NULL_VECTOR;
            std::ostringstream out; out << distance;
            result.pairs.push_back({{x,y},out.str(),parallel});
        }
    }
public:
    NearBoundaryPairs result;
    NearSelfSearch(const Tree& tree,const Eigen::MatrixXi& faces,double threshold)
        :tree_(tree),faces_(faces),limit_(Q(threshold)*Q(threshold)) {
        for (std::size_t x=0;x<tree.triangles.size();++x) visit(x,0);
        std::sort(result.pairs.begin(),result.pairs.end(),[](const auto& a,const auto& b) {
            return a.faces<b.faces;
        });
    }
};
} // namespace

struct ExactBoundaryDistanceIndex::Data {
    Tree tree;
    std::string inputIdentity;
    Data(const Eigen::MatrixXd& vertices,const Eigen::MatrixXi& faces)
        :tree(vertices,faces),inputIdentity(meshIdentity(vertices,faces)) {}
};

const std::string& ExactBoundaryDistanceIndex::inputIdentity() const {
    if(!data_)throw std::logic_error("empty exact boundary distance index");
    return data_->inputIdentity;
}

bool ExactBoundaryDistanceIndex::matches(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces) const {
    return data_&&data_->inputIdentity==meshIdentity(vertices,faces);
}

std::shared_ptr<const ExactBoundaryDistanceIndex> makeExactBoundaryDistanceIndex(
    const Eigen::MatrixXd& vertices,const Eigen::MatrixXi& faces) {
    if(vertices.cols()!=3||faces.cols()!=3||vertices.rows()==0||faces.rows()==0||!vertices.allFinite())
        throw std::invalid_argument("exact boundary index requires finite nonempty triangle mesh data");
    for(Eigen::Index i=0;i<faces.rows();++i)for(int j=0;j<3;++j)
        if(faces(i,j)<0||faces(i,j)>=vertices.rows())
            throw std::invalid_argument("exact boundary index face index out of range");
    auto data=std::make_shared<const ExactBoundaryDistanceIndex::Data>(vertices,faces);
    for(const auto& triangle:data->tree.triangles)
        if(triangle.is_degenerate())throw std::invalid_argument("exact boundary index contains a degenerate triangle");
    return std::shared_ptr<const ExactBoundaryDistanceIndex>(
        new ExactBoundaryDistanceIndex(std::move(data)));
}

BoundaryDistance exactBoundaryDistance(const ExactBoundaryDistanceIndex& a,
    const ExactBoundaryDistanceIndex& b,double minimumClearance) {
    if(!std::isfinite(minimumClearance)||minimumClearance<0)
        throw std::invalid_argument("minimum clearance must be finite and nonnegative");
    if(!a.data_||!b.data_)throw std::invalid_argument("empty exact boundary distance index");
    const Search search(a.data_->tree,b.data_->tree);
    if(!search.best||*search.best<Q(0))throw std::runtime_error("incomplete exact boundary distance");
    std::ostringstream out;out<<*search.best;const Q threshold(minimumClearance);
    return {out.str(),search.witness,*search.best>=threshold*threshold,*search.best==Q(0)};
}

BoundaryDistance exactBoundaryDistance(const Eigen::MatrixXd& av,const Eigen::MatrixXi& af,
    const Eigen::MatrixXd& bv,const Eigen::MatrixXi& bf,double minimumClearance) {
    const auto a=makeExactBoundaryDistanceIndex(av,af),b=makeExactBoundaryDistanceIndex(bv,bf);
    return exactBoundaryDistance(*a,*b,minimumClearance);
}

BoundaryDistance exactBoundaryDistance(
    const std::vector<std::array<std::string,3>>& av,const Eigen::MatrixXi& af,
    const std::vector<std::array<std::string,3>>& bv,const Eigen::MatrixXi& bf,
    double minimumClearance) {
    if(!std::isfinite(minimumClearance)||minimumClearance<0)
        throw std::invalid_argument("minimum clearance must be finite and nonnegative");
    const Tree a(av,af),b(bv,bf);
    const Search search(a,b);
    if(!search.best||*search.best<Q(0))
        throw std::runtime_error("incomplete exact rational boundary distance");
    std::ostringstream out;out<<*search.best;const Q threshold(minimumClearance);
    return {out.str(),search.witness,*search.best>=threshold*threshold,*search.best==Q(0)};
}

NearBoundaryPairs exactNonincidentNearPairs(const Eigen::MatrixXd& v,
    const Eigen::MatrixXi& f,double threshold) {
    if (!std::isfinite(threshold) || threshold<=0 || v.cols()!=3 || f.cols()!=3 ||
        v.rows()==0 || f.rows()==0 || !v.allFinite())
        throw std::invalid_argument("finite nonempty triangles and positive threshold required");
    for (Eigen::Index i=0;i<f.rows();++i) for (int j=0;j<3;++j)
        if (f(i,j)<0 || f(i,j)>=v.rows()) throw std::invalid_argument("invalid face index");
    const Tree tree(v,f);
    for (const auto& t:tree.triangles)
        if (t.is_degenerate()) throw std::invalid_argument("degenerate triangle");
    NearSelfSearch search(tree,f,threshold);
    return std::move(search.result);
}
} // namespace verification::detail
