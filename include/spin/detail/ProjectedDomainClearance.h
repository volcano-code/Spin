#pragma once

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h>
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <iterator>
#include <map>
#include <set>
#include <vector>

namespace spin::detail {

// Exact predicates on the represented vertices and axis. Orthogonal
// projection is rational (no normalization/square roots). Dropping one
// coordinate is used only for incidence; distances use the 3-D plane metric.
class ProjectedDomainClearance {
    using Q=CGAL::Gmpq;
    using K=CGAL::Cartesian<Q>;
    using P=K::Point_2;
    using Polygon=CGAL::Polygon_2<K>;
    std::array<Q,3> n_;
    Q n2_,margin2_;
    Q domainHeightMaximum_;
    int drop_=0;
    std::vector<P> hull_;
    bool valid_=false;

    P project(const Eigen::Vector3d& x) const {
        std::array<Q,3> v={Q(x.x()),Q(x.y()),Q(x.z())};
        Q dot=0;for(int i=0;i<3;++i)dot+=v[i]*n_[i];
        std::array<Q,2> p;int j=0;
        for(int i=0;i<3;++i)if(i!=drop_)p[j++]=v[i]-n_[i]*dot/n2_;
        return P(p[0],p[1]);
    }
    K::Point_3 lift(const P& p) const {
        std::array<Q,3> x;int j=0;Q dot=0;
        for(int i=0;i<3;++i)if(i!=drop_) {
            x[i]=j++==0?p.x():p.y();dot+=x[i]*n_[i];
        }
        x[drop_]=-dot/n_[drop_];return {x[0],x[1],x[2]};
    }
    bool boundariesSeparated(const std::vector<P>& a,const std::vector<P>& b) const {
        for(std::size_t i=0;i<a.size();++i)for(std::size_t j=0;j<b.size();++j) {
            const P& a1=a[i];const P& a2=a[(i+1)%a.size()];
            const P& b1=b[j];const P& b2=b[(j+1)%b.size()];
            if(CGAL::do_intersect(K::Segment_2(a1,a2),K::Segment_2(b1,b2)))return false;
            if(CGAL::squared_distance(K::Segment_3(lift(a1),lift(a2)),
                                     K::Segment_3(lift(b1),lift(b2)))<=margin2_)return false;
        }
        return true;
    }
    std::vector<P> convex(const std::vector<Eigen::Vector3d>& vertices) const {
        std::vector<P> points,result;
        for(const auto& v:vertices) {
            if(!v.allFinite())return {};
            points.push_back(project(v));
        }
        CGAL::convex_hull_2(points.begin(),points.end(),std::back_inserter(result));
        return result;
    }
public:
    ProjectedDomainClearance(const Eigen::MatrixXd& domain,const Eigen::Vector3d& axis,double margin) {
        if(domain.cols()!=3 || domain.rows()==0 || !domain.allFinite() ||
           !axis.allFinite() || !std::isfinite(margin) || margin<0)return;
        Eigen::Index drop=0;axis.cwiseAbs().maxCoeff(&drop);drop_=static_cast<int>(drop);
        n2_=0;for(int i=0;i<3;++i){n_[i]=Q(axis[i]);n2_+=n_[i]*n_[i];}
        if(n2_==0)return;
        margin2_=Q(margin)*Q(margin);
        std::vector<Eigen::Vector3d> vertices;
        for(int i=0;i<domain.rows();++i) {
            vertices.push_back(domain.row(i));Q h=0;
            for(int j=0;j<3;++j)h+=n_[j]*Q(domain(i,j));
            if(i==0 || h>domainHeightMaximum_)domainHeightMaximum_=h;
        }
        hull_=convex(vertices);valid_=!hull_.empty();
    }
    bool inside(const std::vector<Eigen::Vector3d>& rim) const {
        if(!valid_ || rim.size()<3)return false;
        std::vector<P> boundary;
        for(const auto& x:rim){if(!x.allFinite())return false;boundary.push_back(project(x));}
        Polygon polygon(boundary.begin(),boundary.end());
        if(!polygon.is_simple() || polygon.area()==0)return false;
        for(const auto& p:hull_)if(polygon.bounded_side(p)!=CGAL::ON_BOUNDED_SIDE)return false;
        // Vertices inside a concave polygon do not imply its convex hull is
        // inside. Exclude every edge crossing, contact and sub-margin gap.
        return boundariesSeparated(boundary,hull_);
    }
    bool separated(const std::vector<Eigen::Vector3d>& vertices) const {
        if(!valid_ || vertices.empty())return false;
        const auto other=convex(vertices);if(other.empty())return false;
        if(hull_.size()>=3) {
            Polygon domain(hull_.begin(),hull_.end());
            for(const auto& p:other)if(domain.bounded_side(p)!=CGAL::ON_UNBOUNDED_SIDE)return false;
        }
        if(other.size()>=3) {
            Polygon polygon(other.begin(),other.end());
            if(polygon.bounded_side(hull_.front())!=CGAL::ON_UNBOUNDED_SIDE)return false;
        }
        return boundariesSeparated(other,hull_);
    }
    // Sufficient certificate for an extruded far graph: only its portion
    // above the design projection needs to exceed the design's axial maximum.
    // A simple enclosing boundary, coherently oriented nondegenerate projected
    // triangles, and cancelling interior edges establish degree one over that
    // boundary. Thus the clipped triangles cover the design without folds.
    bool graphAbove(const std::vector<Eigen::Vector3d>& vertices,
                    const std::vector<std::array<int,3>>& faces,
                    const std::vector<int>& rim) const {
        if(!valid_ || hull_.size()<3 || faces.empty() || rim.size()<3)return false;
        const auto validIndex=[&](int i){return i>=0 && std::size_t(i)<vertices.size();};
        std::vector<Eigen::Vector3d> boundary;
        std::set<std::pair<int,int>> rimEdges;
        for(std::size_t i=0;i<rim.size();++i) {
            int a=rim[i],b=rim[(i+1)%rim.size()];
            if(!validIndex(a) || !validIndex(b) || a==b)return false;
            boundary.push_back(vertices[a]);
            if(!rimEdges.emplace(std::min(a,b),std::max(a,b)).second)return false;
        }
        if(!inside(boundary))return false;
        using Node=std::pair<P,Q>;
        std::vector<Node> nodes;
        for(const auto& v:vertices) {
            if(!v.allFinite())return false;
            Q h=0;for(int j=0;j<3;++j)h+=n_[j]*Q(v[j]);
            nodes.emplace_back(project(v),h);
        }
        const auto cross=[](const P& a,const P& b,const P& c) {
            return (b.x()-a.x())*(c.y()-a.y())-(b.y()-a.y())*(c.x()-a.x());
        };
        std::map<std::pair<int,int>,std::pair<int,int>> edges;
        int orientation=0;
        for(const auto& face:faces) {
            for(int id:face)if(!validIndex(id))return false;
            Q a=cross(nodes[face[0]].first,nodes[face[1]].first,nodes[face[2]].first);
            if(a==0)return false;
            int sign=a>0?1:-1;
            if(orientation && sign!=orientation)return false;
            orientation=sign;
            std::vector<Node> polygon;
            for(int j=0;j<3;++j) {
                int x=face[j],y=face[(j+1)%3];polygon.push_back(nodes[x]);
                auto& edge=edges[{std::min(x,y),std::max(x,y)}];
                ++edge.first;edge.second+=x<y?1:-1;
            }
            // Clip in the exact orthogonal projection, interpolating the
            // original represented axial height with the same rational t.
            for(std::size_t i=0;i<hull_.size() && !polygon.empty();++i) {
                std::vector<Node> clipped;
                const P& p=hull_[i];const P& q=hull_[(i+1)%hull_.size()];
                for(std::size_t j=0;j<polygon.size();++j) {
                    const auto& x=polygon[j];const auto& y=polygon[(j+1)%polygon.size()];
                    Q sx=cross(p,q,x.first),sy=cross(p,q,y.first);
                    if(sx>=0)clipped.push_back(x);
                    if((sx<0 && sy>0) || (sx>0 && sy<0)) {
                        Q t=sx/(sx-sy);
                        clipped.emplace_back(P(x.first.x()+t*(y.first.x()-x.first.x()),
                            x.first.y()+t*(y.first.y()-x.first.y())),x.second+t*(y.second-x.second));
                    }
                }
                polygon=std::move(clipped);
            }
            for(const auto& p:polygon) {
                Q gap=p.second-domainHeightMaximum_;
                if(gap<=0 || gap*gap<=margin2_*n2_)return false;
            }
        }
        std::set<std::pair<int,int>> actualRim;
        for(const auto& e:edges) {
            if(e.second.first==1)actualRim.insert(e.first);
            else if(e.second.first!=2 || e.second.second!=0)return false;
        }
        return actualRim==rimEdges;
    }
};
} // namespace spin::detail
