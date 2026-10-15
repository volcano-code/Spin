#pragma once
#include <CGAL/number_utils.h>
#include <array>
#include <vector>
#include <algorithm>
#include <limits>

namespace verification::detail {
// Exact projections converted to enclosing intervals. Strict separation is
// sufficient only; uncertain intervals and touching triangles use the caller's
// exact predicate. No floating-point vertex perturbation is performed.
template<class Point> class ConservativeDiagonalSlabs {
    static constexpr std::array<std::array<int,3>,10> directions{{
        {{1,1,0}},{{1,-1,0}},{{1,0,1}},{{1,0,-1}},{{0,1,1}},
        {{0,1,-1}},{{1,1,1}},{{1,1,-1}},{{1,-1,1}},{{1,-1,-1}}
    }};
    using Intervals=std::array<std::pair<double,double>,10>;
    std::vector<Intervals> projections;
public:
    struct Box { std::array<double,10> low, high; };
    explicit ConservativeDiagonalSlabs(const std::vector<Point>& points) {
        projections.reserve(points.size());
        for(const auto& p:points) {
            Intervals values;
            for(std::size_t k=0;k<directions.size();++k) {
                const auto& d=directions[k];
                values[k]=CGAL::to_interval(p.x()*d[0]+p.y()*d[1]+p.z()*d[2]);
            }
            projections.push_back(values);
        }
    }
    Box box(const std::array<int,3>& face) const {
        Box result;
        result.low.fill(std::numeric_limits<double>::infinity());
        result.high.fill(-std::numeric_limits<double>::infinity());
        for(int index:face) {
            const auto& values=projections.at(static_cast<std::size_t>(index));
            for(std::size_t k=0;k<directions.size();++k) {
                result.low[k]=std::min(result.low[k],values[k].first);
                result.high[k]=std::max(result.high[k],values[k].second);
            }
        }
        return result;
    }
    static bool separated(const Box& a,const Box& b) {
        for(std::size_t k=0;k<directions.size();++k)
            if(a.high[k]<b.low[k] || b.high[k]<a.low[k]) return true;
        return false;
    }
};
}
