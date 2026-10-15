#pragma once
#include <algorithm>
#include <array>
#include <climits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace verification::detail {

struct ExactRankBox3 {
    std::array<int,3> low,high;
};

inline bool exactRankBoxesOverlap(const ExactRankBox3& a,const ExactRankBox3& b) {
    for(int axis=0;axis<3;++axis)
        if(a.high[axis]<b.low[axis] || b.high[axis]<a.low[axis])return false;
    return true;
}

// Coordinate compression of a finite EXACT point tape, local to one audit.
// Equal values share a rank; strict order is preserved without conversion,
// quantization or epsilon. Thus closed AABB relations are exactly unchanged.
// The consumer still uses the original points for triangle predicates.
template<class Point>
class ExactCoordinateRanks {
    std::vector<std::array<int,3>> ranks_;
public:
    explicit ExactCoordinateRanks(const std::vector<Point>& points) {
        if(points.empty() || points.size()>static_cast<std::size_t>(INT_MAX))
            throw std::invalid_argument("exact rank point count is empty or exceeds index capacity");
        ranks_.resize(points.size());
        std::vector<int> order(points.size());
        std::iota(order.begin(),order.end(),0);
        for(int axis=0;axis<3;++axis) {
            std::sort(order.begin(),order.end(),[&](int a,int b) {
                if(points[a][axis]!=points[b][axis])return points[a][axis]<points[b][axis];
                return a<b;
            });
            int value=0;
            ranks_[order.front()][axis]=value;
            for(std::size_t i=1;i<order.size();++i) {
                if(points[order[i]][axis]!=points[order[i-1]][axis])++value;
                ranks_[order[i]][axis]=value;
            }
        }
    }
    int rank(int vertex,int axis) const {
        if(axis<0 || axis>=3 || vertex<0 || static_cast<std::size_t>(vertex)>=ranks_.size())
            throw std::out_of_range("exact coordinate rank index out of range");
        return ranks_[vertex][axis];
    }
    ExactRankBox3 box(const std::array<int,3>& vertices) const {
        ExactRankBox3 result;
        for(int axis=0;axis<3;++axis) {
            result.low[axis]=result.high[axis]=rank(vertices[0],axis);
            for(int i=1;i<3;++i) {
                const int value=rank(vertices[i],axis);
                result.low[axis]=std::min(result.low[axis],value);
                result.high[axis]=std::max(result.high[axis],value);
            }
        }
        return result;
    }
};

} // namespace verification::detail
