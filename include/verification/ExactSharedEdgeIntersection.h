#pragma once
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/number_utils.h>
#include <stdexcept>

namespace verification::detail {

// AB is an indexed common edge, and both triangles must be nondegenerate.
// Distinct supporting planes intersect on the line AB; each triangle meets
// that line exactly on segment AB, so there can be no excess intersection.
// In a common plane, the two triangles have excess intersection precisely
// when C and D lie on the same side of line AB. An exact nonzero component
// of their parallel cross products decides that side without constructing
// rational triangle-triangle intersection vertices. No tolerance or rounding.
template<class Kernel>
bool sharedEdgeHasExcessIntersection(const typename Kernel::Point_3& a,
                                    const typename Kernel::Point_3& b,
                                    const typename Kernel::Point_3& c,
                                    const typename Kernel::Point_3& d) {
    if(CGAL::orientation(a,b,c,d)!=CGAL::COPLANAR) return false;
    const auto edge=b-a;
    const auto first=CGAL::cross_product(edge,c-a);
    const auto second=CGAL::cross_product(edge,d-a);
    if(first==CGAL::NULL_VECTOR || second==CGAL::NULL_VECTOR)
        throw std::invalid_argument("shared-edge predicate requires nondegenerate triangles");
    for(int axis=0;axis<3;++axis) {
        const auto side=CGAL::sign(first[axis]);
        if(side!=CGAL::ZERO) return side==CGAL::sign(second[axis]);
    }
    throw std::logic_error("nonzero exact normal has no nonzero component");
}

} // namespace verification::detail
