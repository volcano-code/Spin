#pragma once
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/number_utils.h>
#include <stdexcept>

namespace verification::detail {

// r and p,a,b are coplanar; p,a,b must be noncollinear. These exact relative
// orientations test both closed sides of the positive cone span(a-p,b-p).
template<class Kernel>
bool inPositivePlanarCone(const typename Kernel::Point_3& p,
                         const typename Kernel::Point_3& a,
                         const typename Kernel::Point_3& b,
                         const typename Kernel::Point_3& r) {
    return CGAL::coplanar_orientation(p,a,b,r)!=CGAL::NEGATIVE &&
           CGAL::coplanar_orientation(p,b,a,r)!=CGAL::NEGATIVE;
}

// The triangles share the INDEXED vertex p and are already certified
// nondegenerate by the caller. A nonzero intersection of their positive cones
// is equivalent to extra triangle contact: scaling the common ray small enough
// lies in both finite triangles. No geometric epsilon, vertex perturbation or
// nearest-neighbor pairing is used, and all constructions remain exact.
template<class Kernel>
bool sharedVertexHasExcessIntersection(const typename Kernel::Point_3& p,
                                      const typename Kernel::Point_3& a,
                                      const typename Kernel::Point_3& b,
                                      const typename Kernel::Point_3& c,
                                      const typename Kernel::Point_3& d) {
    const auto sideC=CGAL::orientation(p,a,b,c);
    const auto sideD=CGAL::orientation(p,a,b,d);
    if(sideC==CGAL::COPLANAR && sideD==CGAL::COPLANAR) {
        // Two noncollinear closed planar cones intersect away from their apex
        // iff a boundary ray of either cone belongs to the other.
        return inPositivePlanarCone<Kernel>(p,a,b,c) ||
               inPositivePlanarCone<Kernel>(p,a,b,d) ||
               inPositivePlanarCone<Kernel>(p,c,d,a) ||
               inPositivePlanarCone<Kernel>(p,c,d,b);
    }
    if(sideC==sideD)return false; // Both are strictly on one side of plane PAB.
    if(sideC==CGAL::COPLANAR)return inPositivePlanarCone<Kernel>(p,a,b,c);
    if(sideD==CGAL::COPLANAR)return inPositivePlanarCone<Kernel>(p,a,b,d);
    const auto sideA=CGAL::orientation(p,c,d,a);
    const auto sideB=CGAL::orientation(p,c,d,b);
    if(sideA==sideB && sideA!=CGAL::COPLANAR)return false;
    if(sideA==CGAL::COPLANAR)return inPositivePlanarCone<Kernel>(p,c,d,a);
    if(sideB==CGAL::COPLANAR)return inPositivePlanarCone<Kernel>(p,c,d,b);
    // Opposite plane sides: a POSITIVE combination of C-P and D-P lies in
    // plane PAB. Absolute signed tetrahedron volumes are exact weights, so
    // no division or triangle-triangle intersection polygon is needed.
    const auto weightC=CGAL::abs(CGAL::volume(p,a,b,d));
    const auto weightD=CGAL::abs(CGAL::volume(p,a,b,c));
    const auto ray=weightC*(c-p)+weightD*(d-p);
    if(ray==CGAL::NULL_VECTOR)
        throw std::invalid_argument("shared-vertex predicate requires nondegenerate triangles");
    return inPositivePlanarCone<Kernel>(p,a,b,p+ray);
}

} // namespace verification::detail
