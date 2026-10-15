#pragma once

namespace spin { namespace detail {
/** Shared formula ordering for legacy and fixed-precision objective paths.
 * Scalar supplies second-order arithmetic; no eigenvectors are selected.
 */
template<class Scalar> struct CentralInertia10 {
    Scalar i11, i22, i33, i12, i13, i23;
};
template<class Moments>
auto centralInertia10(const Moments& r) -> CentralInertia10<typename Moments::value_type> {
    CentralInertia10<typename Moments::value_type> result;
    result.i11=r[5]+r[6]-(r[2]*r[2]+r[3]*r[3])/r[0];
    result.i22=r[4]+r[6]-(r[1]*r[1]+r[3]*r[3])/r[0];
    result.i33=r[4]+r[5]-(r[1]*r[1]+r[2]*r[2])/r[0];
    result.i12=-r[9]+(r[1]*r[2])/r[0];
    result.i13=-r[8]+(r[1]*r[3])/r[0];
    result.i23=-r[7]+(r[2]*r[3])/r[0];
    return result;
}
template<class Scalar>
Scalar transverseInertiaRatio(const CentralInertia10<Scalar>& inertia) {
    const Scalar numerator=inertia.i11*inertia.i11+inertia.i22*inertia.i22+
                           2.0*(inertia.i12*inertia.i12);
    return numerator/(inertia.i33*inertia.i33);
}
}} // namespace spin::detail
