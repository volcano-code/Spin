#include <CGAL/basic.h>
#if !CGAL_USE_CORE
#error "Rank-two algebraic backend requires the pinned CORE configuration"
#endif
#include <CGAL/CORE_BigInt.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool ok,const char* message) {
    if(!ok)throw std::runtime_error(message);
}
}
int main() {
    try {
        using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
        using Arrangement=CGAL::Arrangement_2<Traits>;
        using Polynomial=Traits::Polynomial_2;
        Traits traits;
        const auto curve=traits.construct_curve_2_object();
        const auto x=CGAL::shift(Polynomial(1),1,0);
        const auto y=CGAL::shift(Polynomial(1),1,1);
        const auto s=x*x+y*y;
        Arrangement disk(&traits);
        CGAL::insert(disk,curve(s-1));
        require(disk.number_of_faces()==2,"hidden disk not recovered");
        Arrangement rings(&traits);
        CGAL::insert(rings,curve(4*s*s-12*s+7));
        require(rings.number_of_faces()==3,"hidden annulus not recovered");
        Arrangement isolated(&traits);
        CGAL::insert(isolated,curve(s));
        require(isolated.number_of_vertices()==1 && isolated.number_of_edges()==0
                && isolated.number_of_faces()==1,"isolated point lost");
        CGAL::insert(disk,curve(x));
        require(disk.number_of_faces()==4,"line/circle incidences lost");
        std::cout<<"PASS backend: disk=2 annulus=3 isolated=1 bisected=4\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n';return 1;
    }
}
