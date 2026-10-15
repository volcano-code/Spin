#include "verification/ExactSharedVertexIntersection.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#ifdef TEST_FILTERED_SHARED_VERTEX_KERNEL
#include <CGAL/Filtered_kernel.h>
#endif
#include <CGAL/intersections.h>
#include <array>
#include <iostream>
#include <random>
#include <stdexcept>

#ifdef TEST_FILTERED_SHARED_VERTEX_KERNEL
using Kernel=CGAL::Filtered_kernel<CGAL::Simple_cartesian<CGAL::Gmpq>>;
#else
using Kernel=CGAL::Simple_cartesian<CGAL::Gmpq>;
#endif
using Point=Kernel::Point_3;
using FT=Kernel::FT;

static void require(bool value,const char* reason){if(!value)throw std::runtime_error(reason);}

static bool reference(const Point& p,const Point& a,const Point& b,const Point& c,const Point& d) {
    const auto value=CGAL::intersection(Kernel::Triangle_3(p,a,b),Kernel::Triangle_3(p,c,d));
    if(!value)return false;
    if(const auto* x=boost::get<Point>(&*value))return *x!=p;
    if(const auto* x=boost::get<Kernel::Segment_3>(&*value))return x->source()!=p || x->target()!=p;
    if(const auto* x=boost::get<Kernel::Triangle_3>(&*value))
        return x->vertex(0)!=p || x->vertex(1)!=p || x->vertex(2)!=p;
    if(const auto* x=boost::get<std::vector<Point>>(&*value)) {
        for(const auto& vertex:*x)if(vertex!=p)return true;
        return false;
    }
    throw std::runtime_error("unrecognized exact reference variant");
}

static void check(const Point& p,const Point& a,const Point& b,const Point& c,const Point& d,bool expected) {
    require(reference(p,a,b,c,d)==expected,"literal expected intersection disagrees with independent CGAL");
    for(int swapA=0;swapA<2;++swapA)for(int swapB=0;swapB<2;++swapB) {
        const auto& first=swapA?b:a;const auto& second=swapA?a:b;
        const auto& third=swapB?d:c;const auto& fourth=swapB?c:d;
        require(verification::detail::sharedVertexHasExcessIntersection<Kernel>(p,first,second,third,fourth)==expected,
                "shared-vertex predicate wrong or orientation dependent");
        require(verification::detail::sharedVertexHasExcessIntersection<Kernel>(p,third,fourth,first,second)==expected,
                "triangle exchange changed the result");
    }
}

int main() {
    try {
        const Point p(0,0,0),a(1,0,0),b(0,1,0);
        check(p,a,b,Point(-1,0,0),Point(0,-1,0),false);
        check(p,a,b,Point(1,0,0),Point(1,1,0),true);
        check(p,a,b,Point(0,-1,0),Point(1,0,0),true);
        check(p,a,b,Point(10,10,0),Point(-10,10,0),true);
        check(p,a,b,Point(0,0,1),Point(1,1,1),false);
        check(p,a,b,Point(1,0,0),Point(0,-1,1),true);
        check(p,a,b,Point(-1,0,0),Point(0,-1,1),false);
        check(p,a,b,Point(-1,1,-1),Point(-1,1,1),false);
        check(p,a,b,Point(1,1,-1),Point(1,1,1),true);
        CGAL::Gmpz denominator(1);for(int i=0;i<200;++i)denominator*=2;
        const FT delta(CGAL::Gmpz(1),denominator);
        check(p,a,b,Point(1,1,-1),Point(-1,-1,FT(1)+delta),true);
        check(p,a,b,Point(-1,-1,-1),Point(1,1,FT(1)+delta),false);
        const Point shift(FT(denominator),-FT(denominator),FT(denominator));
        const auto translation=shift-p;
        check(shift,a+translation,b+translation,Point(1,1,-1)+translation,
              Point(-1,-1,FT(1)+delta)+translation,true);
        CGAL::Gmpz extreme(denominator);for(int i=200;i<1200;++i)extreme*=2;
        const FT scale(CGAL::Gmpz(1),extreme);
        const Point offset(FT(extreme),-FT(extreme),FT(extreme));
        const auto shifted=[&](const Point& value){return offset+scale*(value-p);};
        check(offset,shifted(a),shifted(b),shifted(Point(1,1,-1)),
              shifted(Point(-1,-1,FT(1)+delta)),true);
        std::mt19937 generator(0);std::uniform_int_distribution<int> integer(-4,4);
        int verified=0;
        while(verified<256) {
            std::array<Point,4> vertices;
            for(auto& vertex:vertices)vertex=Point(integer(generator),integer(generator),integer(generator));
            if(CGAL::collinear(p,vertices[0],vertices[1]) || CGAL::collinear(p,vertices[2],vertices[3]))continue;
            check(p,vertices[0],vertices[1],vertices[2],vertices[3],
                reference(p,vertices[0],vertices[1],vertices[2],vertices[3]));
            ++verified;
        }
        std::cout<<"PASS 13 literal, sub-binary64/extreme/translated, 256 randomized, permutation and exchange cases\n";
        return 0;
    } catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
