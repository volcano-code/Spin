#include "verification/ConservativeDiagonalSlabs.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
#include <iostream>
#include <stdexcept>

using K=CGAL::Simple_cartesian<CGAL::Gmpq>;
using P=K::Point_3;
using Slabs=verification::detail::ConservativeDiagonalSlabs<P>;
static void require(bool b,const char* reason){if(!b)throw std::runtime_error(reason);}

int main(){
    try {
        // Ordinary coordinate AABBs overlap; x+y proves strict separation.
        std::vector<P> p={{0,0,0},{2,0,0},{0,2,0},
            {2,2,0},{2,CGAL::Gmpq(1,2),0},{CGAL::Gmpq(1,2),2,0},
            {3,0,0},{2,1,0}};
        Slabs slabs(p);
        require(Slabs::separated(slabs.box({0,1,2}),slabs.box({3,4,5})),"diagonal separator not used");
        require(!Slabs::separated(slabs.box({0,1,2}),slabs.box({1,6,7})),"shared point discarded");
        CGAL::Gmpz huge(1);for(int i=0;i<1200;++i)huge*=2;
        const CGAL::Gmpq tiny(CGAL::Gmpz(1),huge),big(huge);
        p.push_back(P(1,tiny,0));p.push_back(P(1,-tiny,0));p.push_back(P(big,0,0));
        Slabs wide(p);
        // Each returned separator must agree with independent exact triangles,
        // including nonrepresentable doubles and contacts on the same plane.
        std::vector<std::array<int,3>> faces;
        for(int a=0;a<int(p.size());++a)for(int b=a+1;b<int(p.size());++b)
            for(int c=b+1;c<int(p.size());++c)
                if(!K::Triangle_3(p[a],p[b],p[c]).is_degenerate())faces.push_back({a,b,c});
        std::size_t compared=0;
        for(const auto& a:faces)for(const auto& b:faces){
            if(wide.separated(wide.box(a),wide.box(b)))
                require(!CGAL::do_intersect(K::Triangle_3(p[a[0]],p[a[1]],p[a[2]]),
                    K::Triangle_3(p[b[0]],p[b[1]],p[b[2]])),"false negative exact intersection");
            ++compared;
        }
        std::cout<<"PASS diagonal separation, closed contact, extreme rational values, "<<compared<<" exact pair comparisons\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
