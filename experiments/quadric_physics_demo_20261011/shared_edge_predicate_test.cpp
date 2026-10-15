// The oracle constructs the complete intersection; the candidate may only use
// exact predicates. Both triangles are nondegenerate and share indexed edge AB.
#if __has_include("verification/ExactSharedEdgeIntersection.h")
#include "verification/ExactSharedEdgeIntersection.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
#include <iostream>
#include <stdexcept>

using K=CGAL::Simple_cartesian<CGAL::Gmpq>;
using P=K::Point_3;
using Q=CGAL::Gmpq;

bool oracle(const P& a,const P& b,const P& c,const P& d) {
    const auto hit=CGAL::intersection(K::Triangle_3(a,b,c),K::Triangle_3(b,a,d));
    if(!hit) throw std::runtime_error("shared edge was lost by oracle");
    const K::Segment_3 allowed(a,b);
    if(const auto* p=boost::get<P>(&*hit)) return !allowed.has_on(*p);
    if(const auto* s=boost::get<K::Segment_3>(&*hit)) return !allowed.has_on(s->source()) || !allowed.has_on(s->target());
    if(const auto* t=boost::get<K::Triangle_3>(&*hit))
        return !allowed.has_on(t->vertex(0)) || !allowed.has_on(t->vertex(1)) || !allowed.has_on(t->vertex(2));
    if(const auto* ring=boost::get<std::vector<P>>(&*hit)) {
        for(const auto& p:*ring) if(!allowed.has_on(p)) return true;
        return false;
    }
    throw std::runtime_error("unknown oracle intersection");
}

int main() {
    unsigned count=0;
    auto check=[&](P a,P b,P c,P d) {
        for(int reverse=0;reverse<2;++reverse) {
            if(verification::detail::sharedEdgeHasExcessIntersection<K>(a,b,c,d)!=oracle(a,b,c,d))
                throw std::runtime_error("shared-edge predicate disagrees with full exact construction");
            std::swap(a,b); ++count;
        }
    };
    const P a(0,0,0),b(1,0,0),c(0,1,0);
    check(a,b,c,P(0,-1,0));
    check(a,b,c,P(Q(1,3),Q(1,7),0));
    check(a,b,c,P(0,0,1));
    check(a,b,c,P(0,1,Q("1/100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001")));
    const P x(Q(1,3),Q(-2,7),Q(4,5)), y(Q(2,7),Q(3,5),Q(-7,11));
    for(int i=1;i<=80;++i) {
        const P z(Q(i%13,17),Q((i*7)%19,23),Q((i*11)%29,31));
        if(CGAL::collinear(x,y,z)) continue;
        const auto edge=y-x, side=z-x;
        check(x,y,z,x+edge*Q(i%5,7)+side*Q(1,3));
        check(x,y,z,x+edge*Q(i%5,7)-side*Q(1,3));
        const P w(Q((i*13)%23,29),Q((i*3)%17,19),Q((i*5)%31,37));
        if(!CGAL::collinear(x,y,w)) check(x,y,z,w);
    }
    bool rejected=false;
    try { verification::detail::sharedEdgeHasExcessIntersection<K>(a,b,P(2,0,0),c); }
    catch(const std::invalid_argument&) { rejected=true; }
    if(!rejected) throw std::runtime_error("degenerate triangle accepted");
    std::cout<<"full exact oracle comparisons passed: "<<count<<"\n";
}
#else
#include <iostream>
int main() { std::cerr<<"shared-edge exact predicate is missing\n"; return 1; }
#endif
