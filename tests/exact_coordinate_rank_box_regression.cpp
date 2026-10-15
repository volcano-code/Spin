#include "verification/ExactCoordinateRankBox.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <stdexcept>

using Kernel=CGAL::Simple_cartesian<CGAL::Gmpq>;
using Point=Kernel::Point_3;
using FT=Kernel::FT;
using Ranks=verification::detail::ExactCoordinateRanks<Point>;

static void require(bool value,const char* reason) {
    if(!value)throw std::runtime_error(reason);
}

static bool referenceOverlap(const std::vector<Point>& points,
                             const std::array<int,3>& a,const std::array<int,3>& b) {
    for(int axis=0;axis<3;++axis) {
        FT al=points[a[0]][axis],ah=al,bl=points[b[0]][axis],bh=bl;
        for(int i=1;i<3;++i) {
            al=std::min(al,points[a[i]][axis]);ah=std::max(ah,points[a[i]][axis]);
            bl=std::min(bl,points[b[i]][axis]);bh=std::max(bh,points[b[i]][axis]);
        }
        if(ah<bl || bh<al)return false;
    }
    return true;
}

int main() {
    try {
        CGAL::Gmpz large(1);
        for(int i=0;i<1200;++i)large*=2;
        const FT huge(large),tiny(CGAL::Gmpz(1),large);
        const std::vector<Point> points={{-huge,0,0},{FT(1)-tiny,tiny,-tiny},
            {1,0,0},{FT(1)+tiny,-tiny,tiny},{huge,0,0}};
        Ranks rank(points);
        for(int i=0;i<5;++i)require(rank.rank(i,0)==i,"exact x ordering changed");
        require(rank.rank(0,1)==rank.rank(2,1) && rank.rank(2,1)==rank.rank(4,1),
                "equal coordinates did not share a rank");
        require(rank.rank(3,1)<rank.rank(0,1) && rank.rank(0,1)<rank.rank(1,1),
                "sub-binary64 differences were rounded");
        const auto a=rank.box({0,1,2}),b=rank.box({3,4,3}),c=rank.box({2,3,4});
        require(!verification::detail::exactRankBoxesOverlap(a,b),"exactly separated boxes overlap");
        require(verification::detail::exactRankBoxesOverlap(a,c),"closed-box contact was discarded");
        std::vector<std::array<int,3>> faces;
        for(int i=0;i<5;++i)for(int j=i+1;j<5;++j)for(int k=j+1;k<5;++k)faces.push_back({i,j,k});
        for(const auto& first:faces)for(const auto& second:faces)
            require(verification::detail::exactRankBoxesOverlap(rank.box(first),rank.box(second))==
                    referenceOverlap(points,first,second),"rank box differs from exact rational comparison");
        bool badAxis=false,badIndex=false,empty=false;
        try{rank.rank(0,3);}catch(const std::out_of_range&){badAxis=true;}
        try{rank.box({0,-1,2});}catch(const std::out_of_range&){badIndex=true;}
        try{Ranks absent(std::vector<Point>{});}catch(const std::invalid_argument&){empty=true;}
        require(badAxis && badIndex && empty,"invalid rank input was silently admitted");
        std::cout<<"PASS exact ordering, closed contact, 100 paired boxes, invalid inputs\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
