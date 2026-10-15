#include "spin/detail/RankTwoAlgebraic.h"
#include "spin/detail/RankTwoObservedCartesian.h"
#include "spin/RankTwoSource.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
namespace {
using C=spin::detail::CertifiedAlgebraicReal;
using K=spin::detail::RankTwoObservedCartesian;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
// Actual recorded predicates from native constraint23. Not generated to suit
// the fix: coordinates retain original certificate identities/ordinals.
const std::array<const char*,4> coordinates{
"COORD 878addad331b75d72e8a1f526e41493a1a52b24bbe76164ecdef2e662ae120b3 12 127/128 257/256 1135621175606454381287569,0,-14728846771971196802324436,0,74954940150452821641802396,0,-197452685239399446510108672,0,299280892452829089062977536,0,-270436611218602258379833344,0,143332618753540348060893184,0,-40921666757306910599282688,0,4835703278458516698824704,",
"COORD 931d9ec8beaf418d03fd71c40cd0c74dd753cf570c2ea21622fbaa269c34083b 4 -283/256 -35/32 -12874169408340480301974244978533596475861281250,0,-91168941603779598820500560882597313353932800000,0,546013254309289904776596612111908076081468560625,0,1984992374872406092032450540618137206007427072000,0,-10504451722949725691120665572436393876469789491200,0,-1957439090330625147554285068792997420581351587840,0,61741327006717256209125965580491379602636287246336,0,-89038870478104082700840788380085684092062085939200,0,1349385625118595957407604883329233887370609164288,0,91974502412032932135200017332146101839658022338560,0,-67693584794564218068027373507699812388395376705536,0,-2923247626894165840699849756390373102387975946240,0,23508443317885195048289161471750176310060386353152,0,-10340409533654679338267281848661592263092999290880,0,1461501637330902918203684832716283019655932542976,",
"COORD dc78579ebe28e146ad85a4b4dfc79c627a7b82eb9ca15eaebf39d389f15fb181 12 255/256 129/128 13350437130625,0,-181622512362900,0,960335237004508,0,-2608267047736320,0,4055178823139328,0,-3745189503959040,0,2022985430990848,0,-587139209232384,0,70368744177664,",
"COORD 3e3702e81f542eea2bb8ecc50af783ad795cddd78b4efe72f0d554f22497e2d6 4 -141/128 -279/256 -1824048070005489726793122,0,-24032546954099816984497152,0,82633556318631280787198673,0,587068578421954164242523648,0,-2051040853128531878032932864,0,-2234261763674422808017895424,0,15904530455565734118827229184,0,-19023942250183401660487827456,0,-3049909545894960402733727744,0,22278273587850391722022404096,0,-14360263475982990752480231424,0,-1613591011341927480771477504,0,5461209053313498425250021376,0,-2263109134318585815049961472,0,309485009821345068724781056,"
};
C read(const char* row){
    std::istringstream input(row);std::string tag,polynomial;
    spin::detail::AlgebraicRootCertificate c;
    input>>tag>>c.rootIdentity>>c.rootIndex>>c.lower>>c.upper>>polynomial;
    require(bool(input) && tag=="COORD","invalid captured root record");
    std::istringstream coefficients(polynomial);std::string coefficient;
    while(std::getline(coefficients,coefficient,','))c.coefficients.push_back(coefficient);
    return C::reload(c);
}
}
int main(int argc,char** argv){try{
    require(argc==2,"predicate replay evidence path required");
    spin::detail::requireExactAlgebraicEnvironment();
    const auto px=read(coordinates[0]),py=read(coordinates[1]),qx=read(coordinates[2]),qy=read(coordinates[3]);
    const K::Point_2 p(px.expression(),py.expression()),q(qx.expression(),qy.expression());
    const auto orientation=K().orientation_2_object();
    std::cerr<<"REPLAY_BEGIN ORIENTATION2 1429 702 702\n";
    require(orientation(p,K::Point_2(q),K::Point_2(q))==CGAL::COLLINEAR,"actual duplicated point did not yield exact zero");
    std::cerr<<"REPLAY_END ORIENTATION2 1429 702 702 0\n";
    for(const auto& triplet:std::array<std::array<K::Point_2,3>,6>{{{p,q,q},{q,p,q},{q,q,p},{p,p,q},{p,q,p},{q,p,p}}})
        require(orientation(triplet[0],triplet[1],triplet[2])==CGAL::COLLINEAR,"duplicate-point permutation changed sign");
    CORE::BigInt denominator(1);denominator<<=100;
    const CORE::Expr delta(CORE::BigRat(CORE::BigInt(1),denominator));
    const K::Point_2 above(q.x(),q.y()+delta),below(q.x(),q.y()-delta);
    require(orientation(p,q,above)==CGAL::LEFT_TURN && orientation(p,q,below)==CGAL::RIGHT_TURN,"distinct near point was incorrectly welded to structural zero");
    spin::detail::requireExactAlgebraicEnvironment();
    std::ofstream out(argv[1]);require(bool(out),"cannot save replay evidence");
    out<<"rank-two-actual-finite-predicate-replay-v1\nactual_constraint 23\nactual_points 1429 702 702\n";
    for(const auto row:coordinates)out<<row<<'\n';
    out<<"orientation 0\nduplicate_permutations 6\nnear_distinct_offset 1/1267650600228229401496703205376\nnear_distinct_signs 1 -1\n";
    out.flush();require(bool(out),"incomplete replay evidence");std::cout<<"PASS actual high-degree predicate replay and nonzero near-point controls\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
