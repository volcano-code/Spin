#include "verification/ExactMaterialRegion.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace verification;
void require(bool value,const char* message) {
    if (!value) throw std::runtime_error(message);
}
template<class F> void rejects(F f) {
    bool caught=false;
    try { f(); } catch (const std::invalid_argument&) { caught=true; }
    require(caught,"invalid argument must throw");
}
int main() {
    try {
        Eigen::MatrixXd v(4,3); v << 0,0,0,4,0,0,0,3,0,0,0,2;
        Eigen::MatrixXi f(4,3); f << 0,2,1,0,1,3,0,3,2,1,2,3;
        const auto a=makeExactMaterialRegion(v,f), e=emptyExactMaterialRegion();
        const auto p=checkExactMaterialPartition(a,{e,a});
        require(p.completed && p.completePartition==true && p.unionEqualsWhole==true &&
                p.noPositiveVolumeOverlap==true,"explicit empty plus self partitions");
        require(p.partSubsetWhole==std::vector<bool>({true,true}),"inclusive subsets computed");
        const auto duplicate=checkExactMaterialPartition(a,{a,a});
        require(duplicate.completed && duplicate.unionEqualsWhole==true && duplicate.completePartition==false &&
                duplicate.positiveVolumeOverlapPairs==std::vector<std::array<std::size_t,2>>({{0,1}}),
                "duplicate solid fails partition despite union equality");
        const auto missing=checkExactMaterialPartition(a,{e});
        require(missing.missingWholeContent==true && missing.excessUnionContent==false,"missing material");
        rejects([&]{checkExactMaterialPartition(a,{});});
        const auto empty=checkExactBoundaryClearance(a,e,1.);
        require(empty.completed && empty.applicable==false && empty.state=="NOT_APPLICABLE" &&
                !empty.squaredDistanceExact && !empty.witnessFaces && !empty.atLeastMinimum && !empty.contact,
                "empty distance has absent facts");
        const auto same=checkExactBoundaryClearance(a,a,0.);
        require(same.completed && same.applicable==true && same.contact==true && same.atLeastMinimum==true &&
                same.squaredDistanceExact=="0/1" && same.witnessFaces==std::array<std::size_t,2>({0,0}),
                "contact has exact zero and lexicographic witness; zero threshold passes");
        require(checkExactBoundaryClearance(a,a,1.).atLeastMinimum==false,"positive threshold fails contact");
        for (double threshold : {-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
            rejects([&]{checkExactBoundaryClearance(a,a,threshold);});
        const auto invalid=makeExactMaterialRegion(v,Eigen::MatrixXi(f.rowwise().reverse()));
        rejects([&]{checkExactMaterialPartition(a,{invalid});});
        rejects([&]{checkExactBoundaryClearance(invalid,a,0.);});
        v.setZero(); f.setZero();
        require(checkExactMaterialPartition(a,{e,a}).partitionDigest==p.partitionDigest &&
                checkExactBoundaryClearance(a,a,0.).clearanceDigest==same.clearanceDigest,
                "immutable operands and repeatable evidence");
        std::cout << "exact partition/clearance native checks passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
