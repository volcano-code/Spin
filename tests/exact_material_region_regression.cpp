#include "verification/ExactMaterialRegion.h"
#include <iostream>
#include <stdexcept>
#include <type_traits>

using namespace verification;
void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}
int main() {
    static_assert(!std::is_default_constructible<ExactMaterialRegion>::value,
                  "callers cannot construct uncomputed material regions");
    try {
        Eigen::MatrixXd v(4,3); v << 0,0,0,1,0,0,0,2,0,0,0,3;
        Eigen::MatrixXi f(4,3); f << 0,2,1,0,1,3,0,3,2,1,2,3;
        const auto a=makeExactMaterialRegion(v,f);
        require(a.facts().valid,a.facts().error.c_str());
        require(a.facts().materialComponents==1 && a.facts().boundedVoidComponents==0,"tetra topology");
        require(a.boundaryFacts().componentSignedVolumesExact==std::vector<std::string>{"1/1"},"original signed volume retained");
        const auto same=relateExactMaterialRegions(a,a);
        require(same.completed && same.equal==true && same.positiveVolumeOverlap==true,"nonempty self relation");
        const auto empty=emptyExactMaterialRegion();
        const auto algebra=relateExactMaterialRegions(empty,a);
        require(algebra.completed && algebra.aSubsetB==true && algebra.bSubsetA==false &&
                algebra.intersectionNonempty==false,"empty algebra");
        Eigen::MatrixXi reversed=f.rowwise().reverse();
        const auto invalid=makeExactMaterialRegion(v,reversed);
        require(!invalid.facts().valid && !invalid.facts().materialComponents.has_value(),"negative root invalid/unknown");
        bool rejected=false;
        try { relateExactMaterialRegions(invalid,empty); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"invalid operand must never yield clearance");
        rejected=false;
        try { makeExactMaterialRegion(Eigen::MatrixXd(0,3),Eigen::MatrixXi(0,3)); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"malformed empty not algebraic empty");
        // Coplanar shared face: reflection through z=0, then fix winding.
        Eigen::MatrixXd below=v; below.col(2)*=-1;
        const auto contact=relateExactMaterialRegions(a,makeExactMaterialRegion(below,reversed));
        require(contact.completed && contact.intersectionNonempty==true &&
                contact.positiveVolumeOverlap==false,"face-only contact retains boundary");
        v.setZero(); f.setZero();
        require(relateExactMaterialRegions(a,a).equal==true,"operand owns original immutable geometry");
        std::cout << "exact material native checks passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
