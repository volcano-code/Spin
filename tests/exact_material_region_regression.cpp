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
        // A certified single shell should not perform a second expensive
        // regularity overlay after union with empty. Its ownership transfer
        // must keep the live Nef charge until the last shared owner disappears.
        const auto singleBudget=std::make_shared<ExactNefBudget>(100000);
        std::optional<ExactMaterialRegion> retained;
        {
            const auto single=makeExactMaterialRegionWithBudget(v,f,singleBudget);
            require(single.facts().valid,single.facts().error.c_str());
            const auto& operations=singleBudget->facts().observationsByOperation;
            require(operations.at("symmetric difference")==1,
                    "single component repeats its exact regularity proof");
            require(operations.count("union")==0,"single component still unions with empty");
            require(singleBudget->facts().currentComplexity>0,"live solid lost its resource charge");
            retained.emplace(single);
        }
        require(singleBudget->facts().currentComplexity>0,"shared owner lost its resource charge");
        retained.reset();
        require(singleBudget->facts().currentComplexity==0,"transferred solid leaks its resource charge");
        const auto tinyBudget=std::make_shared<ExactNefBudget>(1);
        {
            const auto limited=makeExactMaterialRegionWithBudget(v,f,tinyBudget);
            require(!limited.facts().valid && tinyBudget->facts().exceeded,
                    "resource rejection was bypassed by single-component path");
            require(!limited.facts().materialComponents.has_value(),"rejected region acquired topology facts");
        }
        require(tinyBudget->facts().currentComplexity==0,"failed construction leaks its resource charge");
        const auto a=makeExactMaterialRegion(v,f);
        require(a.facts().valid,a.facts().error.c_str());
        require(a.facts().materialComponents==1 && a.facts().boundedVoidComponents==0,"tetra topology");
        require(a.boundaryFacts().componentSignedVolumesExact==std::vector<std::string>{"1/1"},"original signed volume retained");
        const auto same=relateExactMaterialRegions(a,a);
        require(same.completed && same.equal==true && same.positiveVolumeOverlap==true,"nonempty self relation");
        const auto selfOverlap=checkExactPositiveVolumeOverlap(a,a);
        require(selfOverlap.completed && selfOverlap.positiveVolumeOverlap==true,
                "narrow overlap predicate detects a positive-volume self-overlap");
        Eigen::MatrixXd nestedVertices=v;
        const Eigen::Vector3d centroid(.25,.5,.75);
        for(Eigen::Index row=0;row<v.rows();++row)
            nestedVertices.row(row)=centroid.transpose()+.25*(v.row(row)-centroid.transpose());
        const auto nested=makeExactMaterialRegion(nestedVertices,f);
        const auto nestedRelation=relateExactMaterialRegions(a,nested);
        const auto nestedOverlap=checkExactPositiveVolumeOverlap(a,nested);
        require(nestedRelation.completed && nestedRelation.positiveVolumeOverlap==true &&
                nestedOverlap.completed && nestedOverlap.positiveVolumeOverlap==true,
                "narrow overlap predicate detects strict containment as positive-volume overlap");
        const auto generalBudget=std::make_shared<ExactNefBudget>(100000);
        const auto narrowBudget=std::make_shared<ExactNefBudget>(100000);
        require(relateExactMaterialRegionsWithBudget(a,nested,generalBudget).positiveVolumeOverlap==true &&
                checkExactPositiveVolumeOverlapWithBudget(a,nested,narrowBudget).positiveVolumeOverlap==true,
                "general and narrow exact overlap paths agree");
        require(narrowBudget->facts().observationsByOperation.count("relation difference A-B")==0 &&
                narrowBudget->facts().observationsByOperation.count("relation symmetric difference")==0,
                "overlap-only query avoids unrelated containment and symmetric-difference Booleans");
        const auto empty=emptyExactMaterialRegion();
        const auto algebra=relateExactMaterialRegions(empty,a);
        require(algebra.completed && algebra.aSubsetB==true && algebra.bSubsetA==false &&
                algebra.intersectionNonempty==false,"empty algebra");
        const auto emptyOverlap=checkExactPositiveVolumeOverlap(empty,a);
        require(emptyOverlap.completed && emptyOverlap.positiveVolumeOverlap==false,
                "narrow overlap predicate treats explicit empty as non-overlapping");
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
        const auto contactRegion=makeExactMaterialRegion(below,reversed);
        const auto contactOverlap=checkExactPositiveVolumeOverlap(a,contactRegion);
        require(contactOverlap.completed && contactOverlap.positiveVolumeOverlap==false,
                "narrow overlap predicate distinguishes face contact from positive-volume overlap");
        require(contactOverlap.algorithm!="" && contactOverlap.operationDigest.size()==64,
                "narrow overlap evidence has an independent operation identity");
        v.setZero(); f.setZero();
        require(relateExactMaterialRegions(a,a).equal==true,"operand owns original immutable geometry");
        std::cout << "exact material native checks passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
