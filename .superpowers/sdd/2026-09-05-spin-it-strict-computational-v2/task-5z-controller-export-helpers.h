#pragma once
#include "spin2014/SourceCageMaterialRegions.h"
#include "task-5w-controller-export-helpers.h"

template<class T> void optionalJson(const std::optional<T>& v) {
    if(v)std::cout<<*v;else std::cout<<"null";
}
void optionalJson(const std::optional<std::string>& v) {
    if(v)std::cout<<std::quoted(*v);else std::cout<<"null";
}
void regionJson(const spin2014::SourceCageRegionEvidence& e) {
    const auto& r=e.region;
    std::cout<<"{\"input\":"<<std::quoted(r.inputDigest)<<",\"digest\":"<<std::quoted(r.regionDigest)
        <<",\"algorithm\":"<<std::quoted(r.algorithm)<<",\"backend\":"<<std::quoted(r.backend)
        <<",\"semantics\":"<<std::quoted(r.semantics)<<",\"completed\":"<<r.completed
        <<",\"valid\":"<<r.valid<<",\"empty\":"<<r.explicitEmpty<<",\"material_components\":";
    optionalJson(r.materialComponents);std::cout<<",\"void_components\":";optionalJson(r.boundedVoidComponents);
    std::cout<<",\"parents\":";arrayJson(r.parents);std::cout<<",\"depths\":";arrayJson(r.depths);
    std::cout<<",\"signs\":";arrayJson(r.originalSigns);
    std::cout<<",\"boundary_input\":";
    if(e.boundary)std::cout<<std::quoted(e.boundary->inputDigest);else std::cout<<"null";
    std::cout<<'}';
}
void partitionJson(const verification::ExactMaterialPartitionFacts& p) {
    std::cout<<"{\"whole\":"<<std::quoted(p.wholeDigest)<<",\"digest\":"<<std::quoted(p.partitionDigest)
        <<",\"completed\":"<<p.completed<<",\"whole_empty\":"<<p.wholeExplicitEmpty<<",\"parts\":[";
    for(std::size_t i=0;i<p.partDigests.size();++i){if(i)std::cout<<',';std::cout<<std::quoted(p.partDigests[i]);}
    std::cout<<"],\"parts_empty\":";arrayJson(p.partExplicitEmpty);
    std::cout<<",\"complete_partition\":";optionalJson(p.completePartition);
    std::cout<<",\"equal_union\":";optionalJson(p.unionEqualsWhole);
    std::cout<<",\"no_overlap\":";optionalJson(p.noPositiveVolumeOverlap);
    std::cout<<",\"missing\":";optionalJson(p.missingWholeContent);
    std::cout<<",\"excess\":";optionalJson(p.excessUnionContent);
    std::cout<<",\"subsets\":";if(p.partSubsetWhole)arrayJson(*p.partSubsetWhole);else std::cout<<"null";
    std::cout<<'}';
}
void regionsJson(const spin2014::SourceCageMaterialRegionsAudit& audit,const spin2014::SourceCageBinaryMaterialAudit& binary) {
    const auto& f=audit.facts();const auto& c=f.outerClearance;
    const auto product=audit.validatedRegions();
    std::cout<<"{\"outcome\":"<<std::quoted(spin2014::sourceCageMaterialRegionsOutcomeToken(f.outcome))
        <<",\"published\":"<<bool(product)<<",\"same_binary\":"<<(product&&product->binaryMaterial()==binary.validatedMaterial())
        <<",\"admitted\":"<<f.admitted<<",\"identities\":"<<f.identitiesChecked
        <<",\"design_member\":"<<f.representedDesignMembershipEstablished
        <<",\"minimum\":"<<f.settings.minimumOuterWall<<",\"nef_limit\":"<<f.settings.maximumNefComplexity
        <<",\"full_thickness\":"<<std::quoted(f.fullManufacturingThickness)
        <<",\"scientific_hash\":"<<std::quoted(f.scientificSha256)
        <<",\"binary_input_hash\":"<<std::quoted(f.inputSha256)
        <<",\"partition_hash\":"<<std::quoted(f.partitionSha256)
        <<",\"source_hash\":"<<std::quoted(f.sourceSha256)
        <<",\"evaluation_hash\":"<<std::quoted(f.evaluationSha256)
        <<",\"outer\":";regionJson(f.outer);std::cout<<",\"carved\":";regionJson(f.carved);
    std::cout<<",\"material\":";regionJson(f.material);
    std::cout<<",\"material_partition\":";partitionJson(f.materialPartition);
    std::cout<<",\"selected_partition\":";partitionJson(f.selectedPartition);
    std::cout<<",\"selected_count\":"<<f.selectedCount<<",\"selected\":[";
    for(std::size_t i=0;i<f.selected.size();++i) {
        if(i)std::cout<<',';const auto& s=f.selected[i];
        std::cout<<"{\"primitive\":"<<s.primitiveId<<",\"owner\":"<<s.betaOwner<<",\"beta\":"<<s.binaryBeta
            <<",\"level\":"<<s.sourceKey.level<<",\"key\":";arrayJson(s.sourceKey.coordinates);
        std::cout<<",\"mesh_hash\":"<<std::quoted(s.meshSha256)<<",\"region\":";regionJson(s.region);
        std::cout<<",\"subset_completed\":"<<s.subsetOuter.completed<<",\"subset\":";optionalJson(s.subsetOuter.aSubsetB);
        std::cout<<",\"subset_a\":"<<std::quoted(s.subsetOuter.aDigest)<<",\"subset_b\":"<<std::quoted(s.subsetOuter.bDigest)<<'}';
    }
    std::cout<<"],\"clearance\":{\"completed\":"<<c.completed<<",\"state\":"<<std::quoted(c.state)
        <<",\"a\":"<<std::quoted(c.aDigest)<<",\"b\":"<<std::quoted(c.bDigest)
        <<",\"a_input\":"<<std::quoted(c.aInputDigest)<<",\"b_input\":"<<std::quoted(c.bInputDigest)
        <<",\"minimum_exact\":"<<std::quoted(c.minimumClearanceExact)<<",\"applicable\":";
    optionalJson(c.applicable);std::cout<<",\"contact\":";optionalJson(c.contact);
    std::cout<<",\"at_least\":";optionalJson(c.atLeastMinimum);
    std::cout<<",\"squared_distance\":";optionalJson(c.squaredDistanceExact);
    std::cout<<",\"witness\":";if(c.witnessFaces)arrayJson(*c.witnessFaces);else std::cout<<"null";
    std::cout<<"},\"nef_peak\":"<<f.nefBudget.peakComplexity<<",\"nef_exceeded\":"<<f.nefBudget.exceeded<<'}';
}
