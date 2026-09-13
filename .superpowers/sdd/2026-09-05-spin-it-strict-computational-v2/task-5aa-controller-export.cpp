// Fixed existing synthetic fixture helpers; do not execute its regression main.
#define main task5z_unused_regression_main
#include "spin_it_2014_source_cage_material_regions_regression.cpp"
#undef main
#include "spin2014/SourceCageComponentClearance.h"
#include "task-5z-controller-export-helpers.h"

template<class T> void optionalArrayJson(const std::optional<T>& a) {
    if(a)arrayJson(*a);else std::cout<<"null";
}
void meshJson(const SourceCageBinaryMesh& m) {
    std::cout<<"{\"vertices\":";matrixJson(m.vertices);
    std::cout<<",\"faces\":";matrixJson(m.triangles);std::cout<<'}';
}
void dump(const char* id,const SourceCageMaterialRegionsAudit& input,
          const SourceCageBinaryMaterial& binary,const SourceCageComponentClearanceSettings& settings) {
    const auto audit=validateSourceCageComponentClearance(input,settings);
    const auto& f=audit.facts();const auto product=audit.validatedClearance();
    std::cout<<"{\"id\":"<<std::quoted(id)<<",\"scope\":\"KNOWN_SOLUTION_VALIDATION\""
        <<",\"sealed\":false,\"paper_reproduction\":false,\"outcome\":"<<int(f.outcome)
        <<",\"token\":"<<std::quoted(sourceCageComponentClearanceOutcomeToken(f.outcome))
        <<",\"published\":"<<bool(product)<<",\"same_input\":"<<(product&& &product->regionsAudit().facts()==&input.facts())
        <<",\"same_binary\":"<<(product&&product->regions()->binaryMaterial().get()==&binary)
        <<",\"admitted\":"<<f.admitted<<",\"identities\":"<<f.identitiesChecked<<",\"coverage\":"<<f.coverageChecked
        <<",\"metric_scope\":"<<std::quoted(f.scope)<<",\"full_thickness\":"<<std::quoted(f.fullManufacturingThickness)
        <<",\"state\":"<<std::quoted(f.state)<<",\"minimum\":"<<settings.minimumComponentSeparation
        <<",\"nef_limit\":"<<settings.maximumNefComplexity<<",\"nef_peak\":"<<f.nefBudget.peakComplexity
        <<",\"nef_exceeded\":"<<f.nefBudget.exceeded<<",\"input_hash\":"<<std::quoted(f.inputSha256)
        <<",\"upstream_hash\":"<<std::quoted(input.facts().scientificSha256)
        <<",\"upstream_validated\":"<<bool(input.validatedRegions())
        <<",\"material_input\":"<<std::quoted(f.materialInputSha256)
        <<",\"scientific_hash\":"<<std::quoted(f.scientificSha256)
        <<",\"component_count\":"<<f.componentCount<<",\"pair_count\":"<<f.plannedPairCount
        <<",\"face_components\":";arrayJson(f.sourceFaceComponents);
    std::cout<<",\"material\":";meshJson(binary.materialBoundary());
    std::cout<<",\"source\":";dumpSource(binary.functional(),binary.binaryState());
    std::cout<<",\"components\":[";
    for(std::size_t i=0;i<f.components.size();++i) {
        if(i)std::cout<<',';const auto& c=f.components[i];
        std::cout<<"{\"id\":"<<c.componentId<<",\"vertices\":";arrayJson(c.sourceVertexIds);
        std::cout<<",\"faces\":";arrayJson(c.sourceFaceIds);
        std::cout<<",\"sign\":"<<c.originalSign<<",\"parent\":"<<c.parent<<",\"depth\":"<<c.depth
            <<",\"reversed\":"<<c.checkerReversed<<",\"original_input\":"<<std::quoted(c.originalInputSha256)
            <<",\"checker_input\":"<<std::quoted(c.checkerInputSha256)<<",\"region\":";regionJson(c.region);std::cout<<'}';
    }
    std::cout<<"],\"pairs\":[";
    for(std::size_t i=0;i<f.pairs.size();++i) {
        if(i)std::cout<<',';const auto& p=f.pairs[i];const auto& c=p.clearance;
        std::cout<<"{\"components\":";arrayJson(p.componentIds);
        std::cout<<",\"completed\":"<<c.completed<<",\"applicable\":";optionalJson(c.applicable);
        std::cout<<",\"contact\":";optionalJson(c.contact);std::cout<<",\"at_least\":";optionalJson(c.atLeastMinimum);
        std::cout<<",\"minimum_exact\":"<<std::quoted(c.minimumClearanceExact)
            <<",\"a_input\":"<<std::quoted(c.aInputDigest)<<",\"b_input\":"<<std::quoted(c.bInputDigest)
            <<",\"a\":"<<std::quoted(c.aDigest)<<",\"b\":"<<std::quoted(c.bDigest)
            <<",\"distance\":";optionalJson(c.squaredDistanceExact);
        std::cout<<",\"local_witness\":";optionalArrayJson(c.witnessFaces);
        std::cout<<",\"parent_witness\":";optionalArrayJson(p.parentWitnessFaces);std::cout<<'}';
    }
    std::cout<<"],\"applicable\":";optionalJson(f.applicable);
    std::cout<<",\"contact\":";optionalJson(f.contact);std::cout<<",\"at_least\":";optionalJson(f.atLeastMinimum);
    std::cout<<",\"distance\":";optionalJson(f.minimumSquaredDistanceExact);
    std::cout<<",\"minimum_pair\":";optionalJson(f.minimumPairIndex);
    std::cout<<",\"minimum_components\":";optionalArrayJson(f.minimumComponentIds);
    std::cout<<",\"minimum_local_witness\":";optionalArrayJson(f.minimumLocalWitnessFaces);
    std::cout<<",\"minimum_parent_witness\":";optionalArrayJson(f.minimumParentWitnessFaces);
    std::cout<<"}\n";
}
int main(){try {
    std::cout<<std::setprecision(17)<<std::boolalpha;
    for(int i=0;i<3;++i) {
        const auto f=fixture(4,i==2);auto x=state(f,i==1?1.:0.);
        if(i==2) {
            const auto& leaves=f.source().boundary().binding().domain().leaves();
            for(std::size_t j=0;j<leaves.size();++j)
                if(leaves[j].coordinates[0]==6||leaves[j].coordinates[0]==9)x[j]=1.;
        }
        const auto b=validateSourceCageBinaryMaterial(endpoint(f,x));
        require(bool(b.validatedMaterial()),"actual validated binary fixture");
        SourceCageMaterialRegionsSettings rs;rs.minimumOuterWall=.05;
        const auto r=validateSourceCageMaterialRegions(b,rs);require(bool(r.validatedRegions()),"actual exact region fixture");
        SourceCageComponentClearanceSettings cs;cs.minimumComponentSeparation=.05;
        dump(i==0?"empty":i==1?"one_cavity":"two_cavities",r,*b.validatedMaterial(),cs);
        if(i==2) {cs.minimumComponentSeparation=.2;dump("two_cavities_rejected",r,*b.validatedMaterial(),cs);}
        if(i==1) {
            cs.minimumComponentSeparation=.2;dump("threshold_rejected",r,*b.validatedMaterial(),cs);
            cs.minimumComponentSeparation=.05;cs.maximumNefComplexity=1;
            dump("resource_limited",r,*b.validatedMaterial(),cs);
            cs.maximumNefComplexity=1000000;rs.minimumOuterWall=.2;
            const auto rejected=validateSourceCageMaterialRegions(b,rs);
            require(!rejected.validatedRegions(),"real upstream rejected");
            dump("not_admitted",rejected,*b.validatedMaterial(),cs);
        }
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
