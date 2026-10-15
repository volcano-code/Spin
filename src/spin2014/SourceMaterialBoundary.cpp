#include "spin2014/SourceMaterialBoundary.h"
#include "SourceAdaptiveSpinEndpoint.h"
#include "verification/Sha256.h"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace spin2014 {
namespace {
std::string digest(std::initializer_list<std::string> fields) {
    verification::Sha256 hash;
    const unsigned char zero=0;
    for(const auto& field:fields) {
        hash.update(reinterpret_cast<const unsigned char*>(field.data()),field.size());
        hash.update(&zero,1);
    }
    return hash.finish();
}
std::string scalar(double value) {
    std::ostringstream out;out.imbue(std::locale::classic());out<<std::hexfloat<<value;
    return out.str();
}
std::string productDigest(const SourceMaterialBoundaryFacts& f) {
    std::ostringstream evidence;
    evidence.imbue(std::locale::classic());evidence<<std::hexfloat;
    for(int i=0;i<10;++i)
        evidence<<f.sourceSubtractionMoments.values()(i)<<','<<f.integratedMoments.values()(i)
            <<','<<f.boundaryMaterialMoments.values()(i)<<','<<f.interiorMaterialMoments.values()(i)
            <<','<<f.normalizedMomentResiduals(i)<<';';
    evidence<<f.normalizedMomentError<<';'<<f.minimumCentralInertiaEigenvalue<<';'
        <<f.materialComponentsExported<<';'<<std::hexfloat<<f.boundaryDensity<<','<<f.interiorDensity<<';';
    for(const auto* r:{&f.outerRegion,&f.designRegion,&f.materialRegion,&f.carvedRegion,
                       &f.boundaryMaterialRegion,&f.interiorMaterialRegion}) {
        evidence<<r->regionDigest<<';'<<r->completed<<r->valid<<r->explicitEmpty<<';';
        if(r->materialComponents) evidence<<*r->materialComponents;
        evidence<<',';if(r->boundedVoidComponents) evidence<<*r->boundedVoidComponents;
        evidence<<';';for(int sign:r->originalSigns) evidence<<sign<<',';
        evidence<<';';for(int depth:r->depths) evidence<<depth<<',';
    }
    evidence<<f.clearance.state<<';'<<f.clearance.completed<<';';
    if(f.clearance.squaredDistanceExact) evidence<<*f.clearance.squaredDistanceExact;
    if(f.clearance.witnessFaces) evidence<<';'<<(*f.clearance.witnessFaces)[0]<<','<<(*f.clearance.witnessFaces)[1];
    const auto optionalBool=[&](const std::optional<bool>& value) {
        evidence<<';';if(value) evidence<<*value;else evidence<<"UNKNOWN";
    };
    optionalBool(f.clearance.applicable);optionalBool(f.clearance.contact);optionalBool(f.clearance.atLeastMinimum);
    optionalBool(f.partition.completePartition);optionalBool(f.dualMaterialPartition.completePartition);
    optionalBool(f.carvedSubsetDesign.aSubsetB);
    return digest({f.algorithm,f.inputDigest,f.status,f.diagnostic,f.unionFacts.outputDigest,
        f.outerMeshHash,f.materialMeshHash,f.partition.partitionDigest,f.carvedSubsetDesign.relationDigest,
        f.clearance.clearanceDigest,evidence.str(),verification::hashExactNefBudgetFacts(f.resourceBudget)});
}
void requireRegion(const verification::ExactMaterialRegion& r,const char* label) {
    if(!r.facts().completed || !r.facts().valid)
        throw std::runtime_error(std::string(label)+": "+r.facts().error);
}
void differenceMesh(const Eigen::MatrixXd& wholeVertices,const Eigen::MatrixXi& wholeFaces,
    const Eigen::MatrixXd& subtractVertices,const Eigen::MatrixXi& subtractFaces,
    Eigen::MatrixXd& vertices,Eigen::MatrixXi& faces) {
    vertices.resize(wholeVertices.rows()+subtractVertices.rows(),3);
    faces.resize(wholeFaces.rows()+subtractFaces.rows(),3);
    vertices<<wholeVertices,subtractVertices;
    faces.topRows(wholeFaces.rows())=wholeFaces;
    faces.bottomRows(subtractFaces.rows())=
        (subtractFaces.rowwise().reverse().array()+wholeVertices.rows()).matrix();
}
void geometricMoments(const SourceAdaptiveDomain& domain,const Eigen::VectorXd& binary,
    const Eigen::MatrixXd& boundaryVertices,const Eigen::MatrixXi& boundaryFaces,
    const Eigen::MatrixXd& interiorVertices,const Eigen::MatrixXi& interiorFaces,
    const Eigen::MatrixXd& compositeVertices,const Eigen::MatrixXi& compositeFaces,
    double boundaryDensity,double interiorDensity,bool splitMaterials,SourceMaterialBoundaryFacts& facts) {
    const auto& source=domain.source();
    facts.sourceSubtractionMoments=(source.outerSurfaceMoments()-source.designSurfaceMoments()).scaledDensity(boundaryDensity)
        +source.designSurfaceMoments().scaledDensity(interiorDensity)
        -domain.removedMoments(binary).scaledDensity(interiorDensity);
    // Independently integrate both exported material boundaries; do not call
    // the optimizer's cell-moment columns on this validation path.
    MassMoments boundaryMoments,interiorMoments,moments;
    if(splitMaterials) {
        if(boundaryFaces.rows())boundaryMoments=computeSurfaceMassMoments(boundaryVertices,boundaryFaces,boundaryDensity);
        if(interiorFaces.rows())interiorMoments=computeSurfaceMassMoments(interiorVertices,interiorFaces,interiorDensity);
        facts.boundaryMaterialMoments=physics::RawMassMoments10::fromMassMoments(boundaryMoments);
        facts.interiorMaterialMoments=physics::RawMassMoments10::fromMassMoments(interiorMoments);
        moments=boundaryMoments+interiorMoments;
    } else {
        moments=computeSurfaceMassMoments(compositeVertices,compositeFaces,boundaryDensity);
    }
    facts.integratedMoments=physics::RawMassMoments10::fromMassMoments(moments);
    if(!moments.isFinite() || !facts.boundaryMaterialMoments.isFinite() || !facts.interiorMaterialMoments.isFinite() ||
       !facts.sourceSubtractionMoments.isFinite() || moments.mass<=0)
        throw std::runtime_error("nonfinite or nonpositive final material moments");
    const double side=domain.source().config().rootSide;
    const double densityScale=std::max(boundaryDensity,interiorDensity);
    const double scales[]={densityScale*std::pow(side,3),densityScale*std::pow(side,4),densityScale*std::pow(side,5)};
    for(double scale:scales)
        if(!std::isfinite(scale) || scale<=0) throw std::runtime_error("invalid fixed root moment scale");
    for(int i=0;i<10;++i)
        facts.normalizedMomentResiduals(i)=std::abs(facts.integratedMoments.values()(i)-facts.sourceSubtractionMoments.values()(i))/scales[i==0?0:(i<4?1:2)];
    facts.normalizedMomentError=facts.normalizedMomentResiduals.maxCoeff();
    const auto inertia=moments.inertiaAboutCenter();
    if(!inertia.allFinite()) throw std::runtime_error("nonfinite final central inertia");
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(inertia);
    if(eigen.info()!=Eigen::Success) throw std::runtime_error("final central inertia eigensolve failed");
    facts.minimumCentralInertiaEigenvalue=eigen.eigenvalues().minCoeff();
    if(!std::isfinite(facts.normalizedMomentError) || facts.normalizedMomentError>1e-10 ||
            facts.minimumCentralInertiaEigenvalue < -128*std::numeric_limits<double>::epsilon()*scales[2])
        throw std::runtime_error("final mesh moment/central-inertia gate failed");
}
}

void validateSourceMaterialBoundarySettings(const SourceMaterialBoundarySettings& s) {
    verification::validateExactCuboidUnionSettings(s.unionSettings);
    if(!std::isfinite(s.minimumOuterWall) || s.minimumOuterWall<0)
        throw std::invalid_argument("minimum outer wall must be finite nonnegative binary64");
}
std::string hashSourceMaterialBoundarySettings(const SourceMaterialBoundarySettings& s) {
    validateSourceMaterialBoundarySettings(s);
    return digest({"source-material-boundary-settings-v1","IMPLEMENTATION_DEFAULT",
        verification::hashExactCuboidUnionSettings(s.unionSettings),scalar(s.minimumOuterWall),
        "fixed-root-powers-3-4-5;moment-1e-10;inertia-128eps;strict-outer-noncontact"});
}

SourceMaterialBoundaryResult makeSourceMaterialBoundary(const SourceAdaptiveDomain& domain,
    const Eigen::VectorXd& binary,double density,const SourceMaterialBoundarySettings& settings) {
    return makeSourceMaterialBoundaryImpl(domain,binary,density,density,false,settings);
}

SourceMaterialBoundaryResult makeSourceMaterialBoundary(const SourceAdaptiveDomain& domain,
    const Eigen::VectorXd& binary,double boundaryDensity,double interiorDensity,
    const SourceMaterialBoundarySettings& settings) {
    return makeSourceMaterialBoundaryImpl(domain,binary,boundaryDensity,interiorDensity,true,settings);
}

SourceMaterialBoundaryResult makeSourceMaterialBoundaryImpl(const SourceAdaptiveDomain& domain,
    const Eigen::VectorXd& binary,double boundaryDensity,double interiorDensity,bool splitMaterials,
    const SourceMaterialBoundarySettings& settings) {
    validateSourceMaterialBoundarySettings(settings);
    if(!std::isfinite(boundaryDensity) || boundaryDensity<=0 || !std::isfinite(interiorDensity) || interiorDensity<=0 ||
            binary.size()!=static_cast<Eigen::Index>(domain.leaves().size()) ||
            !binary.allFinite() || ((binary.array()!=0)&&(binary.array()!=1)).any())
        throw std::invalid_argument("source material requires positive material densities and finite binary source beta");
    SourceMaterialBoundaryResult result;
    auto& facts=result.facts_;
    facts.density=boundaryDensity==interiorDensity?boundaryDensity:0.;
    facts.boundaryDensity=boundaryDensity;facts.interiorDensity=interiorDensity;facts.minimumOuterWall=settings.minimumOuterWall;
    facts.sourceHash=domain.source().partitionSha256();facts.partitionHash=domain.partitionHash();
    std::string binaryText=std::to_string(binary.size())+":";
    for(Eigen::Index i=0;i<binary.size();++i) binaryText+=binary(i)==0?'0':'1';
    facts.binaryHash=digest({"source-material-binary-v1",facts.sourceHash,facts.partitionHash,binaryText});
    facts.inputDigest=digest({facts.algorithm,facts.sourceHash,facts.partitionHash,facts.binaryHash,
        scalar(boundaryDensity),scalar(interiorDensity),hashSourceMaterialBoundarySettings(settings)});
    const auto budget=std::make_shared<verification::ExactNefBudget>(settings.unionSettings.maximumIntermediateComplexity);
    const auto finishFacts=[&] {
        facts.resourceBudget=budget->facts();
        if(facts.resourceBudget.exceeded) facts.status="RESOURCE_LIMIT_UNKNOWN";
        facts.productDigest=productDigest(facts);
    };
    try {
        domain.validateCoverage();
        const auto output=sourceBinaryCuboids(domain,binary,settings.unionSettings.maximumInputCuboids);
        std::vector<verification::ExactCuboid> boxes;
        boxes.reserve(output.cuboids.size());
        for(const auto& box:output.cuboids) boxes.push_back({box.lower,box.upper});
        auto united=verification::makeExactCuboidUnionWithBudget(boxes,settings.unionSettings,budget);
        facts.unionFacts=united.facts();
        if(!united.region()) {
            facts.status=united.facts().status;facts.diagnostic=united.facts().diagnostic;
            finishFacts();return result;
        }
        SourceMaterialBoundary product;
        const auto& source=domain.source();
        product.outerVertices_=source.outerVertices();product.outerFaces_=source.outerFaces();
        product.cavityVertices_=united.vertices();product.cavityFaces_=united.faces();
        const auto outer=verification::makeExactMaterialRegionWithBudget(product.outerVertices_,product.outerFaces_,budget);
        requireRegion(outer,"original outer");
        const auto design=verification::makeExactMaterialRegionWithBudget(source.designVertices(),source.designFaces(),budget);
        requireRegion(design,"original design");
        facts.outerRegion=outer.facts();facts.designRegion=design.facts();facts.carvedRegion=united.region()->facts();
        facts.outerMeshHash=outer.facts().inputDigest;
        facts.carvedSubsetDesign=verification::relateExactMaterialRegionsWithBudget(*united.region(),design,budget);
        if(!facts.carvedSubsetDesign.completed || facts.carvedSubsetDesign.aSubsetB!=true)
            throw std::runtime_error("newly carved union is not an established subset of original design");
        const auto designInOuter=verification::relateExactMaterialRegionsWithBudget(design,outer,budget);
        if(!designInOuter.completed || designInOuter.aSubsetB!=true)
            throw std::runtime_error("design domain is not established as a subset of outer material");
        facts.clearance=verification::checkExactBoundaryClearance(outer,*united.region(),settings.minimumOuterWall);
        if(!facts.clearance.completed || (!boxes.empty() &&
                (facts.clearance.applicable!=true || facts.clearance.contact!=false || facts.clearance.atLeastMinimum!=true))) {
            facts.status="OUTER_CLEARANCE_REJECTED";facts.diagnostic="strict original-boundary separation/threshold not established";
            finishFacts();return result;
        }
        const auto vertexCount=product.outerVertices_.rows()+product.cavityVertices_.rows();
        const auto faceCount=product.outerFaces_.rows()+product.cavityFaces_.rows();
        if(static_cast<std::size_t>(vertexCount)>settings.unionSettings.maximumOutputVertices ||
                static_cast<std::size_t>(faceCount)>settings.unionSettings.maximumOutputTriangles)
            throw std::length_error("RESOURCE_LIMIT: final material mesh output");
        product.materialVertices_.resize(vertexCount,3);product.materialFaces_.resize(faceCount,3);
        product.materialVertices_<<product.outerVertices_,product.cavityVertices_;
        product.materialFaces_.topRows(product.outerFaces_.rows())=product.outerFaces_;
        // Complement each already SIGNED union triangle: outward union roots
        // become cavities, inward void shells become outward retained islands.
        product.materialFaces_.bottomRows(product.cavityFaces_.rows())=
            (product.cavityFaces_.rowwise().reverse().array()+product.outerVertices_.rows()).matrix();
        const auto emptyRegion=verification::emptyExactMaterialRegionWithBudget(budget);
        const auto shellEqualsEmpty=designInOuter.equal==true;
        if(!splitMaterials) {
            product.boundaryMaterialVertices_.resize(0,3);product.boundaryMaterialFaces_.resize(0,3);
            product.interiorMaterialVertices_.resize(0,3);product.interiorMaterialFaces_.resize(0,3);
        } else if(shellEqualsEmpty) {
            product.boundaryMaterialVertices_.resize(0,3);product.boundaryMaterialFaces_.resize(0,3);
            facts.boundaryMaterialRegion=emptyRegion.facts();
        } else {
            differenceMesh(product.outerVertices_,product.outerFaces_,source.designVertices(),source.designFaces(),
                product.boundaryMaterialVertices_,product.boundaryMaterialFaces_);
            const auto boundaryRegion=verification::makeExactMaterialRegionWithBudget(
                product.boundaryMaterialVertices_,product.boundaryMaterialFaces_,budget);
            requireRegion(boundaryRegion,"fixed boundary material");facts.boundaryMaterialRegion=boundaryRegion.facts();
        }
        if(!splitMaterials) {
            // The legacy one-material path intentionally does not claim
            // per-material component surfaces when regions touch at interfaces.
        } else if(facts.carvedSubsetDesign.equal==true) {
            product.interiorMaterialVertices_.resize(0,3);product.interiorMaterialFaces_.resize(0,3);
            facts.interiorMaterialRegion=emptyRegion.facts();
        } else if(united.region()->facts().explicitEmpty) {
            product.interiorMaterialVertices_=source.designVertices();product.interiorMaterialFaces_=source.designFaces();
            facts.interiorMaterialRegion=design.facts();
        } else {
            differenceMesh(source.designVertices(),source.designFaces(),product.cavityVertices_,product.cavityFaces_,
                product.interiorMaterialVertices_,product.interiorMaterialFaces_);
            const auto interiorRegion=verification::makeExactMaterialRegionWithBudget(
                product.interiorMaterialVertices_,product.interiorMaterialFaces_,budget);
            requireRegion(interiorRegion,"retained interior material");facts.interiorMaterialRegion=interiorRegion.facts();
        }
        const auto material=verification::makeExactMaterialRegionWithBudget(product.materialVertices_,product.materialFaces_,budget);
        facts.materialRegion=material.facts();facts.materialBoundary=material.boundaryFacts();
        requireRegion(material,"final signed material");
        facts.materialMeshHash=material.facts().inputDigest;
        facts.partition=verification::checkExactMaterialPartitionWithBudget(outer,{material,*united.region()},budget);
        if(!facts.partition.completed || facts.partition.completePartition!=true)
            throw std::runtime_error("exact inclusive final-material/new-carve partition not established");
        if(splitMaterials) {
            const auto shellRegion=shellEqualsEmpty?emptyRegion:
                verification::makeExactMaterialRegionWithBudget(product.boundaryMaterialVertices_,product.boundaryMaterialFaces_,budget);
            const auto interiorRegion=facts.carvedSubsetDesign.equal==true?emptyRegion:
                (united.region()->facts().explicitEmpty?design:
                 verification::makeExactMaterialRegionWithBudget(product.interiorMaterialVertices_,product.interiorMaterialFaces_,budget));
            facts.dualMaterialPartition=verification::checkExactMaterialPartitionWithBudget(
                outer,{shellRegion,interiorRegion,*united.region()},budget);
            if(!facts.dualMaterialPartition.completed || facts.dualMaterialPartition.completePartition!=true)
                throw std::runtime_error("exact shell/interior/carve material partition not established");
            facts.materialComponentsExported=true;
        }
        geometricMoments(domain,binary,product.boundaryMaterialVertices_,product.boundaryMaterialFaces_,
            product.interiorMaterialVertices_,product.interiorMaterialFaces_,product.materialVertices_,product.materialFaces_,
            boundaryDensity,interiorDensity,splitMaterials,facts);
        facts.status="MATERIAL_BOUNDARY_VALIDATED";
        facts.diagnostic="native exact material geometry and independent moments only; no solver/stability/manufacturing authority";
        finishFacts();
        product.facts_=facts;
        result.product_=std::shared_ptr<const SourceMaterialBoundary>(new SourceMaterialBoundary(std::move(product)));
    } catch(const std::length_error& e) {facts.status="RESOURCE_LIMIT_UNKNOWN";facts.diagnostic=e.what();}
    catch(const std::bad_alloc&) {facts.status="RESOURCE_LIMIT_UNKNOWN";facts.diagnostic="source material allocation failure";}
    catch(const std::exception& e) {facts.status="GEOMETRY_UNKNOWN";facts.diagnostic=e.what();}
    catch(...) {facts.status="GEOMETRY_UNKNOWN";facts.diagnostic="unknown source material operation failure";}
    if(!result.product_) finishFacts();
    return result;
}

SourceAdaptiveSpinMaterialResult solveSourceAdaptiveSpinMaterial(
    const SourceAdaptiveSpinRequest& request,const SourceMaterialBoundarySettings& settings) {
    std::shared_ptr<const SourceAdaptiveCheckpoint> checkpoint;
    return solveSourceAdaptiveSpinMaterialStages(request,settings,nullptr,checkpoint);
}

SourceAdaptiveSpinMaterialResult solveSourceAdaptiveSpinMaterialStages(
    const SourceAdaptiveSpinRequest& request,const SourceMaterialBoundarySettings& settings,
    const SourceAdaptiveStageControl* control,std::shared_ptr<const SourceAdaptiveCheckpoint>& checkpoint) {
    const auto settingsHash=hashSourceMaterialBoundarySettings(settings);
    SourceAdaptiveSpinMaterialResult result;
    detail::SourceAdaptiveSpinEndpoint endpoint;
    result.solve_=detail::solveSourceAdaptiveSpinEndpoint(request,&endpoint,control,&checkpoint);
    if(result.solve_.status==SourceAdaptiveSpinStatus::BinaryMomentValidated) {
        if(!endpoint.domain) throw std::logic_error("successful source solve missing private terminal domain");
        if(request.boundaryDensity)
            result.boundary_=makeSourceMaterialBoundary(*endpoint.domain,endpoint.binary,
                request.resolvedBoundaryDensity(),request.resolvedInteriorDensity(),settings);
        else
            result.boundary_=makeSourceMaterialBoundary(*endpoint.domain,endpoint.binary,request.density,settings);
    }
    result.outputDigest_=digest({"source-adaptive-spin-material-wrapper-v1",settingsHash,
        result.solve_.inputHash,result.solve_.traceHash,result.solve_.outputHash,
        result.solve_.finalPartitionHash,result.boundary_.facts().inputDigest,
        result.boundary_.facts().productDigest,result.boundary_.facts().status});
    return result;
}
}
