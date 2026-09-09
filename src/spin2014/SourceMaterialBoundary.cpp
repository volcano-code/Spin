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
            <<','<<f.normalizedMomentResiduals(i)<<';';
    evidence<<f.normalizedMomentError<<';'<<f.minimumCentralInertiaEigenvalue<<';';
    for(const auto* r:{&f.outerRegion,&f.designRegion,&f.materialRegion,&f.carvedRegion}) {
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
    optionalBool(f.partition.completePartition);optionalBool(f.carvedSubsetDesign.aSubsetB);
    return digest({f.algorithm,f.inputDigest,f.status,f.diagnostic,f.unionFacts.outputDigest,
        f.outerMeshHash,f.materialMeshHash,f.partition.partitionDigest,f.carvedSubsetDesign.relationDigest,
        f.clearance.clearanceDigest,evidence.str()});
}
void requireRegion(const verification::ExactMaterialRegion& r,const char* label) {
    if(!r.facts().completed || !r.facts().valid)
        throw std::runtime_error(std::string(label)+": "+r.facts().error);
}
void geometricMoments(const SourceAdaptiveDomain& domain,const Eigen::VectorXd& binary,
    const Eigen::MatrixXd& vertices,const Eigen::MatrixXi& faces,double density,
    SourceMaterialBoundaryFacts& facts) {
    facts.sourceSubtractionMoments=domain.source().outerSurfaceMoments().scaledDensity(density)
        -domain.removedMoments(binary).scaledDensity(density);
    // Independent signed tetrahedra on the ACTUAL final exported triangle mesh.
    // Does not call Algorithm1/oracle columns for this geometric integration.
    const auto moments=computeSurfaceMassMoments(vertices,faces,density);
    facts.integratedMoments=physics::RawMassMoments10::fromMassMoments(moments);
    if(!moments.isFinite() || !facts.sourceSubtractionMoments.isFinite() || moments.mass<=0)
        throw std::runtime_error("nonfinite or nonpositive final material moments");
    const double side=domain.source().config().rootSide;
    const double scales[]={density*std::pow(side,3),density*std::pow(side,4),density*std::pow(side,5)};
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
    validateSourceMaterialBoundarySettings(settings);
    if(!std::isfinite(density) || density<=0 || binary.size()!=static_cast<Eigen::Index>(domain.leaves().size()) ||
            !binary.allFinite() || ((binary.array()!=0)&&(binary.array()!=1)).any())
        throw std::invalid_argument("source material requires positive density and finite binary source beta");
    SourceMaterialBoundaryResult result;
    auto& facts=result.facts_;
    facts.density=density;facts.minimumOuterWall=settings.minimumOuterWall;
    facts.sourceHash=domain.source().partitionSha256();facts.partitionHash=domain.partitionHash();
    std::string binaryText=std::to_string(binary.size())+":";
    for(Eigen::Index i=0;i<binary.size();++i) binaryText+=binary(i)==0?'0':'1';
    facts.binaryHash=digest({"source-material-binary-v1",facts.sourceHash,facts.partitionHash,binaryText});
    facts.inputDigest=digest({facts.algorithm,facts.sourceHash,facts.partitionHash,facts.binaryHash,
        scalar(density),hashSourceMaterialBoundarySettings(settings)});
    try {
        domain.validateCoverage();
        const auto output=sourceBinaryCuboids(domain,binary,settings.unionSettings.maximumInputCuboids);
        std::vector<verification::ExactCuboid> boxes;
        boxes.reserve(output.cuboids.size());
        for(const auto& box:output.cuboids) boxes.push_back({box.lower,box.upper});
        auto united=verification::makeExactCuboidUnion(boxes,settings.unionSettings);
        facts.unionFacts=united.facts();
        if(!united.region()) {
            facts.status=united.facts().status;facts.diagnostic=united.facts().diagnostic;
            facts.productDigest=productDigest(facts);return result;
        }
        SourceMaterialBoundary product;
        const auto& source=domain.source();
        product.outerVertices_=source.outerVertices();product.outerFaces_=source.outerFaces();
        product.cavityVertices_=united.vertices();product.cavityFaces_=united.faces();
        const auto outer=verification::makeExactMaterialRegion(product.outerVertices_,product.outerFaces_);
        const auto design=verification::makeExactMaterialRegion(source.designVertices(),source.designFaces());
        requireRegion(outer,"original outer");requireRegion(design,"original design");
        facts.outerRegion=outer.facts();facts.designRegion=design.facts();facts.carvedRegion=united.region()->facts();
        facts.outerMeshHash=outer.facts().inputDigest;
        facts.carvedSubsetDesign=verification::relateExactMaterialRegions(*united.region(),design);
        if(!facts.carvedSubsetDesign.completed || facts.carvedSubsetDesign.aSubsetB!=true)
            throw std::runtime_error("newly carved union is not an established subset of original design");
        facts.clearance=verification::checkExactBoundaryClearance(outer,*united.region(),settings.minimumOuterWall);
        if(!facts.clearance.completed || (!boxes.empty() &&
                (facts.clearance.applicable!=true || facts.clearance.contact!=false || facts.clearance.atLeastMinimum!=true))) {
            facts.status="OUTER_CLEARANCE_REJECTED";facts.diagnostic="strict original-boundary separation/threshold not established";
            facts.productDigest=productDigest(facts);return result;
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
        const auto material=verification::makeExactMaterialRegion(product.materialVertices_,product.materialFaces_);
        facts.materialRegion=material.facts();facts.materialBoundary=material.boundaryFacts();
        requireRegion(material,"final signed material");
        facts.materialMeshHash=material.facts().inputDigest;
        facts.partition=verification::checkExactMaterialPartition(outer,{material,*united.region()});
        if(!facts.partition.completed || facts.partition.completePartition!=true)
            throw std::runtime_error("exact inclusive final-material/new-carve partition not established");
        geometricMoments(domain,binary,product.materialVertices_,product.materialFaces_,density,facts);
        facts.status="MATERIAL_BOUNDARY_VALIDATED";
        facts.diagnostic="native exact material geometry and independent moments only; no solver/stability/manufacturing authority";
        facts.productDigest=productDigest(facts);
        product.facts_=facts;
        result.product_=std::shared_ptr<const SourceMaterialBoundary>(new SourceMaterialBoundary(std::move(product)));
    } catch(const std::length_error& e) {facts.status="RESOURCE_LIMIT_UNKNOWN";facts.diagnostic=e.what();}
    catch(const std::bad_alloc&) {facts.status="RESOURCE_LIMIT_UNKNOWN";facts.diagnostic="source material allocation failure";}
    catch(const std::exception& e) {facts.status="GEOMETRY_UNKNOWN";facts.diagnostic=e.what();}
    catch(...) {facts.status="GEOMETRY_UNKNOWN";facts.diagnostic="unknown source material operation failure";}
    if(!result.product_) facts.productDigest=productDigest(facts);
    return result;
}

SourceAdaptiveSpinMaterialResult solveSourceAdaptiveSpinMaterial(
    const SourceAdaptiveSpinRequest& request,const SourceMaterialBoundarySettings& settings) {
    const auto settingsHash=hashSourceMaterialBoundarySettings(settings);
    SourceAdaptiveSpinMaterialResult result;
    detail::SourceAdaptiveSpinEndpoint endpoint;
    result.solve_=detail::solveSourceAdaptiveSpinEndpoint(request,&endpoint);
    if(result.solve_.status==SourceAdaptiveSpinStatus::BinaryMomentValidated) {
        if(!endpoint.domain) throw std::logic_error("successful source solve missing private terminal domain");
        result.boundary_=makeSourceMaterialBoundary(*endpoint.domain,endpoint.binary,request.density,settings);
    }
    result.outputDigest_=digest({"source-adaptive-spin-material-wrapper-v1",settingsHash,
        result.solve_.inputHash,result.solve_.traceHash,result.solve_.outputHash,
        result.solve_.finalPartitionHash,result.boundary_.facts().inputDigest,
        result.boundary_.facts().productDigest,result.boundary_.facts().status});
    return result;
}
}
