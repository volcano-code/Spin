#include "spin2014/SourceCageBinaryMaterial.h"
#include "verification/Sha256.h"
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>

namespace spin2014 {
namespace {
using Outcome=SourceCageBinaryMaterialOutcome;
using Ten=SourceCageBinaryMaterialFacts::Ten;
[[noreturn]] void resource(const char* s) {throw std::length_error(s);}
std::size_t add(std::size_t a,std::size_t b) {if(a>SIZE_MAX-b)resource("binary addition overflow");return a+b;}
std::size_t mul(std::size_t a,std::size_t b) {if(b&&a>SIZE_MAX/b)resource("binary multiplication overflow");return a*b;}
class Hash {
    verification::Sha256 h_;
public:
    void integer(std::uint64_t n) {unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=n>>(8*i);h_.update(b,8);}
    void string(const std::string& s) {integer(s.size());h_.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());}
    void real(double d) {std::uint64_t n;std::memcpy(&n,&d,8);integer(n);}
    template<class T> void matrix(const Eigen::MatrixBase<T>& m) {integer(m.rows());integer(m.cols());for(Eigen::Index i=0;i<m.rows();++i)for(Eigen::Index j=0;j<m.cols();++j)real(m(i,j));}
    void values(const Eigen::VectorXd& x) {integer(x.size());for(Eigen::Index i=0;i<x.size();++i)real(x[i]);}
    std::string finish() {return h_.finish();}
};
// These domains match the adaptive audit's immutable terminal identity chain.
std::string stateHash(const Eigen::VectorXd& x) {Hash h;h.string("source-adaptive-full-state-v1");h.values(x);return h.finish();}
std::string functionalHash(const SourceCagePaperFunctional& f) {
    Hash h;h.string("source-adaptive-functional-v1");
    for(const auto* s:{&f.configSha256(),&f.upstreamSha256(),&f.source().configSha256(),
        &f.source().upstreamSha256(),&f.source().boundary().configSha256(),&f.source().boundary().outputSha256(),
        &f.source().boundary().binding().configSha256(),&f.source().boundary().binding().outputSha256(),
        &f.source().boundary().binding().domain().partitionHash()})h.string(*s);
    return h.finish();
}
std::string configHash(const SourceCageBinaryMaterialSettings& s) {
    Hash h;h.string("actual-deformed-binary-material-v1 fixed-input-moment-gate-1e-10 normalized-constraints inertia-psd-128eps beta-ge-half-hollowed no-seal");
    h.real(s.constraintTolerance);h.real(s.geometry.minimumOuterDomainClearance);
    const auto& g=s.geometry;
    for(auto n:{s.maximumPoints,s.maximumTriangles,s.maximumOwnedBytes,s.maximumExactDiagnosticBytes,
        g.maximumPoints,g.maximumTriangles,g.maximumPrimitives,g.maximumPairs,g.maximumRecords,g.maximumCheckerBytes,
        g.maximumNefComplexity,g.maximumConnectedShells,g.maximumExactStringBytes})h.integer(n);
    return h.finish();
}
void terminal(SourceCageBinaryMaterialFacts& f,Outcome o,const char* s) {
    f.outcome=o;const auto n=std::min(f.reason.size()-1,std::strlen(s));
    std::memcpy(f.reason.data(),s,n);f.reason[n]=0;
}
void hashMesh(Hash& h,const SourceCageBinaryMesh& m) {h.matrix(m.vertices);h.matrix(m.triangles);}
void hashBoundary(Hash& h,const verification::ExactBoundaryFacts& f) {
    for(const auto* s:{&f.algorithm,&f.backend,&f.inputDigest,&f.orientation,&f.intersectionError})h.string(*s);
    for(bool b:{f.boundaryValid,f.closed,f.edgeManifold,f.vertexManifold,f.consistentlyOriented,
        f.componentOrientationsValid,f.intersectionCompleted,f.intersectionFree})h.integer(b);
    h.integer(f.exactPairTests);
    for(const auto* v:{&f.unreferencedVertices,&f.degenerateFaces,&f.nonManifoldVertices,&f.faceComponents,&f.componentVolumeSigns}) {
        h.integer(v->size());for(auto n:*v)h.integer(n);
    }
    for(const auto* v:{&f.boundaryEdges,&f.nonManifoldEdges,&f.inconsistentEdges,&f.intersectingPairs}) {
        h.integer(v->size());for(const auto& pair:*v)for(auto n:pair)h.integer(n);
    }
    h.integer(f.componentSignedVolumesExact.size());for(const auto& s:f.componentSignedVolumesExact)h.string(s);
}
void finish(SourceCageBinaryMaterialFacts& f,bool candidate,bool validated) {
    Hash h;h.string(f.semantics);h.string(f.roundingConvention);h.integer(unsigned(f.outcome));h.string(f.reason.data());
    for(const auto* s:{&f.inputSha256,&f.functionalSha256,&f.partitionSha256,&f.configSha256,&f.outputSha256,
        &f.relaxedStateSha256,&f.binaryStateSha256,&f.binaryFunctionalOutputSha256})h.string(*s);
    for(bool b:{f.endpointAdmitted,f.emptyCavity,f.momentsAvailable,f.constraintsAvailable,f.momentGatePassed,
        f.constraintGatePassed,f.positiveMass,f.centralInertiaPsd,candidate,validated})h.integer(b);
    h.integer(f.selectedPrimitives);h.integer(f.plannedOwnedBytes);h.integer(f.plannedExactDiagnosticBytes);
    h.real(f.rootSide);h.real(f.density);h.real(f.normalizedMomentTolerance);h.real(f.normalizedInertiaTolerance);
    h.matrix(f.stateDelta);h.matrix(f.worldIntegratedMoments);h.matrix(f.spinIntegratedMoments);
    h.matrix(f.binaryFunctionalWorldMoments);h.matrix(f.normalizedMomentResiduals);h.matrix(f.rawConstraints);h.matrix(f.normalizedConstraints);
    h.matrix(f.centralInertia);h.matrix(f.centralInertiaEigenvalues);
    h.integer(bool(f.geometry));if(f.geometry)h.string(f.geometry->facts().scientificSha256);
    h.integer(bool(f.exactBoundary));if(f.exactBoundary)hashBoundary(h,f.exactBoundary->facts());
    f.scientificSha256=h.finish();
}
void plan(const SourceCagePaperFunctional& f,SourceCageBinaryMaterialFacts& out) {
    const auto& s=out.settings;
    for(auto n:{s.maximumPoints,s.maximumTriangles})if(!n||n>std::size_t(INT_MAX/3))resource("invalid binary mesh count cap");
    for(auto n:{s.maximumOwnedBytes,s.maximumExactDiagnosticBytes})if(!n||n>std::size_t(PTRDIFF_MAX))resource("invalid binary byte cap");
    mul(s.maximumTriangles,s.maximumTriangles);mul(s.maximumPoints,3*sizeof(double));
    const auto& b=f.source().boundary();
    const auto points=std::size_t(b.points().rows()),outerPoints=b.outerPointCount();
    const auto triangles=add(b.triangles().size(),std::size_t(b.binding().source().outerFaces().rows()));
    const auto materialPoints=add(points,outerPoints); // Conservative: includes unused internal points.
    if(points>s.maximumPoints||materialPoints>s.maximumPoints||triangles>s.maximumTriangles)
        resource("binary source/emitted mesh count cap before allocation");
    // Bound all mesh storage, coefficient/index scratch, state snapshots, delta,
    // value-only evaluation/decomposition and full one-row BBW scratch.
    auto bytes=add(16384,mul(materialPoints,6*sizeof(double)+2*sizeof(std::size_t)));
    bytes=add(bytes,mul(triangles,6*sizeof(int)+sizeof(std::int64_t)));
    bytes=add(bytes,mul(add(mul(f.stateCount(),6),add(mul(f.handleCount(),2),mul(f.constraintCount(),4))),sizeof(double)));
    bytes=add(bytes,mul(points,3*sizeof(double)));
    if(f.shapeFunctional())bytes=add(bytes,mul(outerPoints,3*sizeof(double)));
    out.plannedOwnedBytes=bytes;
    // ExactBoundaryAudit has no internal budget argument. Reserve a conservative
    // worst case for every retained face pair, edge/vertex diagnostic and every
    // possible component volume string BEFORE calling it. Geometric vector
    // capacity growth is covered by four payloads; native rational scratch is
    // opaque and remains under the caller's external process guard.
    const auto pairs=triangles?mul(triangles/2,triangles-1+triangles%2):0;
    out.plannedExactDiagnosticBytes=mul(4,add(add(mul(materialPoints,16),mul(triangles,128)),
        add(mul(pairs,sizeof(std::array<int,2>)),mul(triangles,8192+sizeof(std::string)))));
    if(bytes>s.maximumOwnedBytes||out.plannedExactDiagnosticBytes>s.maximumExactDiagnosticBytes)
        resource("binary known buffer/exact diagnostic cap before allocation");
}
bool identity(const SourceCageAdaptiveSolverAudit& audit,const SourceCageAdaptiveSolvedEndpoint& endpoint,
    const SourceCageBinaryMaterialFacts& out) {
    const auto& facts=audit.facts();if(facts.passes.empty())return false;
    const auto& last=facts.passes.back();const auto& g=endpoint.solve().facts();
    return last.source&&last.solve&&last.rebind&&last.rebind->checkpoint()&&
        &endpoint.functional()==last.source.get()&&&endpoint.solve()==last.solve.get()&&
        last.oldFunctionalSha256==out.functionalSha256&&last.newFunctionalSha256==out.functionalSha256&&
        last.oldPartitionSha256==out.partitionSha256&&last.newPartitionSha256==out.partitionSha256&&
        last.solvedStateSha256==out.relaxedStateSha256&&last.newStateSha256==out.relaxedStateSha256&&
        last.guardedScientificSha256==g.scientificSha256&&
        g.source.numericalResultAvailable&&g.source.numerical.status==SlqpStatus::Converged&&
        g.source.finalEvaluationAvailable&&g.finalStatus==SourceCageGuardedFinalStatus::ChecksPassed&&g.finalRecheckCount==1;
}
// Signed tetrahedra (0,a,b,c): E(lambda_i)=1/4,
// E(lambda_i lambda_j)=(1+delta_ij)/20. This does not call Algorithm 1,
// the source functional's integrator, or any shared moment transform routine.
Ten integrate(const SourceCageBinaryMesh& m,double density) {
    Eigen::Matrix<long double,10,1> total=Eigen::Matrix<long double,10,1>::Zero();
    for(Eigen::Index i=0;i<m.triangles.rows();++i) {
        const Eigen::Matrix<long double,3,1> a=m.vertices.row(m.triangles(i,0)).transpose().cast<long double>();
        const Eigen::Matrix<long double,3,1> b=m.vertices.row(m.triangles(i,1)).transpose().cast<long double>();
        const Eigen::Matrix<long double,3,1> c=m.vertices.row(m.triangles(i,2)).transpose().cast<long double>();
        const long double v=density*a.dot(b.cross(c))/6.L;
        const Eigen::Matrix<long double,3,1> sum=a+b+c;
        const Eigen::Matrix<long double,3,3> q=v*(sum*sum.transpose()+a*a.transpose()+b*b.transpose()+c*c.transpose())/20.L;
        total[0]+=v;total.template segment<3>(1)+=v*sum/4.L;
        total[4]+=q(0,0);total[5]+=q(1,1);total[6]+=q(2,2);
        total[7]+=q(1,2);total[8]+=q(0,2);total[9]+=q(0,1);
    }
    return total.cast<double>();
}
Eigen::Matrix3d second(const Ten& m) {Eigen::Matrix3d q;q<<m[4],m[9],m[8],m[9],m[5],m[7],m[8],m[7],m[6];return q;}
Ten transform(const Ten& m,const SourceCagePaperFrame& frame) {
    const Eigen::Vector3d l=m.segment<3>(1),o=frame.origin;
    const Eigen::Vector3d first=frame.rotation*(l-m[0]*o);
    const Eigen::Matrix3d q=frame.rotation*(second(m)-o*l.transpose()-l*o.transpose()+m[0]*o*o.transpose())*frame.rotation.transpose();
    Ten result;result<<m[0],first,q(0,0),q(1,1),q(2,2),q(1,2),q(0,2),q(0,1);return result;
}
void numerical(const SourceCagePaperFunctional& f,const Eigen::VectorXd& binary,
    const SourceCageBinaryMesh& mesh,const SourceCagePaperEvaluation& evaluated,SourceCageBinaryMaterialFacts& out) {
    out.worldIntegratedMoments=integrate(mesh,out.density);
    out.binaryFunctionalWorldMoments=evaluated.world.finalMoments.values();
    out.spinIntegratedMoments=transform(out.worldIntegratedMoments,f.frame());
    out.momentsAvailable=true;
    const int degrees[10]={0,1,1,1,2,2,2,2,2,2};
    for(int i=0;i<10;++i) {
        const double scale=out.density*std::pow(out.rootSide,3+degrees[i]);
        if(!std::isfinite(scale)||scale<=0.)throw std::overflow_error("nonfinite/nonpositive fixed input normalization scale");
        out.normalizedMomentResiduals[i]=(out.worldIntegratedMoments[i]-out.binaryFunctionalWorldMoments[i])/scale;
    }
    out.momentGatePassed=out.worldIntegratedMoments.allFinite()&&out.spinIntegratedMoments.allFinite()&&
        out.binaryFunctionalWorldMoments.allFinite()&&out.normalizedMomentResiduals.allFinite()&&
        out.normalizedMomentResiduals.cwiseAbs().maxCoeff()<=out.normalizedMomentTolerance;
    const auto& m=out.spinIntegratedMoments;out.positiveMass=std::isfinite(m[0])&&m[0]>0.;
    if(!out.positiveMass)return;
    const Eigen::Vector3d first=m.segment<3>(1);
    const Eigen::Matrix3d covariance=second(m)-first*first.transpose()/m[0];
    out.centralInertia=covariance.trace()*Eigen::Matrix3d::Identity()-covariance;
    const double firstScale=out.density*std::pow(out.rootSide,4),secondScale=out.density*std::pow(out.rootSide,5);
    if(out.centralInertia.allFinite()) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(.5*(out.centralInertia+out.centralInertia.transpose()));
        if(eigen.info()==Eigen::Success) {
            out.centralInertiaEigenvalues=eigen.eigenvalues();
            out.centralInertiaPsd=out.centralInertiaEigenvalues.allFinite()&&
                out.centralInertiaEigenvalues.minCoeff()/secondScale>=-out.normalizedInertiaTolerance;
        }
    }
    out.rawConstraints.resize(f.constraintCount());out.rawConstraints[0]=m[1];out.rawConstraints[1]=m[2];
    if(f.hasAngle()) {
        Eigen::Index k=2;if(f.config().objective==Spin2014Objective::YoYo)out.rawConstraints[k++]=m[3];
        out.rawConstraints[k++]=m[8];out.rawConstraints[k++]=m[7];
        const double c=std::cos(binary[f.angleOffset()]),s=std::sin(binary[f.angleOffset()]);
        out.rawConstraints[k]=c*s*(m[4]-m[5])+(c*c-s*s)*m[9];
    }
    const Eigen::Index firstCount=f.config().objective==Spin2014Objective::YoYo?3:2;
    out.normalizedConstraints.resize(out.rawConstraints.size());
    for(Eigen::Index i=0;i<out.rawConstraints.size();++i)out.normalizedConstraints[i]=out.rawConstraints[i]/(i<firstCount?firstScale:secondScale);
    out.constraintsAvailable=true;out.constraintGatePassed=out.rawConstraints.allFinite()&&out.normalizedConstraints.allFinite()&&
        out.normalizedConstraints.cwiseAbs().maxCoeff()<=out.settings.constraintTolerance;
}
}
struct SourceCageBinaryMaterial::Data {
    explicit Data(const SourceCagePaperFunctional& f):functional(f) {}
    SourceCagePaperFunctional functional;
    Eigen::VectorXd relaxed,binary;
    SourceCageBinaryMesh outer,carved,material;
};
const SourceCageBinaryMaterial::Data& SourceCageBinaryMaterial::data() const {
    if(!data_)throw std::logic_error("moved-from binary material has no data");return *data_;
}
const SourceCagePaperFunctional& SourceCageBinaryMaterial::functional() const {return data().functional;}
const Eigen::VectorXd& SourceCageBinaryMaterial::relaxedState() const {return data().relaxed;}
const Eigen::VectorXd& SourceCageBinaryMaterial::binaryState() const {return data().binary;}
const SourceCageBinaryMesh& SourceCageBinaryMaterial::outerBoundary() const {return data().outer;}
const SourceCageBinaryMesh& SourceCageBinaryMaterial::carvedBoundary() const {return data().carved;}
const SourceCageBinaryMesh& SourceCageBinaryMaterial::materialBoundary() const {return data().material;}
struct SourceCageBinaryMaterialAudit::Data {
    SourceCageBinaryMaterialFacts facts;
    std::shared_ptr<const SourceCageBinaryMaterial> candidate,validated;
};
const SourceCageBinaryMaterialFacts& SourceCageBinaryMaterialAudit::facts() const {
    if(!data_)throw std::logic_error("moved-from binary material audit has no evidence");return data_->facts;
}
std::shared_ptr<const SourceCageBinaryMaterial> SourceCageBinaryMaterialAudit::candidate() const {(void)facts();return data_->candidate;}
std::shared_ptr<const SourceCageBinaryMaterial> SourceCageBinaryMaterialAudit::validatedMaterial() const {(void)facts();return data_->validated;}
const char* sourceCageBinaryMaterialOutcomeToken(Outcome o) {
    switch(o) {case Outcome::Validated:return "BINARY_MATERIAL_VALIDATED_PRE_SEAL";
    case Outcome::NotAdmitted:return "NOT_ADMITTED";case Outcome::BinaryValidationFailed:return "BINARY_VALIDATION_FAILED";
    case Outcome::ResourceLimit:return "RESOURCE_LIMIT";default:return "UNKNOWN";}
}
SourceCageBinaryMaterialAudit validateSourceCageBinaryMaterial(const SourceCageAdaptiveSolverAudit& audit,
    const SourceCageBinaryMaterialSettings& settings) {
    if(!std::isfinite(settings.constraintTolerance)||settings.constraintTolerance<=0.)
        throw std::invalid_argument("binary constraint tolerance must be finite positive");
    auto data=std::make_shared<SourceCageBinaryMaterialAudit::Data>();auto& out=data->facts;
    out.settings=settings;out.configSha256=configHash(settings);out.inputSha256=audit.facts().scientificSha256;
    auto run=[&] {
        const auto endpoint=audit.lastSolved();
        if(audit.facts().outcome!=SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete||!endpoint||audit.pendingCheckpoint()) {
            terminal(out,Outcome::NotAdmitted,"requires immutable RelaxedAdaptiveComplete endpoint without pending checkpoint");return;
        }
        const auto& f=endpoint->functional();const auto& b=f.source().boundary();
        out.functionalSha256=functionalHash(f);out.relaxedStateSha256=stateHash(endpoint->state());
        out.partitionSha256=b.binding().domain().partitionHash();
        if(!identity(audit,*endpoint,out)) {terminal(out,Outcome::NotAdmitted,"immutable terminal endpoint identity mismatch");return;}
        out.endpointAdmitted=true;out.rootSide=b.binding().source().config().rootSide;out.density=f.source().density();
        const SourceCageGeometryChecker checker(f,settings.geometry);
        plan(f,out);
        auto material=std::make_shared<SourceCageBinaryMaterial::Data>(f);
        material->relaxed=endpoint->state();material->binary=material->relaxed;
        for(std::size_t i=0;i<f.betaCount();++i)material->binary[i]=material->relaxed[i]>=.5?1.:0.;
        out.stateDelta=material->binary-material->relaxed;out.binaryStateSha256=stateHash(material->binary);
        const auto values=f.evaluateValues(material->binary);out.binaryFunctionalOutputSha256=values.outputSha256;
        material->outer.vertices=values.world.deformedPoints.topRows(b.outerPointCount());
        material->outer.triangles=b.binding().source().outerFaces();
        std::vector<std::int64_t> coefficient(b.triangles().size(),0);
        const auto& primitives=b.binding().primitives();
        for(std::size_t i=0;i<primitives.size();++i)if(material->binary[primitives[i].betaIndex]==1.) {
            ++out.selectedPrimitives;
            for(const auto& ref:b.primitiveBoundaries().at(i)) {
                if(ref.orientation!=1&&ref.orientation!=-1)throw std::logic_error("invalid immutable triangle orientation");
                auto& sum=coefficient.at(ref.triangle);
                if((ref.orientation>0&&sum==INT64_MAX)||(ref.orientation<0&&sum==INT64_MIN))resource("oriented coefficient overflow");
                sum+=ref.orientation;
            }
        }
        std::vector<std::size_t> map(std::size_t(b.points().rows()),SIZE_MAX);
        std::size_t faceCount=0,pointCount=0;
        for(std::size_t i=0;i<coefficient.size();++i)if(coefficient[i]) {
            if(coefficient[i]!=1&&coefficient[i]!=-1)throw std::logic_error("nonunit surviving oriented boundary coefficient");
            ++faceCount;for(auto index:b.triangles()[i])map.at(index)=0;
        }
        for(auto& index:map)if(index!=SIZE_MAX)index=pointCount++;
        out.emptyCavity=out.selectedPrimitives==0;
        if((faceCount==0)!=out.emptyCavity)throw std::logic_error("nonempty selected cavity lacks boundary");
        material->carved.vertices.resize(pointCount,3);material->carved.triangles.resize(faceCount,3);
        for(std::size_t i=0;i<map.size();++i)if(map[i]!=SIZE_MAX)material->carved.vertices.row(map[i])=values.world.deformedPoints.row(i);
        std::size_t row=0;
        for(std::size_t i=0;i<coefficient.size();++i)if(coefficient[i]) {
            auto triangle=b.triangles()[i];if(coefficient[i]<0)std::swap(triangle[1],triangle[2]);
            for(int j=0;j<3;++j)material->carved.triangles(row,j)=int(map[triangle[j]]);++row;
        }
        const auto outerPoints=material->outer.vertices.rows(),outerFaces=material->outer.triangles.rows();
        material->material.vertices.resize(outerPoints+Eigen::Index(pointCount),3);
        material->material.vertices.topRows(outerPoints)=material->outer.vertices;
        material->material.vertices.bottomRows(pointCount)=material->carved.vertices;
        material->material.triangles.resize(outerFaces+Eigen::Index(faceCount),3);
        material->material.triangles.topRows(outerFaces)=material->outer.triangles;
        for(std::size_t i=0;i<faceCount;++i)for(int j=0;j<3;++j)
            material->material.triangles(outerFaces+i,j)=int(outerPoints)+material->carved.triangles(i,j==0?0:3-j);
        Hash h;h.string("actual-deformed-binary-emitted-triangles-v1");
        hashMesh(h,material->outer);hashMesh(h,material->carved);hashMesh(h,material->material);out.outputSha256=h.finish();
        data->candidate=std::shared_ptr<const SourceCageBinaryMaterial>(new SourceCageBinaryMaterial(material));
        numerical(f,material->binary,material->material,values,out);
        out.geometry=std::make_shared<const SourceCageGeometryAudit>(checker.check(material->binary));
        // Always audit final emitted MATERIAL, including a legitimate empty
        // cavity's nonempty outer mesh. No fabricated empty-mesh audit exists.
        out.exactBoundary=std::make_shared<const verification::ExactBoundaryAudit>(verification::auditExactTriangleBoundary(
            material->material.vertices,material->material.triangles,verification::ExactBoundaryOrientation::AnyNonzero));
        const bool numericalFailure=!out.positiveMass||!out.momentGatePassed||!out.constraintGatePassed||!out.centralInertiaPsd;
        if(numericalFailure||out.geometry->facts().outcome==SourceCageGeometryOutcome::RejectedGeometry||
            (!out.exactBoundary->facts().boundaryValid&&out.exactBoundary->facts().intersectionCompleted))
            terminal(out,Outcome::BinaryValidationFailed,"binary numerical or actual geometry/final boundary gate failed; raw evidence retained");
        else if(out.geometry->facts().outcome==SourceCageGeometryOutcome::ResourceLimit)
            terminal(out,Outcome::ResourceLimit,"actual binary geometry checker resource limit");
        else if(out.geometry->facts().outcome!=SourceCageGeometryOutcome::ChecksPassed||!out.exactBoundary->facts().boundaryValid)
            terminal(out,Outcome::Unknown,"actual binary geometry/final boundary audit incomplete");
        else {
            terminal(out,Outcome::Validated,"bounded binary material checks passed; no final exact region partition, thickness, feasibility, spin, physics, seal or paper grant");
            data->validated=data->candidate;
        }
    };
    auto failure=[&](Outcome fallback,const char* reason) {
        // A later exact checker/allocation exception cannot erase an already
        // established numerical failure or turn it into an unknown result.
        const bool established=(out.momentsAvailable&&(!out.momentGatePassed||!out.positiveMass))||
            (out.constraintsAvailable&&(!out.constraintGatePassed||!out.centralInertiaPsd));
        terminal(out,established?Outcome::BinaryValidationFailed:fallback,reason);
    };
    try {run();}
    catch(const std::bad_alloc&) {failure(Outcome::ResourceLimit,"binary allocation failure; available evidence retained");}
    catch(const std::length_error& e) {failure(Outcome::ResourceLimit,e.what());}
    catch(const std::exception& e) {failure(Outcome::Unknown,e.what());}
    finish(out,bool(data->candidate),bool(data->validated));return SourceCageBinaryMaterialAudit(std::move(data));
}
}
