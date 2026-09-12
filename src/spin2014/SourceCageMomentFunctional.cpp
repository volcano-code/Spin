#include "spin2014/SourceCageMomentFunctional.h"
#include "SurfaceTriangleMoments.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>

namespace spin2014 {
namespace {
constexpr auto absent=SourceConformingBoundary::noOwner;
using Embedding=SourceCageGrid::Embedding;
using Triangle=SourceConformingBoundary::Triangle;
[[noreturn]] void resource(const char* reason) {
    throw std::length_error(std::string("RESOURCE_LIMIT: source cage moments ")+reason);
}
std::size_t product(std::size_t a,std::size_t b) {
    if (b && a>std::numeric_limits<std::size_t>::max()/b) resource("size multiplication overflow");
    return a*b;
}
void doubleCapacity(std::size_t values) {
    if (values>std::size_t(std::numeric_limits<Eigen::Index>::max()) ||
        values>std::numeric_limits<std::size_t>::max()/sizeof(double))
        resource("Eigen index or byte capacity");
}
class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {
        unsigned char bytes[8];
        for (int i=0;i<8;++i) bytes[7-i]=static_cast<unsigned char>(x>>(8*i));
        hash_.update(bytes,8);
    }
    void string(const std::string& x) {
        integer(x.size());
        hash_.update(reinterpret_cast<const unsigned char*>(x.data()),x.size());
    }
    void real(double x) {
        static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,"IEEE754 required");
        std::uint64_t bits;
        std::memcpy(&bits,&x,8);
        integer(bits);
    }
    template<class Derived> void matrix(const Eigen::MatrixBase<Derived>& matrix) {
        integer(matrix.rows()); integer(matrix.cols());
        for (Eigen::Index row=0;row<matrix.rows();++row)
            for (Eigen::Index col=0;col<matrix.cols();++col) real(matrix(row,col));
    }
    std::string finish() {return hash_.finish();}
};
void validateConfig(const SourceConformingBoundary& b,double density,
                    const SourceCageMomentFunctionalConfig& c) {
    if (!std::isfinite(density) || density<=0)
        throw std::invalid_argument("source cage moments requires positive finite density");
    if (!c.maximumWeightValues || !c.maximumJacobianValues) resource("zero budget");
    doubleCapacity(c.maximumWeightValues);
    doubleCapacity(c.maximumJacobianValues);
    const auto handles=b.binding().handleCount();
    doubleCapacity(product(handles,4));
    doubleCapacity(product(std::size_t(b.points().rows()),3));
    doubleCapacity(b.binding().domain().leaves().size());
    if (handles>c.maximumWeightValues || handles>b.binding().config().maximumWeightValues ||
        handles>b.binding().bbw().fem().grid().config().maximumInterpolationValues)
        resource("full handle row exceeds weight budget");
    const auto points=std::size_t(b.points().rows());
    if (points>b.config().maximumPoints || points>std::vector<Embedding>().max_size())
        resource("point or embedding budget");
    product(points,sizeof(Embedding));
}
physics::RawMassMoments10 raw(const detail::TriangleMomentValues& values) {
    MassMoments mass;
    mass.mass=values[0];
    mass.firstMoment=values.segment<3>(1);
    mass.secondMoment<<values[4],values[9],values[8],values[9],values[5],values[7],values[8],values[7],values[6];
    return physics::RawMassMoments10::fromMassMoments(mass);
}
}

struct SourceCageMomentFunctional::Data {
    SourceConformingBoundary boundary;
    SourceCageMomentFunctionalConfig config;
    double density;
    std::size_t betaCount,handleCount;
    std::vector<Embedding> embeddings;
    std::string configHash;

    Data(const SourceConformingBoundary& b,double rho,const SourceCageMomentFunctionalConfig& c)
        :boundary(b),config(c),density(rho),betaCount(b.binding().domain().leaves().size()),
         handleCount(b.binding().handleCount()) {
        const auto& grid=boundary.binding().bbw().fem().grid();
        embeddings.reserve(std::size_t(boundary.points().rows()));
        for (Eigen::Index point=0;point<boundary.points().rows();++point)
            embeddings.push_back(grid.embed(boundary.points().row(point).transpose()));
        Digest hash;
        hash.string("source-cage-moment-config-v1");
        hash.string("supplement-algorithms-1-2/fixed-rest-full-bbw/identity-rotations");
        hash.string("original-handle-order:[scale,tx,ty,tz]/native-raw10/hollowed-beta");
        hash.string(boundary.configSha256());
        hash.real(density);
        hash.integer(config.maximumWeightValues);
        hash.integer(config.maximumJacobianValues);
        hash.integer(betaCount);
        hash.integer(handleCount);
        configHash=hash.finish();
    }
    Eigen::Vector3d center(std::size_t handle) const {
        const auto& cage=boundary.binding().bbw().fem().grid().cage();
        return cage.boundaryVertices().row(cage.handleVertexIndices()[handle]).transpose();
    }
    void weights(std::size_t point,Eigen::VectorXd& row) const {
        row.setZero();
        const auto& embedding=embeddings[point];
        const auto& field=boundary.binding().bbw().finalWeights();
        for (int node=0;node<8;++node)
            row+=embedding.coefficients[node]*field.row(embedding.nodeIndices[node]).transpose();
        if (!row.allFinite()) throw std::overflow_error("source cage moments weight arithmetic");
    }
};

SourceCageMomentFunctional::SourceCageMomentFunctional(const SourceConformingBoundary& b,
    double density,const SourceCageMomentFunctionalConfig& c) {
    validateConfig(b,density,c);
    try {data_=std::make_shared<const Data>(b,density,c);}
    catch (const std::bad_alloc&) {resource("retained allocation failure");}
}
const SourceConformingBoundary& SourceCageMomentFunctional::boundary() const {return data_->boundary;}
const SourceCageMomentFunctionalConfig& SourceCageMomentFunctional::config() const {return data_->config;}
double SourceCageMomentFunctional::density() const {return data_->density;}
std::size_t SourceCageMomentFunctional::betaCount() const {return data_->betaCount;}
std::size_t SourceCageMomentFunctional::handleCount() const {return data_->handleCount;}
const std::string& SourceCageMomentFunctional::configSha256() const {return data_->configHash;}
const std::string& SourceCageMomentFunctional::upstreamSha256() const {return boundary().outputSha256();}
Eigen::VectorXd SourceCageMomentFunctional::initialHandleParameters() const {
    try {
        Eigen::VectorXd result=Eigen::VectorXd::Zero(4*handleCount());
        for (std::size_t h=0;h<handleCount();++h) result[4*h]=1.;
        return result;
    } catch (const std::bad_alloc&) {resource("initial state allocation failure");}
}
SourceCageMomentEvaluation SourceCageMomentFunctional::evaluate(const Eigen::VectorXd& beta,
    const Eigen::VectorXd& parameters) const {return run(beta,parameters,true);}
SourceCageMomentEvaluation SourceCageMomentFunctional::evaluateValues(const Eigen::VectorXd& beta,
    const Eigen::VectorXd& parameters) const {return run(beta,parameters,false);}

SourceCageMomentEvaluation SourceCageMomentFunctional::run(const Eigen::VectorXd& beta,
    const Eigen::VectorXd& parameters,bool derivatives) const {
    if (std::size_t(beta.size())!=betaCount() || !beta.allFinite() ||
        (beta.array()<0).any() || (beta.array()>1).any() ||
        std::size_t(parameters.size())!=4*handleCount() || !parameters.allFinite())
        throw std::invalid_argument("source cage moments invalid state dimensions or values");
    for (std::size_t h=0;h<handleCount();++h)
        if (parameters[4*h]<=0) throw std::invalid_argument("source cage moments scale must be positive");
    if (derivatives) {
        const auto betaValues=product(10,betaCount()), handleValues=product(40,handleCount());
        doubleCapacity(betaValues); doubleCapacity(handleValues);
        if (betaValues>config().maximumJacobianValues ||
            handleValues>config().maximumJacobianValues-betaValues)
            resource("combined output Jacobian budget");
    }
    try {
        SourceCageMomentEvaluation result;
        result.hasJacobians=derivatives;
        if (derivatives) {
            result.betaJacobian=Eigen::MatrixXd::Zero(10,betaCount());
            result.handleJacobian=Eigen::MatrixXd::Zero(10,4*handleCount());
        }
        const auto& rest=boundary().points();
        result.deformedPoints.resize(rest.rows(),3);
        Eigen::VectorXd weights(handleCount());
        for (Eigen::Index point=0;point<rest.rows();++point) {
            data_->weights(std::size_t(point),weights);
            const Eigen::Vector3d p=rest.row(point).transpose();
            Eigen::Vector3d moved=Eigen::Vector3d::Zero();
            for (std::size_t h=0;h<handleCount();++h) {
                const Eigen::Vector3d center=data_->center(h);
                const Eigen::Vector3d local=center+parameters[4*h]*(p-center)+parameters.segment<3>(4*h+1);
                if (!local.allFinite()) throw std::overflow_error("source cage moments local deformation overflow");
                moved+=weights[h]*local;
            }
            if (!moved.allFinite()) throw std::overflow_error("source cage moments deformed point overflow");
            result.deformedPoints.row(point)=moved.transpose();
        }
        detail::TriangleMomentValues values=detail::TriangleMomentValues::Zero();
        const auto accumulate=[&](const Triangle& triangle,double coefficient,std::size_t negativeBeta,std::size_t positiveBeta) {
            detail::TriangleMomentJacobian localJacobian;
            const auto local=detail::surfaceTriangleMoments(
                result.deformedPoints.row(triangle[0]).transpose(),
                result.deformedPoints.row(triangle[1]).transpose(),
                result.deformedPoints.row(triangle[2]).transpose(),derivatives ? &localJacobian : nullptr);
            values+=coefficient*local;
            if (!derivatives) return;
            if (negativeBeta!=absent) result.betaJacobian.col(negativeBeta)-=local;
            if (positiveBeta!=absent) result.betaJacobian.col(positiveBeta)+=local;
            for (int vertex=0;vertex<3;++vertex) {
                const auto point=triangle[vertex];
                data_->weights(point,weights);
                const Eigen::Vector3d p=rest.row(point).transpose();
                const Eigen::Matrix<double,10,3> block=coefficient*localJacobian.block<10,3>(0,3*vertex);
                for (std::size_t h=0;h<handleCount();++h) {
                    result.handleJacobian.col(4*h)+=weights[h]*(block*(p-data_->center(h)));
                    result.handleJacobian.block<10,3>(0,4*h+1)+=weights[h]*block;
                }
            }
        };
        const auto& faces=boundary().binding().source().outerFaces();
        for (Eigen::Index face=0;face<faces.rows();++face)
            accumulate({{std::size_t(faces(face,0)),std::size_t(faces(face,1)),std::size_t(faces(face,2))}},1.,absent,absent);
        for (const auto& patch:boundary().patches()) {
            const auto negative=patch.negativeOwner==absent ? absent : boundary().binding().primitives()[patch.negativeOwner].betaIndex;
            const auto positive=patch.positiveOwner==absent ? absent : boundary().binding().primitives()[patch.positiveOwner].betaIndex;
            const double coefficient=(positive==absent ? 0. : beta[positive])-(negative==absent ? 0. : beta[negative]);
            for (std::size_t i=0;i<patch.triangleCount;++i)
                accumulate(boundary().triangles()[patch.triangleBegin+i],coefficient,negative,positive);
        }
        values*=density();
        if (derivatives) {
            result.betaJacobian*=density();
            result.handleJacobian*=density();
        }
        if (!values.allFinite() || !result.betaJacobian.allFinite() || !result.handleJacobian.allFinite())
            throw std::overflow_error("source cage moments accumulated arithmetic overflow");
        if (values[0]<=0) throw std::invalid_argument("source cage moments final mass must be positive; no geometry legality granted");
        result.finalMoments=raw(values);
        Digest hash;
        hash.string("source-cage-moment-output-v1");
        hash.string(configSha256());
        hash.string(upstreamSha256());
        hash.string(boundary().binding().bbw().outputSha256());
        hash.matrix(beta); hash.matrix(parameters);
        hash.integer(derivatives);
        hash.matrix(values); hash.matrix(result.deformedPoints);
        hash.matrix(result.betaJacobian); hash.matrix(result.handleJacobian);
        result.outputSha256=hash.finish();
        return result;
    } catch (const std::bad_alloc&) {resource("evaluation allocation failure");}
}
}
