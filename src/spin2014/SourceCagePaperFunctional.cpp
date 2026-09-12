#include "spin2014/SourceCagePaperFunctional.h"
#include "verification/Sha256.h"
#include <Eigen/LU>
#include <cmath>
#include <cstring>
#include <new>

namespace spin2014 {
namespace {
using Ten=physics::RawMassMoments10::Vector;
[[noreturn]] void resource(const char* reason) {
    throw std::length_error(std::string("RESOURCE_LIMIT: source cage paper ")+reason);
}
std::size_t product(std::size_t a,std::size_t b) {
    if(b && a>std::numeric_limits<std::size_t>::max()/b)resource("size product overflow");
    return a*b;
}
std::size_t sum(std::size_t a,std::size_t b) {
    if(a>std::numeric_limits<std::size_t>::max()-b)resource("size sum overflow");
    return a+b;
}
void capacity(std::size_t n) {
    if(n>std::size_t(std::numeric_limits<Eigen::Index>::max()) ||
       n>std::numeric_limits<std::size_t>::max()/sizeof(double))resource("Eigen index or byte capacity");
}
class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {
        unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=static_cast<unsigned char>(x>>(8*i));hash_.update(b,8);
    }
    void string(const std::string& x) {integer(x.size());hash_.update(reinterpret_cast<const unsigned char*>(x.data()),x.size());}
    void real(double x) {
        static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,"IEEE754 required");
        std::uint64_t b;std::memcpy(&b,&x,8);integer(b);
    }
    template<class Derived> void matrix(const Eigen::MatrixBase<Derived>& m) {
        integer(m.rows());integer(m.cols());
        for(Eigen::Index i=0;i<m.rows();++i)for(Eigen::Index j=0;j<m.cols();++j)real(m(i,j));
    }
    std::string finish() {return hash_.finish();}
};
Eigen::Matrix3d tensor(const Ten& v) {
    Eigen::Matrix3d q;q<<v[4],v[9],v[8],v[9],v[5],v[7],v[8],v[7],v[6];return q;
}
Ten packed(double mass,const Eigen::Vector3d& l,const Eigen::Matrix3d& q) {
    Ten v;v<<mass,l,q(0,0),q(1,1),q(2,2),q(1,2),q(0,2),q(0,1);return v;
}
physics::RawMassMoments10 raw(const Ten& v) {
    MassMoments m;m.mass=v[0];m.firstMoment=v.segment<3>(1);m.secondMoment=tensor(v);
    return physics::RawMassMoments10::fromMassMoments(m);
}
SourceCagePaperFunctional::MomentMap frameMap(const SourceCagePaperFrame& f) {
    const auto& r=f.rotation;const auto& o=f.origin;
    if(!r.allFinite() || !o.allFinite())throw std::invalid_argument("source cage paper frame must be finite");
    const Eigen::Matrix3d residual=r*r.transpose()-Eigen::Matrix3d::Identity();
    const double determinant=r.determinant();
    if(!residual.allFinite() || !std::isfinite(determinant) || residual.cwiseAbs().maxCoeff()>1e-10 ||
       std::abs(determinant-1)>1e-10)throw std::invalid_argument("source cage paper frame must be proper orthonormal");
    SourceCagePaperFunctional::MomentMap map;
    // Apply the affine raw-moment identity to each basis vector. Mixed tensor
    // entries occur twice in Q, but are stored once in the native ten-vector.
    for(int col=0;col<10;++col) {
        const Ten basis=Ten::Unit(col);const Eigen::Vector3d l=basis.segment<3>(1);
        const Eigen::Vector3d first=r*(l-basis[0]*o);
        const Eigen::Matrix3d second=r*(tensor(basis)-o*l.transpose()-l*o.transpose()+basis[0]*(o*o.transpose()))*r.transpose();
        map.col(col)=packed(basis[0],first,second);
    }
    if(!map.allFinite())throw std::overflow_error("source cage paper frame moment map overflow");
    return map;
}
void validateBetaPhi(const Eigen::VectorXd& beta,double phi,std::size_t b,bool angle) {
    if(std::size_t(beta.size())!=b || !beta.allFinite() || (beta.array()<0).any() || (beta.array()>1).any() ||
       !std::isfinite(phi) || (!angle && phi!=0.))throw std::invalid_argument("source cage paper invalid beta or phi");
}
}

struct SourceCagePaperFunctional::Data {
    SourceCageMomentFunctional source;
    SourceCagePaperFunctionalConfig config;
    SourceCagePaperFrame frame;
    Spin2014GivensProblem paper;
    MomentMap map;
    std::size_t handlesOffset,n,k,derivativeValues;
    std::string configHash;
    Data(const SourceCageMomentFunctional& s,const SourceCagePaperFunctionalConfig& c,const SourceCagePaperFrame& f)
        :source(s),config(c),frame(f),paper(c.objective,c.inertiaWeight,c.centerWeight),map(frameMap(f)) {
        if(!std::isfinite(c.fillWeight) || c.fillWeight<0)throw std::invalid_argument("source cage paper invalid fill weight");
        if(!c.maximumJointDerivativeValues)resource("zero derivative budget");
        capacity(c.maximumJointDerivativeValues);
        handlesOffset=sum(source.betaCount(),paper.hasAngle()?1:0);
        n=sum(handlesOffset,product(4,source.handleCount()));capacity(n);capacity(source.betaCount());
        k=paper.hasAngle()?(c.objective==Spin2014Objective::YoYo?6:5):2;
        derivativeValues=sum(product(11+k,n),source.betaCount());capacity(derivativeValues);
        Digest hash;hash.string("source-cage-paper-config-v1");
        hash.string("PAPER_GIVENS/fixed-frame/IMPLEMENTATION_DEFAULT:source-unit-adjacency-fill");
        hash.string(source.configSha256());hash.string(source.upstreamSha256());
        hash.string(spin2014ObjectiveToken(c.objective));
        hash.real(c.inertiaWeight);hash.real(c.centerWeight);hash.real(c.fillWeight);
        hash.integer(c.maximumJointDerivativeValues);hash.matrix(f.origin);hash.matrix(f.rotation);
        hash.integer(n);hash.integer(k);configHash=hash.finish();
    }
};
SourceCagePaperFunctional::SourceCagePaperFunctional(const SourceCageMomentFunctional& source,
    const SourceCagePaperFunctionalConfig& config,const SourceCagePaperFrame& frame) {
    try {data_=std::make_shared<const Data>(source,config,frame);}
    catch(const std::bad_alloc&) {resource("retained allocation failure");}
}
const SourceCageMomentFunctional& SourceCagePaperFunctional::source() const {return data_->source;}
const SourceCagePaperFunctionalConfig& SourceCagePaperFunctional::config() const {return data_->config;}
const SourceCagePaperFrame& SourceCagePaperFunctional::frame() const {return data_->frame;}
const SourceCagePaperFunctional::MomentMap& SourceCagePaperFunctional::worldToSpinMomentMap() const {return data_->map;}
std::size_t SourceCagePaperFunctional::betaCount() const {return source().betaCount();}
std::size_t SourceCagePaperFunctional::handleCount() const {return source().handleCount();}
bool SourceCagePaperFunctional::hasAngle() const {return data_->paper.hasAngle();}
std::size_t SourceCagePaperFunctional::angleOffset() const {return hasAngle()?betaCount():std::numeric_limits<std::size_t>::max();}
std::size_t SourceCagePaperFunctional::handleOffset() const {return data_->handlesOffset;}
std::size_t SourceCagePaperFunctional::stateCount() const {return data_->n;}
std::size_t SourceCagePaperFunctional::constraintCount() const {return data_->k;}
const std::string& SourceCagePaperFunctional::configSha256() const {return data_->configHash;}
const std::string& SourceCagePaperFunctional::upstreamSha256() const {return source().configSha256();}
Eigen::VectorXd SourceCagePaperFunctional::initialState(const Eigen::VectorXd& beta,double phi) const {
    validateBetaPhi(beta,phi,betaCount(),hasAngle());
    try {
        Eigen::VectorXd x(stateCount());x.head(betaCount())=beta;
        if(hasAngle())x[angleOffset()]=phi;
        x.tail(4*handleCount())=source().initialHandleParameters();return x;
    } catch(const std::bad_alloc&) {resource("initial state allocation failure");}
}
SourceCagePaperEvaluation SourceCagePaperFunctional::evaluate(const Eigen::VectorXd& x) const {return runState(x,true);}
SourceCagePaperEvaluation SourceCagePaperFunctional::evaluateValues(const Eigen::VectorXd& x) const {return runState(x,false);}
SourceCagePaperEvaluation SourceCagePaperFunctional::evaluate(const Eigen::VectorXd& b,double phi,const Eigen::VectorXd& p) const {return run(b,phi,p,true);}
SourceCagePaperEvaluation SourceCagePaperFunctional::evaluateValues(const Eigen::VectorXd& b,double phi,const Eigen::VectorXd& p) const {return run(b,phi,p,false);}
SourceCagePaperEvaluation SourceCagePaperFunctional::runState(const Eigen::VectorXd& x,bool derivatives) const {
    if(std::size_t(x.size())!=stateCount() || !x.allFinite())throw std::invalid_argument("source cage paper invalid joint state");
    try {return run(x.head(betaCount()),hasAngle()?x[angleOffset()]:0.,x.tail(4*handleCount()),derivatives);}
    catch(const std::bad_alloc&) {resource("split state allocation failure");}
}
SourceCagePaperEvaluation SourceCagePaperFunctional::run(const Eigen::VectorXd& beta,double phi,
    const Eigen::VectorXd& parameters,bool derivatives) const {
    validateBetaPhi(beta,phi,betaCount(),hasAngle());
    if(std::size_t(parameters.size())!=4*handleCount() || !parameters.allFinite())
        throw std::invalid_argument("source cage paper invalid handle parameters");
    for(std::size_t h=0;h<handleCount();++h)if(parameters[4*h]<=0)
        throw std::invalid_argument("source cage paper handle scale must be positive");
    if(derivatives && data_->derivativeValues>config().maximumJointDerivativeValues)
        resource("combined joint outputs and fill scratch budget");
    try {
        SourceCagePaperEvaluation result;result.hasJacobians=derivatives;
        result.world=derivatives?source().evaluate(beta,parameters):source().evaluateValues(beta,parameters);
        const Ten moments=data_->map*result.world.finalMoments.values();
        if(!moments.allFinite())throw std::overflow_error("source cage paper spin moments overflow");
        result.spinMoments=raw(moments);
        Spin2014GivensProblem::Vector11 paperState;paperState<<moments,phi;
        const auto paper=data_->paper.evaluate(paperState);
        result.paperObjective=paper.objective;result.constraints=paper.constraints;
        const Eigen::VectorXd laplacianBeta=source().boundary().binding().domain().uniformLaplacian()*beta;
        result.fillPenalty=.5*config().fillWeight*beta.dot(laplacianBeta);
        result.totalObjective=result.paperObjective+result.fillPenalty;
        if(!laplacianBeta.allFinite() || !std::isfinite(result.fillPenalty) || !std::isfinite(result.totalObjective))
            throw std::overflow_error("source cage paper fill or objective overflow");
        if(derivatives) {
            result.spinMomentJacobian=Eigen::MatrixXd::Zero(10,stateCount());
            result.spinMomentJacobian.leftCols(betaCount()).noalias()=data_->map*result.world.betaJacobian;
            result.spinMomentJacobian.rightCols(4*handleCount()).noalias()=data_->map*result.world.handleJacobian;
            result.jointGradient.resize(stateCount());result.jointConstraintJacobian.resize(constraintCount(),stateCount());
            // Column-wise contraction keeps the chain-rule workspace fixed; no
            // 11 x N temporary or second copy of a global joint matrix.
            for(std::size_t col=0;col<stateCount();++col) {
                result.jointGradient[col]=paper.objectiveGradient.head<10>().dot(result.spinMomentJacobian.col(col));
                for(std::size_t row=0;row<constraintCount();++row)
                    result.jointConstraintJacobian(row,col)=paper.constraintJacobian.row(row).head<10>().dot(result.spinMomentJacobian.col(col));
            }
            for(std::size_t i=0;i<betaCount();++i)result.jointGradient[i]+=config().fillWeight*laplacianBeta[i];
            if(hasAngle()) {
                result.jointGradient[angleOffset()]+=paper.objectiveGradient[10];
                result.jointConstraintJacobian.col(angleOffset())+=paper.constraintJacobian.col(10);
            }
            if(!result.spinMomentJacobian.allFinite() || !result.jointGradient.allFinite() || !result.jointConstraintJacobian.allFinite())
                throw std::overflow_error("source cage paper joint derivative overflow");
        }
        Digest hash;hash.string("source-cage-paper-output-v1");hash.string(configSha256());
        hash.string(result.world.outputSha256);hash.matrix(beta);hash.real(phi);hash.matrix(parameters);hash.integer(derivatives);
        hash.matrix(moments);hash.real(result.paperObjective);hash.real(result.fillPenalty);hash.real(result.totalObjective);
        hash.matrix(result.constraints);hash.matrix(result.jointGradient);hash.matrix(result.jointConstraintJacobian);hash.matrix(result.spinMomentJacobian);
        result.outputSha256=hash.finish();return result;
    } catch(const std::bad_alloc&) {resource("evaluation allocation failure");}
}
}
