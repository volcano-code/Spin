#include "spin2014/SourceCageShapeFunctional.h"
#include "verification/Sha256.h"
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>

namespace spin2014 {
namespace {
[[noreturn]] void resource(const char* reason) {throw std::length_error(std::string("RESOURCE_LIMIT: source cage shape ")+reason);}
void capacity(std::size_t n) {
    if(n>std::size_t(std::numeric_limits<Eigen::Index>::max()) || n>std::numeric_limits<std::size_t>::max()/sizeof(double))resource("index or byte capacity");
}
std::size_t timesFour(std::size_t n) {if(n>std::numeric_limits<std::size_t>::max()/4)resource("handle product overflow");return n*4;}
class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=static_cast<unsigned char>(x>>(8*i));hash_.update(b,8);}
    void string(const std::string& s) {integer(s.size());hash_.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());}
    void real(double x) {std::uint64_t b;std::memcpy(&b,&x,8);integer(b);}
    template<class T> void matrix(const Eigen::MatrixBase<T>& m) {integer(m.rows());integer(m.cols());for(Eigen::Index i=0;i<m.rows();++i)for(Eigen::Index j=0;j<m.cols();++j)real(m(i,j));}
    std::string finish() {return hash_.finish();}
};
}
void SourceCageShapeFunctional::validateConfig(const SourceCageShapeFunctionalConfig& c) {
    LinearizedRotationLaplacian::validateConfig(c.operatorConfig);
    if(!std::isfinite(c.shapeWeight) || !std::isfinite(c.displacementWeight) || c.shapeWeight<0 || c.displacementWeight<0)
        throw std::invalid_argument("source cage shape weights must be finite nonnegative");
    if(!c.maximumHandleGradientValues || !c.maximumWeightValues)resource("zero budget");
    capacity(c.maximumHandleGradientValues);capacity(c.maximumWeightValues);
}
struct SourceCageShapeFunctional::Data {
    SourceCageBinding binding;
    SourceCageShapeFunctionalConfig config;
    LinearizedRotationLaplacian laplacian;
    std::string configHash;
    Data(const SourceCageBinding& b,const SourceCageShapeFunctionalConfig& c)
        :binding(b),config(c),laplacian(b.source().outerVertices(),b.source().outerFaces(),c.operatorConfig) {
        capacity(timesFour(b.handleCount()));
        if(laplacian.vertexCount()!=b.outerPointCount())throw std::invalid_argument("source cage shape outer count mismatch");
        Digest hash;hash.string("source-cage-shape-config-v1");
        hash.string("actual-source-outer/fixed-linearized-rotation/all-original-handles-centered-scale/full-actual-BBW");
        hash.string(b.configSha256());hash.string(b.outputSha256());hash.string(b.bbw().outputSha256());
        hash.string(laplacian.operatorSha256());hash.real(c.shapeWeight);hash.real(c.displacementWeight);
        hash.integer(c.maximumHandleGradientValues);hash.integer(c.maximumWeightValues);configHash=hash.finish();
    }
};
SourceCageShapeFunctional::SourceCageShapeFunctional(const SourceCageBinding& b,const SourceCageShapeFunctionalConfig& c) {
    validateConfig(c);
    try {data_=std::make_shared<const Data>(b,c);}catch(const std::bad_alloc&) {resource("construction allocation failure");}
}
const SourceCageBinding& SourceCageShapeFunctional::binding() const {return data_->binding;}
const SourceCageShapeFunctionalConfig& SourceCageShapeFunctional::config() const {return data_->config;}
const LinearizedRotationLaplacian& SourceCageShapeFunctional::laplacian() const {return data_->laplacian;}
const std::string& SourceCageShapeFunctional::configSha256() const {return data_->configHash;}
const std::string& SourceCageShapeFunctional::upstreamSha256() const {return binding().outputSha256();}
SourceCageShapeEvaluation SourceCageShapeFunctional::evaluate(Eigen::Ref<const Eigen::MatrixXd> v) const {return run(v,true);}
SourceCageShapeEvaluation SourceCageShapeFunctional::evaluateValues(Eigen::Ref<const Eigen::MatrixXd> v) const {return run(v,false);}
SourceCageShapeEvaluation SourceCageShapeFunctional::run(Eigen::Ref<const Eigen::MatrixXd> v,bool derivatives) const {
    const auto handles=binding().handleCount();
    if(derivatives && (timesFour(handles)>config().maximumHandleGradientValues || handles>config().maximumWeightValues))resource("all-handle gradient or weight-row cap");
    try {
        SourceCageShapeEvaluation r;
        // Move assignment of dynamic Eigen arrays: no second V-sized copy.
        static_cast<LinearizedRotationLaplacianEvaluation&>(r)=derivatives?
            laplacian().evaluate(v,config().shapeWeight,config().displacementWeight):
            laplacian().evaluateValues(v,config().shapeWeight,config().displacementWeight);
        if(derivatives) {
            r.handleGradient=Eigen::VectorXd::Zero(timesFour(handles));Eigen::VectorXd row(handles);
            const auto& field=binding().bbw().finalWeights();const auto& cage=binding().bbw().fem().grid().cage();
            for(std::size_t i=0;i<binding().outerPointCount();++i) {
                row.setZero();const auto& embedding=binding().embeddings()[i];
                for(std::size_t h=0;h<handles;++h)for(int node=0;node<8;++node)
                    row[h]+=embedding.coefficients[node]*field(embedding.nodeIndices[node],h);
                if(!row.allFinite())throw std::overflow_error("source cage shape BBW interpolation overflow");
                const Eigen::Vector3d rest=laplacian().restVertices().row(i),gradient=r.vertexGradient.row(i);
                for(std::size_t h=0;h<handles;++h) {
                    const Eigen::Vector3d center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]);
                    r.handleGradient[4*h]+=row[h]*gradient.dot(rest-center);
                    for(int axis=0;axis<3;++axis)r.handleGradient[4*h+1+axis]+=row[h]*gradient[axis];
                }
            }
            if(!r.handleGradient.allFinite())throw std::overflow_error("source cage shape handle gradient overflow");
        }
        Digest hash;hash.string("source-cage-shape-output-v1");hash.string(configSha256());hash.string(r.outputSha256);
        hash.matrix(r.handleGradient);r.outputSha256=hash.finish();return r;
    } catch(const std::bad_alloc&) {resource("evaluation allocation failure");}
}
}
