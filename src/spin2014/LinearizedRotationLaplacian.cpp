#include "spin2014/LinearizedRotationLaplacian.h"
#include "verification/Sha256.h"
#include <Eigen/SVD>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>

namespace spin2014 {
namespace {
[[noreturn]] void resource(const char* reason) {throw std::length_error(std::string("RESOURCE_LIMIT: linearized rotation ")+reason);}
std::size_t add(std::size_t a,std::size_t b) {if(a>std::numeric_limits<std::size_t>::max()-b)resource("sum overflow");return a+b;}
std::size_t mul(std::size_t a,std::size_t b) {if(b && a>std::numeric_limits<std::size_t>::max()/b)resource("product overflow");return a*b;}
std::size_t doubles(std::size_t bytes) {return add(bytes,sizeof(double)-1)/sizeof(double);}
void capacity(std::size_t n,std::size_t width=sizeof(double)) {
    if(n>std::size_t(std::numeric_limits<Eigen::Index>::max()) || n>std::numeric_limits<std::size_t>::max()/width)resource("index or byte capacity");
}
void workspace(std::size_t required,const LinearizedRotationLaplacianConfig& c) {
    if(required>c.maximumConstructionWorkspaceDoubles)resource("construction workspace cap");
}
class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=static_cast<unsigned char>(x>>(8*i));hash_.update(b,8);}
    void string(const std::string& s) {integer(s.size());hash_.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());}
    void real(double x) {std::uint64_t b;static_assert(sizeof(b)==sizeof(x) && std::numeric_limits<double>::is_iec559,"IEEE754 required");std::memcpy(&b,&x,8);integer(b);}
    template<class T> void matrix(const Eigen::MatrixBase<T>& m) {integer(m.rows());integer(m.cols());for(Eigen::Index i=0;i<m.rows();++i)for(Eigen::Index j=0;j<m.cols();++j)real(m(i,j));}
    std::string finish() {return hash_.finish();}
};
using Block=Eigen::Matrix<double,3,7>;
Block linearBlock(const Eigen::Vector3d& p,bool translation) {
    Block b;b<<p.x(),0,p.z(),-p.y(),0,0,0,p.y(),-p.z(),0,p.x(),0,0,0,p.z(),p.y(),-p.x(),0,0,0,0;
    if(translation)b.rightCols<3>().setIdentity();return b;
}
void weights(double shape,double displacement) {
    if(!std::isfinite(shape) || !std::isfinite(displacement) || shape<0 || displacement<0)
        throw std::invalid_argument("linearized rotation weights must be finite and nonnegative");
}
}
void LinearizedRotationLaplacian::validateConfig(const LinearizedRotationLaplacianConfig& c) {
    for(std::size_t value:{c.maximumVertices,c.maximumOneRingVertices,c.maximumSparseValues,c.maximumConstructionWorkspaceDoubles,c.maximumEvaluationWorkspaceDoubles}) {
        if(!value)resource("zero budget");capacity(value);
    }
    if(!std::isfinite(c.relativeRankTolerance) || c.relativeRankTolerance<=0 || c.relativeRankTolerance>=1)
        throw std::invalid_argument("linearized rotation rank tolerance must lie in (0,1)");
}
struct LinearizedRotationLaplacian::Data {
    Eigen::MatrixXd rest;
    Eigen::MatrixXi faces;
    LinearizedRotationLaplacianConfig config;
    std::vector<std::size_t> offsets,columns;
    std::vector<double> values;
    std::size_t maxRing=0,constructionDoubles=0;
    std::string configHash,operatorHash;
    Data(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f,const LinearizedRotationLaplacianConfig& c):config(c) {
        validateConfig(c);
        if(v.cols()!=3 || v.rows()<3 || f.cols()!=3 || f.rows()<1 || !v.allFinite())
            throw std::invalid_argument("linearized rotation malformed or nonfinite rest mesh");
        const auto count=std::size_t(v.rows()),faceCount=std::size_t(f.rows());
        if(count>c.maximumVertices)resource("vertex cap");
        const auto coordinates=mul(count,3),faceIndices=mul(faceCount,3),edgeCapacity=mul(faceCount,6);
        capacity(coordinates);capacity(faceIndices,sizeof(int));
        using Edge=std::pair<std::size_t,std::size_t>;
        capacity(edgeCapacity,sizeof(Edge));
        if(edgeCapacity>std::vector<Edge>().max_size())resource("edge vector capacity");
        const auto geometryBytes=add(mul(coordinates,sizeof(double)),mul(faceIndices,sizeof(int)));
        const auto edgeBytes=mul(edgeCapacity,sizeof(Edge));
        workspace(doubles(add(geometryBytes,edgeBytes)),c);
        // Reject face indices and degenerate triangles before edge allocation.
        for(Eigen::Index i=0;i<f.rows();++i) {
            for(int a=0;a<3;++a)if(f(i,a)<0 || std::size_t(f(i,a))>=count)
                throw std::invalid_argument("linearized rotation face index out of range");
            if(f(i,0)==f(i,1) || f(i,0)==f(i,2) || f(i,1)==f(i,2))
                throw std::invalid_argument("linearized rotation repeated face vertex");
            const Eigen::Vector3d e=v.row(f(i,1))-v.row(f(i,0)),g=v.row(f(i,2))-v.row(f(i,0));
            if(!e.allFinite() || !g.allFinite())throw std::overflow_error("linearized rotation triangle differences overflow");
            const double s=std::max(e.cwiseAbs().maxCoeff(),g.cwiseAbs().maxCoeff());
            if(s==0 || (e/s).cross(g/s).cwiseAbs().maxCoeff()==0)
                throw std::invalid_argument("linearized rotation zero-area triangle");
        }
        std::vector<Edge> edges;edges.reserve(edgeCapacity);
        for(Eigen::Index i=0;i<f.rows();++i)for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            if(a!=b)edges.emplace_back(std::size_t(f(i,a)),std::size_t(f(i,b)));
        std::sort(edges.begin(),edges.end());edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
        std::size_t cursor=0;
        for(std::size_t i=0;i<count;++i) {
            const auto start=cursor;while(cursor<edges.size() && edges[cursor].first==i)++cursor;
            const auto degree=cursor-start;
            if(!degree)throw std::invalid_argument("linearized rotation isolated vertex");
            if(degree>c.maximumOneRingVertices)resource("one-ring cap");
            maxRing=std::max(maxRing,degree);
        }
        const auto blocks=add(count,edges.size()),sparseValues=mul(blocks,9);
        capacity(blocks,sizeof(std::size_t));capacity(add(count,1),sizeof(std::size_t));capacity(sparseValues);
        if(blocks>columns.max_size() || add(count,1)>offsets.max_size() || sparseValues>values.max_size())resource("sparse vector capacity");
        if(sparseValues>c.maximumSparseValues)resource("sparse contribution cap");
        const auto indexBytes=mul(add(add(count,1),blocks),sizeof(std::size_t));
        const auto retainedBytes=add(geometryBytes,add(indexBytes,mul(sparseValues,sizeof(double))));
        const auto rows=mul(3,add(maxRing,1));capacity(mul(rows,7));
        // A, thin U, QR preconditioner/scaled work and packing conservatively
        // bounded by eight m-by-7 buffers. Never an m-by-m U or a global GEMM.
        const auto localDoubles=add(mul(mul(rows,7),8),2048);
        constructionDoubles=add(doubles(add(retainedBytes,edgeBytes)),localDoubles);
        workspace(constructionDoubles,c);
        rest=v;faces=f;offsets.resize(add(count,1));columns.resize(blocks);values.resize(sparseValues);
        cursor=0;std::size_t output=0;
        for(std::size_t i=0;i<count;++i) {
            const auto start=cursor;while(cursor<edges.size() && edges[cursor].first==i)++cursor;
            const auto degree=cursor-start,local=add(degree,1);offsets[i]=output;
            columns[output]=i;for(std::size_t j=0;j<degree;++j)columns[output+1+j]=edges[start+j].second;
            double scale=0;
            for(std::size_t j=1;j<local;++j) {
                const Eigen::Vector3d difference=rest.row(columns[output+j])-rest.row(i);
                if(!difference.allFinite())throw std::overflow_error("linearized rotation rest difference overflow");
                scale=std::max(scale,difference.stableNorm());
            }
            if(!std::isfinite(scale))throw std::overflow_error("linearized rotation local scale overflow");
            if(scale<=0)throw std::invalid_argument("linearized rotation zero local scale");
            Eigen::MatrixXd a(3*local,7);Eigen::Vector3d delta=Eigen::Vector3d::Zero();
            for(std::size_t j=0;j<local;++j) {
                const Eigen::Vector3d q=(rest.row(columns[output+j])-rest.row(i)).transpose()/scale;
                a.block<3,7>(3*j,0)=linearBlock(q,true);
                if(j)delta-=q/double(degree);
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);
            if(!svd.singularValues().allFinite() || !svd.matrixU().allFinite() || !svd.matrixV().allFinite())
                throw std::overflow_error("linearized rotation nonfinite local SVD");
            if(svd.singularValues().size()!=7 || svd.singularValues()[6]<=c.relativeRankTolerance*svd.singularValues()[0])
                throw std::invalid_argument("linearized rotation rank-deficient or ill-conditioned local fit");
            const Block d=linearBlock(delta,false);Block dv;
            // Explicit seven-direction contractions avoid dynamic packing.
            for(int row=0;row<3;++row)for(int k=0;k<7;++k) {
                double value=0;for(int t=0;t<7;++t)value+=d(row,t)*svd.matrixV()(t,k);
                dv(row,k)=value/svd.singularValues()[k];
            }
            for(std::size_t j=0;j<local;++j)for(int row=0;row<3;++row)for(int col=0;col<3;++col) {
                double value=0;for(int k=0;k<7;++k)value+=dv(row,k)*svd.matrixU()(3*j+col,k);
                if(row==col)value+=j?1./double(degree):-1.;
                if(!std::isfinite(value))throw std::overflow_error("linearized rotation coefficient overflow");
                values[9*(output+j)+3*row+col]=value;
            }
            output+=local;
        }
        offsets[count]=output;
        Digest hash;hash.string("linearized-rotation-laplacian-config-v1");
        hash.string("Sorkine2004:Eq2,5,6,8-12/uniform-one-ring/centered-isotropic-rest/thin-SVD/full-rank-7");
        hash.string("IMPLEMENTATION_DEFAULT:uniform-weights-reconstruction-and-resource-rank-settings");
        hash.integer(c.maximumVertices);hash.integer(c.maximumOneRingVertices);hash.integer(c.maximumSparseValues);
        hash.integer(c.maximumConstructionWorkspaceDoubles);hash.integer(c.maximumEvaluationWorkspaceDoubles);hash.real(c.relativeRankTolerance);
        configHash=hash.finish();Digest op;op.string("linearized-rotation-laplacian-operator-v1/D*Aplus-L/xyz-block-CSR");
        op.string(configHash);op.matrix(rest);op.matrix(faces);op.integer(maxRing);op.integer(constructionDoubles);
        for(auto offset:offsets)op.integer(offset);for(auto column:columns)op.integer(column);for(double value:values)op.real(value);
        operatorHash=op.finish();
    }
};
LinearizedRotationLaplacian::LinearizedRotationLaplacian(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f,const LinearizedRotationLaplacianConfig& c) {
    try {data_=std::make_shared<const Data>(v,f,c);}catch(const std::bad_alloc&) {resource("construction allocation failure");}
}
const Eigen::MatrixXd& LinearizedRotationLaplacian::restVertices() const {return data_->rest;}
const Eigen::MatrixXi& LinearizedRotationLaplacian::faces() const {return data_->faces;}
const LinearizedRotationLaplacianConfig& LinearizedRotationLaplacian::config() const {return data_->config;}
std::size_t LinearizedRotationLaplacian::vertexCount() const {return restVertices().rows();}
std::size_t LinearizedRotationLaplacian::coefficientCount() const {return data_->columns.size();}
std::size_t LinearizedRotationLaplacian::sparseValueCount() const {return data_->values.size();}
std::size_t LinearizedRotationLaplacian::maximumOneRingSize() const {return data_->maxRing;}
std::size_t LinearizedRotationLaplacian::requiredConstructionWorkspaceDoubles() const {return data_->constructionDoubles;}
std::size_t LinearizedRotationLaplacian::requiredEvaluationWorkspaceDoubles(bool derivatives) const {return mul(vertexCount(),derivatives?6:3);}
const std::vector<std::size_t>& LinearizedRotationLaplacian::rowOffsets() const {return data_->offsets;}
const std::vector<std::size_t>& LinearizedRotationLaplacian::columnVertices() const {return data_->columns;}
const std::vector<double>& LinearizedRotationLaplacian::coefficients() const {return data_->values;}
const std::string& LinearizedRotationLaplacian::configSha256() const {return data_->configHash;}
const std::string& LinearizedRotationLaplacian::operatorSha256() const {return data_->operatorHash;}
LinearizedRotationLaplacianEvaluation LinearizedRotationLaplacian::evaluate(Eigen::Ref<const Eigen::MatrixXd> v,double s,double d) const {return run(v,s,d,true);}
LinearizedRotationLaplacianEvaluation LinearizedRotationLaplacian::evaluateValues(Eigen::Ref<const Eigen::MatrixXd> v,double s,double d) const {return run(v,s,d,false);}
LinearizedRotationLaplacianEvaluation LinearizedRotationLaplacian::run(Eigen::Ref<const Eigen::MatrixXd> v,double sw,double dw,bool derivatives) const {
    weights(sw,dw);
    if(v.cols()!=3 || std::size_t(v.rows())!=vertexCount() || !v.allFinite())throw std::invalid_argument("linearized rotation invalid deformed vertices");
    const auto required=requiredEvaluationWorkspaceDoubles(derivatives);capacity(required);
    if(required>config().maximumEvaluationWorkspaceDoubles)resource("evaluation workspace cap");
    try {
        LinearizedRotationLaplacianEvaluation r;r.hasGradient=derivatives;r.residual=Eigen::MatrixXd::Zero(v.rows(),3);
        for(std::size_t i=0;i<vertexCount();++i)for(std::size_t b=data_->offsets[i];b<data_->offsets[i+1];++b)
            for(int row=0;row<3;++row)for(int col=0;col<3;++col)
                r.residual(i,row)+=data_->values[9*b+3*row+col]*v(data_->columns[b],col);
        if(!r.residual.allFinite())throw std::overflow_error("linearized rotation residual overflow");
        if(derivatives)r.vertexGradient=Eigen::MatrixXd::Zero(v.rows(),3);
        for(std::size_t i=0;i<vertexCount();++i)for(int axis=0;axis<3;++axis) {
            const double delta=v(i,axis)-data_->rest(i,axis),residual=r.residual(i,axis);
            if(!std::isfinite(delta))throw std::overflow_error("linearized rotation displacement overflow");
            if(sw)r.shapeEnergy+=sw*residual*residual;
            if(dw)r.displacementEnergy+=dw*delta*delta;
            if(derivatives)r.vertexGradient(i,axis)=2*dw*delta;
        }
        if(derivatives && sw)for(std::size_t i=0;i<vertexCount();++i)for(std::size_t b=data_->offsets[i];b<data_->offsets[i+1];++b)
            for(int row=0;row<3;++row)for(int col=0;col<3;++col)
                r.vertexGradient(data_->columns[b],col)+=2*sw*data_->values[9*b+3*row+col]*r.residual(i,row);
        if(!std::isfinite(r.shapeEnergy) || !std::isfinite(r.displacementEnergy) || !r.vertexGradient.allFinite())
            throw std::overflow_error("linearized rotation energy or gradient overflow");
        Digest hash;hash.string("linearized-rotation-laplacian-output-v1");hash.string(operatorSha256());hash.matrix(v);hash.real(sw);hash.real(dw);
        hash.integer(derivatives);hash.real(r.shapeEnergy);hash.real(r.displacementEnergy);hash.matrix(r.residual);hash.matrix(r.vertexGradient);
        r.outputSha256=hash.finish();return r;
    } catch(const std::bad_alloc&) {resource("evaluation allocation failure");}
}
}
