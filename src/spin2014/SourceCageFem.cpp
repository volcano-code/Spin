#include "spin2014/SourceCageFem.h"
#include "verification/Sha256.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>

namespace spin2014 {
namespace {
using Tet=SourceCageFem::Tetrahedron;
using Entry=SourceCageFem::BoundaryEntry;
using Row=SourceCageFem::BoundaryRow;
using Key=SourceCageGrid::Key;
using Pair=std::pair<int,int>;
using Assembly=std::map<Pair,long double>;
using Rows=std::map<int,Row>;
using Vec=Eigen::Matrix<long double,3,1>;
constexpr int split[6][4]={{0,1,3,7},{0,5,1,7},{0,3,2,7},
                          {0,2,6,7},{0,4,5,7},{0,6,4,7}};
class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {
        unsigned char b[8]; for(int i=0;i<8;++i)b[7-i]=static_cast<unsigned char>(x>>(8*i)); hash_.update(b,8);
    }
    void string(const std::string& s) {integer(s.size()); hash_.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());}
    void real(double x) {
        static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,"IEEE754 binary64 required");
        std::uint64_t b; std::memcpy(&b,&x,8); integer(b);
    }
    std::string finish() {return hash_.finish();}
};
double represented(long double x,const char* context) {
    const double d=static_cast<double>(x);
    if(!std::isfinite(x) || !std::isfinite(d) || (x!=0 && d==0))
        throw std::overflow_error(context);
    return d;
}
void validate(const SourceCageGrid& g,const SourceCageFemConfig& c) {
    if(c.massRule!=SourceCageFemMassRule::RESTRICTED_GRID_VORONOI && c.massRule!=SourceCageFemMassRule::BARYCENTRIC_TET)
        throw std::invalid_argument("source cage FEM invalid mass rule");
    const auto imax=std::size_t(std::numeric_limits<Eigen::Index>::max());
    const auto bytes=std::numeric_limits<std::size_t>::max();
    if(!c.maximumTetrahedra || c.maximumTetrahedra>std::vector<Tet>().max_size() || c.maximumTetrahedra>std::vector<double>().max_size() ||
       !c.maximumAssemblyEntries || c.maximumAssemblyEntries>bytes/sizeof(long double) ||
       !c.maximumStiffnessEntries || c.maximumStiffnessEntries>Assembly().max_size() || c.maximumStiffnessEntries>std::size_t(std::numeric_limits<int>::max()) ||
       !c.maximumBoundaryEntries || c.maximumBoundaryEntries>Rows().max_size() || c.maximumBoundaryEntries>std::vector<Row>().max_size() || c.maximumBoundaryEntries>std::vector<Entry>().max_size() ||
       !c.maximumEvaluationValues || c.maximumEvaluationValues>imax || c.maximumEvaluationValues>bytes/sizeof(double))
        throw std::length_error("RESOURCE_LIMIT: invalid FEM tetrahedron/assembly/sparse/boundary/evaluation budget");
    const auto n=g.nodeKeys().size();
    if(n>=std::size_t(std::numeric_limits<int>::max()) || n>std::vector<long double>().max_size() || n>imax || n>bytes/sizeof(double) || n+1>bytes/sizeof(int))
        throw std::length_error("RESOURCE_LIMIT: FEM nodal storage/index overflow");
    if(g.cells().size()>c.maximumTetrahedra/6 || g.cells().size()>c.maximumAssemblyEntries/96)
        throw std::length_error("RESOURCE_LIMIT: FEM tetrahedron or assembly budget exceeded before allocation");
}
int findNode(const SourceCageGrid& g,const Key& key) {
    const auto it=std::lower_bound(g.nodeKeys().begin(),g.nodeKeys().end(),key);
    if(it==g.nodeKeys().end() || *it!=key) throw std::invalid_argument("source cage FEM missing actual boundary node");
    return static_cast<int>(it-g.nodeKeys().begin());
}
bool consistent(const std::vector<Entry>& a,const std::array<Entry,3>& b,int count) {
    std::size_t i=0; int j=0;
    while(i<a.size() || j<count) {
        double difference;
        if(j==count || (i<a.size() && a[i].handleIndex<b[j].handleIndex)) difference=a[i++].weight;
        else if(i==a.size() || b[j].handleIndex<a[i].handleIndex) difference=b[j++].weight;
        else {difference=a[i++].weight-b[j++].weight;}
        if(std::abs(difference)>SourceCageFem::boundaryValidationTolerance)return false;
    }
    return true;
}
// Only integer face-box nodes are visited. Barycentric side tests then use
// actual represented positions, including rounded unequal fine intervals.
std::vector<Row> boundaryRows(const SourceCageGrid& g,const SourceCageFemConfig& config) {
    const auto& cage=g.cage(); const auto& faces=cage.boundaryFaces();
    const auto factor=std::uint32_t(1)<<(g.config().fineLevel-cage.config().cageLevel);
    std::vector<int> handles(cage.boundaryVertexKeys().size(),-1);
    for(std::size_t h=0;h<cage.handleVertexIndices().size();++h) handles[cage.handleVertexIndices()[h]]=static_cast<int>(h);
    Rows rows; std::size_t entries=0;
    for(Eigen::Index t=0;t<faces.rows();++t) {
        Key lo=cage.boundaryVertexKeys()[faces(t,0)],hi=lo;
        std::array<Vec,3> p;
        for(int j=0;j<3;++j) {
            p[j]=cage.boundaryVertices().row(faces(t,j)).transpose().cast<long double>();
            const auto& key=cage.boundaryVertexKeys()[faces(t,j)];
            for(int a=0;a<3;++a) {lo[a]=std::min(lo[a],key[a]); hi[a]=std::max(hi[a],key[a]);}
            if(handles[faces(t,j)]<0) throw std::invalid_argument("source cage FEM triangle corner is not a handle");
        }
        int axis=-1; for(int a=0;a<3;++a)if(lo[a]==hi[a]) {if(axis!=-1)throw std::invalid_argument("source cage FEM degenerate triangle"); axis=a;}
        if(axis<0) throw std::invalid_argument("source cage FEM triangle is not an axis-aligned cage face");
        const int u=(axis+1)%3,v=(axis+2)%3;
        const Vec a=p[1]-p[0],b=p[2]-p[0];
        const long double det=a[u]*b[v]-a[v]*b[u];
        if(!std::isfinite(det) || det==0)throw std::overflow_error("source cage FEM degenerate or unsupported triangle area");
        const long double extent=std::max(a.cwiseAbs().maxCoeff(),b.cwiseAbs().maxCoeff());
        for(int d=0;d<3;++d) {lo[d]*=factor; hi[d]*=factor;}
        for(std::uint32_t i=lo[u];i<=hi[u];++i) for(std::uint32_t j=lo[v];j<=hi[v];++j) {
            Key key=lo; key[u]=i; key[v]=j; const int n=findNode(g,key);
            const Vec point=g.nodePositions().row(n).transpose().cast<long double>();
            const std::array<Vec,3> offset{p[0]-point,p[1]-point,p[2]-point};
            std::array<long double,3> w;
            // Evaluate all three opposite-edge areas directly. Computing the
            // first weight as 1-w1-w2 can turn an exact axis-aligned edge zero
            // negative through cancellation, excluding it from both faces.
            // Shared edges use reversed identical operands; no tolerance,
            // clamping, normalization or change to represented points occurs.
            for(int k=0;k<3;++k) {
                const Vec& first=offset[(k+1)%3];
                const Vec& second=offset[(k+2)%3];
                w[k]=(first[u]*second[v]-first[v]*second[u])/det;
                if(!std::isfinite(w[k]))
                    throw std::overflow_error("source cage FEM nonfinite barycentric coordinates");
            }
            if(w[0]<0 || w[1]<0 || w[2]<0)continue;
            std::array<Entry,3> row{}; int count=0; long double sum=0; Vec residual=Vec::Zero();
            for(int k=0;k<3;++k) if(w[k]!=0) {
                const double weight=represented(w[k],"source cage FEM boundary weight unrepresentable");
                row[count++]={handles[faces(t,k)],weight}; sum+=weight; residual+=static_cast<long double>(weight)*(p[k]-point);
            }
            if(std::abs(sum-1)>SourceCageFem::boundaryValidationTolerance || residual.cwiseAbs().maxCoeff()>SourceCageFem::boundaryValidationTolerance*extent)
                throw std::overflow_error("source cage FEM boundary partition or reconstruction validation failed");
            std::sort(row.begin(),row.begin()+count,[](const Entry& x,const Entry& y){return x.handleIndex<y.handleIndex;});
            const auto old=rows.find(n);
            if(old!=rows.end()) {
                if(!consistent(old->second.entries,row,count))throw std::overflow_error("source cage FEM conflicting duplicate boundary rows");
                continue;
            }
            if(std::size_t(count)>config.maximumBoundaryEntries-entries)
                throw std::length_error("RESOURCE_LIMIT: FEM boundary entry budget exceeded before row allocation");
            entries+=count; rows.emplace(n,Row{n,std::vector<Entry>(row.begin(),row.begin()+count)});
        }
    }
    // Independent completeness pass over every exposed actual fine-cell face.
    // No global lattice or node-by-triangle search, and no interior row added.
    const auto side=std::uint32_t(1)<<g.config().fineLevel;
    for(std::size_t c=0;c<g.cells().size();++c) for(int a=0;a<3;++a) for(int upper=0;upper<2;++upper) {
        Key adjacent=g.cells()[c]; bool inside=upper?adjacent[a]+1<side:adjacent[a]>0;
        if(inside) {if(upper)++adjacent[a]; else --adjacent[a]; inside=std::binary_search(g.cells().begin(),g.cells().end(),adjacent);}
        if(inside)continue;
        for(int b=0;b<8;++b)if(((b>>a)&1)==upper && rows.find(g.cellNodes()[c][b])==rows.end())
            throw std::overflow_error("source cage FEM incomplete boundary constraints");
    }
    std::vector<Row> output; output.reserve(rows.size());
    for(auto& r:rows)output.push_back(std::move(r.second));
    return output;
}
} // namespace

struct SourceCageFem::Data {
    SourceCageGrid grid;
    SourceCageFemConfig config;
    std::vector<Tet> tets;
    std::vector<double> volumes;
    SparseMatrix k;
    Eigen::VectorXd mass;
    std::vector<Row> boundary;
    std::string configHash,outputHash;
    Data(const SourceCageGrid& g,const SourceCageFemConfig& c):grid(g),config(c) {
        validate(grid,config);
        tets.reserve(6*g.cells().size()); volumes.reserve(6*g.cells().size());
        std::vector<long double> m(g.nodeKeys().size(),0); Assembly coefficients;
        for(const auto& cell:g.cellNodes()) {
            if(config.massRule==SourceCageFemMassRule::RESTRICTED_GRID_VORONOI) {
                const Vec d=g.nodePositions().row(cell[7]).transpose().cast<long double>()-g.nodePositions().row(cell[0]).transpose().cast<long double>();
                const long double dual=d.prod()/8;
                if(!std::isfinite(dual) || dual<=0)throw std::overflow_error("source cage FEM unsupported restricted Voronoi volume");
                represented(dual,"source cage FEM restricted Voronoi contribution unrepresentable");
                for(int n:cell)m[n]+=dual;
            }
            for(const auto& local:split) {
                Tet tet; for(int a=0;a<4;++a)tet[a]=cell[local[a]];
                const Vec origin=g.nodePositions().row(tet[0]).transpose().cast<long double>();
                const Vec a=g.nodePositions().row(tet[1]).transpose().cast<long double>()-origin;
                const Vec b=g.nodePositions().row(tet[2]).transpose().cast<long double>()-origin;
                const Vec c=g.nodePositions().row(tet[3]).transpose().cast<long double>()-origin;
                const long double det=a.dot(b.cross(c)),volume=det/6;
                if(!std::isfinite(volume) || volume<=0)throw std::overflow_error("source cage FEM nonpositive or unsupported tetrahedral volume");
                const double v=represented(volume,"source cage FEM tetrahedral volume unrepresentable");
                std::array<Vec,4> grad; grad[1]=b.cross(c)/det; grad[2]=c.cross(a)/det; grad[3]=a.cross(b)/det; grad[0]=-grad[1]-grad[2]-grad[3];
                for(const auto& gradient:grad)if(!gradient.allFinite())throw std::overflow_error("source cage FEM gradient overflow");
                tets.push_back(tet); volumes.push_back(v);
                for(int i=0;i<4;++i) {
                    if(config.massRule==SourceCageFemMassRule::BARYCENTRIC_TET) {
                        represented(volume/4,"source cage FEM barycentric contribution unrepresentable"); m[tet[i]]+=volume/4;
                    }
                    for(int j=0;j<4;++j) {
                        const long double value=volume*grad[i].dot(grad[j]);
                        represented(value,"source cage FEM stiffness contribution unrepresentable");
                        if(value==0)continue;
                        const Pair key{tet[i],tet[j]}; auto found=coefficients.find(key);
                        if(found==coefficients.end()) {
                            if(coefficients.size()>=config.maximumStiffnessEntries)throw std::length_error("RESOURCE_LIMIT: FEM unique stiffness entry budget exceeded before allocation");
                            found=coefficients.emplace(key,0).first;
                        }
                        found->second+=value;
                        if(!std::isfinite(found->second))throw std::overflow_error("source cage FEM stiffness accumulation overflow");
                    }
                }
            }
        }
        mass.resize(static_cast<Eigen::Index>(m.size()));
        for(std::size_t i=0;i<m.size();++i) {
            mass(i)=represented(m[i],"source cage FEM mass unrepresentable");
            if(mass(i)<=0)throw std::overflow_error("source cage FEM nonpositive nodal mass");
        }
        k.resize(mass.size(),mass.size()); k.reserve(static_cast<int>(coefficients.size()));
        for(const auto& item:coefficients)if(item.second!=0) k.insert(item.first.first,item.first.second)=represented(item.second,"source cage FEM assembled stiffness unrepresentable");
        k.makeCompressed();
        // Complete nodal cage-face constraints are assembled independently of K.
        boundary=boundaryRows(grid,config);
        Digest cfg; cfg.string("source-cage-fem-Freudenthal-0-7-positive-IMPLEMENTATION_DEFAULT-v1");
        cfg.string("RESTRICTED_GRID_VORONOI-Vbrick/8-or-BARYCENTRIC_TET-Vtet/4-Q=KtMinvK-E=half-sum-xQx-v1");
        cfg.string("NODAL-boundary-three-direct-opposite-edge-areas-v2");
        cfg.string(grid.configSha256()); cfg.integer(static_cast<unsigned>(config.massRule));
        cfg.integer(config.maximumTetrahedra); cfg.integer(config.maximumAssemblyEntries); cfg.integer(config.maximumStiffnessEntries);
        cfg.integer(config.maximumBoundaryEntries); cfg.integer(config.maximumEvaluationValues); cfg.real(boundaryValidationTolerance); configHash=cfg.finish();
        Digest out; out.string("source-cage-fem-actual-tets-volumes-row-major-K-M-sorted-sparse-NODAL-boundary-v1");
        out.string(grid.sourceSha256()); out.string(grid.configSha256()); out.string(grid.outputSha256()); out.string(configHash);
        out.integer(tets.size()); for(std::size_t t=0;t<tets.size();++t){for(int n:tets[t])out.integer(n); out.real(volumes[t]);}
        out.integer(mass.size()); for(Eigen::Index n=0;n<mass.size();++n)out.real(mass(n));
        out.integer(k.rows()); out.integer(k.cols()); out.integer(k.nonZeros());
        for(int i=0;i<k.outerSize();++i)for(SparseMatrix::InnerIterator it(k,i);it;++it){out.integer(it.row());out.integer(it.col());out.real(it.value());}
        out.integer(boundary.size()); for(const auto& row:boundary){out.integer(row.nodeIndex);out.integer(row.entries.size());for(const auto& e:row.entries){out.integer(e.handleIndex);out.real(e.weight);}}
        outputHash=out.finish();
    }
};
SourceCageFem::SourceCageFem(const SourceCageGrid& g,const SourceCageFemConfig& c):data_(std::make_shared<const Data>(g,c)) {}
const SourceCageGrid& SourceCageFem::grid() const{return data_->grid;}
const SourceCageFemConfig& SourceCageFem::config() const{return data_->config;}
const std::vector<Tet>& SourceCageFem::tetrahedra() const{return data_->tets;}
const std::vector<double>& SourceCageFem::tetrahedronVolumes() const{return data_->volumes;}
const SourceCageFem::SparseMatrix& SourceCageFem::stiffness() const{return data_->k;}
const Eigen::VectorXd& SourceCageFem::lumpedMass() const{return data_->mass;}
const std::vector<Row>& SourceCageFem::boundaryConstraints() const{return data_->boundary;}
const std::string& SourceCageFem::sourceSha256() const{return grid().sourceSha256();}
const std::string& SourceCageFem::gridConfigSha256() const{return grid().configSha256();}
const std::string& SourceCageFem::gridOutputSha256() const{return grid().outputSha256();}
const std::string& SourceCageFem::configSha256() const{return data_->configHash;}
const std::string& SourceCageFem::outputSha256() const{return data_->outputHash;}
Eigen::MatrixXd SourceCageFem::massInverseStiffness(const Eigen::MatrixXd& x) const {
    if(x.rows()!=lumpedMass().size() || x.cols()<=0 || !x.allFinite())throw std::invalid_argument("source cage FEM requires finite Nnodes x K field, K>0");
    const auto columns=static_cast<std::size_t>(x.cols()),rows=static_cast<std::size_t>(x.rows());
    if(rows>config().maximumEvaluationValues/2/columns)throw std::length_error("RESOURCE_LIMIT: FEM two-field evaluation workspace exceeded before allocation");
    Eigen::MatrixXd result(x.rows(),x.cols());
    for(int i=0;i<stiffness().outerSize();++i)for(Eigen::Index c=0;c<x.cols();++c) {
        long double value=0;
        for(SparseMatrix::InnerIterator it(stiffness(),i);it;++it)value+=static_cast<long double>(it.value())*x(it.col(),c);
        result(i,c)=represented(value/lumpedMass()(i),"source cage FEM M^-1 K field unrepresentable");
    }
    return result;
}
Eigen::MatrixXd SourceCageFem::applyQ(const Eigen::MatrixXd& x) const {
    const auto intermediate=massInverseStiffness(x);
    Eigen::MatrixXd result(x.rows(),x.cols());
    // K was assembled symmetrically, with bit-identical paired coefficients;
    // its row access is therefore also the row access of K^T.
    for(int i=0;i<stiffness().outerSize();++i)for(Eigen::Index c=0;c<x.cols();++c) {
        long double value=0;
        for(SparseMatrix::InnerIterator it(stiffness(),i);it;++it)value+=static_cast<long double>(it.value())*intermediate(it.col(),c);
        result(i,c)=represented(value,"source cage FEM Q field unrepresentable");
    }
    return result;
}
double SourceCageFem::energy(const Eigen::MatrixXd& x) const {
    const auto intermediate=massInverseStiffness(x); long double value=0;
    for(Eigen::Index c=0;c<x.cols();++c)for(Eigen::Index i=0;i<x.rows();++i)
        value+=static_cast<long double>(intermediate(i,c))*intermediate(i,c)*lumpedMass()(i)/2;
    return represented(value,"source cage FEM energy unrepresentable");
}
} // namespace spin2014
