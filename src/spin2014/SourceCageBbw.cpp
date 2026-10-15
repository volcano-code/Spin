#include "spin2014/SourceCageBbw.h"
#include "verification/Sha256.h"
#include <Highs.h>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <map>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace spin2014 {
namespace {
using Protocol=SourceCageBbwProtocol;
using Terminal=SourceCageBbwTerminal;
using Assessment=SourceCageBbwAssessment;
using Pair=std::pair<int,int>;
using Assembly=std::map<Pair,long double>;
constexpr std::size_t diagnosticCapacity=65536;
constexpr const char* backendRevision="7df0786de3088c832297e5ed821db236d8fab281";

class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {
        unsigned char bytes[8];
        for(int i=0;i<8;++i) bytes[7-i]=static_cast<unsigned char>(x>>(8*i));
        hash_.update(bytes,8);
    }
    void text(const std::string& s) {
        integer(s.size());
        hash_.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());
    }
    void real(double x) {
        std::uint64_t bits;
        static_assert(sizeof(bits)==sizeof(x),"binary64 required");
        std::memcpy(&bits,&x,sizeof(bits));
        integer(bits);
    }
    void field(const Eigen::MatrixXd& w) {
        integer(w.rows()); integer(w.cols());
        for(Eigen::Index i=0;i<w.size();++i) real(w.data()[i]);
    }
    std::string finish() {return hash_.finish();}
};
std::size_t multiply(std::size_t a,std::size_t b) {
    if(b && a>std::numeric_limits<std::size_t>::max()/b)
        throw std::length_error("RESOURCE_LIMIT: BBW count/byte multiplication overflow");
    return a*b;
}
std::size_t add(std::size_t a,std::size_t b) {
    if(a>std::numeric_limits<std::size_t>::max()-b)
        throw std::length_error("RESOURCE_LIMIT: BBW count/byte addition overflow");
    return a+b;
}
double represented(long double x) {
    const double result=static_cast<double>(x);
    if(!std::isfinite(x) || !std::isfinite(result) || (x!=0 && result==0))
        throw std::overflow_error("BBW unrepresentable FEM/Q/reduced arithmetic");
    return result;
}
void validateConfig(const SourceCageBbwConfig& c) {
    if(c.protocol!=Protocol::COUPLED_PARTITION && c.protocol!=Protocol::INDEPENDENT_NORMALIZED)
        throw std::invalid_argument("BBW invalid protocol");
    for(double tolerance:{c.primalTolerance,c.optimalityTolerance,c.backendTolerance})
        if(!std::isfinite(tolerance) || tolerance<=0 || tolerance>=.5)
            throw std::invalid_argument("BBW invalid tolerance");
    if(c.backendTolerance<1e-10 || c.backendTolerance>=c.primalTolerance || c.backendTolerance>=c.optimalityTolerance)
        throw std::invalid_argument("BBW backend tolerance must be supported and stricter than acceptance");
    if(!std::isfinite(c.timeLimitSeconds) || c.timeLimitSeconds<=0 || c.iterationLimit<0)
        throw std::invalid_argument("BBW invalid time or iteration limit");
    const auto indexMax=std::size_t(std::numeric_limits<int>::max());
    const auto byteMax=std::size_t(std::numeric_limits<std::ptrdiff_t>::max());
    if(!c.maximumQPairContributions || c.maximumQPairContributions>byteMax/sizeof(long double) ||
       !c.maximumQNonzeros || c.maximumQNonzeros>indexMax || c.maximumQNonzeros>Assembly().max_size() ||
       !c.maximumVariables || c.maximumVariables>indexMax ||
       !c.maximumDenseValues || c.maximumDenseValues>byteMax/sizeof(double) ||
       !c.maximumBytes || c.maximumBytes>byteMax)
        throw std::length_error("RESOURCE_LIMIT: BBW invalid count/index/byte budget");
}
struct Layout {
    std::vector<int> freeNodes, freeIndex;
    std::size_t denseValues=0, variables=0, bytes=0;
};
Layout layout(const SourceCageFem& fem,const SourceCageBbwConfig& c) {
    const std::size_t n=fem.grid().nodeKeys().size();
    const std::size_t h=fem.grid().handleNodeIndices().size();
    const std::size_t nf=n-fem.boundaryConstraints().size();
    Layout l;
    l.variables=multiply(nf,h);
    if(n>=std::size_t(std::numeric_limits<int>::max()) || h>=std::size_t(std::numeric_limits<int>::max()) || l.variables>c.maximumVariables)
        throw std::length_error("RESOURCE_LIMIT: BBW coupled variables/index budget");
    // Conservative live dense bound: raw/final/boundary/linear/assessment,
    // FEM temporaries and solver-facing arrays. Never an N x N dense matrix.
    l.denseValues=multiply(multiply(n,h),16);
    if(l.denseValues>c.maximumDenseValues || multiply(multiply(n,h),2)>fem.config().maximumEvaluationValues)
        throw std::length_error("RESOURCE_LIMIT: BBW nodal/FEM evaluation values");
    l.bytes=add(multiply(l.denseValues,sizeof(double)),multiply(n,sizeof(long double)+3*sizeof(int)+16));
    l.bytes=add(l.bytes,multiply(h,2*diagnosticCapacity+1024+sizeof(SourceCageBbwAttempt)));
    if(l.bytes>c.maximumBytes)
        throw std::length_error("RESOURCE_LIMIT: BBW dense/trace byte budget");
    l.freeIndex.assign(n,0);
    for(const auto& row:fem.boundaryConstraints()) l.freeIndex[row.nodeIndex]=-1;
    l.freeNodes.reserve(nf);
    for(std::size_t i=0;i<n;++i) if(l.freeIndex[i]!=-1) {
        l.freeIndex[i]=static_cast<int>(l.freeNodes.size());
        l.freeNodes.push_back(static_cast<int>(i));
    }
    return l;
}
double inputScale(const SourceCageFem& fem) {
    std::vector<long double> diagonal(fem.grid().nodeKeys().size(),0);
    const auto& k=fem.stiffness();
    for(int row=0;row<k.rows();++row) {
        for(SourceCageFem::SparseMatrix::InnerIterator a(k,row);a;++a) {
            const long double term=static_cast<long double>(a.value())*a.value()/fem.lumpedMass()(row);
            represented(term);
            diagonal[a.col()]+=term;
        }
    }
    const double scale=represented(*std::max_element(diagonal.begin(),diagonal.end()));
    if(scale<=0) throw std::overflow_error("BBW input Q has no positive finite diagonal scale");
    return scale;
}
Eigen::MatrixXd boundaryField(const SourceCageFem& fem) {
    Eigen::MatrixXd w=Eigen::MatrixXd::Zero(fem.grid().nodeKeys().size(),fem.grid().handleNodeIndices().size());
    for(const auto& row:fem.boundaryConstraints())
        for(const auto& entry:row.entries) w(row.nodeIndex,entry.handleIndex)=entry.weight;
    return w;
}
Assessment assess(const SourceCageFem& fem,const Eigen::MatrixXd& w,
                  const SourceCageBbwConfig& c,const Layout& l,double scale,int onlyHandle=-1) {
    Assessment a;
    const int handles=static_cast<int>(fem.grid().handleNodeIndices().size());
    if(w.rows()!=static_cast<Eigen::Index>(fem.grid().nodeKeys().size()) ||
       w.cols()!=(onlyHandle<0?handles:1) || !w.allFinite()) {
        a.reason="invalid candidate shape or nonfinite value";
        return a;
    }
    a.energyScale=scale;
    a.energy=fem.energy(w);
    Eigen::MatrixXd gradient=fem.applyQ(w);
    for(Eigen::Index i=0;i<gradient.size();++i)
        gradient.data()[i]=represented(static_cast<long double>(gradient.data()[i])/scale);
    a.boundsResidual=std::max({0.,-w.minCoeff(),w.maxCoeff()-1});
    if(onlyHandle<0) a.partitionResidual=(w.rowwise().sum().array()-1).abs().maxCoeff();
    for(const auto& row:fem.boundaryConstraints()) {
        for(int j=0;j<w.cols();++j) {
            const int h=onlyHandle<0?j:onlyHandle;
            double expected=0;
            for(const auto& entry:row.entries) if(entry.handleIndex==h) expected=entry.weight;
            a.boundaryResidual=std::max(a.boundaryResidual,std::abs(w(row.nodeIndex,j)-expected));
        }
    }
    for(int node:l.freeNodes) {
        double lambda=0;
        if(c.protocol==Protocol::COUPLED_PARTITION) {
            long double sum=0;
            int count=0;
            double lower=-std::numeric_limits<double>::infinity();
            double upper=std::numeric_limits<double>::infinity();
            for(int h=0;h<w.cols();++h) {
                const double x=w(node,h),g=gradient(node,h);
                if(x<=c.primalTolerance) lower=std::max(lower,-g);
                else if(x>=1-c.primalTolerance) upper=std::min(upper,-g);
                else {sum-=g; ++count;}
            }
            if(count) lambda=represented(sum/count);
            else if(std::isfinite(lower) && std::isfinite(upper)) lambda=represented((static_cast<long double>(lower)+upper)/2);
            else if(std::isfinite(lower)) lambda=lower;
            else if(std::isfinite(upper)) lambda=upper;
        }
        for(int h=0;h<w.cols();++h) {
            const double x=w(node,h);
            const double r=represented(static_cast<long double>(gradient(node,h))+lambda);
            double lowerDual=0,upperDual=0;
            if(x<=c.primalTolerance) {
                lowerDual=std::max(0.,r);
                a.dualResidual=std::max(a.dualResidual,std::max(0.,-r));
            } else if(x>=1-c.primalTolerance) {
                upperDual=std::max(0.,-r);
                a.dualResidual=std::max(a.dualResidual,std::max(0.,r));
            }
            a.stationarityResidual=std::max(a.stationarityResidual,std::abs(r-lowerDual+upperDual));
            a.complementarityResidual=std::max({a.complementarityResidual,std::abs(lowerDual*x),std::abs(upperDual*(1-x))});
            // Assessment-only projection; never changes the candidate field.
            const double projected=std::max(0.,std::min(1.,x-r));
            a.projectedGradientResidual=std::max(a.projectedGradientResidual,std::abs(x-projected));
        }
    }
    a.unscaledStationarityResidual=represented(static_cast<long double>(a.stationarityResidual)*scale);
    a.unscaledDualResidual=represented(static_cast<long double>(a.dualResidual)*scale);
    a.unscaledComplementarityResidual=represented(static_cast<long double>(a.complementarityResidual)*scale);
    a.unscaledProjectedGradientResidual=represented(static_cast<long double>(a.projectedGradientResidual)*scale);
    a.available=true;
    a.primalPassed=std::max(a.boundsResidual,a.boundaryResidual)<=c.primalTolerance &&
        (c.protocol==Protocol::INDEPENDENT_NORMALIZED || a.partitionResidual<=c.primalTolerance);
    a.optimalityPassed=a.primalPassed && std::max({a.stationarityResidual,a.dualResidual,
        a.complementarityResidual,a.projectedGradientResidual})<=c.optimalityTolerance;
    a.reason=a.optimalityPassed?"original FEM KKT passed":"original FEM primal/KKT rejected";
    return a;
}
void hashAssessment(Digest& d,const Assessment& a) {
    d.integer(a.available); d.integer(a.primalPassed); d.integer(a.optimalityPassed); d.text(a.reason);
    for(double v:{a.energy,a.energyScale,a.boundaryResidual,a.boundsResidual,a.partitionResidual,
                  a.stationarityResidual,a.dualResidual,a.complementarityResidual,a.projectedGradientResidual,
                  a.unscaledStationarityResidual,a.unscaledDualResidual,a.unscaledComplementarityResidual,a.unscaledProjectedGradientResidual}) d.real(v);
}
std::string settingsDescription(const SourceCageBbwConfig& c) {
    std::ostringstream s;
    s.imbue(std::locale::classic());
    s<<std::setprecision(17)<<"HiGHS revision="<<backendRevision
        <<";solver=qpasm;threads=1;parallel=off;random_seed=0;presolve=off"
        <<";qp_regularization_value=0;small_matrix_value=1e-12;output_flag=true"
        <<";log_to_console=false;log_dev_level=0;timeless_log=true"
        <<";primal_feasibility_tolerance="<<c.backendTolerance
        <<";dual_feasibility_tolerance="<<c.backendTolerance
        <<";optimality_tolerance="<<c.backendTolerance
        <<";qp_iteration_limit="<<c.iterationLimit
        <<";qp_nullspace_limit="<<c.maximumVariables
        <<";time_limit_per_call="<<c.timeLimitSeconds
        <<";diagnostic_byte_cap=65536;trace=backend-calls-statuses-counts"
        <<";independent-active-set-polish=exact-reduced-fem-hessian;apply-only-if-independent-kkt-passes"
        <<";independent-redundant-row=0=0;no-implicit-retry";
    return s.str();
}
std::string configDigest(const SourceCageBbwConfig& c) {
    Digest d;
    d.text("SourceCageBbw-v1"); d.text(backendRevision);
    d.text(settingsDescription(c));
    d.integer(static_cast<int>(c.protocol));
    for(double v:{c.primalTolerance,c.optimalityTolerance,c.backendTolerance,c.timeLimitSeconds}) d.real(v);
    d.integer(c.iterationLimit);
    for(auto v:{c.maximumQPairContributions,c.maximumQNonzeros,c.maximumVariables,c.maximumDenseValues,c.maximumBytes}) d.integer(v);
    return d.finish();
}
Assembly assemble(const SourceCageFem& fem,const SourceCageBbwConfig& c,Layout& l,
                  std::size_t& pairs) {
    const auto& k=fem.stiffness();
    // Count every triangular row-outer-product contribution before allocation.
    for(int row=0;row<k.rows();++row) {
        std::size_t count=0;
        for(SourceCageFem::SparseMatrix::InnerIterator a(k,row);a;++a) ++count;
        pairs=add(pairs,multiply(count,add(count,1))/2);
        if(pairs>c.maximumQPairContributions)
            throw std::length_error("RESOURCE_LIMIT: BBW Q pair-contribution budget");
    }
    Assembly q;
    for(int row=0;row<k.rows();++row) {
        for(SourceCageFem::SparseMatrix::InnerIterator a(k,row);a;++a) {
            for(SourceCageFem::SparseMatrix::InnerIterator b(k,row);b;++b) {
                if(a.col()>b.col()) continue;
                const long double term=static_cast<long double>(a.value())*b.value()/fem.lumpedMass()(row);
                represented(term);
                const Pair key{a.col(),b.col()};
                auto entry=q.find(key);
                if(entry==q.end()) {
                    if(q.size()>=c.maximumQNonzeros)
                        throw std::length_error("RESOURCE_LIMIT: BBW unique Q entry budget");
                    // 128 conservatively bounds a map node on this pinned
                    // platform; reduced/backend sparse copies are added later.
                    const auto bytes=add(l.bytes,multiply(add(q.size(),1),128));
                    if(bytes>c.maximumBytes)
                        throw std::length_error("RESOURCE_LIMIT: BBW sparse assembly byte budget");
                    entry=q.emplace(key,0).first;
                }
                entry->second+=term;
                represented(entry->second);
            }
        }
    }
    l.bytes=add(l.bytes,multiply(q.size(),128));
    return q;
}
double unixSeconds() {
    return std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
}
void ensure(HighsStatus status,const char* context) {
    if(status!=HighsStatus::kOk) throw std::runtime_error(context);
}
void configure(Highs& highs,const SourceCageBbwConfig& c,SourceCageBbwAttempt& attempt) {
    ensure(highs.setOptionValue("output_flag",true),"HiGHS output setting rejected");
    ensure(highs.setOptionValue("log_to_console",false),"HiGHS log option rejected");
    ensure(highs.setOptionValue("log_dev_level",0),"HiGHS developer log setting rejected");
    ensure(highs.setOptionValue("timeless_log",true),"HiGHS timeless option rejected");
    ensure(highs.setOptionValue("solver","qpasm"),"HiGHS explicit QP active-set solver unavailable");
    ensure(highs.setOptionValue("threads",1),"HiGHS thread setting rejected");
    ensure(highs.setOptionValue("parallel","off"),"HiGHS parallel setting rejected");
    ensure(highs.setOptionValue("random_seed",0),"HiGHS seed setting rejected");
    ensure(highs.setOptionValue("presolve","off"),"HiGHS presolve setting rejected");
    ensure(highs.setOptionValue("qp_regularization_value",0.),"HiGHS zero QP regularization unavailable");
    ensure(highs.setOptionValue("small_matrix_value",1e-12),"HiGHS matrix threshold rejected");
    ensure(highs.setOptionValue("qp_iteration_limit",c.iterationLimit),"HiGHS QP iteration limit rejected");
    ensure(highs.setOptionValue("qp_nullspace_limit",static_cast<int>(c.maximumVariables)),"HiGHS nullspace limit rejected");
    ensure(highs.setOptionValue("time_limit",c.timeLimitSeconds),"HiGHS time limit rejected");
    for(const char* option:{"primal_feasibility_tolerance","dual_feasibility_tolerance","optimality_tolerance"})
        ensure(highs.setOptionValue(option,c.backendTolerance),"HiGHS stricter tolerance unavailable");
    HighsCallbackFunctionType callback=[](int,const std::string& message,const HighsCallbackOutput*,HighsCallbackInput*,void* data) {
        auto& a=*static_cast<SourceCageBbwAttempt*>(data);
        if(message.size()>diagnosticCapacity-a.diagnostics.size()) a.diagnosticsTruncated=true;
        if(a.diagnostics.size()<diagnosticCapacity)
            a.diagnostics.append(message,0,diagnosticCapacity-a.diagnostics.size());
    };
    ensure(highs.setCallback(callback,&attempt),"HiGHS diagnostic callback rejected");
    ensure(highs.startCallback(HighsCallbackType::kCallbackLogging),"HiGHS logging callback unavailable");
}
bool refineIndependentColumn(const SourceCageFem& fem,const Assembly& q,const Layout& layout,
                             const SourceCageBbwConfig& config,double scale,int handle,
                             Eigen::VectorXd& column,SourceCageBbwAttempt& attempt) {
    attempt.activeSetRefinementAttempted=true;
    const Eigen::MatrixXd initial=column;
    const auto before=assess(fem,initial,config,layout,scale,handle);
    attempt.backendStationarityResidual=before.stationarityResidual;
    attempt.refinedStationarityResidual=before.stationarityResidual;
    if(before.optimalityPassed) return false;
    const std::size_t nodeCount=fem.grid().nodeKeys().size();
    std::vector<bool> boundary(nodeCount,false);
    for(const auto& row:fem.boundaryConstraints()) boundary[row.nodeIndex]=true;
    Eigen::VectorXd fixed=Eigen::VectorXd::Zero(static_cast<Eigen::Index>(nodeCount));
    for(const auto& row:fem.boundaryConstraints())
        for(const auto& entry:row.entries)
            if(entry.handleIndex==handle) fixed(row.nodeIndex)=entry.weight;
    std::vector<int> freeIndex(nodeCount,-1),freeNodes;
    freeNodes.reserve(layout.freeNodes.size());
    for(std::size_t node=0;node<nodeCount;++node) {
        if(boundary[node]) continue;
        if(column(static_cast<Eigen::Index>(node))<=config.primalTolerance) fixed(node)=0.;
        else if(column(static_cast<Eigen::Index>(node))>=1.-config.primalTolerance) fixed(node)=1.;
        else {
            freeIndex[node]=static_cast<int>(freeNodes.size());
            freeNodes.push_back(static_cast<int>(node));
        }
    }
    if(freeNodes.empty()) {
        attempt.refinedStationarityResidual=before.stationarityResidual;
        return false;
    }
    const long double refinementPeakBytes=static_cast<long double>(layout.bytes)+
        static_cast<long double>(q.size())*2.*sizeof(Eigen::Triplet<double>)+
        static_cast<long double>(freeNodes.size())*freeNodes.size()*4.*sizeof(double);
    if(refinementPeakBytes>static_cast<long double>(config.maximumBytes)) return false;
    using Sparse=Eigen::SparseMatrix<double>;
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(q.size()*2);
    Eigen::VectorXd rhs=Eigen::VectorXd::Zero(static_cast<Eigen::Index>(freeNodes.size()));
    for(const auto& item:q) {
        const int i=item.first.first,j=item.first.second;
        const double value=represented(item.second/scale);
        const int fi=freeIndex[static_cast<std::size_t>(i)],fj=freeIndex[static_cast<std::size_t>(j)];
        if(fi>=0 && fj>=0) {
            entries.emplace_back(fi,fj,value);
            if(i!=j) entries.emplace_back(fj,fi,value);
        } else if(fi>=0) rhs(fi)-=value*fixed(j);
        else if(fj>=0 && i!=j) rhs(fj)-=value*fixed(i);
    }
    Sparse reduced(static_cast<int>(freeNodes.size()),static_cast<int>(freeNodes.size()));
    reduced.setFromTriplets(entries.begin(),entries.end());
    reduced.makeCompressed();
    Eigen::SimplicialLDLT<Sparse> factor;
    factor.compute(reduced);
    if(factor.info()!=Eigen::Success || !factor.vectorD().allFinite() ||
       (factor.vectorD().array()<=0.).any()) {
        attempt.refinedStationarityResidual=before.stationarityResidual;
        return false;
    }
    const Eigen::VectorXd solution=factor.solve(rhs);
    if(factor.info()!=Eigen::Success || !solution.allFinite()) {
        attempt.refinedStationarityResidual=before.stationarityResidual;
        return false;
    }
    Eigen::VectorXd refined=fixed;
    for(std::size_t i=0;i<freeNodes.size();++i)
        refined(freeNodes[i])=solution(static_cast<Eigen::Index>(i));
    const auto after=assess(fem,Eigen::MatrixXd(refined),config,layout,scale,handle);
    attempt.refinedStationarityResidual=after.stationarityResidual;
    if(!after.primalPassed || !after.optimalityPassed) return false;
    column=std::move(refined);
    attempt.activeSetRefinementApplied=true;
    return true;
}
}

struct SourceCageBbw::Data {
    SourceCageFem fem;
    SourceCageBbwConfig config;
    Terminal terminal=Terminal::BACKEND_ERROR;
    std::string reason="run not completed",configHash,outputHash,settings;
    Eigen::MatrixXd raw,final;
    Assessment rawAssessment,finalAssessment;
    std::vector<SourceCageBbwAttempt> attempts;
    double normalizationChange=0,scale=0;
    std::size_t pairs=0,nonzeros=0,variables=0,bytes=0;
    Data(const SourceCageFem& f,const SourceCageBbwConfig& c)
        :fem(f),config(c),configHash(configDigest(c)),settings(settingsDescription(c)) {}
    void solve();
    void digest();
};

void SourceCageBbw::Data::solve() {
    validateConfig(config);
    auto l=layout(fem,config);
    variables=l.variables;
    scale=inputScale(fem);
    auto q=assemble(fem,config,l,pairs);
    nonzeros=q.size();
    const int nf=static_cast<int>(l.freeNodes.size());
    const int handles=static_cast<int>(fem.grid().handleNodeIndices().size());
    const bool coupled=config.protocol==Protocol::COUPLED_PARTITION;
    const int blocks=coupled?handles:1;
    std::size_t reducedEntries=0;
    for(const auto& item:q)
        if(l.freeIndex[item.first.first]>=0 && l.freeIndex[item.first.second]>=0 && item.second!=0) ++reducedEntries;
    const auto hessianEntries=multiply(reducedEntries,blocks);
    const auto backendVariables=multiply(nf,blocks);
    if(hessianEntries>std::size_t(std::numeric_limits<int>::max()))
        throw std::length_error("RESOURCE_LIMIT: BBW replicated Hessian index overflow");
    // Our model plus the HiGHS model copy. Internal backend allocator peak
    // remains covered by the external process limit, not exact accounting.
    l.bytes=add(l.bytes,multiply(hessianEntries,2*(sizeof(double)+sizeof(HighsInt))));
    l.bytes=add(l.bytes,multiply(add(backendVariables,1),8*(sizeof(double)+sizeof(HighsInt))));
    if(l.bytes>config.maximumBytes)
        throw std::length_error("RESOURCE_LIMIT: BBW reduced/backend sparse byte budget");
    bytes=l.bytes;
    raw=boundaryField(fem);
    const Eigen::MatrixXd linear=fem.applyQ(raw);
    HighsModel model;
    auto& lp=model.lp_;
    lp.num_col_=static_cast<HighsInt>(backendVariables);
    // Explicit from the outset, never a retry: the pinned backend's zero-row
    // bounded startup uses an incomplete sparse Cholesky and a 1e-4 cutoff.
    // A redundant 0=0 row preserves the exact independent box QP and routes
    // all columns through its general feasible-start/active-set algorithm.
    lp.num_row_=coupled?nf:1;
    lp.col_lower_.assign(backendVariables,0);
    lp.col_upper_.assign(backendVariables,1);
    lp.col_cost_.resize(backendVariables);
    lp.row_lower_.assign(lp.num_row_,coupled?1:0);
    lp.row_upper_.assign(lp.num_row_,coupled?1:0);
    lp.a_matrix_.format_=MatrixFormat::kColwise;
    lp.a_matrix_.start_.clear();
    lp.a_matrix_.start_.reserve(backendVariables+1);
    lp.a_matrix_.index_.reserve(coupled?backendVariables:0);
    lp.a_matrix_.value_.reserve(coupled?backendVariables:0);
    for(int j=0;j<lp.num_col_;++j) {
        lp.a_matrix_.start_.push_back(lp.a_matrix_.value_.size());
        if(coupled) {
            lp.a_matrix_.index_.push_back(j%nf);
            lp.a_matrix_.value_.push_back(1);
        }
    }
    lp.a_matrix_.start_.push_back(lp.a_matrix_.value_.size());
    auto& hs=model.hessian_;
    hs.dim_=lp.num_col_;
    hs.format_=HessianFormat::kTriangular;
    hs.start_.clear();
    hs.start_.reserve(backendVariables+1);
    hs.index_.reserve(hessianEntries);
    hs.value_.reserve(hessianEntries);
    for(int block=0;block<blocks;++block) {
        for(int col=0;col<nf;++col) {
            hs.start_.push_back(hs.value_.size());
            const int node=l.freeNodes[col];
            auto item=q.lower_bound({node,node});
            while(item!=q.end() && item->first.first==node) {
                const int row=l.freeIndex[item->first.second];
                if(row>=0 && item->second!=0) {
                    hs.index_.push_back(block*nf+row);
                    hs.value_.push_back(represented(item->second/scale));
                }
                ++item;
            }
        }
    }
    hs.start_.push_back(hs.value_.size());
    if(nf==0) {
        rawAssessment=assess(fem,raw,config,l,scale);
        final=raw;
        auto coupledConfig=config; coupledConfig.protocol=Protocol::COUPLED_PARTITION;
        finalAssessment=assess(fem,final,coupledConfig,l,scale);
        terminal=rawAssessment.optimalityPassed && finalAssessment.optimalityPassed?
            Terminal::VALIDATED:Terminal::KKT_REJECTED;
        reason="no free variables; prescribed field evaluated against original FEM";
        return;
    }
    attempts.reserve(coupled?1:handles);
    for(int call=0;call<(coupled?1:handles);++call) {
        for(int block=0;block<blocks;++block) {
            const int h=coupled?block:call;
            for(int i=0;i<nf;++i)
                lp.col_cost_[block*nf+i]=represented(static_cast<long double>(linear(l.freeNodes[i],h))/scale);
        }
        attempts.emplace_back();
        auto& attempt=attempts.back();
        attempt.diagnostics.reserve(diagnosticCapacity);
        attempt.handle=coupled?-1:call;
        attempt.startUnixSeconds=unixSeconds();
        const auto start=std::chrono::steady_clock::now();
        try {
            Highs highs;
            configure(highs,config,attempt);
            const auto passStatus=highs.passModel(model);
            attempt.passModelStatus=static_cast<int>(passStatus);
            ensure(passStatus,"HiGHS rejected exact reduced QP");
            attempt.runInvoked=true;
            const auto status=highs.run();
            attempt.returnStatus=static_cast<int>(status);
            const auto modelStatus=highs.getModelStatus();
            attempt.modelStatus=static_cast<int>(modelStatus);
            attempt.terminalReason=highs.modelStatusToString(modelStatus);
            attempt.qpIterations=highs.getInfo().qp_iteration_count;
            const auto& solution=highs.getSolution();
            attempt.candidateAvailable=solution.value_valid && solution.col_value.size()==backendVariables;
            if(attempt.candidateAvailable) {
                for(double x:solution.col_value) if(!std::isfinite(x)) attempt.candidateAvailable=false;
            }
            if(attempt.candidateAvailable) {
                for(int block=0;block<blocks;++block) {
                    const int h=coupled?block:call;
                    for(int i=0;i<nf;++i) raw(l.freeNodes[i],h)=solution.col_value[block*nf+i];
                }
                if(!coupled) {
                    Eigen::VectorXd candidate=raw.col(call);
                    refineIndependentColumn(fem,q,l,config,scale,call,candidate,attempt);
                    if(attempt.activeSetRefinementApplied) raw.col(call)=candidate;
                }
                attempt.assessment=coupled?assess(fem,raw,config,l,scale):
                    assess(fem,Eigen::MatrixXd(raw.col(call)),config,l,scale,call);
            }
            if(modelStatus==HighsModelStatus::kIterationLimit) terminal=Terminal::ITERATION_LIMIT;
            else if(modelStatus==HighsModelStatus::kTimeLimit) terminal=Terminal::TIME_LIMIT;
            else if(status!=HighsStatus::kOk || modelStatus!=HighsModelStatus::kOptimal || !attempt.candidateAvailable)
                terminal=Terminal::BACKEND_ERROR;
            else if(!attempt.assessment.optimalityPassed) terminal=Terminal::KKT_REJECTED;
            else terminal=Terminal::VALIDATED;
        } catch(const std::bad_alloc&) {
            attempt.terminalReason="RESOURCE_LIMIT: backend call allocation failed";
            terminal=Terminal::RESOURCE_LIMIT;
        } catch(const std::length_error& e) {
            attempt.terminalReason=e.what();
            terminal=Terminal::RESOURCE_LIMIT;
        } catch(const std::overflow_error& e) {
            attempt.terminalReason=e.what();
            terminal=Terminal::NUMERICAL_ERROR;
        } catch(const std::exception& e) {
            attempt.terminalReason=e.what();
            terminal=Terminal::BACKEND_ERROR;
        }
        attempt.endUnixSeconds=unixSeconds();
        attempt.elapsedSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(terminal!=Terminal::VALIDATED) {
            reason="UNKNOWN feasibility; backend call "+std::to_string(call)+": "+attempt.terminalReason;
            if(terminal==Terminal::KKT_REJECTED) reason+="; original FEM KKT rejected";
            return;
        }
    }
    rawAssessment=assess(fem,raw,config,l,scale);
    if(!rawAssessment.optimalityPassed) {
        terminal=Terminal::KKT_REJECTED;
        reason="complete raw field failed original FEM KKT; UNKNOWN feasibility";
        return;
    }
    final=raw;
    if(!coupled) {
        for(int i=0;i<final.rows();++i) {
            const double sum=final.row(i).sum();
            if(!std::isfinite(sum) || sum<=0) {
                final.resize(0,0);
                terminal=Terminal::NORMALIZATION_REJECTED;
                reason="nonpositive/nonfinite normalization denominator; UNKNOWN feasibility";
                return;
            }
            final.row(i)/=sum;
        }
        if(!final.allFinite()) {
            final.resize(0,0);
            terminal=Terminal::NORMALIZATION_REJECTED;
            reason="nonfinite normalized field; UNKNOWN feasibility";
            return;
        }
        normalizationChange=(final-raw).cwiseAbs().maxCoeff();
    }
    auto coupledConfig=config;
    coupledConfig.protocol=Protocol::COUPLED_PARTITION;
    finalAssessment=assess(fem,final,coupledConfig,l,scale);
    if(!finalAssessment.primalPassed || (coupled && !finalAssessment.optimalityPassed)) {
        terminal=Terminal::KKT_REJECTED;
        reason="final field rejected; UNKNOWN feasibility";
        final.resize(0,0);
        return;
    }
    terminal=Terminal::VALIDATED;
    reason=coupled?"coupled field validated by original FEM KKT":
        "all raw column QPs validated; explicit normalization validated for primal constraints; final coupled KKT assessed separately";
}
void SourceCageBbw::Data::digest() {
    Digest d;
    d.text(configHash); d.text(fem.outputSha256()); d.integer(static_cast<int>(terminal)); d.text(reason);
    d.real(scale); d.real(normalizationChange);
    d.integer(pairs); d.integer(nonzeros); d.integer(variables); d.integer(bytes);
    d.field(raw); d.field(final);
    hashAssessment(d,rawAssessment); hashAssessment(d,finalAssessment);
    d.integer(attempts.size());
    for(const auto& a:attempts) {
        d.integer(a.handle); d.integer(a.passModelStatus); d.integer(a.runInvoked);
        d.integer(a.returnStatus); d.integer(a.modelStatus);
        d.integer(a.qpIterations); d.integer(a.candidateAvailable); d.text(a.terminalReason);
        d.integer(a.activeSetRefinementAttempted); d.integer(a.activeSetRefinementApplied);
        d.real(a.backendStationarityResidual); d.real(a.refinedStationarityResidual);
        hashAssessment(d,a.assessment);
    }
    outputHash=d.finish();
}
SourceCageBbw::SourceCageBbw(const SourceCageFem& fem,const SourceCageBbwConfig& config) {
    auto data=std::make_shared<Data>(fem,config);
    try {data->solve();}
    catch(const std::invalid_argument& e) {data->terminal=Terminal::INVALID_CONFIG; data->reason=e.what();}
    catch(const std::length_error& e) {data->terminal=Terminal::RESOURCE_LIMIT; data->reason=e.what();}
    catch(const std::bad_alloc&) {data->terminal=Terminal::RESOURCE_LIMIT; data->reason="RESOURCE_LIMIT: allocation failed; UNKNOWN feasibility";}
    catch(const std::overflow_error& e) {data->terminal=Terminal::NUMERICAL_ERROR; data->reason=e.what();}
    catch(const std::exception& e) {data->terminal=Terminal::BACKEND_ERROR; data->reason=e.what();}
    if(data->terminal!=Terminal::VALIDATED) data->final.resize(0,0);
    data->digest();
    data_=std::move(data);
}
Assessment assessSourceCageBbw(const SourceCageFem& fem,const Eigen::MatrixXd& w,const SourceCageBbwConfig& c) {
    try {
        validateConfig(c);
        const auto l=layout(fem,c);
        return assess(fem,w,c,l,inputScale(fem));
    } catch(const std::exception& e) {
        Assessment a;
        a.reason=e.what();
        return a;
    }
}
bool SourceCageBbw::success() const {return data_->terminal==Terminal::VALIDATED;}
SourceCageBbwTerminal SourceCageBbw::terminal() const {return data_->terminal;}
const std::string& SourceCageBbw::terminalReason() const {return data_->reason;}
const SourceCageBbwConfig& SourceCageBbw::config() const {return data_->config;}
const std::string& SourceCageBbw::backendSettings() const {return data_->settings;}
const SourceCageFem& SourceCageBbw::fem() const {return data_->fem;}
const Eigen::MatrixXd& SourceCageBbw::rawWeights() const {return data_->raw;}
const Eigen::MatrixXd& SourceCageBbw::finalWeights() const {return data_->final;}
const Assessment& SourceCageBbw::rawAssessment() const {return data_->rawAssessment;}
const Assessment& SourceCageBbw::finalAssessment() const {return data_->finalAssessment;}
const std::vector<SourceCageBbwAttempt>& SourceCageBbw::attempts() const {return data_->attempts;}
double SourceCageBbw::maximumNormalizationChange() const {return data_->normalizationChange;}
double SourceCageBbw::energyScale() const {return data_->scale;}
std::size_t SourceCageBbw::qPairContributions() const {return data_->pairs;}
std::size_t SourceCageBbw::qNonzeros() const {return data_->nonzeros;}
std::size_t SourceCageBbw::variablesPerCoupledProblem() const {return data_->variables;}
std::size_t SourceCageBbw::checkedStorageBytes() const {return data_->bytes;}
const std::string& SourceCageBbw::configSha256() const {return data_->configHash;}
const std::string& SourceCageBbw::outputSha256() const {return data_->outputHash;}
}
