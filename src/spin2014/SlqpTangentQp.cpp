#include "spin2014/SlqpTangentQp.h"
#include "SchurFactorCache.h"
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <set>

namespace spin2014 {
namespace {
double infinityNorm(const Eigen::VectorXd& v) {
    return v.size() ? v.cwiseAbs().maxCoeff() : 0.;
}
bool validHistory(const LbfgsHistory& h, Eigen::Index n) {
    if (h.capacity < 0 || h.s.size() > static_cast<std::size_t>(h.capacity) ||
        h.s.size() != h.y.size() || h.s.size() != h.inverseCurvature.size()) return false;
    for (std::size_t i = 0; i < h.s.size(); ++i) {
        if (h.s[i].size() != n || h.y[i].size() != n ||
            !h.s[i].allFinite() || !h.y[i].allFinite()) return false;
        const double sy = h.s[i].dot(h.y[i]);
        const double normProduct = h.s[i].norm() * h.y[i].norm();
        const double rho = h.inverseCurvature[i];
        if (!std::isfinite(sy) || !std::isfinite(normProduct) ||
            sy <= 1e-12 * normProduct || !std::isfinite(rho) || rho <= 0. ||
            !std::isfinite(sy * rho) || std::abs(sy * rho - 1.) > 1e-12) return false;
    }
    return true;
}

// Independent direct BFGS recursion. Store rank-two updates, never B itself.
struct HessianAction {
    double scale = 1.;
    std::vector<Eigen::VectorXd> bs;
    std::vector<double> sbs;
    const LbfgsHistory& history;

    explicit HessianAction(const LbfgsHistory& h) : history(h) {}
    Eigen::VectorXd apply(const Eigen::VectorXd& v) const {
        Eigen::VectorXd result = v / scale;
        for (std::size_t i = 0; i < bs.size(); ++i) {
            const double sy = history.s[i].dot(history.y[i]);
            result += history.y[i] * (history.y[i].dot(v) / sy) -
                      bs[i] * (bs[i].dot(v) / sbs[i]);
        }
        return result;
    }
    bool initialize() {
        if (!history.s.empty()) {
            const double yy = history.y.back().squaredNorm();
            if (!std::isfinite(yy) || yy <= 0.) return false;
            scale = std::clamp(history.s.back().dot(history.y.back()) / yy, 1e-6, 1e6);
        }
        for (const auto& s : history.s) {
            Eigen::VectorXd v = apply(s);
            const double curvature = s.dot(v);
            if (!v.allFinite() || !std::isfinite(curvature) || curvature <= 0.) return false;
            bs.push_back(std::move(v));
            sbs.push_back(curvature);
        }
        return true;
    }
};
}

static TangentQpResult solveQp(const Eigen::VectorXd& x,
                             const Eigen::VectorXd& g,
                             const Eigen::MatrixXd& a,
                             const Eigen::VectorXd& c,
                             Eigen::Index bounded,
                             const LbfgsHistory& history,
                             const TangentQpSettings& settings, bool finiteBox) {
    TangentQpResult result;
    result.attempted = true;
    result.status = finiteBox ? "BOX_QP_UNKNOWN" : "QP_UNKNOWN";
    auto record = [&](const char* action, int changed = -1) {
        result.trace.push_back({result.workingSet, action, changed,
                               result.equalityResidual, result.stationarity});
    };
    try {
        const Eigen::Index n = x.size(), m = c.size();
        if (n <= 0 || g.size() != n || a.cols() != n || a.rows() != m ||
            bounded < 0 || bounded > n || bounded > std::numeric_limits<int>::max()/2 ||
            !x.allFinite() || !g.allFinite() || !a.allFinite() || !c.allFinite() ||
            settings.maximumFaces < 0 || settings.maximumIterations < 0 ||
            !std::isfinite(settings.activeTolerance) || settings.activeTolerance < 0. ||
            settings.activeTolerance >= .5 || !std::isfinite(settings.residualTolerance) ||
            settings.residualTolerance <= 0. || !validHistory(history, n)) {
            record("INVALID_INPUT_OR_HISTORY");
            return result;
        }
        // Count before allocating face vectors or entering working-set solves.
        int faceCount = 0;
        for (Eigen::Index i = 0; i < bounded; ++i) {
            if (x[i] < 0. || x[i] > 1.) {
                record("INVALID_BOX_PREFIX");
                return result;
            }
            if (!finiteBox && (x[i] <= settings.activeTolerance || x[i] >= 1. - settings.activeTolerance)) {
                if (faceCount == settings.maximumFaces) {
                    record("FACE_LIMIT");
                    return result;
                }
                ++faceCount;
            }
        }
        if (settings.maximumIterations == 0) {
            record("ITERATION_LIMIT");
            return result;
        }
        if (m > n || faceCount > n - m) {
            record("RANK_DEFICIENT");
            return result;
        }
        std::vector<int> faces, signs;
        faces.reserve(faceCount);
        signs.reserve(faceCount);
        for (Eigen::Index i = 0; i < bounded; ++i) {
            if (finiteBox) {
                faces.push_back(static_cast<int>(2*i)); signs.push_back(-1);
                faces.push_back(static_cast<int>(2*i+1)); signs.push_back(1);
            } else if (x[i] <= settings.activeTolerance) {
                faces.push_back(static_cast<int>(i)); signs.push_back(-1);
            } else if (x[i] >= 1. - settings.activeTolerance) {
                faces.push_back(static_cast<int>(i)); signs.push_back(1);
            }
        }
        if (!finiteBox) result.workingSet = faces;
        HessianAction hessian(history);
        const Eigen::VectorXd hg = history.applyInverse(g);
        if (!hg.allFinite() || !hessian.initialize()) {
            record("INVALID_METRIC_ACTION");
            return result;
        }
        std::set<std::vector<int>> seen;
        for (int iteration = 0; iteration < settings.maximumIterations; ++iteration) {
            // Until this set has been solved, no previous set's values apply.
            result.direction.resize(0);
            result.multipliers.resize(0);
            result.equalityResidual = 0.;
            result.stationarity = 0.;
            if (!seen.insert(result.workingSet).second) {
                record("WORKING_SET_CYCLE");
                return result;
            }
            const Eigen::Index rows = m + static_cast<Eigen::Index>(result.workingSet.size());
            if (rows > n) { record("RANK_DEFICIENT"); return result; }
            Eigen::MatrixXd constraints = Eigen::MatrixXd::Zero(rows, n);
            constraints.topRows(m) = a;
            for (std::size_t j = 0; j < result.workingSet.size(); ++j)
                constraints(m + static_cast<Eigen::Index>(j), finiteBox ? result.workingSet[j]/2 : result.workingSet[j]) = 1.;
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(rows);
            rhs.head(m) = c;
            if (finiteBox) for (std::size_t j=0;j<result.workingSet.size();++j) {
                const int face=result.workingSet[j];
                rhs[m+j]=x[face/2]-(face%2);
            }
            Eigen::MatrixXd hc(n, rows);
            for (Eigen::Index j = 0; j < rows; ++j)
                hc.col(j) = history.applyInverse(constraints.row(j).transpose());
            if (!hc.allFinite()) { record("NONFINITE_METRIC_ACTION"); return result; }
            if (rows) {
                const Eigen::MatrixXd schur = constraints * hc;
                const Eigen::VectorXd schurRhs = rhs - constraints * hg;
                if (!schur.allFinite() || !schurRhs.allFinite()) {
                    record("NONFINITE_SCHUR_SYSTEM"); return result;
                }
                Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition(schur);
                decomposition.setThreshold(1e-12);
                if (decomposition.rank() != rows) { record("RANK_DEFICIENT"); return result; }
                result.multipliers = decomposition.solve(schurRhs);
            }
            result.direction = -hg;
            if (rows) result.direction -= hc * result.multipliers;
            // Only equality-fixed coordinates get an exact-zero representation.
            // Free coordinates are never clamped or otherwise repaired.
            for (int index : result.workingSet) {
                if (finiteBox) result.direction[index/2]=(index%2)-x[index/2];
                else result.direction[index] = 0.;
            }
            Eigen::VectorXd equality(rows);
            Eigen::VectorXd stationarity = g + hessian.apply(result.direction);
            // Older supported Eigen versions dereference null data pointers in
            // some zero-row products even when the mathematical result is empty.
            if (rows) {
                equality = constraints * result.direction + rhs;
                stationarity += constraints.transpose() * result.multipliers;
            }
            result.equalityResidual = infinityNorm(equality);
            result.stationarity = infinityNorm(stationarity);
            if (!result.direction.allFinite() || !result.multipliers.allFinite() ||
                !equality.allFinite() || !stationarity.allFinite() ||
                !std::isfinite(result.equalityResidual) || !std::isfinite(result.stationarity) ||
                result.equalityResidual > settings.residualTolerance ||
                result.stationarity > settings.residualTolerance) {
                record("RESIDUAL_CHECK_FAILED"); return result;
            }
            int add = -1, release = -1;
            const double boundTolerance=finiteBox ? std::min(settings.activeTolerance,settings.residualTolerance) : settings.activeTolerance;
            double worstPrimal = boundTolerance, worstDual = boundTolerance;
            // Ascending traversal and strict comparison break equal violations
            // deterministically in favor of the smaller coordinate index.
            for (std::size_t k = 0; k < faces.size(); ++k) {
                const auto found = std::lower_bound(result.workingSet.begin(), result.workingSet.end(), faces[k]);
                if (found == result.workingSet.end() || *found != faces[k]) {
                    const int index=finiteBox ? faces[k]/2 : faces[k];
                    const double outward = finiteBox ? signs[k]*(result.direction[index]+x[index]-(faces[k]%2)) : signs[k] * result.direction[index];
                    // At an exact starting boundary even a sub-tolerance
                    // outward component implies a zero feasible line-search
                    // step. Enforce that face by re-solving, never by clamping.
                    const bool zeroStepFace=finiteBox && signs[k]*result.direction[index]>0. &&
                        x[index]==double(faces[k]%2);
                    if (outward > worstPrimal || (zeroStepFace && add<0)) {
                        worstPrimal = outward; add = faces[k];
                    }
                } else {
                    const Eigen::Index offset = found - result.workingSet.begin();
                    const double wrongSign = -signs[k] * result.multipliers[m + offset];
                    if (wrongSign > worstDual) { worstDual = wrongSign; release = faces[k]; }
                }
            }
            if (add >= 0) {
                if (finiteBox && result.workingSet.size()>=static_cast<std::size_t>(settings.maximumFaces)) {
                    record("FACE_LIMIT");return result;
                }
                record("ADD_FACE", add);
                result.workingSet.insert(std::lower_bound(result.workingSet.begin(), result.workingSet.end(), add), add);
            } else if (release >= 0) {
                record("RELEASE_FACE", release);
                result.workingSet.erase(std::lower_bound(result.workingSet.begin(), result.workingSet.end(), release));
            } else {
                record("CHECKED");
                result.checked = true;
                result.status = finiteBox ? "BOX_QP_FIRST_ORDER_CHECKED" : "QP_FIRST_ORDER_CHECKED";
                return result;
            }
        }
        // The last trace carries the unsolved next set. Clear stale solve data
        // so multipliers cannot appear aligned with that different set.
        result.direction.resize(0);
        result.multipliers.resize(0);
        result.equalityResidual = 0.;
        result.stationarity = 0.;
        record("ITERATION_LIMIT");
    } catch (const std::bad_alloc&) {
        result.checked = false;
        result.status = finiteBox ? "BOX_QP_UNKNOWN" : "QP_UNKNOWN";
        // Allocating an additional trace record may fail for the same reason.
    }
    return result;
}
TangentQpResult solveTangentQp(const Eigen::VectorXd& x,const Eigen::VectorXd& g,
    const Eigen::MatrixXd& a,const Eigen::VectorXd& c,Eigen::Index bounded,
    const LbfgsHistory& history,const TangentQpSettings& settings) {
    return solveQp(x,g,a,c,bounded,history,settings,false);
}
TangentQpResult solveFiniteBoxQp(const Eigen::VectorXd& x,const Eigen::VectorXd& g,
    const Eigen::MatrixXd& a,const Eigen::VectorXd& c,Eigen::Index bounded,
    const LbfgsHistory& history,const TangentQpSettings& settings) {
    return solveQp(x,g,a,c,bounded,history,settings,true);
}

TangentQpResult solvePrimalBoxQp(const Eigen::VectorXd& x,const Eigen::VectorXd& g,
    const Eigen::MatrixXd& a,const Eigen::VectorXd& c,Eigen::Index bounded,
    const LbfgsHistory& history,const TangentQpSettings& settings) {
    TangentQpResult r;r.attempted=true;r.status="BOX_QP_UNKNOWN";
    auto record=[&](const std::string& action,int index=-1){
        r.trace.push_back({r.workingSet,action,index,r.equalityResidual,r.stationarity});
    };
    try {
        const Eigen::Index n=x.size(),m=c.size();
        if(n<=0 || bounded<0 || bounded>n || bounded>std::numeric_limits<int>::max()/2 ||
           g.size()!=n || a.cols()!=n || a.rows()!=m || m>n ||
           !x.allFinite() || !g.allFinite() || !a.allFinite() || !c.allFinite() ||
           settings.maximumFaces < -1 || settings.maximumIterations<=0 ||
           !std::isfinite(settings.activeTolerance) || settings.activeTolerance<0. || settings.activeTolerance>=.5 ||
           !std::isfinite(settings.residualTolerance) || settings.residualTolerance<=0. || !validHistory(history,n)) {
            record("INVALID_INPUT_OR_HISTORY");return r;
        }
        for(Eigen::Index i=0;i<bounded;++i)if(x[i]<0. || x[i]>1.){record("INVALID_BOX_PREFIX");return r;}
        // A feasible displacement cannot be at both ends of a unit-width box.
        // Resolve per call, including after an adaptive partition changes size.
        const int maximumFaces=settings.maximumFaces==-1 ? static_cast<int>(bounded) : settings.maximumFaces;
        const double boundTolerance=std::min(settings.activeTolerance,settings.residualTolerance);
        auto primal=[&](const Eigen::VectorXd& d){
            double error=0.;
            for(Eigen::Index i=0;i<bounded;++i)error=std::max({error,-x[i]-d[i],d[i]-(1.-x[i])});
            return error;
        };
        Eigen::VectorXd z=Eigen::VectorXd::Zero(n);
        int used=0;
        if(infinityNorm(c)>settings.residualTolerance) {
            // Feasibility-only minimum-norm problem, not the objective step.
            auto phaseSettings=settings;phaseSettings.maximumFaces=maximumFaces;
            const auto seed=solveFiniteBoxQp(x,Eigen::VectorXd::Zero(n),a,c,bounded,{},phaseSettings);
            for(const auto& t:seed.trace){auto copy=t;copy.action="PHASE_I_"+copy.action;r.trace.push_back(std::move(copy));}
            used=static_cast<int>(seed.trace.size());
            if(!seed.checked){record("PHASE_I_UNKNOWN");return r;}
            z=seed.direction;
        }
        if(!z.allFinite() || primal(z)>boundTolerance || (m && infinityNorm(a*z+c)>settings.residualTolerance)){
            record("PHASE_I_PRIMAL_CHECK_FAILED");return r;
        }
        HessianAction B(history);const Eigen::VectorXd hg=history.applyInverse(g);
        if(!B.initialize() || !hg.allFinite()){record("INVALID_METRIC_ACTION");return r;}
        // H and A stay fixed within this QP. Retain only current working faces,
        // not every visited face (which could grow into an n*n cache).
        detail::SchurFactorCache factorCache;
        Eigen::MatrixXd equalityColumns, equalitySchur;
        bool equalityColumnsReady=false;
        struct MetricColumn { Eigen::VectorXd value, equalityProduct; };
        std::map<int,MetricColumn> faceColumns;
        Eigen::VectorXd equalityHg=Eigen::VectorXd::Zero(m);
        if(m)equalityHg.noalias()=a*hg;
        for(;used<settings.maximumIterations;++used) {
            r.direction.resize(0);r.multipliers.resize(0);r.equalityResidual=0.;r.stationarity=0.;
            const Eigen::Index rows=m+r.workingSet.size();
            if(rows>n){record("RANK_DEFICIENT");return r;}
            // C=[A; e_i^T] stays implicit: its box rows contain exactly one
            // unit entry. Multiplying those rows densely costs O(n*k*k) work
            // on zeros and allocates another n*(m+k) matrix each iteration.
            Eigen::VectorXd rhs(rows);rhs.head(m)=c;
            for(std::size_t j=0;j<r.workingSet.size();++j){const int f=r.workingSet[j];rhs[m+j]=x[f/2]-f%2;}
            if(!equalityColumnsReady) {
                equalityColumns.resize(n,m);
                for(Eigen::Index j=0;j<m;++j)equalityColumns.col(j)=history.applyInverse(a.row(j).transpose());
                equalitySchur.resize(m,m);
                if(m)equalitySchur.noalias()=a*equalityColumns;
                if(!equalityColumns.allFinite()) {
                    record("INVALID_METRIC_ACTION");return r;
                }
                equalityColumnsReady=true;
            }
            for(auto it=faceColumns.begin();it!=faceColumns.end();) {
                if(!std::binary_search(r.workingSet.begin(),r.workingSet.end(),it->first))it=faceColumns.erase(it);
                else ++it;
            }
            for(std::size_t j=0;j<r.workingSet.size();++j) {
                const int f=r.workingSet[j];auto found=faceColumns.find(f);
                if(found==faceColumns.end()) {
                    MetricColumn column;
                    column.value=history.applyInverse(Eigen::VectorXd::Unit(n,f/2));
                    column.equalityProduct.resize(m);
                    if(m)column.equalityProduct.noalias()=a*column.value;
                    if(!column.value.allFinite()) {
                        record("INVALID_METRIC_ACTION");return r;
                    }
                    faceColumns.emplace(f,std::move(column));
                }
            }
            bool fastSolve=false;
            Eigen::MatrixXd scaled;
            Eigen::VectorXd scale,scaledRhs;
            auto referenceSolve=[&]() {
                factorCache.clear();
                Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> qr(scaled);qr.setThreshold(1e-12);
                if(qr.rank()!=rows){record("RANK_DEFICIENT");return false;}
                r.multipliers=scale.asDiagonal()*qr.solve(scaledRhs);
                return true;
            };
            if(rows) {
                Eigen::MatrixXd S(rows,rows);
                if(m)S.topLeftCorner(m,m)=equalitySchur;
                // Read retained columns directly. H and A are fixed in this
                // QP, so neither a second n*rows materialization nor repeated
                // A*H*C' products are needed. Released faces are still erased.
                for(std::size_t j=0;j<r.workingSet.size();++j) {
                    const auto& column=faceColumns.at(r.workingSet[j]);
                    if(m) {
                        S.col(m+j).head(m)=column.equalityProduct;
                        S.row(m+j).head(m)=equalityColumns.row(r.workingSet[j]/2);
                    }
                    for(std::size_t k=0;k<r.workingSet.size();++k)
                        S(m+k,m+j)=column.value[r.workingSet[k]/2];
                }
                if(!S.allFinite() || (S.diagonal().array()<=0.).any()){record("INVALID_SCHUR_SYSTEM");return r;}
                // Equilibrate only the linear algebra, not the raw acceptance
                // tolerances or the defining objective/constraints.
                scale=S.diagonal().array().sqrt().inverse();
                scaled=scale.asDiagonal()*S*scale.asDiagonal();
                Eigen::VectorXd schurRhs=rhs;
                if(m)schurRhs.head(m)-=equalityHg;
                for(std::size_t j=0;j<r.workingSet.size();++j)
                    schurRhs[m+j]-=hg[r.workingSet[j]/2];
                scaledRhs=scale.asDiagonal()*schurRhs;
                std::vector<int> keys;keys.reserve(rows);
                for(Eigen::Index j=0;j<m;++j)keys.push_back(-1-static_cast<int>(j));
                keys.insert(keys.end(),r.workingSet.begin(),r.workingSet.end());
                Eigen::VectorXd solved;
                fastSolve=factorCache.trySolve(scaled,scaledRhs,keys,solved);
                if(fastSolve)r.multipliers=scale.asDiagonal()*solved;
                else if(!referenceSolve())return r;
            }
            Eigen::VectorXd candidate;
            auto evaluateCandidate=[&]() {
                // By linearity of the same L-BFGS inverse metric used for the
                // cached columns, assemble its covector first and apply H^-1
                // once. This avoids an O(n * working-set-size) sum of cached
                // face columns on every active-set iteration.
                Eigen::VectorXd covector=g;
                if(m)covector.noalias()+=a.transpose()*r.multipliers.head(m);
                for(std::size_t j=0;j<r.workingSet.size();++j)
                    covector[r.workingSet[j]/2]+=r.multipliers[m+j];
                candidate=-history.applyInverse(covector);
            for(int f:r.workingSet)candidate[f/2]=f%2-x[f/2];
            Eigen::VectorXd stationarity=g+B.apply(candidate);
            if(m)stationarity.noalias()+=a.transpose()*r.multipliers.head(m);
            Eigen::VectorXd constraintResidual=rhs;
            if(m)constraintResidual.head(m).noalias()+=a*candidate;
            for(std::size_t j=0;j<r.workingSet.size();++j) {
                const int index=r.workingSet[j]/2;
                stationarity[index]+=r.multipliers[m+j];
                constraintResidual[m+j]+=candidate[index];
            }
            r.equalityResidual=rows?infinityNorm(constraintResidual):0.;r.stationarity=infinityNorm(stationarity);
            return candidate.allFinite() && r.multipliers.allFinite() && stationarity.allFinite() &&
                std::isfinite(r.equalityResidual) && r.equalityResidual<=settings.residualTolerance &&
                std::isfinite(r.stationarity) && r.stationarity<=settings.residualTolerance;
            };
            if(!evaluateCandidate()) {
                // Retry the untouched reference algebra on raw certificate failure.
                if(!fastSolve){record("RESIDUAL_CHECK_FAILED");return r;}
                if(!referenceSolve())return r;
                if(!evaluateCandidate()){record("RESIDUAL_CHECK_FAILED");return r;}
            }
            const Eigen::VectorXd delta=candidate-z;
            double alpha=1.;int blocking=-1;
            for(Eigen::Index i=0;i<bounded;++i) {
                if(delta[i]==0.)continue;
                const int face=2*static_cast<int>(i)+(delta[i]>0.);
                if(std::binary_search(r.workingSet.begin(),r.workingSet.end(),face))continue;
                const double sign=face%2?1.:-1.;
                const double outward=sign*(candidate[i]-(face%2-x[i]));
                // A rounding-sized endpoint overshoot may otherwise add a
                // redundant face to a full-rank working set. Preserve the raw
                // box gate, but enforce exact original endpoints because an
                // outward direction there has zero nonlinear feasible step.
                const bool zeroStep=x[i]==double(face%2) && sign*candidate[i]>0.;
                if(outward<=boundTolerance && !zeroStep)continue;
                const double distance=(face%2-x[i])-z[i];
                const double t=distance/delta[i];
                if(t<alpha){alpha=std::max(0.,t);blocking=face;}
            }
            if(blocking>=0) {
                if(r.workingSet.size()>=static_cast<std::size_t>(maximumFaces)){record("FACE_LIMIT");return r;}
                z+=alpha*delta;
                for(int f:r.workingSet)z[f/2]=f%2-x[f/2];
                z[blocking/2]=blocking%2-x[blocking/2];
                r.equalityResidual=m?infinityNorm(a*z+c):0.;
                if(!z.allFinite() || primal(z)>boundTolerance || r.equalityResidual>settings.residualTolerance){record("PRIMAL_PATH_CHECK_FAILED");return r;}
                record("BLOCKING_FACE",blocking);
                r.workingSet.insert(std::lower_bound(r.workingSet.begin(),r.workingSet.end(),blocking),blocking);
                continue;
            }
            z=candidate;
            int release=-1;double wrong=boundTolerance;
            for(std::size_t j=0;j<r.workingSet.size();++j){const int f=r.workingSet[j];
                const double violation=-(f%2?1.:-1.)*r.multipliers[m+j];
                if(violation>wrong){wrong=violation;release=f;}}
            if(release>=0){record("RELEASE_FACE",release);r.workingSet.erase(std::lower_bound(r.workingSet.begin(),r.workingSet.end(),release));continue;}
            if(primal(z)>boundTolerance){record("PRIMAL_CHECK_FAILED");return r;}
            r.direction=z;r.checked=true;r.status="BOX_QP_FIRST_ORDER_CHECKED";record("CHECKED");return r;
        }
        r.direction.resize(0);r.multipliers.resize(0);r.equalityResidual=0.;r.stationarity=0.;record("ITERATION_LIMIT");
    } catch(const std::bad_alloc&) {r.checked=false;r.status="BOX_QP_UNKNOWN";}
    return r;
}
}
