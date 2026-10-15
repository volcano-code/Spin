#include "spin/FullRankHollowingSolver.h"
#include "spin/DualMonomialCoordinates.h"
#include "manufacturing/ManufacturingGeometry.h"
#include "manufacturing/SlicerValidation.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactMaterialRegion.h"
#include "verification/ExactRationalCoordinates.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <igl/copyleft/cgal/assign_scalar.h>
#endif

namespace {
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now()-start).count();
}
void validate(const spin::SpinDesignRequest& request,
              const spin::QuadricMomentEvaluatorSettings& settings) {
    spin::validateSpinDesignRequest(request, settings.exactBinary64SourceBoundary);
    spin::validateQuadricMomentEvaluatorSettings(settings);
    if(request.settings.equationProfile!=spin::SpinEquationProfile::DualMonomialRawV1)
        throw std::invalid_argument("full-rank hollowing requires DualMonomialRawV1");
    if(request.staticObjective!=spin::StaticObjective::MaximumMass &&
       request.staticObjective!=spin::StaticObjective::LowerCenterOfMass &&
       request.staticObjective!=spin::StaticObjective::ParetoTradeoff)
        throw std::invalid_argument("unsupported full-rank static objective");
}
spin::SpinMomentProblem problemForRequest(const spin::SpinDesignRequest& r) {
    switch(r.mode) {
    case spin::SpinMode::StaticBalance:
        return r.staticObjective==spin::StaticObjective::ParetoTradeoff
            ? spin::SpinMomentProblem::staticBalancePareto(r.staticParetoWeight,spin::staticObjectiveMassScale(r))
            : spin::SpinMomentProblem::staticBalance(r.staticObjective,spin::staticObjectiveMassScale(r));
    case spin::SpinMode::BuoyantStability:
        return spin::SpinMomentProblem::buoyantStability(
            r.buoyancy.liquidDensity*r.buoyancy.displacedLiquidVolume,
            r.buoyancy.centerOfBuoyancyTarget.head<2>());
    case spin::SpinMode::FidgetSpinner: return spin::SpinMomentProblem::fidgetSpinner();
    case spin::SpinMode::YoYo: return spin::SpinMomentProblem::yoYo();
    case spin::SpinMode::SpinningTop: return spin::SpinMomentProblem::spinningTop(r.centerOfMassWeight,r.spinningTopObjective);
    }
    throw std::invalid_argument("unsupported full-rank hollowing mode");
}
template<class T> void appendVector(std::ostream& out,const T& v) {
    out<<v.size()<<':';for(Eigen::Index i=0;i<v.size();++i)out<<v(i)<<';';
}
std::string inputIdentity(const spin::SpinDesignRequest& r,
                          const spin::QuadricMomentEvaluatorSettings& es,
                          const std::vector<spin::QuadricCoefficients10>& guesses) {
    std::ostringstream out;out.imbue(std::locale::classic());out<<std::hexfloat;
    out<<"FullRankHollowingInput-v1;FullRankMaterialDifferential-v3-strict-inactive-selector-margin-1e-10;"
       <<spin::computeQuadricEvaluatorRequestHash(r,es)<<';'
       <<spin::quadricBooleanBackendIdentity()<<';'
       <<spin::spinEquationProfileToken(r.settings.equationProfile)<<';'
       <<spin::spinSolverCoordinateIdentity(r.settings.equationProfile)<<';'
       <<"effective-value-derivative=DiscreteFullRankQuadric;effective-map=FullRankNativeP;"
       <<(guesses.empty()?"default-mode-aware-restarts-v1;":"explicit-restarts;")<<guesses.size()<<';';
    if(es.triangulation.rationalSourceCoordinates)
        out<<"RationalAuthorityFinalSignedTetra-v1;";
    for(const auto& q:guesses)appendVector(out,q.values());
    return manufacturing::computeContentSha256(out.str());
}
std::string topologyIdentity(const spin::TriangleMesh& m) {
    std::ostringstream out;out.imbue(std::locale::classic());
    out<<"FullRankNativeMaterialTopology-v1;"<<m.vertices.rows()<<';'<<m.faces.rows()<<';';
    for(Eigen::Index i=0;i<m.faces.rows();++i)for(int j=0;j<3;++j)out<<m.faces(i,j)<<';';
    return manufacturing::computeContentSha256(out.str());
}
std::string outputIdentity(const spin::FullRankHollowingResult& r) {
    std::ostringstream out;out.imbue(std::locale::classic());out<<std::hexfloat;
    const auto& n=r.numerical;
    out<<"FullRankHollowingOutput-v1;"<<r.inputHash<<';'<<spin::spinSolverStatusToken(n.status)<<';'
       <<n.diagnostic<<';'<<n.eventHash<<';'<<n.traceHash<<';'
       <<spin::spinEquationProfileToken(n.equationProfile)<<';'
       <<n.functionEvaluations<<';'<<n.restartCount<<';'<<n.normalizedEquationResidual<<';'
       <<n.physicalRelativeResidual<<';'<<n.multiplierMu<<';';
    appendVector(out,n.multipliers);
    out<<n.quadric.has_value()<<';';if(n.quadric)appendVector(out,n.quadric->values());
    out<<n.rawEquationResidual.has_value()<<';';if(n.rawEquationResidual)out<<*n.rawEquationResidual<<';';
    out<<r.finalGeometryVerified<<';'<<r.verificationDiagnostic<<';'<<r.finalEvaluation.has_value()<<';';
    if(r.finalEvaluation) {
        out<<r.finalEvaluation->geometryHash<<';'<<r.finalEvaluation->scientificHash<<';';
        appendVector(out,r.finalEvaluation->totalMomentsTarget.values());
    }
    out<<r.independentTotalMoments.has_value()<<';';
    if(r.independentTotalMoments)appendVector(out,r.independentTotalMoments->values());
    out<<r.independentlyRecomputedRawResidual.has_value()<<';';
    if(r.independentlyRecomputedRawResidual)out<<*r.independentlyRecomputedRawResidual<<';';
    out<<r.theoremOneNumericalAudit.has_value()<<';';
    if(r.theoremOneNumericalAudit) {
        const auto& audit=*r.theoremOneNumericalAudit;
        out<<spin::fasterTheoremConditionStatusToken(audit.status)<<';'
           <<audit.unitNormError<<';'<<audit.constraintResidual<<';'
           <<audit.stationarityResidual<<';'<<audit.multiplierMu<<';'
           <<audit.diagnostic<<';'<<audit.proofScope<<';';
    }
    return manufacturing::computeContentSha256(out.str());
}
void requireGate(bool condition,const char* diagnostic) {
    if(!condition)throw std::runtime_error(diagnostic);
}
void verifyFinal(const spin::SpinDesignRequest& r,
                 const spin::QuadricMomentEvaluatorSettings& settings,
                 const spin::MassMomentProblem& problem,
                 spin::FullRankHollowingResult& out) {
    const auto& n=out.numerical;
    requireGate(n.quadric.has_value(),"converged numerical result lacks quadric");
    out.finalEvaluation=spin::evaluateFullRankMaterialMoments(r,*n.quadric,settings);
    const auto& final=*out.finalEvaluation;
    requireGate(final.coefficientJacobian.has_value() && final.coefficientJacobian->allFinite(),
                "final native material derivative unavailable");
    const auto variable=[&]() {
        if(!settings.triangulation.rationalSourceCoordinates)
            return spin::integrateFullRankBoundaryMoments(final.materialBoundaryTarget,r.materialDensity,
                settings.exactBinary64SelfIntersection,settings.exactBinary64MaterialBoundary);
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        const auto exact=spin::integrateFullRankBoundaryMoments(final.materialExactVerticesTarget,
            final.materialBoundaryTarget.faces,r.materialDensity);
        requireGate(exact==final.variableMomentsExactTarget,
            "independent authority moments disagree with same-call exact material moments");
        physics::RawMassMoments10::Vector values;
        for(int i=0;i<10;++i)igl::copyleft::cgal::assign_scalar(
            CGAL::Epeck::FT(verification::parseExactRationalCoordinate(exact.at(i))),values(i));
        MassMoments moments;moments.mass=values(0);
        moments.firstMoment=values.segment<3>(1);
        moments.secondMoment<<values(4),values(9),values(8),values(9),values(5),values(7),
                              values(8),values(7),values(6);
        return physics::RawMassMoments10::fromMassMoments(moments);
#else
        throw std::runtime_error("rational final authority requires exact geometry backend");
#endif
    }();
    const Eigen::Matrix3d rotation=r.targetFrame.axes.transpose();
    out.independentTotalMoments=variable+r.fixedShellMoments.transformed(rotation,-rotation*r.targetFrame.origin);
    const auto& moments=out.independentTotalMoments->values();
    requireGate(moments.allFinite() && std::isfinite(moments.norm()) && moments(0)>0,
                "independent total moments are not finite positive material");
    const double discrepancy=(moments-final.totalMomentsTarget.values()).norm();
    requireGate(std::isfinite(discrepancy) && discrepancy<=1e-10*std::max(1.,final.totalMomentsTarget.values().norm()),
                "independent tetrahedron moments disagree with native total moments");
    const auto evaluated=problem.evaluate(moments);
    requireGate(n.multipliers.size()==evaluated.constraints.size() && n.multipliers.allFinite(),
                "numerical multipliers missing or invalid for final constraints");
    requireGate(std::isfinite(n.multiplierMu) && n.multiplierMu>0,"final multiplier mu must be finite positive");
    // Deliberately reuse numerical multipliers. A new least-squares multiplier
    // solve here could hide stale or inconsistent numerical output.
    const auto s=spin::DualMonomialCoordinates10::fromGeometry(*n.quadric);
    const spin::QuadricCoefficients10::Vector stationarity=evaluated.objectiveGradient-
        evaluated.constraintJacobian.transpose()*n.multipliers+n.multiplierMu*s.values();
    const double raw=std::sqrt(stationarity.squaredNorm()+evaluated.constraints.squaredNorm());
    out.independentlyRecomputedRawResidual=raw;
    requireGate(stationarity.allFinite() && evaluated.constraints.allFinite() && std::isfinite(raw),
                "final raw stationarity/constraint residual is non-finite");
    requireGate(raw<=r.settings.equationTolerance,"independent raw equation residual exceeds original tolerance");
    const double physical=evaluated.constraints.norm()/std::max(1.,moments.norm());
    requireGate(std::isfinite(physical) && physical<=r.settings.physicalRelativeTolerance,
                "independent physical constraint residual exceeds original tolerance");
    requireGate(n.rawEquationResidual && std::isfinite(*n.rawEquationResidual) &&
                std::abs(raw-*n.rawEquationResidual)<=1e-10*std::max(1.,raw),
                "independent raw residual disagrees with numerical result");
    out.theoremOneNumericalAudit=auditFasterTheoremOneNumericalConditions(
        *n.quadric,evaluated,n.multiplierMu,n.multipliers,r.settings.equationTolerance);
    requireGate(
        out.theoremOneNumericalAudit->status==
            spin::FasterTheoremConditionStatus::SatisfiedNumerically,
        "independent Eq. 7 theorem-condition audit failed"
    );
    out.finalGeometryVerified=true;
    out.verificationDiagnostic=settings.triangulation.rationalSourceCoordinates
        ? "strict positive rational authority, independent exact moments and raw-moment KKT, narrow Eq. 7 audit; floating derivative NOT qualified; not a fabrication seal or theorem proof"
        : "strict positive native boundary, independent raw-moment KKT, and narrow Eq. 7 condition audit passed; not a fabrication seal or theorem proof";
}
}
namespace spin {
std::vector<std::string> integrateFullRankBoundaryMoments(
    const ExactBooleanVertexCoordinates& coordinates, const Eigen::MatrixXi& faces, double density) {
    if(!std::isfinite(density) || density<=0)
        throw std::invalid_argument("independent integration density must be finite positive");
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    const auto region=verification::makeExactMaterialRegionFromRationals(coordinates,faces);
    requireGate(region.facts().completed && region.facts().valid && !region.facts().explicitEmpty,
                "independent rational material is not a valid nonempty signed region");
    using Q=CGAL::Gmpq;
    using Q3=std::array<Q,3>;
    using Block=std::array<Q,10>;
    std::vector<Q3> points(coordinates.size());
    for(std::size_t i=0;i<coordinates.size();++i)for(int a=0;a<3;++a)
        points[i][a]=verification::parseExactRationalCoordinate(coordinates[i][a]);
    std::vector<std::optional<Block>> partials;
    for(int f=0;f<faces.rows();++f) {
        const std::array<Q3,3> v{{points.at(faces(f,0)),points.at(faces(f,1)),points.at(faces(f,2))}};
        const Q volume=Q(density)*(v[0][0]*(v[1][1]*v[2][2]-v[1][2]*v[2][1])-
            v[0][1]*(v[1][0]*v[2][2]-v[1][2]*v[2][0])+
            v[0][2]*(v[1][0]*v[2][1]-v[1][1]*v[2][0]))/6;
        Q3 first{};std::array<Q3,3> second{};
        // Independent barycentric E[lambda_i]=1/4 and
        // E[lambda_i lambda_j]=1/10 (i=j), 1/20 otherwise.
        for(int i=0;i<3;++i)for(int a=0;a<3;++a) {
            first[a]+=volume*v[i][a]/4;
            for(int j=0;j<3;++j)for(int b=0;b<3;++b)
                second[a][b]+=volume*v[i][a]*v[j][b]/(i==j?10:20);
        }
        Block block{{volume,first[0],first[1],first[2],second[0][0],second[1][1],
                     second[2][2],second[1][2],second[0][2],second[0][1]}};
        std::size_t level=0;
        while(level<partials.size() && partials[level]) {
            for(int i=0;i<10;++i)block[i]+=(*partials[level])[i];
            partials[level].reset();++level;
        }
        if(level==partials.size())partials.emplace_back(block);
        else partials[level]=block;
    }
    Block moments{};
    for(const auto& block:partials)if(block)for(int i=0;i<10;++i)moments[i]+=(*block)[i];
    requireGate(moments[0]>0,"independent rational integration has nonpositive mass");
    std::vector<std::string> exact;
    for(const auto& value:moments) {
        std::ostringstream out;out.imbue(std::locale::classic());out<<value;
        exact.push_back(out.str());
    }
    return exact;
#else
    (void)coordinates;(void)faces;
    throw std::runtime_error("rational independent integration requires exact geometry backend");
#endif
}
physics::RawMassMoments10 integrateFullRankBoundaryMoments(
    const TriangleMesh& mesh, double density) {
    return integrateFullRankBoundaryMoments(mesh,density,false);
}
physics::RawMassMoments10 integrateFullRankBoundaryMoments(
    const TriangleMesh& mesh, double density, bool exactBinary64SelfIntersection) {
    return integrateFullRankBoundaryMoments(mesh,density,exactBinary64SelfIntersection,false);
}
physics::RawMassMoments10 integrateFullRankBoundaryMoments(
    const TriangleMesh& mesh, double density, bool exactBinary64SelfIntersection,
    bool exactBinary64MaterialBoundary) {
    if(!std::isfinite(density) || density<=0)
        throw std::invalid_argument("independent integration density must be finite positive");
    const auto audit=auditQuadricMaterialBoundary(mesh,exactBinary64SelfIntersection,
        exactBinary64MaterialBoundary);
    requireGate(audit.strictBoundaryValid() && std::isfinite(audit.signedVolume) && audit.signedVolume>0,
                "independent material boundary is not strict closed positive native geometry");
    using L3=Eigen::Matrix<long double,3,1>;
    long double mass=0;L3 first=L3::Zero();
    Eigen::Matrix<long double,3,3> second=Eigen::Matrix<long double,3,3>::Zero();
    // For each oriented tetrahedron (0,a,b,c), barycentric integrals are
    // E[lambda_i]=1/4, E[lambda_i^2]=1/10 and E[lambda_i lambda_j]=1/20.
    // Long-double accumulation is independent of native surface integration.
    for(Eigen::Index f=0;f<mesh.faces.rows();++f) {
        L3 v[3];for(int i=0;i<3;++i)v[i]=mesh.vertices.row(mesh.faces(f,i)).transpose().cast<long double>();
        const long double weightedVolume=static_cast<long double>(density)*v[0].dot(v[1].cross(v[2]))/6;
        mass+=weightedVolume;
        for(int i=0;i<3;++i) {
            first+=weightedVolume*v[i]/4;
            for(int j=0;j<3;++j)second+=weightedVolume*(v[i]*v[j].transpose())/(i==j?10:20);
        }
    }
    MassMoments moments;moments.mass=static_cast<double>(mass);
    moments.firstMoment=first.cast<double>();moments.secondMoment=second.cast<double>();
    requireGate(moments.isFinite() && moments.mass>0,"independent signed-tetrahedron integration failed");
    return physics::RawMassMoments10::fromMassMoments(moments);
}
SpinMomentEvaluator makeFullRankMomentEvaluator(
    const SpinDesignRequest& request, const QuadricMomentEvaluatorSettings& settings) {
    validate(request,settings);
    (void)problemForRequest(request);
    return [request,settings](const QuadricCoefficients10& q) {
        const auto start=std::chrono::steady_clock::now();
        const auto value=evaluateFullRankMaterialMoments(request,q,settings);
        SpinMomentMapEvaluation map;
        map.rawMoments=value.totalMomentsTarget.values();
        map.derivativeDefined=value.coefficientJacobian.has_value();
        if(value.coefficientJacobian) map.momentJacobian=*value.coefficientJacobian;
        map.derivativeDiagnostic=value.diagnostic;
        map.geometryHash=value.geometryHash;
        map.topologyHash=value.topologyHash.empty()
            ? topologyIdentity(value.materialBoundaryTarget)
            : value.topologyHash;
        map.evaluationMapIdentity=value.scientificHash;
        map.booleanBackendIdentity=value.booleanBackendIdentity;
        map.derivativeIdentity=value.effectiveDerivativeMethod;
        map.residualEvaluationSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        return map;
    };
}
static FullRankHollowingResult solveFullRankWithProblemImpl(
    const SpinDesignRequest& request, const QuadricMomentEvaluatorSettings& settings,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const MassMomentProblem& problem, const std::string& problemIdentity,
    const SpinNewtonEventSink& eventSink) {
    const auto start=Clock::now();
    // Invalid caller configuration is a programming error, not a caught
    // numerical or post-convergence geometry failure.
    const auto evaluator=makeFullRankMomentEvaluator(request,settings);
    FullRankHollowingResult out;out.inputHash=inputIdentity(request,settings,initialGuesses);
    if(!problemIdentity.empty())
        out.inputHash=manufacturing::computeContentSha256(
            "ExplicitMassMomentProblem-v1;"+out.inputHash+';'+problemIdentity);
    const auto guesses=initialGuesses.empty()
        ? deterministicInitialQuadrics(request.mode,request.settings):initialGuesses;
    out.numerical=solveMassMomentSystem(problem,request.settings,evaluator,guesses,eventSink);
    if(out.numerical.status==SpinSolverStatus::Converged) {
        const auto verificationStart=Clock::now();
        try { verifyFinal(request,settings,problem,out); }
        catch(const std::exception& error) { out.verificationDiagnostic=error.what(); }
        out.finalVerificationSeconds=elapsed(verificationStart);
    } else out.verificationDiagnostic="not audited: numerical search did not converge; no final material returned";
    out.outputHash=outputIdentity(out);
    out.totalWallSeconds=elapsed(start);
    return out;
}

FullRankHollowingResult solveFullRankQuadricHollowing(
    const SpinDesignRequest& request, const QuadricMomentEvaluatorSettings& settings,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const SpinNewtonEventSink& eventSink) {
    return solveFullRankWithProblemImpl(request,settings,initialGuesses,
        MassMomentProblem(problemForRequest(request)),{},eventSink);
}

FullRankHollowingResult solveFullRankQuadricHollowing(
    const SpinDesignRequest& request, const MassMomentProblem& problem,
    const std::string& problemIdentity, const QuadricMomentEvaluatorSettings& settings,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const SpinNewtonEventSink& eventSink) {
    if(problemIdentity.empty())
        throw std::invalid_argument("explicit mass moment problem requires a bound scientific identity");
    return solveFullRankWithProblemImpl(request,settings,initialGuesses,problem,problemIdentity,eventSink);
}
}
