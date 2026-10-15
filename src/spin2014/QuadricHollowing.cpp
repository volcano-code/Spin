#include "spin2014/QuadricHollowing.h"
#include "verification/ExactMaterialRegion.h"
#include "manufacturing/SlicerValidation.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactLibiglBoolean.h"
#include "verification/ExactRationalCoordinates.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <igl/copyleft/cgal/assign_scalar.h>
#endif

namespace spin2014
{
spin::FullRankMomentEvaluation exportSpin2014RationalCavity(
    const spin::SpinDesignRequest& request,const spin::FullRankMomentEvaluation& retained)
{
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    const auto require=[](bool ok,const char* message){if(!ok)throw std::runtime_error(message);};
    require(!retained.sourceConstructionIdentity.empty() && !retained.scientificHash.empty(),
            "rational cavity requires the same-call source/material identity");
    const auto independent=spin::integrateCanonicalFullRankBoundaryMoments(
        retained.materialExactVerticesTarget,retained.materialBoundaryTarget.faces,request.materialDensity);
    require(independent==retained.variableMomentsExactTarget,
            "rational cavity retained authority moments or density differ");
    require(retained.designExactVerticesTarget.size()==std::size_t(request.designDomain.vertices.rows()),
            "rational cavity design authority is missing or mismatched");
    using Q=CGAL::Gmpq;
    for(int v=0;v<request.designDomain.vertices.rows();++v)for(int a=0;a<3;++a) {
        Q coordinate=0;
        for(int k=0;k<3;++k)coordinate+=Q(request.targetFrame.axes(k,a))*
            (Q(request.designDomain.vertices(v,k))-Q(request.targetFrame.origin(k)));
        require(coordinate==verification::parseExactRationalCoordinate(retained.designExactVerticesTarget[v][a]),
                "rational cavity request/frame differs from the actual design authority");
    }
    const auto design=verification::makeExactMaterialRegionFromRationals(
        retained.designExactVerticesTarget,request.designDomain.faces);
    const auto material=verification::makeExactMaterialRegionFromRationals(
        retained.materialExactVerticesTarget,retained.materialBoundaryTarget.faces);
    const auto subset=verification::checkExactMaterialSubset(material,design);
    require(subset.completed && subset.aSubsetB.value_or(false),
            "rational cavity retained material is not inside the design authority");
    const auto difference=verification::computeExactLibiglBooleanConstructionFromRationals(
        retained.designExactVerticesTarget,request.designDomain.faces,
        retained.materialExactVerticesTarget,retained.materialBoundaryTarget.faces,
        verification::ExactMaterialBooleanOperation::Difference);
    require(difference.succeeded,"same-retained rational cavity Boolean failed");
    const auto empty=difference.faces.rows()==0;
    const auto cavity=empty?verification::emptyExactMaterialRegion():
        verification::makeExactMaterialRegionFromRationals(difference.vertices,difference.faces);
    const auto partition=verification::checkExactMaterialPartition(design,{material,cavity});
    require(partition.completed && partition.completePartition.value_or(false),
            "same-retained rational cavity does not form a complete disjoint design partition");
    spin::FullRankMomentEvaluation out;
    out.materialExactVerticesTarget=difference.vertices;
    out.materialBoundaryTarget.faces=difference.faces;
    out.materialBirthFaces=difference.birthFaces;
    out.designExactVerticesTarget=retained.designExactVerticesTarget;
    out.variableMomentsExactTarget=empty?std::vector<std::string>(10,"0/1"):
        spin::integrateCanonicalFullRankBoundaryMoments(difference.vertices,difference.faces,request.materialDensity);
    const auto whole=spin::integrateCanonicalFullRankBoundaryMoments(
        retained.designExactVerticesTarget,request.designDomain.faces,request.materialDensity);
    for(int i=0;i<10;++i)require(verification::parseExactRationalCoordinate(whole[i])==
        verification::parseExactRationalCoordinate(independent[i])+
        verification::parseExactRationalCoordinate(out.variableMomentsExactTarget[i]),
        "same-retained rational cavity independent moments do not partition the design");
    out.materialBoundaryTarget.vertices.resize(difference.vertices.size(),3);
    for(int v=0;v<out.materialBoundaryTarget.vertices.rows();++v)for(int a=0;a<3;++a)
        igl::copyleft::cgal::assign_scalar(CGAL::Epeck::FT(
            verification::parseExactRationalCoordinate(difference.vertices[v][a])),
            out.materialBoundaryTarget.vertices(v,a));
    physics::RawMassMoments10::Vector values;
    for(int i=0;i<10;++i)igl::copyleft::cgal::assign_scalar(CGAL::Epeck::FT(
        verification::parseExactRationalCoordinate(out.variableMomentsExactTarget[i])),values(i));
    MassMoments moments;moments.mass=values(0);moments.firstMoment=values.segment<3>(1);
    moments.secondMoment<<values(4),values(9),values(8),values(9),values(5),values(7),
                          values(8),values(7),values(6);
    out.variableMomentsTarget=physics::RawMassMoments10::fromMassMoments(moments);
    out.totalMomentsTarget=out.variableMomentsTarget;
    out.sourceScientificHash=retained.scientificHash;
    std::ostringstream identity;identity.imbue(std::locale::classic());
    identity<<"Spin2014SameRetainedRationalCavity-v1;"<<retained.geometryHash<<';'
        <<retained.sourceConstructionIdentity<<';'<<partition.partitionDigest<<';';
    for(const auto& point:difference.vertices)for(const auto& x:point)identity<<x<<';';
    for(int f=0;f<difference.faces.rows();++f)for(int a=0;a<3;++a)identity<<difference.faces(f,a)<<';';
    out.geometryHash=manufacturing::computeContentSha256(identity.str());
    out.scientificHash=out.geometryHash;out.booleanBackendIdentity=retained.booleanBackendIdentity;
    out.effectiveDerivativeMethod="NOT_DEFINED_SAME_RETAINED_COMPLEMENT_EXPORT";
    out.diagnostic="exact represented DESIGN-minus-SAME-retained and independent partition; no cavity Jacobian, floating qualification or manufacturing seal";
    return out;
#else
    (void)request;(void)retained;
    throw std::runtime_error("rational cavity export requires exact geometry backend");
#endif
}

spin::SpinDesignRequest prepareSpin2014QuadricRequest(
    const spin::SpinDesignRequest& request,const QuadricObjectiveSettings& objective,
    double fixedDensity)
{
    validateSpin2014QuadricObjectiveSettings(objective);
    spin::validateSpinDesignRequest(request);
    if(request.mode!=spin2014ObjectiveMode(objective.objective))
        throw std::invalid_argument("2014 objective and requested physical mode differ");
    if(!std::isfinite(fixedDensity)||fixedDensity<=0)
        throw std::invalid_argument("fixed shell density must be finite positive");
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    const auto outer=verification::makeExactMaterialRegion(request.outerSurface.vertices,request.outerSurface.faces);
    const auto design=verification::makeExactMaterialRegion(request.designDomain.vertices,request.designDomain.faces);
    if(!outer.facts().valid || !design.facts().valid)
        throw std::invalid_argument("2014 quadric input is not an independently valid exact solid");
    const auto subset=verification::checkExactMaterialSubset(design,outer);
    if(!subset.completed || !subset.aSubsetB.has_value())
        throw std::runtime_error("2014 quadric input containment UNKNOWN: "+subset.error);
    if(!*subset.aSubsetB)
        throw std::invalid_argument("design material is outside the unchanged outer solid");
#else
    throw std::runtime_error("2014 quadric input admission requires exact geometry checks");
#endif
    const auto measured=physics::RawMassMoments10::fromMassMoments(
        computeSurfaceDifferenceMassMomentsExtended(
            request.outerSurface.vertices,request.outerSurface.faces,
            request.designDomain.vertices,request.designDomain.faces,fixedDensity));
    if(!measured.isFinite() || measured.mass()<0)
        throw std::invalid_argument("integrated fixed shell is not finite nonnegative material");
    if(request.fixedShellMoments.values().norm()!=0 &&
        (measured.values()-request.fixedShellMoments.values()).norm()>
            1e-10*std::max(1.,measured.values().norm()))
        throw std::invalid_argument("supplied fixed moments disagree with actual fixed shell geometry/density");
    auto owned=request;
    owned.fixedShellMoments=measured;
    spin::validateSpinDesignRequest(owned);
    return owned;
}

QuadricMaterialRun solveSpin2014QuadricHollowing(
    const spin::SpinDesignRequest& request,const QuadricObjectiveSettings& objective,
    double fixedDensity,const spin::QuadricMomentEvaluatorSettings& evaluatorSettings,
    const std::vector<spin::QuadricCoefficients10>& initialGuesses,
    const spin::SpinNewtonEventSink& eventSink)
{
    if(objective.objective==Spin2014Objective::StaticLowering)
        throw std::invalid_argument("static squared-height requests require the dedicated plane-family adapter");
    if(request.settings.equationProfile!=spin::SpinEquationProfile::DualMonomialRawV1 ||
       request.settings.equationTolerance>1e-12 || request.settings.maximumRestarts>8 ||
       request.settings.maximumFunctionEvaluationsPerRestart>60)
        throw std::invalid_argument("2014 quadric profile requires strict raw Newton tolerance and bounded attempts");
    if(evaluatorSettings.derivativeMethod!=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric)
        throw std::invalid_argument("2014 quadric profile requires the actual full-rank discrete moment derivative");
    QuadricMaterialRun out;
    out.objectiveSettings=objective;
    out.objectiveProfile=spin2014QuadricObjectiveProfile();
    out.objectiveIdentity=spin2014QuadricObjectiveHash(objective);
    out.measuredRequest=prepareSpin2014QuadricRequest(request,objective,fixedDensity);
    out.fixedMaterialDensity=fixedDensity;
    const auto problem=makeSpin2014QuadricMomentProblem(objective);
    out.solve=spin::solveCanonicalFullRankQuadricHollowing(
        out.measuredRequest,problem,out.objectiveIdentity,evaluatorSettings,initialGuesses,eventSink);
    if(out.solve.numerical.status==spin::SpinSolverStatus::Converged && out.solve.finalGeometryVerified) {
        try {
            if(evaluatorSettings.triangulation.rationalSourceCoordinates) {
                out.cavityEvaluation=exportSpin2014RationalCavity(
                    out.measuredRequest,*out.solve.finalEvaluation);
            } else {
                const auto complement=spin::QuadricCoefficients10::fromUnnormalized(
                    -out.solve.numerical.quadric->values());
                out.cavityEvaluation=spin::evaluateCanonicalFullRankMaterialMoments(
                    out.measuredRequest,complement,evaluatorSettings);
            }
            out.materialExportComplete=true;
            out.materialExportDiagnostic=evaluatorSettings.triangulation.rationalSourceCoordinates
                ? "rational DESIGN-minus-SAME-retained, exact partition and independent moments; floating derivative NOT qualified; sealing still required"
                : "retained material and complement from the same quadric; full partition/sealing validation still required";
        } catch(const std::exception& error) {
            out.materialExportDiagnostic=error.what();
            // Keep the retained-material audit and its bound hash untouched.
            // The absent complement/export flag independently blocks sealing.
        }
    } else out.materialExportDiagnostic="not exported: no verified converged terminal material";
    return out;
}
} // namespace spin2014
