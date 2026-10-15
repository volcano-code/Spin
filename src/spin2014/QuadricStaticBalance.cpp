#include "spin2014/QuadricStaticBalance.h"
#include "manufacturing/SlicerValidation.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace spin2014 {
StaticSquaredHeightVerification verifySpin2014StaticPlane(
    const Eigen::Matrix<double,10,1>& moments,const Eigen::Vector4d& plane,
    double length,const QuadricObjectiveSettings& objective,double constraintScale,
    double tolerance)
{
    validateSpin2014QuadricObjectiveSettings(objective);
    if(objective.objective!=Spin2014Objective::StaticLowering || !moments.allFinite() ||
       moments[0]<=0 || moments[3]<0 || !plane.allFinite() ||
       std::abs(plane.norm()-1)>1e-12 || !std::isfinite(length) || length<=0 ||
       !std::isfinite(constraintScale) || constraintScale<=0 ||
       !std::isfinite(tolerance) || tolerance<=0 || tolerance>1e-10)
        throw std::invalid_argument("static squared-height requires positive finite mass, nonnegative height, normalized plane and strict fixed scales");
    StaticSquaredHeightVerification out;
    out.evaluation=makeSpin2014QuadricMomentProblem(objective).evaluate(moments);
    // Use the original objective's transformed gradient, not an independently
    // retyped formula or a scaled old plane solver acceptance boolean.
    const auto& gradient=out.evaluation.objectiveGradient;
    if(objective.centerWeight==0 || moments[3]==0) {
        out.degenerateObjective=true;
        out.diagnostic="ZERO_WEIGHT_OR_ZERO_HEIGHT: zero objective sensitivity; ordinary positive-multiplier condition is not applicable";
        return out;
    }
    if(plane[3]>=0) {
        out.diagnostic="squared-height candidate does not retain the below-plane halfspace";
        return out;
    }
    Eigen::Vector4d physical=plane;physical.tail<3>()/=length;
    out.multiplier=-gradient[3]/physical[3];
    // Static equalities are the declared first moments. Their derivative is
    // taken from the same 2014 objective adapter, including unit conversion.
    for(int axis=0;axis<2;++axis) {
        const double equalityDerivative=out.evaluation.constraintJacobian(axis,axis+1);
        out.equalityMultipliers[axis]=(-out.multiplier*physical[axis+1]-gradient[axis+1])/
                                      equalityDerivative;
    }
    Eigen::Vector4d polynomialResidual=gradient.head<4>();
    polynomialResidual+=out.evaluation.constraintJacobian.leftCols<4>().transpose()*out.equalityMultipliers;
    polynomialResidual+=out.multiplier*physical;
    out.stationarityResidual=polynomialResidual.lpNorm<Eigen::Infinity>();
    out.squaredHeightConditionsVerified=std::isfinite(out.multiplier) && out.multiplier>0 &&
        out.equalityMultipliers.allFinite() && std::isfinite(out.stationarityResidual) &&
        out.stationarityResidual<=tolerance &&
        out.evaluation.constraints.lpNorm<Eigen::Infinity>()/constraintScale<=tolerance;
    out.diagnostic=out.squaredHeightConditionsVerified?
        "fresh 2014 squared-height derivatives, equalities and positive plane multiplier numerically verified":
        "fresh squared-height stationarity or horizontal constraints failed; no infeasibility claim";
    return out;
}

QuadricStaticMaterialRun solveSpin2014StaticQuadricHollowing(
    const spin::SpinDesignRequest& request,const QuadricObjectiveSettings& objective,
    double fixedDensity,const std::optional<Eigen::Vector4d>& initial,
    const spin::PlaneStaticBalanceSettings& settings)
{
    if(objective.objective!=Spin2014Objective::StaticLowering || settings.constraintTolerance>1e-10)
        throw std::invalid_argument("static plane adapter requires the original static squared-height objective and strict tolerance");
    QuadricStaticMaterialRun out;
    out.objectiveSettings=objective;out.objectiveProfile=spin2014QuadricObjectiveProfile();
    out.objectiveIdentity=spin2014QuadricObjectiveHash(objective);
    out.measuredRequest=prepareSpin2014QuadricRequest(request,objective,fixedDensity);
    out.fixedMaterialDensity=fixedDensity;
    Eigen::Vector3d lo=Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());
    Eigen::Vector3d hi=-lo;
    for(Eigen::Index i=0;i<request.outerSurface.vertices.rows();++i) {
        const Eigen::Vector3d x=request.targetFrame.axes.transpose()*(
            request.outerSurface.vertices.row(i).transpose()-request.targetFrame.origin);
        if(x.z()<0)throw std::invalid_argument("outer material is below its declared support plane");
        lo=lo.cwiseMin(x);hi=hi.cwiseMax(x);
    }
    out.targetLength=(hi-lo).maxCoeff();
    if(initial) {
        out.initialCoefficients=*initial;
        out.initialization="EXPLICIT_DIAGNOSTIC_PLANE";
    } else {
        const auto fullModel=out.measuredRequest.fixedShellMoments+
            physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(
                request.designDomain.vertices,request.designDomain.faces,request.materialDensity));
        const auto target=fullModel.transformed(request.targetFrame.axes.transpose(),
            -request.targetFrame.axes.transpose()*request.targetFrame.origin);
        if(!target.isFinite() || target.mass()<=0 || target.values()[3]<0)
            throw std::invalid_argument("static input full material requires finite positive mass and nonnegative height");
        out.initialCoefficients=Eigen::Vector4d(
            target.values()[3]/target.mass()/out.targetLength,0,0,-1).normalized();
        out.initialization="INPUT_FULL_MATERIAL_HEIGHT_V1";
    }
    auto candidateRequest=out.measuredRequest;
    candidateRequest.staticObjective=spin::StaticObjective::LowerCenterOfMass;
    candidateRequest.staticParetoWeight=0;
    candidateRequest.staticObjectiveMeasure=spin::StaticObjectiveMeasure::DensityWeightedMass;
    out.candidate=spin::solveStaticBalanceByPlaneFamily(candidateRequest,out.initialCoefficients,settings);
    const auto& moments=out.candidate.evaluation.totalMomentsTarget.values();
    if(moments.allFinite() && moments[0]>0) {
        const double declaredConstraintScale=out.candidate.constraintScale*
            objective.solverMassUnitInDeclaredUnits*objective.solverLengthUnitInDeclaredUnits;
        out.squaredHeight=verifySpin2014StaticPlane(moments,out.candidate.coefficients,
            out.targetLength,objective,declaredConstraintScale,settings.constraintTolerance);
    } else out.squaredHeight.diagnostic="candidate has no finite positive terminal material; squared-height verification unavailable";
    // candidate.paperOptimalityConditionsVerified and boundaryOptimality are
    // intentionally not read. A strict full region needs its own fresh
    // squared-objective negative-sensitivity sign witness, not the interface
    // equation for an arbitrary accepted plane which does not cut the region.
    if(out.candidate.evaluation.materialRegion==spin::PlaneMaterialRegion::StrictFull) {
        out.squaredHeight.squaredHeightConditionsVerified=false;
        auto& checked=out.squaredHeight;
        if(!checked.degenerateObjective && checked.evaluation.constraints.size()==2 &&
           checked.evaluation.constraints.lpNorm<Eigen::Infinity>()/
               (out.candidate.constraintScale*objective.solverMassUnitInDeclaredUnits*
                objective.solverLengthUnitInDeclaredUnits)<=settings.constraintTolerance) {
            const Eigen::Vector4d negative=-checked.evaluation.objectiveGradient.head<4>();
            double minimum=std::numeric_limits<double>::infinity(),magnitude=0;
            for(Eigen::Index i=0;i<request.designDomain.vertices.rows();++i) {
                const Eigen::Vector3d x=request.targetFrame.axes.transpose()*(
                    request.designDomain.vertices.row(i).transpose()-request.targetFrame.origin);
                minimum=std::min(minimum,negative[0]+negative.tail<3>().dot(x));
                magnitude=std::max(magnitude,std::abs(negative[0])+
                    negative.tail<3>().cwiseAbs().dot(x.cwiseAbs()));
            }
            Eigen::Vector4d witness=negative;witness.tail<3>()*=out.targetLength;
            const double multiplier=witness.norm();
            if(negative.allFinite() && std::isfinite(minimum) && std::isfinite(magnitude) &&
               std::isfinite(multiplier) && multiplier>0 &&
               minimum>64*std::numeric_limits<double>::epsilon()*magnitude) {
                checked.boundaryWitnessVerified=true;
                checked.boundaryWitnessCoefficients=witness/multiplier;
                checked.boundaryMultiplier=multiplier;
                checked.boundaryGradientMinimum=minimum;
                checked.squaredHeightConditionsVerified=true;
            }
        }
        checked.diagnostic=checked.boundaryWitnessVerified?
            "fresh squared-height negative-sensitivity sign witness with zero equality multipliers numerically verified on the strict full design domain":
            "strict full domain failed fresh squared-height sign witness; no inherited height-objective pass or infeasibility claim";
    }
    return out;
}
} // namespace spin2014
