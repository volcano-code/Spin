#include "spin/PlaneStaticBalanceSolver.h"

#include "spin/SpinDesignRequest.h"
#include "physics/MassProperties.h"

#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace spin {
namespace {

void validateSettings(const PlaneStaticBalanceSettings& settings) {
    if (settings.maximumIterations < 0 || settings.maximumLineSearchSteps <= 0 ||
        !std::isfinite(settings.constraintTolerance) || settings.constraintTolerance <= 0.0 ||
        !std::isfinite(settings.backtrackingFactor) || settings.backtrackingFactor <= 0.0 ||
        settings.backtrackingFactor >= 1.0 || !std::isfinite(settings.minimumStep) ||
        settings.minimumStep <= 0.0)
        throw std::invalid_argument("invalid static plane solver settings");
}

double fixedConstraintScale(const SpinDesignRequest& request) {
    const double infinity=std::numeric_limits<double>::infinity();
    Eigen::Vector3d lower=Eigen::Vector3d::Constant(infinity);
    Eigen::Vector3d upper=Eigen::Vector3d::Constant(-infinity);
    for(Eigen::Index i=0;i<request.outerSurface.vertices.rows();++i) {
        const Eigen::Vector3d world=request.outerSurface.vertices.row(i).transpose();
        const Eigen::Vector3d target=request.targetFrame.axes.transpose()*(world-request.targetFrame.origin);
        lower=lower.cwiseMin(target);
        upper=upper.cwiseMax(target);
    }
    const double length=(upper-lower).maxCoeff();
    const double scale=request.materialDensity*std::pow(length,4);
    if(!std::isfinite(length) || length<=0.0 || !std::isfinite(scale) || scale<=0.0)
        throw std::invalid_argument("static plane fixed-input first-moment scale is not representable");
    return scale;
}

struct EvaluatedPlane {
    Eigen::Vector2d lambda = Eigen::Vector2d::Zero();
    PlaneMomentEvaluation evaluation;
    Eigen::Vector2d constraints = Eigen::Vector2d::Zero();
    Eigen::Vector2d normalizedConstraints = Eigen::Vector2d::Zero();
    double residual = std::numeric_limits<double>::infinity();
    bool succeeded = false;
};

struct EvaluatedPlaneFamily {
    std::optional<PlaneStaticBalanceResult::FullDomainOptimality> boundaryOptimality;
    Eigen::Vector4d coefficients=Eigen::Vector4d::Zero();
    PlaneMomentEvaluation evaluation;
    Eigen::Vector3d residual=Eigen::Vector3d::Constant(
        std::numeric_limits<double>::infinity());
    Eigen::Matrix<double,3,4> jacobian=Eigen::Matrix<double,3,4>::Zero();
    Eigen::Vector2d constraints=Eigen::Vector2d::Zero();
    Eigen::Vector2d normalizedConstraints=Eigen::Vector2d::Zero();
    double stationarityResidual=std::numeric_limits<double>::infinity();
};

double targetLength(const SpinDesignRequest& request) {
    const double infinity=std::numeric_limits<double>::infinity();
    Eigen::Vector3d lower=Eigen::Vector3d::Constant(infinity);
    Eigen::Vector3d upper=Eigen::Vector3d::Constant(-infinity);
    for(Eigen::Index i=0;i<request.outerSurface.vertices.rows();++i) {
        const Eigen::Vector3d world=request.outerSurface.vertices.row(i).transpose();
        const Eigen::Vector3d target=request.targetFrame.axes.transpose()*(world-request.targetFrame.origin);
        lower=lower.cwiseMin(target);upper=upper.cwiseMax(target);
    }
    const double length=(upper-lower).maxCoeff();
    if(!std::isfinite(length)||length<=0.0)
        throw std::invalid_argument("plane-family characteristic length is not representable");
    return length;
}

EvaluatedPlaneFamily evaluateFamily(const SpinDesignRequest& request,
                                    const Eigen::Vector4d& coefficients,
                                    const PlaneStaticBalanceSettings& settings,
                                    double length,double constraintScale) {
    Eigen::Vector4d physical=coefficients;
    physical.tail<3>()/=length;
    EvaluatedPlaneFamily state;state.coefficients=coefficients;
    state.evaluation=evaluatePlaneMaterialMoments(request,physical,settings.planeEvaluation);
    const auto& moments=state.evaluation.totalMomentsTarget.values();
    state.constraints<<moments[1],moments[2];
    state.normalizedConstraints=state.constraints/constraintScale;
    const double mass=moments[0];
    const double heightMoment=moments[3];
    const double centerHeight=heightMoment/mass;
    if(!std::isfinite(mass)||mass<=0.0||!std::isfinite(heightMoment)||!std::isfinite(centerHeight))
        throw std::domain_error("plane-family KKT requires positive finite total mass");
    const double kappa=request.staticObjective==StaticObjective::LowerCenterOfMass
        ? 0.0 : request.staticParetoWeight;
    state.residual.head<2>()=state.normalizedConstraints;
    const double objectiveMassDerivative=kappa*staticObjectiveMassScale(request)+(1.0-kappa)*centerHeight/mass;
    const double objectiveHeightDerivative=-(1.0-kappa)/mass;
    state.residual[2]=coefficients[0]*objectiveHeightDerivative-
        (coefficients[3]/length)*objectiveMassDerivative;
    state.stationarityResidual=std::abs(state.residual[2]);
    // At a strictly full domain there is no moving interface. The interior
    // plane equation need not hold for the particular accepted plane. A
    // strictly positive negative-gradient polynomial defines the SAME full
    // region and satisfies Eq. 7 with positive multiplier. This uses the
    // explicitly selected mass or verified uniform-density volume objective.
    if(state.evaluation.materialRegion==PlaneMaterialRegion::StrictFull &&
       state.normalizedConstraints.norm()<=settings.constraintTolerance) {
        double minimum=std::numeric_limits<double>::infinity(), magnitude=0.0;
        for(Eigen::Index i=0;i<request.designDomain.vertices.rows();++i) {
            const Eigen::Vector3d x=request.targetFrame.axes.transpose()*(
                request.designDomain.vertices.row(i).transpose()-request.targetFrame.origin);
            const double term=objectiveHeightDerivative*x.z();
            minimum=std::min(minimum,objectiveMassDerivative+term);
            magnitude=std::max(magnitude,std::abs(objectiveMassDerivative)+std::abs(term));
        }
        Eigen::Vector4d witness(objectiveMassDerivative,0.0,0.0,length*objectiveHeightDerivative);
        if(witness.allFinite()&&std::isfinite(minimum)&&std::isfinite(magnitude)&&
           witness.norm()>0.0&&std::isfinite(witness.norm())&&
           minimum>64*std::numeric_limits<double>::epsilon()*magnitude) {
            PlaneStaticBalanceResult::FullDomainOptimality proof;
            proof.coefficients=witness.normalized();proof.minimumGradientPolynomial=minimum;
            state.boundaryOptimality=proof;
        }
    }
    if(!state.evaluation.coefficientJacobian)return state;

    Eigen::Matrix4d coefficientMap=Eigen::Matrix4d::Identity();
    coefficientMap.diagonal().tail<3>().setConstant(1.0/length);
    const Eigen::Matrix<double,10,4> momentJacobian=(*state.evaluation.coefficientJacobian)*coefficientMap;
    state.jacobian.row(0)=momentJacobian.row(1)/constraintScale;
    state.jacobian.row(1)=momentJacobian.row(2)/constraintScale;
    const Eigen::RowVector4d massDerivative=momentJacobian.row(0);
    const Eigen::RowVector4d heightMomentDerivative=momentJacobian.row(3);
    const Eigen::RowVector4d dObjectiveMassDerivative=(1.0-kappa)*(
        heightMomentDerivative/(mass*mass)-2.0*heightMoment*massDerivative/(mass*mass*mass));
    const Eigen::RowVector4d dObjectiveHeightDerivative=(1.0-kappa)*massDerivative/(mass*mass);
    state.jacobian.row(2)=Eigen::RowVector4d::Unit(0)*objectiveHeightDerivative+
        coefficients[0]*dObjectiveHeightDerivative-
        Eigen::RowVector4d::Unit(3)*(objectiveMassDerivative/length)-
        (coefficients[3]/length)*dObjectiveMassDerivative;
    return state;
}

void appendFamily(PlaneStaticBalanceResult& result,const EvaluatedPlaneFamily& current,
                  int iteration,double step,bool accepted,int evaluation) {
    PlaneStaticBalanceTraceRow row;row.evaluation=evaluation;row.iteration=iteration;
    row.coefficients=current.coefficients;row.constraints=current.constraints;
    row.normalizedConstraints=current.normalizedConstraints;
    row.normalizedStationarityResidual=current.stationarityResidual;
    row.residualNorm=current.residual.allFinite()?current.residual.norm():-1.0;
    row.variableMass=current.evaluation.variableMomentsTarget.mass();row.step=step;
    row.accepted=accepted;row.evaluationSucceeded=true;
    row.derivativeAvailable=current.evaluation.coefficientJacobian.has_value();
    row.diagnostic=current.evaluation.diagnostic;result.trace.push_back(std::move(row));
}

void retainFamily(PlaneStaticBalanceResult& result,const EvaluatedPlaneFamily& current) {
    result.coefficients=current.coefficients;result.constraints=current.constraints;
    result.normalizedConstraints=current.normalizedConstraints;
    result.stationarityResidual=current.stationarityResidual;
    result.variableMoments=current.evaluation.variableMomentsTarget;
    result.evaluation=current.evaluation;
    result.boundaryOptimality=current.boundaryOptimality;
}

EvaluatedPlane evaluate(const SpinDesignRequest& request,
                        const Eigen::Vector2d& lambda,
                        const PlaneStaticBalanceSettings& settings,
                        double constraintScale) {
    Eigen::Vector4d coefficients(1.0, lambda.x(), lambda.y(), 0.0);
    EvaluatedPlane result;
    result.lambda = lambda;
    result.evaluation = evaluatePlaneMaterialMoments(request, coefficients,
                                                       settings.planeEvaluation);
    const auto& total = result.evaluation.totalMomentsTarget.values();
    result.constraints << total[1], total[2];
    result.normalizedConstraints=result.constraints/constraintScale;
    result.residual = result.normalizedConstraints.norm();
    result.succeeded = true;
    return result;
}

Eigen::Matrix2d reducedJacobian(const PlaneMomentEvaluation& evaluation) {
    if (!evaluation.coefficientJacobian)
        throw std::runtime_error("plane-family differential unavailable");
    Eigen::Matrix2d result;
    result << (*evaluation.coefficientJacobian)(1,1), (*evaluation.coefficientJacobian)(1,2),
              (*evaluation.coefficientJacobian)(2,1), (*evaluation.coefficientJacobian)(2,2);
    return result;
}

void append(PlaneStaticBalanceResult& result, const EvaluatedPlane& current,
            int iteration, double step, bool accepted, int evaluation) {
    PlaneStaticBalanceTraceRow row;
    row.evaluation = evaluation;
    row.iteration = iteration;
    row.lambda = current.lambda;
    row.constraints = current.constraints;
    row.normalizedConstraints=current.normalizedConstraints;
    row.residualNorm = current.succeeded ? current.residual : -1.0;
    row.variableMass = current.succeeded ? current.evaluation.variableMomentsTarget.mass() : 0.0;
    row.step = step;
    row.accepted = accepted;
    row.evaluationSucceeded = current.succeeded;
    row.derivativeAvailable = current.succeeded && current.evaluation.coefficientJacobian.has_value();
    row.diagnostic = current.evaluation.diagnostic;
    result.trace.push_back(std::move(row));
}

void retain(PlaneStaticBalanceResult& result, const EvaluatedPlane& current) {
    result.lambda = current.lambda;
    result.constraints = current.constraints;
    result.normalizedConstraints=current.normalizedConstraints;
    result.variableMoments = current.evaluation.variableMomentsTarget;
    result.evaluation = current.evaluation;
}

} // namespace

const char* planeStaticBalanceStatusToken(PlaneStaticBalanceStatus status) noexcept {
    switch(status) {
    case PlaneStaticBalanceStatus::Converged:return "CONVERGED";
    case PlaneStaticBalanceStatus::IterationLimit:return "ITERATION_LIMIT";
    case PlaneStaticBalanceStatus::DerivativeUnavailable:return "DERIVATIVE_UNAVAILABLE";
    case PlaneStaticBalanceStatus::SingularReducedSystem:return "SINGULAR_REDUCED_SYSTEM";
    case PlaneStaticBalanceStatus::LineSearchExhausted:return "LINE_SEARCH_EXHAUSTED";
    case PlaneStaticBalanceStatus::GeometryEvaluationFailed:return "GEOMETRY_EVALUATION_FAILED";
    }
    return "UNKNOWN";
}

PlaneStaticBalanceResult solveStaticBalanceByPlane(
    const SpinDesignRequest& request, const Eigen::Vector2d& initialLambda,
    const PlaneStaticBalanceSettings& settings) {
    validateSpinDesignRequest(request);
    validateSettings(settings);
    if(request.mode!=SpinMode::StaticBalance ||
       request.staticObjective!=StaticObjective::MaximumMass)
        throw std::invalid_argument("plane solver implements static maximum-mass only");
    if(!initialLambda.allFinite())
        throw std::invalid_argument("initial plane multipliers must be finite");
    const double constraintScale=fixedConstraintScale(request);

    PlaneStaticBalanceResult result;
    result.constraintScale=constraintScale;
    EvaluatedPlane current;
    try {
        current=evaluate(request,initialLambda,settings,constraintScale);
    } catch(const std::exception& error) {
        result.status=PlaneStaticBalanceStatus::GeometryEvaluationFailed;
        result.diagnostic=error.what();
        result.lambda=initialLambda;
        EvaluatedPlane failed;failed.lambda=initialLambda;
        append(result,failed,0,0.0,false,0);
        result.trace.back().diagnostic=result.diagnostic;
        return result;
    }
    retain(result,current);
    append(result,current,0,0.0,false,0);
    if(current.residual<=settings.constraintTolerance) {
        result.status=PlaneStaticBalanceStatus::Converged;
        result.paperOptimalityConditionsVerified=true;
        result.diagnostic="vertical-plane threshold satisfies both static constraints";
        return result;
    }

    int evaluationCount=1;
    for(int iteration=0;iteration<settings.maximumIterations;++iteration) {
        Eigen::Matrix2d jacobian;
        try { jacobian=reducedJacobian(current.evaluation)/constraintScale; }
        catch(const std::exception& error) {
            result.status=PlaneStaticBalanceStatus::DerivativeUnavailable;
            result.diagnostic=error.what();
            return result;
        }
        if(!jacobian.allFinite()) {
            result.status=PlaneStaticBalanceStatus::DerivativeUnavailable;
            result.diagnostic="nonfinite reduced plane Jacobian";
            return result;
        }
        Eigen::FullPivLU<Eigen::Matrix2d> factor(jacobian);
        if(factor.rank()!=2 || factor.rcond()<1e-13) {
            result.status=PlaneStaticBalanceStatus::SingularReducedSystem;
            result.diagnostic="static two-constraint plane Jacobian is singular or ill-conditioned";
            return result;
        }
        const Eigen::Vector2d direction=factor.solve(-current.normalizedConstraints);
        if(!direction.allFinite() || direction.norm()==0.0) {
            result.status=PlaneStaticBalanceStatus::SingularReducedSystem;
            result.diagnostic="reduced Newton direction is zero or nonfinite";
            return result;
        }

        bool accepted=false;
        double step=1.0;
        for(int line=0;line<settings.maximumLineSearchSteps && step>=settings.minimumStep;++line) {
            const Eigen::Vector2d candidateLambda=current.lambda+step*direction;
            try {
                auto candidate=evaluate(request,candidateLambda,settings,constraintScale);
                ++evaluationCount;
                const bool improves=candidate.residual<current.residual;
                const bool reaches=candidate.residual<=settings.constraintTolerance;
                if(improves && (candidate.evaluation.coefficientJacobian || reaches)) {
                    current=std::move(candidate);
                    ++result.acceptedSteps;
                    result.iterations=iteration+1;
                    append(result,current,iteration+1,step,true,evaluationCount-1);
                    retain(result,current);
                    accepted=true;
                    if(reaches) {
                        result.status=PlaneStaticBalanceStatus::Converged;
                        result.paperOptimalityConditionsVerified=true;
                        result.diagnostic="vertical-plane threshold satisfies both static constraints";
                        return result;
                    }
                    break;
                }
                append(result,candidate,iteration+1,step,false,evaluationCount-1);
            } catch(const std::exception& error) {
                EvaluatedPlane failed;
                failed.lambda=candidateLambda;
                ++evaluationCount;
                result.diagnostic=error.what();
                append(result,failed,iteration+1,step,false,evaluationCount-1);
                result.trace.back().diagnostic=result.diagnostic;
            }
            step*=settings.backtrackingFactor;
        }
        if(!accepted) {
            result.status=PlaneStaticBalanceStatus::LineSearchExhausted;
            if(result.diagnostic.empty())
                result.diagnostic="backtracking could not reduce static constraint residual";
            return result;
        }
    }
    result.status=PlaneStaticBalanceStatus::IterationLimit;
    result.diagnostic="static plane Newton iteration limit reached; no infeasibility claim";
    return result;
}

PlaneStaticBalanceResult solveStaticBalanceByPlaneFamily(
    const SpinDesignRequest& request,const Eigen::Vector4d& initialCoefficients,
    const PlaneStaticBalanceSettings& settings) {
    validateSpinDesignRequest(request);validateSettings(settings);
    if(request.mode!=SpinMode::StaticBalance)
        throw std::invalid_argument("plane-family solver is restricted to static balance");
    if(!initialCoefficients.allFinite()||initialCoefficients.norm()<=0.0)
        throw std::invalid_argument("plane-family initial coefficients must be finite and nonzero");
    const double length=targetLength(request);
    const double constraintScale=fixedConstraintScale(request);

    const bool maximumMassObjective =
        request.staticObjective == StaticObjective::MaximumMass ||
        (request.staticObjective == StaticObjective::ParetoTradeoff &&
         request.staticParetoWeight == 1.0);
    const double coefficientScale = initialCoefficients.norm();
    const Eigen::Vector4d normalizedInitial =
        initialCoefficients / coefficientScale;
    if (maximumMassObjective && normalizedInitial[0] > 0.0 &&
        normalizedInitial.tail<3>().norm() <= 1e-14) {
        PlaneStaticBalanceResult fullDomain;
        fullDomain.constraintScale = constraintScale;
        fullDomain.coefficients = normalizedInitial;
        const Eigen::Matrix3d modelToTarget = request.targetFrame.axes.transpose();
        const Eigen::Vector3d targetTranslation =
            -modelToTarget * request.targetFrame.origin;
        TriangleMesh boundaryTarget = request.designDomain;
        for (Eigen::Index vertex = 0; vertex < boundaryTarget.vertices.rows(); ++vertex) {
            const Eigen::Vector3d world =
                request.designDomain.vertices.row(vertex).transpose();
            boundaryTarget.vertices.row(vertex) =
                (modelToTarget * world + targetTranslation).transpose();
        }
        const MassMoments variable = computeSurfaceMassMoments(
            boundaryTarget.vertices, boundaryTarget.faces, request.materialDensity);
        if (!variable.isFinite() || variable.mass <= 0.0)
            throw std::runtime_error(
                "constant positive plane produced invalid full-domain moments");
        const physics::RawMassMoments10 variableRaw =
            physics::RawMassMoments10::fromMassMoments(variable);
        const physics::RawMassMoments10 fixedTarget =
            request.fixedShellMoments.transformed(modelToTarget, targetTranslation);
        const physics::RawMassMoments10 total = variableRaw + fixedTarget;
        const auto& totalValues = total.values();
        fullDomain.constraints << totalValues[1], totalValues[2];
        fullDomain.normalizedConstraints =
            fullDomain.constraints / constraintScale;
        fullDomain.variableMoments = variableRaw;
        fullDomain.evaluation.variableMaterialBoundaryTarget =
            std::move(boundaryTarget);
        fullDomain.evaluation.variableMomentsTarget = variableRaw;
        fullDomain.evaluation.totalMomentsTarget = total;
        fullDomain.evaluation.normalizedCoefficients = normalizedInitial;
        fullDomain.diagnostic =
            "constant positive plane retains the full design domain; its mass "
            "is the attainable upper bound";
        if (fullDomain.normalizedConstraints.norm() <=
            settings.constraintTolerance) {
            fullDomain.status = PlaneStaticBalanceStatus::Converged;
            fullDomain.paperOptimalityConditionsVerified = true;
            fullDomain.evaluation.diagnostic = fullDomain.diagnostic;
        } else {
            fullDomain.status = PlaneStaticBalanceStatus::DerivativeUnavailable;
            fullDomain.diagnostic +=
                ", but the full-domain horizontal balance constraints are not satisfied "
                "and a constant plane has no cut differential";
            fullDomain.evaluation.diagnostic = fullDomain.diagnostic;
        }
        PlaneStaticBalanceTraceRow row;
        row.coefficients = normalizedInitial;
        row.constraints = fullDomain.constraints;
        row.normalizedConstraints = fullDomain.normalizedConstraints;
        row.residualNorm = fullDomain.normalizedConstraints.norm();
        row.variableMass = variable.mass;
        row.evaluationSucceeded = true;
        row.diagnostic = fullDomain.diagnostic;
        fullDomain.trace.push_back(std::move(row));
        return fullDomain;
    }

    if(request.staticObjective==StaticObjective::MaximumMass ||
       (request.staticObjective==StaticObjective::ParetoTradeoff&&request.staticParetoWeight==1.0)) {
        if(initialCoefficients[0]<=0.0||std::abs(initialCoefficients[3])>1e-12)
            throw std::invalid_argument("maximum-mass plane initialization must have positive intercept and vertical interface");
        const Eigen::Vector2d lambda(initialCoefficients[1]/(initialCoefficients[0]*length),
                                     initialCoefficients[2]/(initialCoefficients[0]*length));
        auto maximumMassRequest=request;
        maximumMassRequest.staticObjective=StaticObjective::MaximumMass;
        maximumMassRequest.staticParetoWeight=0.0;
        auto result=solveStaticBalanceByPlane(maximumMassRequest,lambda,settings);
        result.coefficients<<1.0,result.lambda.x()*length,result.lambda.y()*length,0.0;
        result.coefficients.normalize();
        return result;
    }

    Eigen::Vector4d currentCoefficients=initialCoefficients.normalized();
    if(currentCoefficients[3]>=-1e-12)
        throw std::invalid_argument("lowering/Pareto plane initialization must select the below-plane halfspace");

    PlaneStaticBalanceResult result;result.constraintScale=constraintScale;
    EvaluatedPlaneFamily current;
    try {current=evaluateFamily(request,currentCoefficients,settings,length,constraintScale);}
    catch(const std::exception& error) {
        result.status=PlaneStaticBalanceStatus::GeometryEvaluationFailed;
        result.diagnostic=error.what();result.coefficients=currentCoefficients;
        PlaneStaticBalanceTraceRow failed;
        failed.coefficients=currentCoefficients;failed.residualNorm=-1.0;
        failed.diagnostic=result.diagnostic;
        result.trace.push_back(std::move(failed));return result;
    }
    retainFamily(result,current);appendFamily(result,current,0,0.0,false,0);
    const double tolerance=settings.constraintTolerance;
    auto stationary=[&](const EvaluatedPlaneFamily& value) {
        if(value.evaluation.materialRegion==PlaneMaterialRegion::StrictFull)
            return value.boundaryOptimality.has_value();
        return value.normalizedConstraints.norm()<=tolerance&&
            (value.boundaryOptimality.has_value() ||
             (value.stationarityResidual<=tolerance&&value.coefficients[3]<-1e-12));
    };
    if(stationary(current)) {
        result.status=PlaneStaticBalanceStatus::Converged;
        result.paperOptimalityConditionsVerified=true;
        result.diagnostic=current.boundaryOptimality?
            "strict full domain satisfies numerical negative-gradient sign optimality":
            "static plane satisfies horizontal balance and the paper objective KKT relation";
        return result;
    }

    int evaluationCount=1;
    for(int iteration=0;iteration<settings.maximumIterations;++iteration) {
        if(!current.evaluation.coefficientJacobian) {
            result.status=PlaneStaticBalanceStatus::DerivativeUnavailable;
            result.diagnostic="static plane-family differential unavailable at current geometry";
            return result;
        }
        Eigen::Matrix4d bordered;bordered.topRows<3>()=current.jacobian;
        bordered.row(3)=current.coefficients.transpose();
        Eigen::FullPivLU<Eigen::Matrix4d> factor(bordered);
        if(factor.rank()!=4||factor.rcond()<1e-13) {
            result.status=PlaneStaticBalanceStatus::SingularReducedSystem;
            result.diagnostic="three KKT equations on the normalized plane sphere are singular";
            return result;
        }
        Eigen::Vector4d rhs;rhs.head<3>()=-current.residual;rhs[3]=0.0;
        const Eigen::Vector4d direction=factor.solve(rhs);
        if(!direction.allFinite()||direction.norm()==0.0) {
            result.status=PlaneStaticBalanceStatus::SingularReducedSystem;
            result.diagnostic="plane-family tangent Newton direction is zero or nonfinite";
            return result;
        }

        bool accepted=false;double step=1.0;
        for(int line=0;line<settings.maximumLineSearchSteps&&step>=settings.minimumStep;++line) {
            Eigen::Vector4d candidateCoefficients=current.coefficients+step*direction;
            if(!candidateCoefficients.allFinite()||candidateCoefficients.norm()<=0.0) {
                step*=settings.backtrackingFactor;continue;
            }
            candidateCoefficients.normalize();
            if(candidateCoefficients[3]>=-1e-12) {
                PlaneStaticBalanceTraceRow row;row.evaluation=evaluationCount++;row.iteration=iteration+1;
                row.coefficients=candidateCoefficients;row.step=step;row.evaluationSucceeded=false;
                row.diagnostic="trial crossed the paper's below-plane orientation branch";
                result.trace.push_back(std::move(row));step*=settings.backtrackingFactor;continue;
            }
            try {
                auto candidate=evaluateFamily(request,candidateCoefficients,settings,length,constraintScale);
                const bool improves=candidate.residual.norm()<current.residual.norm();
                const bool reaches=stationary(candidate);
                const int thisEvaluation=evaluationCount++;
                if(improves&&(candidate.evaluation.coefficientJacobian||reaches)) {
                    current=std::move(candidate);++result.acceptedSteps;result.iterations=iteration+1;
                    appendFamily(result,current,iteration+1,step,true,thisEvaluation);
                    retainFamily(result,current);accepted=true;
                    if(reaches) {
                        result.status=PlaneStaticBalanceStatus::Converged;
                        result.paperOptimalityConditionsVerified=true;
                        result.diagnostic=current.boundaryOptimality?
                            "strict full domain satisfies numerical negative-gradient sign optimality":
                            "static plane satisfies horizontal balance and the paper objective KKT relation";
                        return result;
                    }
                    break;
                }
                appendFamily(result,candidate,iteration+1,step,false,thisEvaluation);
            } catch(const std::exception& error) {
                ++evaluationCount;PlaneStaticBalanceTraceRow row;row.evaluation=evaluationCount-1;
                row.iteration=iteration+1;row.coefficients=candidateCoefficients;row.step=step;
                row.evaluationSucceeded=false;row.diagnostic=error.what();result.trace.push_back(std::move(row));
                result.diagnostic=error.what();
            }
            step*=settings.backtrackingFactor;
        }
        if(!accepted) {
            result.status=PlaneStaticBalanceStatus::LineSearchExhausted;
            if(result.diagnostic.empty())result.diagnostic="backtracking could not reduce the static plane-family KKT residual";
            return result;
        }
    }
    result.status=PlaneStaticBalanceStatus::IterationLimit;
    result.diagnostic="static plane-family iteration limit reached; no infeasibility claim";
    return result;
}

} // namespace spin
