#pragma once
#include "spin2014/SourceCageBinding.h"
#include "spin2014/LinearizedRotationLaplacian.h"

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCageShapeFunctionalConfig {
    LinearizedRotationLaplacianConfig operatorConfig;
    double shapeWeight=0.,displacementWeight=0.;
    std::size_t maximumHandleGradientValues=1000000;
    std::size_t maximumWeightValues=1000000;
};
struct SourceCageShapeEvaluation : LinearizedRotationLaplacianEvaluation {
    Eigen::VectorXd handleGradient; // All original handles [scale,tx,ty,tz].
};
/** Immutable actual-binding BBW derivative chain on original outer rows only.
 * Input rows are supplied from the actual moment evaluation by the paper
 * functional. This mathematical action does not certify supplied deformations.
 * No supplied weights constructor; uses every cached original rest embedding
 * and the immutable actual final BBW field without clamping/renormalization.
 * Vertex arrays move from operator result; only H row + 4H gradient additional
 * doubles, separately capped. Value-only requires neither row nor gradient.
 */
class SourceCageShapeFunctional {
public:
    SourceCageShapeFunctional(const SourceCageBinding&,const SourceCageShapeFunctionalConfig&);
    static void validateConfig(const SourceCageShapeFunctionalConfig&);
    const SourceCageBinding& binding() const;
    const SourceCageShapeFunctionalConfig& config() const;
    const LinearizedRotationLaplacian& laplacian() const;
    const std::string& configSha256() const;
    const std::string& upstreamSha256() const;
    SourceCageShapeEvaluation evaluate(Eigen::Ref<const Eigen::MatrixXd> deformedOuter) const;
    SourceCageShapeEvaluation evaluateValues(Eigen::Ref<const Eigen::MatrixXd> deformedOuter) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    SourceCageShapeEvaluation run(Eigen::Ref<const Eigen::MatrixXd>,bool) const;
};
}
#endif
