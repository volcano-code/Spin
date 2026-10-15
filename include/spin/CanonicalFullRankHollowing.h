#pragma once
#include "spin/FullRankHollowingSolver.h"
#include "spin/PolyhedralMomentDifferential.h"

namespace spin {
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh&, const TriangleMesh&, const std::vector<SourceVertexJacobian>&,
    const TriangleMesh&, const Eigen::VectorXi&, double, bool, bool, bool,
    const ExactBooleanVertexCoordinates&, const ExactBooleanVertexCoordinates&,
    const ExactBooleanVertexCoordinates&);
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh&, const TriangleMesh&, const std::vector<SourceVertexJacobian>&,
    const TriangleMesh&, const Eigen::VectorXi&, double, bool, bool, bool,
    const ExactBooleanVertexCoordinates&, const ExactBooleanVertexCoordinates&,
    const ExactBooleanVertexCoordinates&,
    const std::function<void(const std::string&)>& diagnosticObserver);
physics::RawMassMoments10 integrateCanonicalFullRankBoundaryMoments(
    const TriangleMesh&, double density);
physics::RawMassMoments10 integrateCanonicalFullRankBoundaryMoments(
    const TriangleMesh&, double density, bool exactBinary64SelfIntersection);
physics::RawMassMoments10 integrateCanonicalFullRankBoundaryMoments(
    const TriangleMesh&, double density, bool exactBinary64SelfIntersection,
    bool exactBinary64MaterialBoundary);
std::vector<std::string> integrateCanonicalFullRankBoundaryMoments(
    const ExactBooleanVertexCoordinates&, const Eigen::MatrixXi&, double density);

QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary, const ExactBooleanVertexCoordinates& exactVertices);
/** Order-only, orientation-preserving canonicalization of the material mesh.
 * No welding, dropping faces, rounding coordinates, tolerance changes or cache.
 * Returned vertexProvenance is indexed in canonical material-vertex order.
 * Same supported geometric/differential domain as the existing evaluator.
 */
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density);
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection);

QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary);
QuadricDiscreteMomentDifferential differentiateCanonicalPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& velocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces, double density,
    bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary);

FullRankMomentEvaluation evaluateCanonicalFullRankMaterialMoments(
    const SpinDesignRequest&, const QuadricCoefficients10&,
    const QuadricMomentEvaluatorSettings& = {});
SpinMomentEvaluator makeCanonicalFullRankMomentEvaluator(
    const SpinDesignRequest&, const QuadricMomentEvaluatorSettings& = {});
FullRankHollowingResult solveCanonicalFullRankQuadricHollowing(
    const SpinDesignRequest&, const QuadricMomentEvaluatorSettings& = {},
    const std::vector<QuadricCoefficients10>& = {}, const SpinNewtonEventSink& = {});
FullRankHollowingResult solveCanonicalFullRankQuadricHollowing(
    const SpinDesignRequest&, const MassMomentProblem&, const std::string& problemIdentity,
    const QuadricMomentEvaluatorSettings& = {},
    const std::vector<QuadricCoefficients10>& = {}, const SpinNewtonEventSink& = {});
}
