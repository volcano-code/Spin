#pragma once
#include "spin/QuadricDiscreteMomentDifferential.h"
#include <array>
#include <functional>

namespace verification { class ExactMaterialRegion; }

namespace spin {
using SourceVertexJacobian = Eigen::Matrix<double,3,10>;
using ExactBooleanVertexCoordinates = std::vector<std::array<std::string,3>>;

using RationalSourceVertexJacobian=std::array<std::array<std::string,10>,3>;
/** Versioned rational transport attachment for the SAME audited Boolean.
 * Exact source velocities, intersection transport and moment flux remain
 * rational; the summary Jacobian is a compatibility view only. The caller
 * binds source velocities to its coefficient/branch identity. This does not
 * certify differentiability across topology changes or manufacture. */
struct RationalPolyhedralMomentDifferential {
    bool defined=false;
    std::string diagnostic;
    QuadricDiscreteMomentDifferential summary;
    std::array<std::array<std::string,10>,10> coefficientJacobian;
    std::vector<RationalSourceVertexJacobian> materialVertexJacobian;
    std::string maximumConstraintResidualExact;
};
RationalPolyhedralMomentDifferential differentiateRationalPolyhedralMaterialMoments(
    const TriangleMesh& design,const TriangleMesh& source,
    const std::vector<RationalSourceVertexJacobian>& sourceVelocity,
    const TriangleMesh& material,const Eigen::VectorXi& birthFaces,double density,
    const ExactBooleanVertexCoordinates& exactMaterial,
    const ExactBooleanVertexCoordinates& exactSource,
    const ExactBooleanVertexCoordinates& exactDesign,
    const std::function<void(const std::string&)>& diagnosticObserver = {},
    const verification::ExactMaterialRegion* materialAuthority = nullptr);

/** Rational-authority variant: all three tapes belong to the SAME Boolean.
 * Display vertices must be their nearest assignments. Source/design audits,
 * facet membership and material integration use the tapes, not the displays.
 */
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh&, const TriangleMesh&,
    const std::vector<SourceVertexJacobian>&, const TriangleMesh&, const Eigen::VectorXi&,
    double, bool, bool, bool, const ExactBooleanVertexCoordinates&,
    const ExactBooleanVertexCoordinates&, const ExactBooleanVertexCoordinates&,
    const std::function<void(const std::string&)>& diagnosticObserver = {});

/** Same calculation with an already completed material-region proof. The
 * proof must match every exact coordinate and indexed face. Source/design,
 * display assignment, provenance and conditioning checks remain unchanged. */
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh&, const TriangleMesh&,
    const std::vector<SourceVertexJacobian>&, const TriangleMesh&, const Eigen::VectorXi&,
    double, bool, bool, bool, const ExactBooleanVertexCoordinates&,
    const ExactBooleanVertexCoordinates&, const ExactBooleanVertexCoordinates&,
    const std::function<void(const std::string&)>&,
    const verification::ExactMaterialRegion*);

/** Fixed-topology material moment differential in an explicit reduced
 * parameter chart. Jacobian rows use the standard ten raw moments; columns
 * correspond exactly to the supplied source-vertex velocity columns.
 */
struct RestrictedPolyhedralMomentDifferential {
    bool defined=false;
    std::string diagnostic;
    Eigen::MatrixXd jacobian;
    std::vector<std::string> vertexProvenance;
    int constructedIntersectionVertexCount=0;
    double maximumConstraintResidual=0;
    double minimumConstraintSingularRatio=1;
};
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh&, const TriangleMesh&,
    const std::vector<Eigen::MatrixXd>&, const TriangleMesh&, const Eigen::VectorXi&,
    double, bool, bool, bool, const ExactBooleanVertexCoordinates&,
    const ExactBooleanVertexCoordinates&, const ExactBooleanVertexCoordinates&,
    const std::function<void(const std::string&)>& diagnosticObserver = {});
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh&, const TriangleMesh&,
    const std::vector<Eigen::MatrixXd>&, const TriangleMesh&, const Eigen::VectorXi&,
    double, bool, bool, bool, const ExactBooleanVertexCoordinates&,
    const ExactBooleanVertexCoordinates&, const ExactBooleanVertexCoordinates&,
    const std::function<void(const std::string&)>&,
    const verification::ExactMaterialRegion*);

/** Local transport of an already audited Boolean result in common coordinates.
 * birthFaces refer to [design faces, source faces].
 * Rows are (mass, lx, ly, lz, Qxx, Qyy, Qzz, Qyz, Qxz, Qxy), with fixed uniform
 * density included. The caller owns the operation and coefficient semantics
 * of sourceVelocity; no normalization
 * is added here. This low-level result is not a geometry/manufacturing seal.
 * Empty material, contact and ill-conditioned source correspondences fail closed.
 */
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& sourceVelocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces,
    double density);
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& sourceVelocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces,
    double density, bool exactBinary64SelfIntersection);

/** Explicit source and design exact boundary proof; material keeps the
 * existing audit and numerical motion/feature conditioning remains required. */
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& sourceVelocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces,
    double density, bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary);
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& sourceVelocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces,
    double density, bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary);

/** Coordinates from the SAME Boolean before binary64 assignment. Every
 * coordinate must round to the supplied material vertex and belong exactly
 * to all its birth facets. Used only for motion-constraint assembly; no
 * geometry, parameter dimension, audit or numerical tolerance is changed.
 * An empty tape preserves legacy assembly. Nonempty invalid tapes fail closed. */
QuadricDiscreteMomentDifferential differentiatePolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<SourceVertexJacobian>& sourceVelocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces,
    double density, bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary, const ExactBooleanVertexCoordinates& exactVertices);

/** Reduced-coordinate counterpart of the fixed ten-dimensional differential.
 * Each source velocity must be 3 x d with the same positive d. The resulting
 * ten raw material moments are differentiated with respect to those d
 * coordinates. No source parameterization or normalization is inferred.
 */
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& sourceVelocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces,
    double density);
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& sourceVelocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces,
    double density, bool exactBinary64SelfIntersection);
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& sourceVelocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces,
    double density, bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary);
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& sourceVelocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces,
    double density, bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary);
RestrictedPolyhedralMomentDifferential differentiateRestrictedPolyhedralMaterialMoments(
    const TriangleMesh& design, const TriangleMesh& source,
    const std::vector<Eigen::MatrixXd>& sourceVelocity,
    const TriangleMesh& material, const Eigen::VectorXi& birthFaces,
    double density, bool exactBinary64SelfIntersection, bool exactBinary64SourceBoundary,
    bool exactBinary64MaterialBoundary, const ExactBooleanVertexCoordinates& exactVertices);
}
