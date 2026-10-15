#include "spin/FullRankMomentEvaluator.h"
#include "manufacturing/SlicerValidation.h"
#include "physics/MassProperties.h"
#ifdef differentiatePolyhedralMaterialMoments
#include "spin/CanonicalFullRankHollowing.h"
#include "spin/QuadricSourceMotion.h"
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace spin {
bool sameOrientedTriangleBoundary(const TriangleMesh& first,
                                 const TriangleMesh& second)
{
    using Point = std::array<double, 3>;
    using Face = std::array<Point, 3>;
    const auto valid = [](const TriangleMesh& mesh) {
        if (mesh.vertices.cols() != 3 || mesh.faces.cols() != 3 ||
            !mesh.vertices.allFinite())
            return false;
        for (Eigen::Index face = 0; face < mesh.faces.rows(); ++face)
            for (int corner = 0; corner < 3; ++corner)
                if (mesh.faces(face, corner) < 0 ||
                    mesh.faces(face, corner) >= mesh.vertices.rows())
                    return false;
        return true;
    };
    if (!valid(first) || !valid(second) || first.faces.rows() != second.faces.rows())
        return false;

    const auto canonicalFaces = [](const TriangleMesh& mesh) {
        std::vector<Face> faces;
        faces.reserve(static_cast<std::size_t>(mesh.faces.rows()));
        for (Eigen::Index row = 0; row < mesh.faces.rows(); ++row)
        {
            Face face;
            for (int corner = 0; corner < 3; ++corner)
            {
                const Eigen::Index vertex = mesh.faces(row, corner);
                for (int axis = 0; axis < 3; ++axis)
                {
                    const double coordinate = mesh.vertices(vertex, axis);
                    face[corner][axis] = coordinate == 0.0 ? 0.0 : coordinate;
                }
            }
            const auto minimum = std::min_element(face.begin(), face.end());
            std::rotate(face.begin(), minimum, face.end());
            faces.push_back(face);
        }
        std::sort(faces.begin(), faces.end());
        return faces;
    };
    return canonicalFaces(first) == canonicalFaces(second);
}

FullRankMomentEvaluation evaluateFullRankMaterialMoments(
 const SpinDesignRequest& request,const QuadricCoefficients10& q,
 const QuadricMomentEvaluatorSettings& settings) {
    // Preserve existing value/audit semantics; never use legacy shape quadrature.
    auto valueSettings=settings;
    valueSettings.derivativeMethod=QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
    const auto native=evaluateQuadricMaterialMoments(request,q,valueSettings);
    FullRankMomentEvaluation out;
    out.variableMomentsTarget=native.variableMomentsTarget;
    out.totalMomentsTarget=native.totalMomentsTarget;
    out.materialBoundaryTarget=native.materialBoundaryTarget;
    out.materialExactVerticesTarget=native.materialExactVerticesTarget;
    out.designExactVerticesTarget=native.designExactVerticesTarget;
    out.sourceExactVerticesTarget=native.quadricTriangulation.exactVertices;
    out.sourceCoordinateLowerTarget=native.quadricTriangulation.coordinateLower;
    out.sourceCoordinateUpperTarget=native.quadricTriangulation.coordinateUpper;
    out.sourceConstructionIdentity=native.quadricTriangulation.rationalConstructionIdentity;
    out.sourceFaces=native.quadricTriangulation.closedBoundary.faces;
    out.materialBirthFaces=native.materialBirthFaces;
    out.variableMomentsExactTarget=native.variableMomentsExactTarget;
    if(settings.triangulation.rationalSourceCoordinates) {
        const auto& tri=native.quadricTriangulation;
        std::ostringstream branch;branch.imbue(std::locale::classic());
        branch<<"RationalProductionBranch-v1;"<<quadricSurfaceTypeToken(tri.type)<<';'
            <<quadricSurfaceSubtypeToken(classifyQuadricSurfaceSubtype(q))<<';'
            <<tri.initialIcosphereSubdivisions<<';'<<tri.effectiveIcosphereSubdivisions<<';'
            <<tri.boundaryLoopCount<<';'<<tri.capFaceCount<<';'<<std::hexfloat
            <<tri.truncationRadius<<';'<<native.quadricSolidIsPositiveRegion<<';'
            <<native.provenanceTopologyHash<<';';
        for(int f=0;f<tri.closedBoundary.faces.rows();++f) {
            for(int k=0;k<3;++k)branch<<tri.closedBoundary.faces(f,k)<<';';
            branch<<tri.faceIsCap(f)<<';';
        }
        for(const auto id:tri.sourceNormalVertexByVertex)branch<<id<<';';
        for(const auto& edge:tri.sourceNormalEdgeByVertex)branch<<edge[0]<<','<<edge[1]<<';';
        out.rationalBranchIdentity=manufacturing::computeContentSha256(branch.str());
    }
    out.geometryHash=native.geometryHash;
    out.booleanBackendIdentity=native.booleanBackendIdentity;
    out.effectiveDerivativeMethod=quadricMomentDerivativeMethodToken(native.derivativeMethod);
    // Keep the value-side topology identity available even when the
    // derivative is rejected.  A subsequent successful differential replaces
    // it with the stronger birth-feature/provenance identity below.
    out.topologyHash=native.topologyHash;
    const auto finish=[&]() {
        std::ostringstream id;id.imbue(std::locale::classic());
        id<<"FullRankMaterialDifferential-v3-strict-inactive-selector-margin-1e-10;"<<computeSpinDesignRequestHash(request,settings.exactBinary64SourceBoundary)<<';'
          <<q.hash()<<';'<<out.booleanBackendIdentity<<';'<<native.derivativeScientificHash<<';'<<out.geometryHash<<';'
          <<out.topologyHash<<';'
          <<out.effectiveDerivativeMethod<<';'
          <<out.sourceScientificHash<<';'<<out.diagnostic<<';'
          <<out.coefficientJacobian.has_value()<<';'<<std::hexfloat;
        if(out.coefficientJacobian)for(int i=0;i<10;++i)for(int j=0;j<10;++j)
            id<<(*out.coefficientJacobian)(i,j)<<';';
        out.scientificHash=manufacturing::computeContentSha256(id.str());
        return out;
    };
    // Replay the same robust selector used by the value map. A volume
    // centroid is only a valid sign witness when it is certified inside the
    // source solid; otherwise the selector uses oriented interface normals.
    // Replacing that fallback with an unconditional centroid sign can report
    // a false branch change on non-convex truncated quadric regions.
    const bool replayedPositiveRegion = quadricSolidInteriorIsPositiveRegion(
        q, native.quadricTriangulation);
    if(replayedPositiveRegion != native.quadricSolidIsPositiveRegion)
        throw std::runtime_error("native material-region selector changed on replay");
    out.topologyHash=!native.provenanceTopologyHash.empty()
        ? native.provenanceTopologyHash
        : native.topologyHash;
    // The source centroid is a candidate sign witness for the native selector.
    // Availability switches at abs(q(centroid))=threshold, NOT at q=0.
    // Recover only the strictly inactive branch with an exact complete-face
    // margin. This is a joint mesh/coefficient selector proof; the matching
    // regular source differential and Boolean motion below remain mandatory.
    // No parameter-neighborhood radius or Jacobian is granted by this alone.
    const auto& sourceBoundary = native.quadricTriangulation.closedBoundary;
    const MassMoments sourceMoments = computeSurfaceMassMoments(
        sourceBoundary.vertices, sourceBoundary.faces, 1.0);
    Eigen::Vector3d centroid;
    if (sourceMoments.tryCenterOfMass(centroid, 1e-15))
    {
        const double value = q.evaluate(centroid);
        const Eigen::Vector3d absolute = centroid.cwiseAbs();
        const double scale = std::abs(q.constant()) +
            q.linearVector().cwiseAbs().dot(absolute) +
            absolute.dot(q.quadraticMatrix().cwiseAbs() * absolute);
        if (std::isfinite(value) &&
            std::abs(value) <= 1e-10 * std::max(1.0, scale) &&
            !quadricSelectorHasStrictInactiveCentroidMargin(
                q,native.quadricTriangulation,replayedPositiveRegion))
        {
            out.diagnostic =
                "native material-region selector is near its sign transition";
            return finish();
        }
    }
#ifdef differentiatePolyhedralMaterialMoments
#ifdef MAKE_IT_STAND_IDENTICAL_TRANSPORT_OVERLAY
    // This opt-in is granted only by a build receipt proving that the native
    // and renamed transport objects compile the SAME source with identical
    // semantic flags. The ordinary coordinate-sorting canonical build must
    // never enable it. Exact-material native evaluation already ran that
    // transport on this same Boolean and source, including rejection checks.
    if (!native.materialExactVerticesTarget.empty() &&
        !settings.triangulation.rationalSourceCoordinates)
    {
        out.sourceScientificHash=native.derivativeScientificHash;
        out.diagnostic=native.derivativeDiagnostic;
        if(native.derivativeDefined)out.coefficientJacobian=native.momentJacobian;
        return finish();
    }
#endif
    if (settings.triangulation.rationalSourceCoordinates)
    {
        // The value evaluator already canonicalized the authority coordinates
        // and differentiated THAT same operation. A binary64-only canonical
        // replay cannot re-identify this rational operand.
        out.sourceScientificHash=native.derivativeScientificHash;
        out.diagnostic=native.derivativeDiagnostic;
        if(native.derivativeDefined)out.coefficientJacobian=native.momentJacobian;
        return finish();
    }
    // The optional canonical-order profile may reorder and re-sum the same
    // derivative vertices. It reuses the native value boundary and birth-face
    // tape; only source velocity is reconstructed, never the Boolean.
    auto design=request.designDomain;
    for(int i=0;i<design.vertices.rows();++i)
        design.vertices.row(i)=request.targetFrame.axes.transpose()*
            (design.vertices.row(i).transpose()-request.targetFrame.origin);
    auto sourceSettings=settings.triangulation;
    sourceSettings.retainFixedNormalProvenance=true;
    sourceSettings.truncationRadius=native.quadricTriangulation.truncationRadius;
    sourceSettings.maximumRadiusDoublings=0;
    if(settings.diagnosticStageObserver)settings.diagnosticStageObserver("canonical-source-motion-begin");
    const auto source=evaluateQuadricSourceMotion(q,design.vertices,sourceSettings);
    if(settings.diagnosticStageObserver)settings.diagnosticStageObserver("canonical-source-motion-end");
    out.sourceScientificHash=source.scientificHash;
    if(!source.defined){out.diagnostic=source.diagnostic;return finish();}
    if(source.triangulation.hash!=native.quadricTriangulation.hash)
        throw std::runtime_error("canonical source motion differs from value-side quadric geometry");
    if(settings.diagnosticStageObserver)settings.diagnosticStageObserver("canonical-material-transport-begin");
    const auto differential=differentiatePolyhedralMaterialMoments(
        design,source.triangulation.closedBoundary,source.velocity,
        native.materialBoundaryTarget,native.materialBirthFaces,request.materialDensity,
        settings.exactBinary64SelfIntersection,settings.exactBinary64SourceBoundary,
        settings.exactBinary64MaterialBoundary,native.materialExactVerticesTarget,
        {},{},settings.diagnosticStageObserver);
    if(settings.diagnosticStageObserver)settings.diagnosticStageObserver("canonical-material-transport-end");
    out.diagnostic=differential.diagnostic;
    if(differential.defined)out.coefficientJacobian=differential.jacobian;
#else
    out.sourceScientificHash=native.derivativeScientificHash;
    out.diagnostic=native.derivativeDiagnostic;
    if(native.derivativeDefined)out.coefficientJacobian=native.momentJacobian;
#endif
    return finish();
}
}
