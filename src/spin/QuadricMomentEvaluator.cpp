#include "spin/QuadricMomentEvaluator.h"
#include "spin/QuadricDiscreteMomentDifferential.h"

#include "manufacturing/ManufacturingGeometry.h"
#include "manufacturing/SlicerValidation.h"
#include "physics/MassProperties.h"

#include <igl/MeshBooleanType.h>
#include <igl/copyleft/cgal/mesh_boolean.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <chrono>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>


namespace
{


spin::TriangleMesh toTarget(
    const spin::TriangleMesh& model,
    const spin::SpinTargetFrame& frame
)
{
    spin::TriangleMesh target = model;
    for (Eigen::Index row = 0; row < model.vertices.rows(); ++row)
    {
        target.vertices.row(row) =
            frame.axes.transpose() *
            (model.vertices.row(row).transpose() - frame.origin);
    }
    return target;
}


spin::TriangleMesh toModel(
    const spin::TriangleMesh& target,
    const spin::SpinTargetFrame& frame
)
{
    spin::TriangleMesh model = target;
    for (Eigen::Index row = 0; row < target.vertices.rows(); ++row)
    {
        model.vertices.row(row) =
            frame.axes * target.vertices.row(row).transpose() + frame.origin;
    }
    return model;
}


std::string meshHash(const spin::TriangleMesh& mesh, const char* version)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << version << ';' << mesh.vertices.rows() << ';'
           << mesh.faces.rows() << ';' << std::hexfloat;
    for (Eigen::Index row = 0; row < mesh.vertices.rows(); ++row)
    {
        stream << mesh.vertices(row, 0) << ';'
               << mesh.vertices(row, 1) << ';'
               << mesh.vertices(row, 2) << ';';
    }
    for (Eigen::Index row = 0; row < mesh.faces.rows(); ++row)
    {
        stream << mesh.faces(row, 0) << ';'
               << mesh.faces(row, 1) << ';'
               << mesh.faces(row, 2) << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}


void canonicalizeBooleanMesh(
    spin::TriangleMesh* mesh,
    Eigen::VectorXi* birthFaces
)
{
    for (Eigen::Index row = 0; row < mesh->vertices.rows(); ++row)
    {
        for (Eigen::Index column = 0; column < mesh->vertices.cols(); ++column)
        {
            if (mesh->vertices(row, column) == 0.0)
            {
                mesh->vertices(row, column) = 0.0;
            }
        }
    }
    std::vector<int> vertexOrder(
        static_cast<std::size_t>(mesh->vertices.rows())
    );
    for (int index = 0; index < mesh->vertices.rows(); ++index)
    {
        vertexOrder[static_cast<std::size_t>(index)] = index;
    }
    std::sort(
        vertexOrder.begin(),
        vertexOrder.end(),
        [&](int left, int right) {
            for (int coordinate = 0; coordinate < 3; ++coordinate)
            {
                if (mesh->vertices(left, coordinate) <
                    mesh->vertices(right, coordinate))
                {
                    return true;
                }
                if (mesh->vertices(right, coordinate) <
                    mesh->vertices(left, coordinate))
                {
                    return false;
                }
            }
            return left < right;
        }
    );
    std::vector<Eigen::Vector3d> uniqueVertices;
    uniqueVertices.reserve(vertexOrder.size());
    std::vector<int> remap(vertexOrder.size());
    for (const int input : vertexOrder)
    {
        const Eigen::Vector3d point = mesh->vertices.row(input);
        if (
            uniqueVertices.empty() ||
            point != uniqueVertices.back()
        )
        {
            uniqueVertices.push_back(point);
        }
        remap[static_cast<std::size_t>(input)] =
            static_cast<int>(uniqueVertices.size()) - 1;
    }

    struct CanonicalFace
    {
        std::array<int, 3> indices;
        int birth = 0;
    };
    std::vector<CanonicalFace> faces;
    faces.reserve(static_cast<std::size_t>(mesh->faces.rows()));
    for (Eigen::Index row = 0; row < mesh->faces.rows(); ++row)
    {
        std::array<int, 3> indices = {{
            remap[static_cast<std::size_t>(mesh->faces(row, 0))],
            remap[static_cast<std::size_t>(mesh->faces(row, 1))],
            remap[static_cast<std::size_t>(mesh->faces(row, 2))]
        }};
        const Eigen::Vector3d& a =
            uniqueVertices[static_cast<std::size_t>(indices[0])];
        const Eigen::Vector3d& b =
            uniqueVertices[static_cast<std::size_t>(indices[1])];
        const Eigen::Vector3d& c =
            uniqueVertices[static_cast<std::size_t>(indices[2])];
        if (
            indices[0] == indices[1] || indices[1] == indices[2] ||
            indices[2] == indices[0] ||
            (b - a).cross(c - a).squaredNorm() == 0.0
        )
        {
            continue;
        }
        const int minimum = static_cast<int>(
            std::min_element(indices.begin(), indices.end()) - indices.begin()
        );
        std::array<int, 3> rotated = {{
            indices[minimum],
            indices[(minimum + 1) % 3],
            indices[(minimum + 2) % 3]
        }};
        faces.push_back({rotated, (*birthFaces)(row)});
    }
    std::sort(
        faces.begin(),
        faces.end(),
        [](const CanonicalFace& left, const CanonicalFace& right) {
            return std::tie(
                left.indices[0],
                left.indices[1],
                left.indices[2],
                left.birth
            ) < std::tie(
                right.indices[0],
                right.indices[1],
                right.indices[2],
                right.birth
            );
        }
    );

    std::vector<unsigned char> referenced(uniqueVertices.size(), 0U);
    for (const CanonicalFace& face : faces)
    {
        for (const int index : face.indices)
        {
            referenced[static_cast<std::size_t>(index)] = 1U;
        }
    }
    std::vector<int> compactRemap(uniqueVertices.size(), -1);
    int compactVertexCount = 0;
    for (std::size_t index = 0; index < uniqueVertices.size(); ++index)
    {
        if (referenced[index] != 0U)
        {
            compactRemap[index] = compactVertexCount++;
        }
    }
    Eigen::MatrixXd sortedVertices(compactVertexCount, 3);
    int destination = 0;
    for (std::size_t index = 0; index < uniqueVertices.size(); ++index)
    {
        if (referenced[index] != 0U)
        {
            sortedVertices.row(destination++) = uniqueVertices[index];
        }
    }
    mesh->vertices = std::move(sortedVertices);
    mesh->faces.resize(static_cast<Eigen::Index>(faces.size()), 3);
    birthFaces->resize(static_cast<Eigen::Index>(faces.size()));
    for (int row = 0; row < static_cast<int>(faces.size()); ++row)
    {
        mesh->faces.row(row) <<
            compactRemap[static_cast<std::size_t>(
                faces[static_cast<std::size_t>(row)].indices[0])],
            compactRemap[static_cast<std::size_t>(
                faces[static_cast<std::size_t>(row)].indices[1])],
            compactRemap[static_cast<std::size_t>(
                faces[static_cast<std::size_t>(row)].indices[2])];
        (*birthFaces)(row) = faces[static_cast<std::size_t>(row)].birth;
    }
}


Eigen::Matrix<double, 10, 1> rawMonomials(const Eigen::Vector3d& point)
{
    Eigen::Matrix<double, 10, 1> values;
    values <<
        1.0,
        point.x(), point.y(), point.z(),
        point.x() * point.x(),
        point.y() * point.y(),
        point.z() * point.z(),
        point.y() * point.z(),
        point.x() * point.z(),
        point.x() * point.y();
    return values;
}


Eigen::Matrix<double, 10, 1> coefficientMonomials(
    const Eigen::Vector3d& point
)
{
    Eigen::Matrix<double, 10, 1> values = rawMonomials(point);
    values.tail<3>() *= 2.0;
    return values;
}


struct QuadraturePoint
{
    std::array<double, 3> barycentric;
    double weight;
};


const std::array<QuadraturePoint, 7>& degreeFiveRule()
{
    static const std::array<QuadraturePoint, 7> rule = {{
        {{{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}}, 0.225},
        {{{0.059715871789770, 0.470142064105115, 0.470142064105115}},
            0.132394152788506},
        {{{0.470142064105115, 0.059715871789770, 0.470142064105115}},
            0.132394152788506},
        {{{0.470142064105115, 0.470142064105115, 0.059715871789770}},
            0.132394152788506},
        {{{0.797426985353087, 0.101286507323456, 0.101286507323456}},
            0.125939180544827},
        {{{0.101286507323456, 0.797426985353087, 0.101286507323456}},
            0.125939180544827},
        {{{0.101286507323456, 0.101286507323456, 0.797426985353087}},
            0.125939180544827}
    }};
    return rule;
}


void accumulateInterfaceDerivative(
    const spin::TriangleMesh& material,
    const Eigen::VectorXi& birthFaces,
    int designFaceCount,
    int quadricInterfaceFaceCount,
    const spin::QuadricCoefficients10& quadric,
    double density,
    double minimumGradient,
    Eigen::Matrix<double, 10, 10>* jacobian,
    int* interfaceTriangleCount,
    bool* derivativeDefined
)
{
    *jacobian = Eigen::Matrix<double, 10, 10>::Zero();
    *interfaceTriangleCount = 0;
    *derivativeDefined = true;
    for (Eigen::Index face = 0; face < material.faces.rows(); ++face)
    {
        const int birth = birthFaces(face);
        if (
            birth < designFaceCount ||
            birth >= designFaceCount + quadricInterfaceFaceCount
        )
        {
            continue;
        }
        ++(*interfaceTriangleCount);
        const Eigen::Vector3d a =
            material.vertices.row(material.faces(face, 0));
        const Eigen::Vector3d b =
            material.vertices.row(material.faces(face, 1));
        const Eigen::Vector3d c =
            material.vertices.row(material.faces(face, 2));
        const double area = 0.5 * (b - a).cross(c - a).norm();
        if (!std::isfinite(area) || area <= 0.0)
        {
            *derivativeDefined = false;
            continue;
        }
        for (const QuadraturePoint& sample : degreeFiveRule())
        {
            const Eigen::Vector3d point =
                sample.barycentric[0] * a +
                sample.barycentric[1] * b +
                sample.barycentric[2] * c;
            const double gradientNorm = quadric.gradient(point).norm();
            if (!std::isfinite(gradientNorm) || gradientNorm <= minimumGradient)
            {
                *derivativeDefined = false;
                continue;
            }
            *jacobian += density * area * sample.weight *
                rawMonomials(point) *
                coefficientMonomials(point).transpose() /
                gradientNorm;
        }
    }
    if (*interfaceTriangleCount == 0 || !jacobian->allFinite())
    {
        *derivativeDefined = false;
    }
    const Eigen::Matrix<double, 10, 10> tangentProjection =
        Eigen::Matrix<double, 10, 10>::Identity() -
        quadric.values() * quadric.values().transpose();
    *jacobian *= tangentProjection;
}


std::string topologyHash(
    const spin::TriangleMesh& material,
    const Eigen::VectorXi& birthFaces
)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << "QuadricMomentTopology-v1;" << material.vertices.rows() << ';'
           << material.faces.rows() << ';';
    for (Eigen::Index face = 0; face < material.faces.rows(); ++face)
    {
        stream << material.faces(face, 0) << ';'
               << material.faces(face, 1) << ';'
               << material.faces(face, 2) << ';'
               << birthFaces(face) << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}


} // namespace


namespace spin
{

const char* quadricMomentDerivativeMethodToken(QuadricMomentDerivativeMethod method) noexcept
{
    switch (method)
    {
    case QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape:
        return "legacy-continuous-implicit-shape-v1";
    case QuadricMomentDerivativeMethod::DiscreteEllipsoid:
        return "discrete-ellipsoid-boolean-transport-v1";
    }
    return "unknown";
}

void validateQuadricMomentEvaluatorSettings(const QuadricMomentEvaluatorSettings& settings)
{
    if (settings.derivativeMethod != QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape &&
        settings.derivativeMethod != QuadricMomentDerivativeMethod::DiscreteEllipsoid)
        throw std::invalid_argument("unknown quadric derivative method");
    if (!std::isfinite(settings.minimumInterfaceGradient) || settings.minimumInterfaceGradient <= 0.0)
        throw std::invalid_argument("minimum interface gradient must be finite and positive");
    const auto& t = settings.triangulation;
    if (t.icosphereSubdivisions < 0 || t.icosphereSubdivisions > 7 ||
        t.paraboloidGridResolution < 1 || t.cylinderAngularResolution < 8 ||
        t.maximumRadiusDoublings < 0 || !std::isfinite(t.truncationRadius) || t.truncationRadius < 0.0)
        throw std::invalid_argument("invalid quadric triangulation settings");
}


std::string computeQuadricEvaluatorRequestHash(
    const SpinDesignRequest& request,
    const QuadricMomentEvaluatorSettings& settings)
{
    validateQuadricMomentEvaluatorSettings(settings);
    const auto& t = settings.triangulation;
    const std::string base = computeSpinDesignRequestHash(request);
    if (settings.derivativeMethod == QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape &&
        settings.minimumInterfaceGradient == 1e-12 && t.icosphereSubdivisions == 3 &&
        t.paraboloidGridResolution == 24 && t.cylinderAngularResolution == 96 &&
        t.truncationRadius == 0.0 && t.maximumRadiusDoublings == 8 &&
        !t.retainFixedNormalProvenance)
        return base;
    std::ostringstream identity;
    identity.imbue(std::locale::classic());
    identity << "QuadricEvaluatorRequest-v1;" << base << ';'
        << static_cast<int>(settings.derivativeMethod) << ';'
        << t.icosphereSubdivisions << ';' << t.paraboloidGridResolution << ';'
        << t.cylinderAngularResolution << ';' << t.maximumRadiusDoublings << ';'
        << t.retainFixedNormalProvenance << ';' << std::hexfloat
        << t.truncationRadius << ';' << settings.minimumInterfaceGradient << ';';
    return manufacturing::computeContentSha256(identity.str());
}


QuadricMomentEvaluation evaluateQuadricMaterialMoments(
    const SpinDesignRequest& request,
    const QuadricCoefficients10& quadricTarget,
    const QuadricMomentEvaluatorSettings& settings
)
{
    const auto evaluationStart = std::chrono::steady_clock::now();
    validateSpinDesignRequest(request);
    validateQuadricMomentEvaluatorSettings(settings);

    const TriangleMesh designTarget =
        toTarget(request.designDomain, request.targetFrame);
    QuadricMomentEvaluation result;
    result.derivativeMethod = settings.derivativeMethod;
    QuadricTriangulationSettings triangulationSettings = settings.triangulation;
    if (settings.derivativeMethod == QuadricMomentDerivativeMethod::DiscreteEllipsoid)
        triangulationSettings.retainFixedNormalProvenance = true;
    result.quadricTriangulation = triangulateQuadric(
        quadricTarget,
        designTarget.vertices,
        triangulationSettings
    );

    const MassMoments quadricSolidMoments = computeSurfaceMassMoments(
        result.quadricTriangulation.closedBoundary.vertices,
        result.quadricTriangulation.closedBoundary.faces,
        1.0
    );
    Eigen::Vector3d quadricInteriorPoint = Eigen::Vector3d::Zero();
    if (!quadricSolidMoments.tryCenterOfMass(quadricInteriorPoint, 1e-15))
    {
        throw std::runtime_error(
            "truncated quadric solid has no positive finite volume"
        );
    }
    result.quadricSolidIsPositiveRegion =
        quadricTarget.evaluate(quadricInteriorPoint) > 0.0;

    Eigen::VectorXi birthFaces;
    const igl::MeshBooleanType operation =
        result.quadricSolidIsPositiveRegion
        ? igl::MESH_BOOLEAN_TYPE_INTERSECT
        : igl::MESH_BOOLEAN_TYPE_MINUS;
    const bool booleanSucceeded = igl::copyleft::cgal::mesh_boolean(
        designTarget.vertices,
        designTarget.faces,
        result.quadricTriangulation.closedBoundary.vertices,
        result.quadricTriangulation.closedBoundary.faces,
        operation,
        result.materialBoundaryTarget.vertices,
        result.materialBoundaryTarget.faces,
        birthFaces
    );
    if (
        !booleanSucceeded || result.materialBoundaryTarget.faces.rows() == 0 ||
        birthFaces.size() != result.materialBoundaryTarget.faces.rows()
    )
    {
        throw std::runtime_error(
            "quadric/design-domain Boolean produced no auditable material"
        );
    }
    canonicalizeBooleanMesh(&result.materialBoundaryTarget, &birthFaces);
    const manufacturing::BoundaryAudit materialAudit =
        manufacturing::auditTriangleBoundary(
            result.materialBoundaryTarget.vertices,
            result.materialBoundaryTarget.faces,
            manufacturing::BoundaryOrientation::AnyNonzero
        );
    if (!materialAudit.strictBoundaryValid() || materialAudit.signedVolume <= 0.0)
    {
        double minimumSquaredDoubleArea =
            std::numeric_limits<double>::infinity();
        for (Eigen::Index face = 0;
             face < result.materialBoundaryTarget.faces.rows(); ++face)
        {
            const Eigen::Vector3d a = result.materialBoundaryTarget.vertices.row(
                result.materialBoundaryTarget.faces(face, 0));
            const Eigen::Vector3d b = result.materialBoundaryTarget.vertices.row(
                result.materialBoundaryTarget.faces(face, 1));
            const Eigen::Vector3d c = result.materialBoundaryTarget.vertices.row(
                result.materialBoundaryTarget.faces(face, 2));
            minimumSquaredDoubleArea = std::min(
                minimumSquaredDoubleArea,
                (b - a).cross(c - a).squaredNorm()
            );
        }
        std::ostringstream diagnostic;
        diagnostic <<
            "quadric material Boolean is not a positive closed manifold: "
            "volume=" << materialAudit.signedVolume <<
            " matrix_shape=" << materialAudit.matrixShapeValid <<
            " finite=" << materialAudit.finiteVertices <<
            " indices=" << materialAudit.faceIndicesValid <<
            " nondegenerate=" << materialAudit.nondegenerateFaces <<
            " degenerate_faces=" << materialAudit.degenerateFaceCount <<
            " min_squared_double_area=" << minimumSquaredDoubleArea <<
            " unreferenced_vertices=" <<
                materialAudit.unreferencedVertexCount <<
            " closed=" << materialAudit.closed <<
            " edge_manifold=" << materialAudit.edgeManifold <<
            " vertex_manifold=" << materialAudit.vertexManifold <<
            " oriented=" << materialAudit.consistentlyOriented <<
            " component_orientation=" <<
                materialAudit.componentOrientationsValid <<
            " intersection_free=" <<
                materialAudit.selfIntersections.intersectionFree <<
            " intersection_completed=" <<
                materialAudit.selfIntersections.completed <<
            " intersection_pairs=" <<
                materialAudit.selfIntersections.intersectingPairCount <<
            " first_faces=" << materialAudit.selfIntersections.firstFace <<
            ',' << materialAudit.selfIntersections.secondFace;
        throw std::runtime_error(diagnostic.str());
    }

    const MassMoments variable = computeSurfaceMassMoments(
        result.materialBoundaryTarget.vertices,
        result.materialBoundaryTarget.faces,
        request.materialDensity
    );
    if (!variable.isFinite() || variable.mass <= 0.0)
    {
        throw std::runtime_error(
            "quadric material has invalid integrated mass moments"
        );
    }
    result.variableMomentsTarget =
        physics::RawMassMoments10::fromMassMoments(variable);
    const Eigen::Matrix3d modelToTarget = request.targetFrame.axes.transpose();
    const Eigen::Vector3d targetTranslation =
        -modelToTarget * request.targetFrame.origin;
    const physics::RawMassMoments10 fixedTarget =
        request.fixedShellMoments.transformed(
            modelToTarget,
            targetTranslation
        );
    result.totalMomentsTarget = result.variableMomentsTarget + fixedTarget;

    const auto residualEnd = std::chrono::steady_clock::now();
    if (settings.derivativeMethod == QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape)
    {
        accumulateInterfaceDerivative(
            result.materialBoundaryTarget,
            birthFaces,
            static_cast<int>(designTarget.faces.rows()),
            result.quadricTriangulation.interfaceFaceCount,
            quadricTarget,
            request.materialDensity,
            settings.minimumInterfaceGradient,
            &result.momentJacobian,
            &result.interfaceTriangleCount,
            &result.derivativeDefined
        );
        result.derivativeDiagnostic = result.derivativeDefined
            ? "legacy continuous implicit-shape quadrature; not the discrete geometry derivative"
            : "legacy interface derivative undefined";
    }
    else
    {
        const auto differential = differentiateDiscreteEllipsoidMoments(
            designTarget, result.quadricTriangulation, quadricTarget,
            result.materialBoundaryTarget, birthFaces, request.materialDensity);
        result.momentJacobian = differential.jacobian;
        result.derivativeDefined = differential.defined;
        result.derivativeDiagnostic = differential.diagnostic;
        result.vertexProvenance = differential.vertexProvenance;
        result.constructedIntersectionVertexCount = differential.constructedIntersectionVertexCount;
        result.maximumConstraintResidual = differential.maximumConstraintResidual;
        result.minimumConstraintSingularRatio = differential.minimumConstraintSingularRatio;
        for (int f = 0; f < birthFaces.size(); ++f)
            if (birthFaces(f) >= designTarget.faces.rows()) ++result.interfaceTriangleCount;
        if (differential.defined)
        {
            std::vector<std::string> vertices = result.vertexProvenance;
            std::sort(vertices.begin(), vertices.end());
            std::vector<std::string> faces;
            for (int f = 0; f < result.materialBoundaryTarget.faces.rows(); ++f)
            {
                std::array<std::string,3> keys;
                for (int k = 0; k < 3; ++k)
                    keys[k] = result.vertexProvenance[result.materialBoundaryTarget.faces(f,k)];
                const int first = static_cast<int>(std::min_element(keys.begin(),keys.end())-keys.begin());
                faces.push_back(keys[first]+";"+keys[(first+1)%3]+";"+keys[(first+2)%3]+";"+std::to_string(birthFaces(f)));
            }
            std::sort(faces.begin(),faces.end());
            std::ostringstream identity;
            identity << "QuadricSourceFeatureTopology-v1;";
            for (const auto& v : vertices) identity << v << '\n';
            for (const auto& f : faces) identity << f << '\n';
            result.provenanceTopologyHash = manufacturing::computeContentSha256(identity.str());
        }
    }
    const auto derivativeEnd = std::chrono::steady_clock::now();
    result.residualEvaluationSeconds =
        std::chrono::duration<double>(residualEnd - evaluationStart).count();
    result.derivativeEvaluationSeconds =
        std::chrono::duration<double>(derivativeEnd - residualEnd).count();
    result.materialBoundary =
        toModel(result.materialBoundaryTarget, request.targetFrame);
    result.geometryHash = meshHash(
        result.materialBoundaryTarget,
        "QuadricMomentGeometry-v1"
    );
    result.topologyHash = topologyHash(
        result.materialBoundaryTarget,
        birthFaces
    );
    std::ostringstream scientific;
    scientific.imbue(std::locale::classic());
    scientific << "QuadricMomentDifferentialScientific-v2;"
        << quadricMomentDerivativeMethodToken(result.derivativeMethod) << ';'
        << result.derivativeDefined << ';' << result.derivativeDiagnostic << ';'
        << computeQuadricEvaluatorRequestHash(request, settings) << ';'
        << result.quadricTriangulation.hash << ';' << result.geometryHash << ';'
        << result.provenanceTopologyHash << ';' << std::hexfloat
        << settings.minimumInterfaceGradient << ';';
    for (int j = 0; j < 10; ++j) scientific << quadricTarget.values()(j) << ';';
    result.derivativeScientificHash = manufacturing::computeContentSha256(scientific.str());
    return result;
}


} // namespace spin
