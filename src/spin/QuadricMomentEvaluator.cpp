#include "spin/QuadricMomentEvaluator.h"

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


QuadricMomentEvaluation evaluateQuadricMaterialMoments(
    const SpinDesignRequest& request,
    const QuadricCoefficients10& quadricTarget,
    const QuadricMomentEvaluatorSettings& settings
)
{
    const auto evaluationStart = std::chrono::steady_clock::now();
    validateSpinDesignRequest(request);
    if (
        !std::isfinite(settings.minimumInterfaceGradient) ||
        settings.minimumInterfaceGradient <= 0.0
    )
    {
        throw std::invalid_argument(
            "minimum interface gradient must be finite and positive"
        );
    }

    const TriangleMesh designTarget =
        toTarget(request.designDomain, request.targetFrame);
    QuadricMomentEvaluation result;
    result.quadricTriangulation = triangulateQuadric(
        quadricTarget,
        designTarget.vertices,
        settings.triangulation
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
    return result;
}


} // namespace spin
