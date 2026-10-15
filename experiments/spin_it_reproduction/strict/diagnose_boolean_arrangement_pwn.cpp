#include "manufacturing/ManufacturingGeometry.h"
#include "spin/QuadricCoefficients.h"
#include "spin/QuadricTriangulator.h"

#include <igl/copyleft/cgal/RemeshSelfIntersectionsParam.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/bfs_orient.h>
#include <igl/piecewise_constant_winding_number.h>
#include <igl/readOBJ.h>
#include <igl/remove_unreferenced.h>
#include <igl/unique_edge_map.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>
#include <Eigen/Core>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using ExactScalar = CGAL::Exact_predicates_exact_constructions_kernel::FT;
using ExactVertices = Eigen::Matrix<ExactScalar, Eigen::Dynamic, 3>;

Eigen::Vector3d exactRowAsDouble(const ExactVertices& vertices, int row) {
    Eigen::Vector3d result;
    for (int column = 0; column < 3; ++column)
        result[column] = CGAL::to_double(vertices(row, column));
    return result;
}

void orientComponentsOutward(const Eigen::MatrixXd& vertices, Eigen::MatrixXi* faces) {
    Eigen::MatrixXi consistentlyOriented;
    Eigen::VectorXi component;
    igl::bfs_orient(*faces, consistentlyOriented, component);
    *faces = std::move(consistentlyOriented);
    for (int id = 0; id <= component.maxCoeff(); ++id) {
        double sixVolume = 0.0;
        for (Eigen::Index face = 0; face < faces->rows(); ++face) {
            if (component[face] != id) continue;
            const Eigen::Vector3d a = vertices.row((*faces)(face, 0));
            const Eigen::Vector3d b = vertices.row((*faces)(face, 1));
            const Eigen::Vector3d c = vertices.row((*faces)(face, 2));
            sixVolume += a.dot(b.cross(c));
        }
        if (sixVolume < 0.0) {
            for (Eigen::Index face = 0; face < faces->rows(); ++face)
                if (component[face] == id) std::swap((*faces)(face, 1), (*faces)(face, 2));
        }
    }
}

std::pair<int, double> gradientAlignment(const spin::TriangleMesh& mesh,
                                        const spin::QuadricCoefficients10& quadric,
                                        int interfaceFaceCount) {
    int negative = 0;
    double minimumCosine = 1.0;
    for (int face = 0; face < interfaceFaceCount; ++face) {
        const Eigen::Vector3d a = mesh.vertices.row(mesh.faces(face, 0));
        const Eigen::Vector3d b = mesh.vertices.row(mesh.faces(face, 1));
        const Eigen::Vector3d c = mesh.vertices.row(mesh.faces(face, 2));
        const Eigen::Vector3d normal = (b - a).cross(c - a).normalized();
        const Eigen::Vector3d desired = -quadric.gradient((a + b + c) / 3.0).normalized();
        const double cosine = normal.dot(desired);
        minimumCosine = std::min(minimumCosine, cosine);
        negative += cosine < 0.0 ? 1 : 0;
    }
    return {negative, minimumCosine};
}

void writeObj(const std::string& path, const spin::TriangleMesh& mesh) {
    std::ofstream output(path, std::ios::out | std::ios::trunc);
    if (!output) throw std::runtime_error("could not open operand export: " + path);
    output << std::setprecision(17);
    for (Eigen::Index row = 0; row < mesh.vertices.rows(); ++row)
        output << "v " << mesh.vertices(row, 0) << ' ' << mesh.vertices(row, 1)
               << ' ' << mesh.vertices(row, 2) << '\n';
    for (Eigen::Index row = 0; row < mesh.faces.rows(); ++row)
        output << "f " << mesh.faces(row, 0) + 1 << ' ' << mesh.faces(row, 1) + 1
               << ' ' << mesh.faces(row, 2) + 1 << '\n';
    if (!output) throw std::runtime_error("failed while writing operand export: " + path);
}

spin::TriangleMesh readTriangleObj(const std::string& path) {
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> faces;
    if (!igl::readOBJ(path, vertices, faces))
        throw std::runtime_error("could not read OBJ: " + path);
    spin::TriangleMesh mesh;
    mesh.vertices.resize(static_cast<Eigen::Index>(vertices.size()), 3);
    mesh.faces.resize(static_cast<Eigen::Index>(faces.size()), 3);
    for (Eigen::Index i = 0; i < mesh.vertices.rows(); ++i) {
        if (vertices[static_cast<std::size_t>(i)].size() != 3)
            throw std::runtime_error("OBJ vertex is not 3D");
        for (int c = 0; c < 3; ++c)
            mesh.vertices(i, c) = vertices[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)];
    }
    for (Eigen::Index i = 0; i < mesh.faces.rows(); ++i) {
        if (faces[static_cast<std::size_t>(i)].size() != 3)
            throw std::runtime_error("diagnostic requires a triangulated OBJ");
        for (int c = 0; c < 3; ++c)
            mesh.faces(i, c) = faces[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)];
    }
    return mesh;
}

bool sameDirectedEdge(const Eigen::MatrixXi& faces, int face, int from, int to) {
    for (int k = 0; k < 3; ++k)
        if (faces(face, k) == from && faces(face, (k + 1) % 3) == to)
            return true;
    return false;
}

struct Witness {
    int edge = -1;
    int from = -1;
    int to = -1;
    int signedIncidence = 0;
    std::vector<int> arrangementFaces;
};

Witness firstPwnWitness(const Eigen::MatrixXi& faces,
                        const Eigen::MatrixXi& uniqueEdges,
                        const std::vector<std::vector<std::size_t>>& edgeToHalfedge) {
    Witness witness;
    for (int edge = 0; edge < uniqueEdges.rows(); ++edge) {
        const int from = uniqueEdges(edge, 0);
        const int to = uniqueEdges(edge, 1);
        int incidence = 0;
        std::vector<int> incidentFaces;
        for (const auto halfedge : edgeToHalfedge[static_cast<std::size_t>(edge)]) {
            const int face = static_cast<int>(halfedge % static_cast<std::size_t>(faces.rows()));
            incidence += sameDirectedEdge(faces, face, from, to) ? 1 : -1;
            incidentFaces.push_back(face);
        }
        if (incidence != 0) {
            witness.edge = edge;
            witness.from = from;
            witness.to = to;
            witness.signedIncidence = incidence;
            witness.arrangementFaces = std::move(incidentFaces);
            return witness;
        }
    }
    return witness;
}

int diagnose(const std::string& designObj, int subdivisions, bool sourceOnly,
             const std::string& exportDirectory, bool alreadyTargetFrame) {
    auto design = readTriangleObj(designObj);
    if (design.vertices.rows() != 98 || design.faces.rows() != 192)
        throw std::runtime_error("input is not the pinned reconstructed Heart design mesh (expected 98V/192F)");

    // Match run_paper_table1.target_domain_points(): the solver evaluates the
    // quadric in the case's target frame, not in the raw OBJ frame.
    constexpr double origin[3] = {-0.5492066981089954, 0.0, -1.616538440729317};
    constexpr double axes[3][3] = {
        {0.9418358198519111, 3.1392484655341483e-16, 0.33607333789498767},
        {-6.106226635438361e-16, 1.0000000000000002, 7.771561172376096e-16},
        {-0.33607333789498767, -9.37167465562917e-16, 0.9418358198519112},
    };
    if (!alreadyTargetFrame) {
        for (Eigen::Index vertex = 0; vertex < design.vertices.rows(); ++vertex) {
            const Eigen::Vector3d point = design.vertices.row(vertex).transpose();
            for (int column = 0; column < 3; ++column) {
                design.vertices(vertex, column) = 0.0;
                for (int row = 0; row < 3; ++row)
                    design.vertices(vertex, column) += axes[row][column] * (point[row] - origin[row]);
            }
        }
    }

    spin::QuadricCoefficients10::Vector coefficients;
    coefficients << -0.5605562977479883, 0.2907335221898638,
        -0.02793678146101007, 0.7298376100915582,
        -0.01413811771179968, -0.0050173448614208305,
        -0.2346716529617043, 0.010029654990118242,
        -0.11140384401058069, -1.5370819321812222e-06;
    const auto quadric = spin::QuadricCoefficients10::fromUnnormalized(coefficients);
    spin::QuadricTriangulationSettings settings;
    settings.icosphereSubdivisions = subdivisions;
    settings.retainFixedNormalProvenance = true;
    const auto sourceTriangulation = spin::triangulateQuadric(
        quadric, design.vertices, settings);
    const auto& source = sourceTriangulation.closedBoundary;
    if (subdivisions == 3 &&
        (source.vertices.rows() != 280 || source.faces.rows() != 552))
        throw std::runtime_error("reconstructed Heart quadric source changed from 280V/552F");

    const auto designAudit = manufacturing::auditTriangleBoundary(
        design.vertices, design.faces, manufacturing::BoundaryOrientation::AnyNonzero);
    const auto sourceAudit = manufacturing::auditTriangleBoundary(
        source.vertices, source.faces, manufacturing::BoundaryOrientation::AnyNonzero);
    Eigen::MatrixXi repairedSourceFaces = source.faces;
    orientComponentsOutward(source.vertices, &repairedSourceFaces);
    const auto repairedSourceAudit = manufacturing::auditTriangleBoundary(
        source.vertices, repairedSourceFaces,
        manufacturing::BoundaryOrientation::OutwardPositive);
    spin::TriangleMesh repairedSource = source;
    repairedSource.faces = repairedSourceFaces;
    const auto rawGradientAlignment = gradientAlignment(
        source, quadric, sourceTriangulation.interfaceFaceCount);
    const auto repairedGradientAlignment = gradientAlignment(
        repairedSource, quadric, sourceTriangulation.interfaceFaceCount);
    const bool designPwn = igl::piecewise_constant_winding_number(design.faces);
    const bool sourcePwn = igl::piecewise_constant_winding_number(source.faces);
    Eigen::MatrixXi sourceEdges, sourceUniqueEdges;
    Eigen::VectorXi sourceEdgeMap;
    std::vector<std::vector<std::size_t>> sourceEdgeToHalfedge;
    igl::unique_edge_map(source.faces, sourceEdges, sourceUniqueEdges,
                         sourceEdgeMap, sourceEdgeToHalfedge);
    const auto sourceWitness = firstPwnWitness(
        source.faces, sourceUniqueEdges, sourceEdgeToHalfedge);

    if (!exportDirectory.empty()) {
        writeObj(exportDirectory + "/design.obj", design);
        writeObj(exportDirectory + "/quadric.obj", source);
    }

    if (sourceOnly) {
        std::cout << std::setprecision(17)
                  << "design_strict=" << designAudit.strictBoundaryValid()
                  << " source_vertices=" << source.vertices.rows()
                  << " source_faces=" << source.faces.rows()
                  << " subdivisions=" << subdivisions
                  << " source_strict=" << sourceAudit.strictBoundaryValid()
                  << " source_pwn=" << sourcePwn
                  << " source_inconsistent_edges=" << sourceAudit.inconsistentEdgeCount
                  << " source_components=" << sourceAudit.faceComponentCount
                  << " bfs_oriented_source_strict=" << repairedSourceAudit.strictBoundaryValid()
                  << " interface_gradient_negative_bfs=" << repairedGradientAlignment.first
                  << '\n';
        if (sourceWitness.edge >= 0) {
            std::cout << "first_source_bad_edge=" << sourceWitness.edge
                      << " endpoints=" << sourceWitness.from << ',' << sourceWitness.to
                      << " signed_incidence=" << sourceWitness.signedIncidence
                      << " q_mid=" << quadric.evaluate((0.5 * (
                          source.vertices.row(sourceWitness.from) +
                          source.vertices.row(sourceWitness.to))).eval())
                      << '\n';
        }
        return sourcePwn == (sourceWitness.edge < 0) ? 0 : 1;
    }

    Eigen::MatrixXd combinedVertices(design.vertices.rows() + source.vertices.rows(), 3);
    combinedVertices.topRows(design.vertices.rows()) = design.vertices;
    combinedVertices.bottomRows(source.vertices.rows()) = source.vertices;
    Eigen::MatrixXi combinedFaces(design.faces.rows() + source.faces.rows(), 3);
    combinedFaces.topRows(design.faces.rows()) = design.faces;
    combinedFaces.bottomRows(source.faces.rows()) =
        source.faces.array() + static_cast<int>(design.vertices.rows());

    igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
    params.stitch_all = true;
    ExactVertices remeshedVertices;
    Eigen::MatrixXi remeshedFaces, intersectingFacePairs;
    Eigen::VectorXi birthFaces, duplicateVertexMap;
    igl::copyleft::cgal::remesh_self_intersections(
        combinedVertices, combinedFaces, params,
        remeshedVertices, remeshedFaces, intersectingFacePairs,
        birthFaces, duplicateVertexMap);
    if (birthFaces.size() != remeshedFaces.rows() ||
        duplicateVertexMap.size() != remeshedVertices.rows())
        throw std::runtime_error("libigl remesh output/provenance shape mismatch");

    for (Eigen::Index i = 0; i < remeshedFaces.size(); ++i)
        remeshedFaces.data()[i] = duplicateVertexMap(remeshedFaces.data()[i]);
    ExactVertices stitchedVertices;
    Eigen::MatrixXi stitchedFaces;
    Eigen::VectorXi unreferencedMap;
    igl::remove_unreferenced(remeshedVertices, remeshedFaces,
                             stitchedVertices, stitchedFaces, unreferencedMap);

    Eigen::MatrixXi edges, uniqueEdges;
    Eigen::VectorXi edgeMap;
    std::vector<std::vector<std::size_t>> edgeToHalfedge;
    igl::unique_edge_map(stitchedFaces, edges, uniqueEdges, edgeMap, edgeToHalfedge);
    const bool arrangementPwn = igl::piecewise_constant_winding_number(
        stitchedFaces, uniqueEdges, edgeToHalfedge);
    const auto witness = firstPwnWitness(stitchedFaces, uniqueEdges, edgeToHalfedge);
    Eigen::MatrixXd booleanVertices;
    Eigen::MatrixXi booleanFaces;
    Eigen::VectorXi booleanBirthFaces;
    const bool booleanSucceeded = igl::copyleft::cgal::mesh_boolean(
        design.vertices, design.faces, source.vertices, source.faces,
        igl::MESH_BOOLEAN_TYPE_MINUS,
        booleanVertices, booleanFaces, booleanBirthFaces);
    Eigen::MatrixXd intersectionVertices;
    Eigen::MatrixXi intersectionFaces;
    Eigen::VectorXi intersectionBirthFaces;
    const bool intersectionSucceeded = igl::copyleft::cgal::mesh_boolean(
        design.vertices, design.faces, source.vertices, source.faces,
        igl::MESH_BOOLEAN_TYPE_INTERSECT,
        intersectionVertices, intersectionFaces, intersectionBirthFaces);
    ExactVertices exactBooleanVertices;
    Eigen::MatrixXi exactBooleanFaces;
    Eigen::VectorXi exactBooleanBirthFaces;
    const bool exactBooleanSucceeded = igl::copyleft::cgal::mesh_boolean(
        design.vertices, design.faces, source.vertices, source.faces,
        igl::MESH_BOOLEAN_TYPE_MINUS,
        exactBooleanVertices, exactBooleanFaces, exactBooleanBirthFaces);

    // Every output face inherits an input triangle. Exact rational plane tests
    // on the exported binary64 vertices quantify loss from constructing an
    // exact intersection point and then returning it through MatrixXd.
    const auto measureOffBirthPlane = [&](const auto& outputVertices,
                                          const Eigen::MatrixXi& outputFaces,
                                          const Eigen::VectorXi& outputBirthFaces) {
        using ExactPoint = Eigen::Matrix<ExactScalar, 3, 1>;
        const auto asExactPoint = [](const auto& vertices, int index) {
            ExactPoint point;
            for (int axis = 0; axis < 3; ++axis)
                point[axis] = ExactScalar(vertices(index, axis));
            return point;
        };
        std::uint64_t incidences = 0U;
        std::uint64_t offPlane = 0U;
        double maximumDistance = 0.0;
        for (Eigen::Index face = 0; face < outputFaces.rows() &&
             face < outputBirthFaces.size(); ++face) {
            const int parent = outputBirthFaces[face];
            const bool fromDesign = parent < design.faces.rows();
            const auto& parentFaces = fromDesign ? design.faces : source.faces;
            const Eigen::MatrixXd& parentVertices = fromDesign
                ? design.vertices : source.vertices;
            const int parentFace = fromDesign
                ? parent
                : parent - static_cast<int>(design.faces.rows());
            if (parentFace < 0 || parentFace >= parentFaces.rows()) continue;
            const ExactPoint a = asExactPoint(parentVertices, parentFaces(parentFace, 0));
            const ExactPoint b = asExactPoint(parentVertices, parentFaces(parentFace, 1));
            const ExactPoint c = asExactPoint(parentVertices, parentFaces(parentFace, 2));
            const ExactPoint normal = (b - a).cross(c - a);
            const double normalLength = std::sqrt(CGAL::to_double(normal.squaredNorm()));
            for (int corner = 0; corner < 3; ++corner) {
                ++incidences;
                const ExactPoint point = asExactPoint(outputVertices, outputFaces(face, corner));
                const ExactScalar residual = normal.dot(point - a);
                if (residual != ExactScalar(0)) {
                    ++offPlane;
                    if (normalLength > 0.0)
                        maximumDistance = std::max(maximumDistance,
                            std::abs(CGAL::to_double(residual)) / normalLength);
                }
            }
        }
        return std::tuple<std::uint64_t, std::uint64_t, double>{
            incidences, offPlane, maximumDistance};
    };
    const auto doublePlaneResidual = measureOffBirthPlane(
        booleanVertices, booleanFaces, booleanBirthFaces);
    const auto exactPlaneResidual = exactBooleanSucceeded
        ? measureOffBirthPlane(exactBooleanVertices, exactBooleanFaces, exactBooleanBirthFaces)
        : std::tuple<std::uint64_t, std::uint64_t, double>{0U, 0U, 0.0};
    Eigen::MatrixXd exactOutputRoundedVertices(exactBooleanVertices.rows(), 3);
    for (Eigen::Index row = 0; row < exactBooleanVertices.rows(); ++row)
        exactOutputRoundedVertices.row(row) = exactRowAsDouble(exactBooleanVertices, static_cast<int>(row));
    const auto exactOutputRoundedAudit = exactBooleanSucceeded && exactBooleanFaces.rows() > 0
        ? manufacturing::auditTriangleBoundary(
            exactOutputRoundedVertices, exactBooleanFaces,
            manufacturing::BoundaryOrientation::AnyNonzero)
        : manufacturing::BoundaryAudit{};
    const auto booleanAudit = booleanSucceeded && booleanFaces.rows() > 0
        ? manufacturing::auditTriangleBoundary(
            booleanVertices, booleanFaces,
            manufacturing::BoundaryOrientation::AnyNonzero)
        : manufacturing::BoundaryAudit{};
    if (!exportDirectory.empty() && booleanSucceeded && booleanFaces.rows() > 0) {
        spin::TriangleMesh libiglMaterial;
        libiglMaterial.vertices = booleanVertices;
        libiglMaterial.faces = booleanFaces;
        writeObj(exportDirectory + "/libigl-material.obj", libiglMaterial);
    }
    if (!exportDirectory.empty() && intersectionSucceeded && intersectionFaces.rows() > 0) {
        spin::TriangleMesh libiglIntersection;
        libiglIntersection.vertices = intersectionVertices;
        libiglIntersection.faces = intersectionFaces;
        writeObj(exportDirectory + "/libigl-intersection.obj", libiglIntersection);
    }

    std::cout << std::setprecision(17)
              << "design_strict=" << designAudit.strictBoundaryValid()
              << " design_pwn=" << designPwn
              << " source_strict=" << sourceAudit.strictBoundaryValid()
              << " source_pwn=" << sourcePwn
              << " source_boundary_edges=" << sourceAudit.boundaryEdgeCount
              << " source_nonmanifold_edges=" << sourceAudit.nonManifoldEdgeCount
              << " source_inconsistent_edges=" << sourceAudit.inconsistentEdgeCount
              << " source_nonmanifold_vertices=" << sourceAudit.nonManifoldVertexCount
              << " source_self_intersection_pairs=" << sourceAudit.selfIntersections.intersectingPairCount
              << " source_components=" << sourceAudit.faceComponentCount
              << " source_signed_volume=" << sourceAudit.signedVolume << '\n'
              << "bfs_oriented_source_strict=" << repairedSourceAudit.strictBoundaryValid()
              << " bfs_oriented_source_pwn=" << igl::piecewise_constant_winding_number(repairedSourceFaces)
              << " bfs_oriented_source_components=" << repairedSourceAudit.faceComponentCount
              << " bfs_oriented_source_signed_volume=" << repairedSourceAudit.signedVolume
              << " interface_gradient_negative_raw=" << rawGradientAlignment.first
              << " interface_gradient_min_cos_raw=" << rawGradientAlignment.second
              << " interface_gradient_negative_bfs=" << repairedGradientAlignment.first
              << " interface_gradient_min_cos_bfs=" << repairedGradientAlignment.second << '\n'
              << "design=98V/192F+source=" << source.vertices.rows() << 'V'
              << '/' << source.faces.rows() << "F subdivisions=" << subdivisions
              << " intersecting_face_pairs=" << intersectingFacePairs.rows()
              << " remeshed=" << stitchedVertices.rows() << "V/" << stitchedFaces.rows() << 'F'
              << " unique_edges=" << uniqueEdges.rows()
              << " arrangement_pwn=" << arrangementPwn << '\n';
    std::cout << "boolean_minus_succeeded=" << booleanSucceeded
              << " boolean_output=" << booleanVertices.rows() << "V/"
              << booleanFaces.rows() << "F birth_faces=" << booleanBirthFaces.size()
              << " output_strict=" << booleanAudit.strictBoundaryValid()
              << " output_pwn=" << (booleanFaces.rows() > 0 &&
                  igl::piecewise_constant_winding_number(booleanFaces))
              << " output_boundary_edges=" << booleanAudit.boundaryEdgeCount
              << " output_nonmanifold_edges=" << booleanAudit.nonManifoldEdgeCount
              << " output_inconsistent_edges=" << booleanAudit.inconsistentEdgeCount
              << " output_degenerate_faces=" << booleanAudit.degenerateFaceCount
              << " output_nonmanifold_vertices=" << booleanAudit.nonManifoldVertexCount
              << " output_unreferenced_vertices=" << booleanAudit.unreferencedVertexCount
              << " output_component_orientations_valid=" << booleanAudit.componentOrientationsValid
              << " output_intersection_check_completed=" << booleanAudit.selfIntersections.completed
              << " output_intersection_free=" << booleanAudit.selfIntersections.intersectionFree
              << " output_self_intersections=" << booleanAudit.selfIntersections.intersectingPairCount
              << " output_components=" << booleanAudit.faceComponentCount
              << " output_signed_volume=" << booleanAudit.signedVolume
              << " intersection_succeeded=" << intersectionSucceeded
              << " intersection_output=" << intersectionVertices.rows() << "V/"
              << intersectionFaces.rows() << "F"
              << " intersection_birth_faces=" << intersectionBirthFaces.size()
              << " double_output_birth_faces=" << booleanBirthFaces.size()
              << " double_output_vertex_incidences=" << std::get<0>(doublePlaneResidual)
              << " double_off_birth_plane_incidences=" << std::get<1>(doublePlaneResidual)
              << " double_max_off_plane_distance=" << std::get<2>(doublePlaneResidual)
              << " exact_output_succeeded=" << exactBooleanSucceeded
              << " exact_output=" << exactBooleanVertices.rows() << "V/"
              << exactBooleanFaces.rows() << "F"
              << " exact_output_birth_faces=" << exactBooleanBirthFaces.size()
              << " exact_output_vertex_incidences=" << std::get<0>(exactPlaneResidual)
              << " exact_off_birth_plane_incidences=" << std::get<1>(exactPlaneResidual)
              << " exact_max_off_plane_distance=" << std::get<2>(exactPlaneResidual)
              << " exact_output_rounded_strict=" << exactOutputRoundedAudit.strictBoundaryValid()
              << " exact_output_rounded_degenerate_faces=" << exactOutputRoundedAudit.degenerateFaceCount
              << " exact_output_rounded_intersection_check_completed="
              << exactOutputRoundedAudit.selfIntersections.completed
              << " exact_output_rounded_intersection_free="
              << exactOutputRoundedAudit.selfIntersections.intersectionFree << '\n';
    int reportedDegenerates = 0;
    for (Eigen::Index face = 0; face < booleanFaces.rows() && reportedDegenerates < 8; ++face) {
        const Eigen::Vector3d a = booleanVertices.row(booleanFaces(face, 0));
        const Eigen::Vector3d b = booleanVertices.row(booleanFaces(face, 1));
        const Eigen::Vector3d c = booleanVertices.row(booleanFaces(face, 2));
        const double squaredDoubleArea = (b - a).cross(c - a).squaredNorm();
        if (squaredDoubleArea < 1e-18) {
            std::cout << "  small_area_face=" << face
                      << " squared_double_area=" << squaredDoubleArea
                      << " birth_face=" << (face < booleanBirthFaces.size()
                            ? booleanBirthFaces[face] : -1)
                      << " vertices=" << booleanFaces(face, 0) << ','
                      << booleanFaces(face, 1) << ',' << booleanFaces(face, 2)
                      << '\n';
            ++reportedDegenerates;
        }
    }
    if (sourceWitness.edge >= 0) {
        const auto from = source.vertices.row(sourceWitness.from).transpose();
        const auto to = source.vertices.row(sourceWitness.to).transpose();
        const Eigen::Vector3d midpoint = 0.5 * (from + to);
        std::cout << "first_source_bad_edge=" << sourceWitness.edge
                  << " endpoints=" << sourceWitness.from << ',' << sourceWitness.to
                  << " signed_incidence=" << sourceWitness.signedIncidence
                  << " p0=" << from.transpose() << " p1=" << to.transpose()
                  << " q_mid=" << quadric.evaluate(midpoint.transpose())
                  << " grad_mid=" << quadric.gradient(midpoint).norm() << '\n';
        for (int face : sourceWitness.arrangementFaces)
            std::cout << "  source_face=" << face
                      << " cap=" << (sourceTriangulation.faceIsCap(face) != 0)
                      << " vertices=" << source.faces(face, 0) << ','
                      << source.faces(face, 1) << ',' << source.faces(face, 2)
                      << " direction="
                      << (sameDirectedEdge(source.faces, face,
                          sourceWitness.from, sourceWitness.to) ? "+" : "-") << '\n';
    }
    if (witness.edge >= 0) {
        const auto from = exactRowAsDouble(stitchedVertices, witness.from);
        const auto to = exactRowAsDouble(stitchedVertices, witness.to);
        std::cout << "first_bad_edge=" << witness.edge
                  << " endpoints=" << witness.from << ',' << witness.to
                  << " signed_incidence=" << witness.signedIncidence
                  << " valence=" << witness.arrangementFaces.size()
                  << " p0=" << from.transpose() << " p1=" << to.transpose() << '\n';
        for (int face : witness.arrangementFaces) {
            const int parent = birthFaces(face);
            const bool fromDesign = parent < design.faces.rows();
            std::cout << "  arrangement_face=" << face
                      << " direction=" << (sameDirectedEdge(stitchedFaces, face,
                          witness.from, witness.to) ? "+" : "-")
                      << " parent=" << (fromDesign ? "design:" : "quadric:")
                      << (fromDesign ? parent : parent - static_cast<int>(design.faces.rows()))
                      << '\n';
        }
    }
    if (!designAudit.strictBoundaryValid() || !designPwn)
        throw std::runtime_error("isolated design boundary failed its strict/PWN audit");
    if (sourcePwn != (sourceWitness.edge < 0))
        throw std::runtime_error("source edge witness disagrees with libigl PWN predicate");
    if (arrangementPwn != (witness.edge < 0))
        throw std::runtime_error("edge witness disagrees with libigl PWN predicate");
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 7) {
            std::cerr << "usage: diagnose-boolean-arrangement-pwn HEART_DESIGN.obj [ICOSPHERE_SUBDIVISIONS] [--source-only] [--already-target-frame] [--export-operands DIR]\n";
            return 2;
        }
        int subdivisions = 3;
        bool sourceOnly = false;
        bool alreadyTargetFrame = false;
        std::string exportDirectory;
        for (int index = 2; index < argc; ++index) {
            if (std::string(argv[index]) == "--source-only") sourceOnly = true;
            else if (std::string(argv[index]) == "--already-target-frame") alreadyTargetFrame = true;
            else if (std::string(argv[index]) == "--export-operands") {
                if (++index >= argc) throw std::invalid_argument("--export-operands requires a directory");
                exportDirectory = argv[index];
            }
            else subdivisions = std::stoi(argv[index]);
        }
        return diagnose(argv[1], subdivisions, sourceOnly, exportDirectory, alreadyTargetFrame);
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
