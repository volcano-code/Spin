#include "verification/ExactTriangleSelfIntersection.h"

#include <Eigen/Core>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct MeshData {
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
};

void appendCube(MeshData& mesh, const Eigen::Vector3d& minimum,
                const Eigen::Vector3d& maximum) {
    const int offset = static_cast<int>(mesh.vertices.rows());
    Eigen::MatrixXd vertices(8, 3);
    vertices << minimum.x(), minimum.y(), minimum.z(),
                maximum.x(), minimum.y(), minimum.z(),
                maximum.x(), maximum.y(), minimum.z(),
                minimum.x(), maximum.y(), minimum.z(),
                minimum.x(), minimum.y(), maximum.z(),
                maximum.x(), minimum.y(), maximum.z(),
                maximum.x(), maximum.y(), maximum.z(),
                minimum.x(), maximum.y(), maximum.z();
    Eigen::MatrixXi faces(12, 3);
    faces << 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7,
             0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5,
             2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7;
    const Eigen::Index oldVertexCount = mesh.vertices.rows();
    const Eigen::Index oldFaceCount = mesh.faces.rows();
    mesh.vertices.conservativeResize(oldVertexCount + vertices.rows(), 3);
    mesh.faces.conservativeResize(oldFaceCount + faces.rows(), 3);
    mesh.vertices.bottomRows(vertices.rows()) = vertices;
    mesh.faces.bottomRows(faces.rows()) = faces.array() + offset;
}

MeshData twoCubes(double offset) {
    MeshData result{Eigen::MatrixXd(0, 3), Eigen::MatrixXi(0, 3)};
    appendCube(result, {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
    appendCube(result, {-1.0 + offset, -1.0 + offset, -1.0 + offset},
               {1.0 + offset, 1.0 + offset, 1.0 + offset});
    return result;
}

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

MeshData readTriangleObj(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open OBJ: " + path);
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector3i> faces;
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream row(line);
        std::string record;
        row >> record;
        if (record == "v") {
            Eigen::Vector3d point;
            if (!(row >> point.x() >> point.y() >> point.z()) || !point.allFinite())
                throw std::runtime_error("invalid OBJ vertex");
            vertices.push_back(point);
        } else if (record == "f") {
            Eigen::Vector3i face;
            std::string token;
            for (int corner = 0; corner < 3; ++corner) {
                if (!(row >> token)) throw std::runtime_error("OBJ face is not a triangle");
                const std::size_t separator = token.find('/');
                const int index = std::stoi(token.substr(0, separator));
                if (index <= 0) throw std::runtime_error("OBJ faces must use positive indices");
                face[corner] = index - 1;
            }
            if (row >> token) throw std::runtime_error("OBJ face has more than three corners");
            faces.push_back(face);
        }
    }
    MeshData mesh{Eigen::MatrixXd(vertices.size(), 3), Eigen::MatrixXi(faces.size(), 3)};
    for (std::size_t i = 0; i < vertices.size(); ++i) mesh.vertices.row(i) = vertices[i];
    for (std::size_t i = 0; i < faces.size(); ++i) mesh.faces.row(i) = faces[i];
    if (mesh.vertices.rows() == 0 || mesh.faces.rows() == 0)
        throw std::runtime_error("OBJ contains no triangle mesh");
    return mesh;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2) {
            const MeshData input = readTriangleObj(argv[1]);
            const auto audit = verification::auditExactTriangleMeshSelfIntersections(
                input.vertices, input.faces);
            std::cout << std::setprecision(17)
                      << "status=" << audit.status
                      << " vertices=" << input.vertices.rows()
                      << " faces=" << audit.checked_face_count
                      << " intersecting_face_pair_count=" << audit.intersecting_face_pairs.size()
                      << " diagnostic=" << audit.diagnostic << '\n';
            for (std::size_t i = 0; i < std::min<std::size_t>(10, audit.intersecting_face_pairs.size()); ++i)
                std::cout << "witness_face_pair=" << audit.intersecting_face_pairs[i].first
                          << ',' << audit.intersecting_face_pairs[i].second << '\n';
            return audit.status == "INTERSECTION_FREE" ? 0 : 2;
        }
        if (argc != 1) throw std::runtime_error("usage: exact-triangle-self-intersection-regression [triangle.obj]");
        const MeshData separated = twoCubes(4.0);
        const auto free = verification::auditExactTriangleMeshSelfIntersections(
            separated.vertices, separated.faces);
        require(free.status == "INTERSECTION_FREE", free.diagnostic);
        require(free.checked_face_count == 24, "did not check all separated-shell faces");
        require(free.intersecting_face_pairs.empty(), "disjoint shells reported an intersection");

        const MeshData crossing = twoCubes(0.5);
        const auto intersecting = verification::auditExactTriangleMeshSelfIntersections(
            crossing.vertices, crossing.faces);
        require(intersecting.status == "INTERSECTIONS_FOUND", intersecting.diagnostic);
        require(!intersecting.intersecting_face_pairs.empty(),
                "crossing closed shells have no exact face-pair witness");

        Eigen::MatrixXi openFaces = separated.faces.topRows(11);
        const auto open = verification::auditExactTriangleMeshSelfIntersections(
            separated.vertices.topRows(8), openFaces);
        require(open.status == "INVALID_OR_OPEN_INPUT",
                "open mesh was incorrectly admitted to the exact shell audit");
        std::cout << "[PASS] exact triangle self-intersection audit\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] exact triangle self-intersection audit: " << error.what() << '\n';
        return 1;
    }
}
