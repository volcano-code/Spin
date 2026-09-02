#include "spin/QuadricTriangulator.h"

#include "manufacturing/SlicerValidation.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>


namespace
{


using Face = std::array<int, 3>;
using Edge = std::pair<int, int>;


struct SphereMesh
{
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Face> faces;
};


Edge undirectedEdge(int first, int second)
{
    return std::minmax(first, second);
}


SphereMesh icosphere(int subdivisions)
{
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    SphereMesh mesh;
    mesh.vertices = {
        {-1.0, phi, 0.0}, {1.0, phi, 0.0},
        {-1.0, -phi, 0.0}, {1.0, -phi, 0.0},
        {0.0, -1.0, phi}, {0.0, 1.0, phi},
        {0.0, -1.0, -phi}, {0.0, 1.0, -phi},
        {phi, 0.0, -1.0}, {phi, 0.0, 1.0},
        {-phi, 0.0, -1.0}, {-phi, 0.0, 1.0}
    };
    for (Eigen::Vector3d& vertex : mesh.vertices)
    {
        vertex.normalize();
    }
    mesh.faces = {
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
        {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
        {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
        {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}
    };

    for (int level = 0; level < subdivisions; ++level)
    {
        std::map<Edge, int> midpointIndices;
        std::vector<Face> refined;
        refined.reserve(mesh.faces.size() * 4);
        const auto midpoint = [&](int first, int second) {
            const Edge edge = undirectedEdge(first, second);
            const auto found = midpointIndices.find(edge);
            if (found != midpointIndices.end())
            {
                return found->second;
            }
            Eigen::Vector3d point =
                (mesh.vertices[first] + mesh.vertices[second]).normalized();
            const int index = static_cast<int>(mesh.vertices.size());
            mesh.vertices.push_back(point);
            midpointIndices.emplace(edge, index);
            return index;
        };
        for (const Face& face : mesh.faces)
        {
            const int ab = midpoint(face[0], face[1]);
            const int bc = midpoint(face[1], face[2]);
            const int ca = midpoint(face[2], face[0]);
            refined.push_back({face[0], ab, ca});
            refined.push_back({face[1], bc, ab});
            refined.push_back({face[2], ca, bc});
            refined.push_back({ab, bc, ca});
        }
        mesh.faces = std::move(refined);
    }
    return mesh;
}


double automaticRadius(const Eigen::MatrixXd& vertices)
{
    if (vertices.rows() == 0 || vertices.cols() != 3 || !vertices.allFinite())
    {
        throw std::invalid_argument(
            "quadric triangulation requires finite design-domain vertices"
        );
    }
    const Eigen::Vector3d minimum = vertices.colwise().minCoeff();
    const Eigen::Vector3d maximum = vertices.colwise().maxCoeff();
    const double diameter = (maximum - minimum).norm();
    double maximumNorm = 0.0;
    for (Eigen::Index row = 0; row < vertices.rows(); ++row)
    {
        maximumNorm = std::max(maximumNorm, vertices.row(row).norm());
    }
    return std::max({1.0, 2.0 * diameter, 2.0 * maximumNorm});
}


spin::QuadricSurfaceType classify(
    const spin::QuadricCoefficients10& quadric,
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>* eigensolver
)
{
    const Eigen::Matrix3d a = quadric.quadraticMatrix();
    eigensolver->compute(a);
    if (eigensolver->info() != Eigen::Success)
    {
        throw std::invalid_argument("quadric eigen-decomposition failed");
    }
    const double scale = std::max(1.0, a.norm());
    const double tolerance = 1e-12 * scale;
    int rank = 0;
    for (int i = 0; i < 3; ++i)
    {
        rank += std::abs(eigensolver->eigenvalues()(i)) > tolerance ? 1 : 0;
    }
    if (rank == 0 && quadric.linearVector().norm() > tolerance)
    {
        return spin::QuadricSurfaceType::Plane;
    }
    if (rank == 2)
    {
        int nullIndex = 0;
        eigensolver->eigenvalues().cwiseAbs().minCoeff(&nullIndex);
        const Eigen::Vector3d nullVector =
            eigensolver->eigenvectors().col(nullIndex);
        if (std::abs(nullVector.dot(quadric.linearVector())) > tolerance)
        {
            return spin::QuadricSurfaceType::Paraboloid;
        }
        std::array<int, 2> nonzero{};
        int cursor = 0;
        for (int index = 0; index < 3; ++index)
        {
            if (index != nullIndex)
            {
                nonzero[cursor++] = index;
            }
        }
        const double first = eigensolver->eigenvalues()(nonzero[0]);
        const double second = eigensolver->eigenvalues()(nonzero[1]);
        if (first * second > 0.0)
        {
            return spin::QuadricSurfaceType::EllipticCylinder;
        }
    }
    if (rank != 3)
    {
        throw std::invalid_argument(
            "unsupported degenerate quadric (only planes/paraboloids allowed)"
        );
    }

    const Eigen::Vector3d eigenvalues = eigensolver->eigenvalues();
    const bool positive = (eigenvalues.array() > tolerance).all();
    const bool negative = (eigenvalues.array() < -tolerance).all();
    if (positive || negative)
    {
        const Eigen::Matrix3d inverse = a.inverse();
        const double shiftedConstant = quadric.constant() -
            0.25 * quadric.linearVector().dot(
                inverse * quadric.linearVector()
            );
        const bool nonempty =
            (positive && shiftedConstant < -tolerance) ||
            (negative && shiftedConstant > tolerance);
        if (!nonempty)
        {
            throw std::invalid_argument("ellipsoid has no real surface");
        }
        return spin::QuadricSurfaceType::Ellipsoid;
    }
    return spin::QuadricSurfaceType::Hyperboloid;
}


bool inverseGaussPoint(
    const spin::QuadricCoefficients10& quadric,
    const Eigen::Matrix3d& inverseA,
    const Eigen::Vector3d& normal,
    Eigen::Vector3d* point
)
{
    const double denominator =
        quadric.linearVector().dot(inverseA * quadric.linearVector()) -
        4.0 * quadric.constant();
    if (std::abs(denominator) <= 1e-15)
    {
        return false;
    }
    const double wSquared =
        4.0 * normal.dot(inverseA * normal) / denominator;
    if (!std::isfinite(wSquared) || wSquared <= 1e-14)
    {
        return false;
    }
    const double w = std::sqrt(wSquared);
    *point = inverseA *
        (normal / w - 0.5 * quadric.linearVector());
    return point->allFinite();
}


void orientInterfaceFace(
    const spin::QuadricCoefficients10& quadric,
    const std::vector<Eigen::Vector3d>& vertices,
    Face* face
)
{
    const Eigen::Vector3d& a = vertices[(*face)[0]];
    const Eigen::Vector3d& b = vertices[(*face)[1]];
    const Eigen::Vector3d& c = vertices[(*face)[2]];
    const Eigen::Vector3d normal = (b - a).cross(c - a);
    const Eigen::Vector3d desired =
        -quadric.gradient((a + b + c) / 3.0);
    if (normal.dot(desired) < 0.0)
    {
        std::swap((*face)[1], (*face)[2]);
    }
}


std::vector<std::vector<int>> boundaryLoops(
    const std::vector<Face>& faces
)
{
    std::map<Edge, int> counts;
    std::map<Edge, std::pair<int, int>> directions;
    for (const Face& face : faces)
    {
        for (int corner = 0; corner < 3; ++corner)
        {
            const int first = face[corner];
            const int second = face[(corner + 1) % 3];
            const Edge key = undirectedEdge(first, second);
            ++counts[key];
            directions[key] = {first, second};
        }
    }
    std::map<int, int> next;
    for (const auto& entry : counts)
    {
        if (entry.second == 1)
        {
            const std::pair<int, int> directed = directions.at(entry.first);
            if (!next.emplace(directed.first, directed.second).second)
            {
                throw std::runtime_error(
                    "quadric truncation boundary is not a simple loop"
                );
            }
        }
    }
    std::vector<std::vector<int>> loops;
    while (!next.empty())
    {
        std::vector<int> loop;
        const int start = next.begin()->first;
        int current = start;
        do
        {
            loop.push_back(current);
            const auto found = next.find(current);
            if (found == next.end())
            {
                throw std::runtime_error("quadric truncation loop is open");
            }
            current = found->second;
            next.erase(found);
        } while (current != start);
        if (loop.size() < 3)
        {
            throw std::runtime_error("quadric truncation loop is degenerate");
        }
        loops.push_back(std::move(loop));
    }
    return loops;
}


int idealNormalEdgeIntersections(
    const SphereMesh& sphere,
    const Eigen::Matrix3d& inverseA
)
{
    std::set<Edge> edges;
    for (const Face& face : sphere.faces)
    {
        edges.insert(undirectedEdge(face[0], face[1]));
        edges.insert(undirectedEdge(face[1], face[2]));
        edges.insert(undirectedEdge(face[2], face[0]));
    }
    int count = 0;
    for (const Edge& edge : edges)
    {
        const Eigen::Vector3d& first = sphere.vertices[edge.first];
        const Eigen::Vector3d difference =
            sphere.vertices[edge.second] - first;
        const double quadratic = difference.dot(inverseA * difference);
        const double linear = 2.0 * first.dot(inverseA * difference);
        const double constant = first.dot(inverseA * first);
        if (std::abs(quadratic) <= 1e-15)
        {
            if (std::abs(linear) > 1e-15)
            {
                const double root = -constant / linear;
                count += root > 0.0 && root < 1.0 ? 1 : 0;
            }
            continue;
        }
        const double discriminant =
            linear * linear - 4.0 * quadratic * constant;
        if (discriminant < 0.0)
        {
            continue;
        }
        const double squareRoot = std::sqrt(std::max(0.0, discriminant));
        const double firstRoot = (-linear - squareRoot) / (2.0 * quadratic);
        const double secondRoot = (-linear + squareRoot) / (2.0 * quadratic);
        count += firstRoot > 0.0 && firstRoot < 1.0 ? 1 : 0;
        if (squareRoot > 1e-15)
        {
            count += secondRoot > 0.0 && secondRoot < 1.0 ? 1 : 0;
        }
    }
    return count;
}


spin::TriangleMesh matrixMesh(
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<Face>& faces
)
{
    spin::TriangleMesh result;
    result.vertices.resize(static_cast<Eigen::Index>(vertices.size()), 3);
    for (std::size_t row = 0; row < vertices.size(); ++row)
    {
        result.vertices.row(static_cast<Eigen::Index>(row)) = vertices[row];
    }
    result.faces.resize(static_cast<Eigen::Index>(faces.size()), 3);
    for (std::size_t row = 0; row < faces.size(); ++row)
    {
        result.faces.row(static_cast<Eigen::Index>(row)) <<
            faces[row][0], faces[row][1], faces[row][2];
    }
    return result;
}


void makePositiveVolume(spin::TriangleMesh* mesh)
{
    double sixVolume = 0.0;
    for (Eigen::Index face = 0; face < mesh->faces.rows(); ++face)
    {
        const Eigen::Vector3d a =
            mesh->vertices.row(mesh->faces(face, 0));
        const Eigen::Vector3d b =
            mesh->vertices.row(mesh->faces(face, 1));
        const Eigen::Vector3d c =
            mesh->vertices.row(mesh->faces(face, 2));
        sixVolume += a.dot(b.cross(c));
    }
    if (sixVolume < 0.0)
    {
        mesh->faces.col(1).swap(mesh->faces.col(2));
    }
}


spin::QuadricTriangulation closeExtrudedInterface(
    spin::QuadricSurfaceType type,
    const std::vector<Eigen::Vector3d>& interfaceVertices,
    const std::vector<Face>& interfaceFaces,
    const Eigen::Vector3d& extrusion,
    double radius
)
{
    const std::vector<std::vector<int>> loops =
        boundaryLoops(interfaceFaces);
    if (loops.size() != 1)
    {
        throw std::runtime_error(
            "plane/paraboloid interface must have one boundary loop"
        );
    }
    std::vector<Eigen::Vector3d> vertices = interfaceVertices;
    const int offset = static_cast<int>(vertices.size());
    for (const Eigen::Vector3d& vertex : interfaceVertices)
    {
        vertices.push_back(vertex + extrusion);
    }
    std::vector<Face> faces = interfaceFaces;
    for (const Face& face : interfaceFaces)
    {
        faces.push_back({
            face[2] + offset,
            face[1] + offset,
            face[0] + offset
        });
    }
    const std::vector<int>& loop = loops.front();
    for (std::size_t edge = 0; edge < loop.size(); ++edge)
    {
        const int first = loop[edge];
        const int second = loop[(edge + 1) % loop.size()];
        faces.push_back({second, first, first + offset});
        faces.push_back({second, first + offset, second + offset});
    }

    spin::QuadricTriangulation result;
    result.type = type;
    result.closedBoundary = matrixMesh(vertices, faces);
    makePositiveVolume(&result.closedBoundary);
    result.interfaceVertexCount = offset;
    result.interfaceFaceCount = static_cast<int>(interfaceFaces.size());
    result.capFaceCount =
        static_cast<int>(faces.size() - interfaceFaces.size());
    result.boundaryLoopCount = 1;
    result.truncationRadius = radius;
    result.faceIsCap = Eigen::VectorXi::Ones(
        static_cast<Eigen::Index>(faces.size())
    );
    result.faceIsCap.head(result.interfaceFaceCount).setZero();
    return result;
}


std::string triangulationHash(const spin::QuadricTriangulation& result)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << "QuadricTriangulation-v1;"
           << spin::quadricSurfaceTypeToken(result.type) << ';'
           << result.interfaceVertexCount << ';'
           << result.interfaceFaceCount << ';'
           << result.capFaceCount << ';'
           << result.boundaryLoopCount << ';'
           << result.idealNormalEdgeIntersectionCount << ';'
           << std::hexfloat << result.truncationRadius << ';';
    for (Eigen::Index row = 0; row < result.closedBoundary.vertices.rows(); ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            stream << result.closedBoundary.vertices(row, column) << ';';
        }
    }
    for (Eigen::Index row = 0; row < result.closedBoundary.faces.rows(); ++row)
    {
        stream << result.closedBoundary.faces(row, 0) << ';'
               << result.closedBoundary.faces(row, 1) << ';'
               << result.closedBoundary.faces(row, 2) << ';'
               << result.faceIsCap(row) << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}


spin::QuadricTriangulation triangulatePlane(
    const spin::QuadricCoefficients10& quadric,
    double radius
)
{
    const Eigen::Vector3d normal = quadric.linearVector();
    const Eigen::Vector3d point =
        -quadric.constant() * normal / normal.squaredNorm();
    const Eigen::Vector3d first = normal.normalized().unitOrthogonal();
    const Eigen::Vector3d second = normal.normalized().cross(first);
    std::vector<Eigen::Vector3d> vertices = {
        point - radius * first - radius * second,
        point + radius * first - radius * second,
        point + radius * first + radius * second,
        point - radius * first + radius * second
    };
    std::vector<Face> faces = {{0, 1, 2}, {0, 2, 3}};
    for (Face& face : faces)
    {
        orientInterfaceFace(quadric, vertices, &face);
    }
    spin::QuadricTriangulation result = closeExtrudedInterface(
        spin::QuadricSurfaceType::Plane,
        vertices,
        faces,
        4.0 * radius * normal.normalized(),
        radius
    );
    result.hash = triangulationHash(result);
    return result;
}


spin::QuadricTriangulation triangulateParaboloid(
    const spin::QuadricCoefficients10& quadric,
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& eigen,
    double radius,
    int resolution
)
{
    int nullIndex = 0;
    eigen.eigenvalues().cwiseAbs().minCoeff(&nullIndex);
    std::array<int, 2> nonzero{};
    int cursor = 0;
    for (int index = 0; index < 3; ++index)
    {
        if (index != nullIndex)
        {
            nonzero[cursor++] = index;
        }
    }
    Eigen::Vector3d n = eigen.eigenvectors().col(nullIndex);
    Eigen::Vector3d u = eigen.eigenvectors().col(nonzero[0]);
    Eigen::Vector3d v = eigen.eigenvectors().col(nonzero[1]);
    double bn = quadric.linearVector().dot(n);
    if (bn < 0.0)
    {
        n = -n;
        bn = -bn;
    }
    const double bu = quadric.linearVector().dot(u);
    const double bv = quadric.linearVector().dot(v);
    const double au = eigen.eigenvalues()(nonzero[0]);
    const double av = eigen.eigenvalues()(nonzero[1]);

    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(static_cast<std::size_t>((resolution + 1) * (resolution + 1)));
    for (int row = 0; row <= resolution; ++row)
    {
        const double t = -radius + 2.0 * radius * row / resolution;
        for (int column = 0; column <= resolution; ++column)
        {
            const double s = -radius + 2.0 * radius * column / resolution;
            const double w = -(
                quadric.constant() + bu * s + bv * t +
                au * s * s + av * t * t
            ) / bn;
            vertices.push_back(u * s + v * t + n * w);
        }
    }
    std::vector<Face> faces;
    faces.reserve(static_cast<std::size_t>(2 * resolution * resolution));
    for (int row = 0; row < resolution; ++row)
    {
        for (int column = 0; column < resolution; ++column)
        {
            const int lower = row * (resolution + 1) + column;
            Face first = {lower, lower + 1, lower + resolution + 2};
            Face second = {lower, lower + resolution + 2, lower + resolution + 1};
            orientInterfaceFace(quadric, vertices, &first);
            orientInterfaceFace(quadric, vertices, &second);
            faces.push_back(first);
            faces.push_back(second);
        }
    }
    spin::QuadricTriangulation result = closeExtrudedInterface(
        spin::QuadricSurfaceType::Paraboloid,
        vertices,
        faces,
        4.0 * radius * n,
        radius
    );
    result.hash = triangulationHash(result);
    return result;
}


spin::QuadricTriangulation triangulateEllipticCylinder(
    const spin::QuadricCoefficients10& quadric,
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& eigen,
    double radius,
    int resolution
)
{
    int nullIndex = 0;
    eigen.eigenvalues().cwiseAbs().minCoeff(&nullIndex);
    std::array<int, 2> nonzero{};
    int cursor = 0;
    for (int index = 0; index < 3; ++index)
    {
        if (index != nullIndex)
        {
            nonzero[cursor++] = index;
        }
    }
    const Eigen::Vector3d n = eigen.eigenvectors().col(nullIndex);
    const Eigen::Vector3d u = eigen.eigenvectors().col(nonzero[0]);
    const Eigen::Vector3d v = eigen.eigenvectors().col(nonzero[1]);
    const double au = eigen.eigenvalues()(nonzero[0]);
    const double av = eigen.eigenvalues()(nonzero[1]);
    const double bu = quadric.linearVector().dot(u);
    const double bv = quadric.linearVector().dot(v);
    const double centerU = -bu / (2.0 * au);
    const double centerV = -bv / (2.0 * av);
    const double shifted = quadric.constant() -
        bu * bu / (4.0 * au) - bv * bv / (4.0 * av);
    const double radiusSquaredU = -shifted / au;
    const double radiusSquaredV = -shifted / av;
    if (radiusSquaredU <= 0.0 || radiusSquaredV <= 0.0)
    {
        throw std::invalid_argument("elliptic cylinder has no real surface");
    }
    const double radiusU = std::sqrt(radiusSquaredU);
    const double radiusV = std::sqrt(radiusSquaredV);

    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(static_cast<std::size_t>(2 * resolution + 2));
    for (int layer = 0; layer < 2; ++layer)
    {
        const double height = layer == 0 ? -2.0 * radius : 2.0 * radius;
        for (int segment = 0; segment < resolution; ++segment)
        {
            const double angle =
                2.0 * std::acos(-1.0) * segment / resolution;
            vertices.push_back(
                u * (centerU + radiusU * std::cos(angle)) +
                v * (centerV + radiusV * std::sin(angle)) +
                n * height
            );
        }
    }
    std::vector<Face> faces;
    faces.reserve(static_cast<std::size_t>(4 * resolution));
    for (int segment = 0; segment < resolution; ++segment)
    {
        const int next = (segment + 1) % resolution;
        Face first = {segment, next, next + resolution};
        Face second = {segment, next + resolution, segment + resolution};
        orientInterfaceFace(quadric, vertices, &first);
        orientInterfaceFace(quadric, vertices, &second);
        faces.push_back(first);
        faces.push_back(second);
    }
    const int interfaceVertexCount = static_cast<int>(vertices.size());
    const int interfaceFaceCount = static_cast<int>(faces.size());
    const int lowerCenter = static_cast<int>(vertices.size());
    vertices.push_back(u * centerU + v * centerV - 2.0 * radius * n);
    const int upperCenter = static_cast<int>(vertices.size());
    vertices.push_back(u * centerU + v * centerV + 2.0 * radius * n);
    for (int segment = 0; segment < resolution; ++segment)
    {
        const int next = (segment + 1) % resolution;
        faces.push_back({next, segment, lowerCenter});
        faces.push_back({
            segment + resolution,
            next + resolution,
            upperCenter
        });
    }

    spin::QuadricTriangulation result;
    result.type = spin::QuadricSurfaceType::EllipticCylinder;
    result.closedBoundary = matrixMesh(vertices, faces);
    makePositiveVolume(&result.closedBoundary);
    result.interfaceVertexCount = interfaceVertexCount;
    result.interfaceFaceCount = interfaceFaceCount;
    result.capFaceCount = static_cast<int>(faces.size()) - interfaceFaceCount;
    result.boundaryLoopCount = 2;
    result.truncationRadius = radius;
    result.faceIsCap = Eigen::VectorXi::Ones(
        static_cast<Eigen::Index>(faces.size())
    );
    result.faceIsCap.head(interfaceFaceCount).setZero();
    result.hash = triangulationHash(result);
    return result;
}


} // namespace


namespace spin
{


const char* quadricSurfaceTypeToken(QuadricSurfaceType type) noexcept
{
    switch (type)
    {
    case QuadricSurfaceType::Ellipsoid:
        return "ELLIPSOID";
    case QuadricSurfaceType::Hyperboloid:
        return "HYPERBOLOID";
    case QuadricSurfaceType::Paraboloid:
        return "PARABOLOID";
    case QuadricSurfaceType::EllipticCylinder:
        return "ELLIPTIC_CYLINDER";
    case QuadricSurfaceType::Plane:
        return "PLANE";
    }
    return "UNKNOWN";
}


QuadricTriangulation triangulateQuadric(
    const QuadricCoefficients10& quadric,
    const Eigen::MatrixXd& designDomainVertices,
    const QuadricTriangulationSettings& settings
)
{
    if (
        settings.icosphereSubdivisions < 0 ||
        settings.icosphereSubdivisions > 7 ||
        settings.paraboloidGridResolution < 1 ||
        settings.cylinderAngularResolution < 8 ||
        settings.maximumRadiusDoublings < 0 ||
        !std::isfinite(settings.truncationRadius) ||
        settings.truncationRadius < 0.0
    )
    {
        throw std::invalid_argument("invalid quadric triangulation settings");
    }
    double radius = settings.truncationRadius;
    if (radius == 0.0)
    {
        radius = automaticRadius(designDomainVertices);
    }
    else
    {
        (void)automaticRadius(designDomainVertices);
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen;
    const QuadricSurfaceType type = classify(quadric, &eigen);
    if (type == QuadricSurfaceType::Plane)
    {
        return triangulatePlane(quadric, radius);
    }
    if (type == QuadricSurfaceType::Paraboloid)
    {
        return triangulateParaboloid(
            quadric,
            eigen,
            radius,
            settings.paraboloidGridResolution
        );
    }
    if (type == QuadricSurfaceType::EllipticCylinder)
    {
        return triangulateEllipticCylinder(
            quadric,
            eigen,
            radius,
            settings.cylinderAngularResolution
        );
    }

    const Eigen::Matrix3d inverseA = quadric.quadraticMatrix().inverse();
    const SphereMesh sphere = icosphere(settings.icosphereSubdivisions);
    if (type == QuadricSurfaceType::Ellipsoid)
    {
        std::vector<Eigen::Vector3d> vertices;
        vertices.reserve(sphere.vertices.size());
        for (const Eigen::Vector3d& normal : sphere.vertices)
        {
            Eigen::Vector3d point;
            if (!inverseGaussPoint(quadric, inverseA, normal, &point))
            {
                throw std::runtime_error(
                    "ellipsoid inverse Gauss map unexpectedly failed"
                );
            }
            vertices.push_back(point);
        }
        std::vector<Face> faces = sphere.faces;
        for (Face& face : faces)
        {
            orientInterfaceFace(quadric, vertices, &face);
        }
        QuadricTriangulation result;
        result.type = type;
        result.closedBoundary = matrixMesh(vertices, faces);
        makePositiveVolume(&result.closedBoundary);
        result.faceIsCap = Eigen::VectorXi::Zero(
            static_cast<Eigen::Index>(faces.size())
        );
        result.interfaceVertexCount = static_cast<int>(vertices.size());
        result.interfaceFaceCount = static_cast<int>(faces.size());
        result.truncationRadius = radius;
        result.hash = triangulationHash(result);
        return result;
    }

    for (int doubling = 0; doubling <= settings.maximumRadiusDoublings;
         ++doubling)
    {
        std::vector<bool> inside(sphere.vertices.size(), false);
        std::vector<Eigen::Vector3d> mapped(sphere.vertices.size());
        for (std::size_t index = 0; index < sphere.vertices.size(); ++index)
        {
            inside[index] = inverseGaussPoint(
                quadric,
                inverseA,
                sphere.vertices[index],
                &mapped[index]
            ) && mapped[index].norm() <= radius;
        }

        std::vector<Eigen::Vector3d> vertices;
        std::vector<Face> faces;
        std::map<int, int> originalVertex;
        std::map<Edge, int> crossingVertex;
        const auto addOriginal = [&](int index) {
            const auto found = originalVertex.find(index);
            if (found != originalVertex.end())
            {
                return found->second;
            }
            const int output = static_cast<int>(vertices.size());
            vertices.push_back(mapped[index]);
            originalVertex.emplace(index, output);
            return output;
        };
        const auto addCrossing = [&](int first, int second) {
            const Edge edge = undirectedEdge(first, second);
            const auto found = crossingVertex.find(edge);
            if (found != crossingVertex.end())
            {
                return found->second;
            }
            int insideIndex = inside[first] ? first : second;
            int outsideIndex = inside[first] ? second : first;
            Eigen::Vector3d insideNormal = sphere.vertices[insideIndex];
            Eigen::Vector3d outsideNormal = sphere.vertices[outsideIndex];
            Eigen::Vector3d insidePoint = mapped[insideIndex];
            for (int iteration = 0; iteration < 80; ++iteration)
            {
                Eigen::Vector3d middleNormal =
                    0.5 * (insideNormal + outsideNormal);
                middleNormal.normalize();
                Eigen::Vector3d middlePoint;
                const bool middleInside = inverseGaussPoint(
                    quadric,
                    inverseA,
                    middleNormal,
                    &middlePoint
                ) && middlePoint.norm() <= radius;
                if (middleInside)
                {
                    insideNormal = middleNormal;
                    insidePoint = middlePoint;
                }
                else
                {
                    outsideNormal = middleNormal;
                }
            }
            const int output = static_cast<int>(vertices.size());
            vertices.push_back(insidePoint);
            crossingVertex.emplace(edge, output);
            return output;
        };

        for (const Face& input : sphere.faces)
        {
            std::vector<int> polygon;
            for (int corner = 0; corner < 3; ++corner)
            {
                const int current = input[corner];
                const int next = input[(corner + 1) % 3];
                if (inside[current])
                {
                    polygon.push_back(addOriginal(current));
                }
                if (inside[current] != inside[next])
                {
                    polygon.push_back(addCrossing(current, next));
                }
            }
            for (std::size_t corner = 1; corner + 1 < polygon.size(); ++corner)
            {
                Face face = {polygon[0], polygon[corner], polygon[corner + 1]};
                orientInterfaceFace(quadric, vertices, &face);
                faces.push_back(face);
            }
        }
        const std::vector<std::vector<int>> loops = boundaryLoops(faces);
        if (loops.size() != 2)
        {
            radius *= 2.0;
            continue;
        }

        const int interfaceVertexCount = static_cast<int>(vertices.size());
        const int interfaceFaceCount = static_cast<int>(faces.size());
        for (const std::vector<int>& loop : loops)
        {
            Eigen::Vector3d center = Eigen::Vector3d::Zero();
            for (int vertex : loop)
            {
                center += vertices[vertex];
            }
            center /= static_cast<double>(loop.size());
            const int centerIndex = static_cast<int>(vertices.size());
            vertices.push_back(center);
            for (std::size_t edge = 0; edge < loop.size(); ++edge)
            {
                faces.push_back({
                    loop[(edge + 1) % loop.size()],
                    loop[edge],
                    centerIndex
                });
            }
        }

        TriangleMesh boundary = matrixMesh(vertices, faces);
        makePositiveVolume(&boundary);

        QuadricTriangulation result;
        result.type = type;
        result.closedBoundary = std::move(boundary);
        result.interfaceVertexCount = interfaceVertexCount;
        result.interfaceFaceCount = interfaceFaceCount;
        result.capFaceCount =
            static_cast<int>(faces.size()) - interfaceFaceCount;
        result.boundaryLoopCount = static_cast<int>(loops.size());
        result.idealNormalEdgeIntersectionCount =
            idealNormalEdgeIntersections(sphere, inverseA);
        result.truncationRadius = radius;
        result.faceIsCap = Eigen::VectorXi::Zero(
            static_cast<Eigen::Index>(faces.size())
        );
        result.faceIsCap.tail(result.capFaceCount).setOnes();
        result.hash = triangulationHash(result);
        return result;
    }
    throw std::runtime_error(
        "hyperboloid truncation did not produce exactly two closed loops"
    );
}


} // namespace spin
