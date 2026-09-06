#pragma once
#include <Eigen/Core>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace verification {

enum class ExactBoundaryOrientation { OutwardPositive = 0, InwardNegative = 1, AnyNonzero = 2 };

struct ExactBoundaryFacts {
    std::string algorithm = "indexed-simplex-exact-rational-bvh-v1";
    std::string backend = "CGAL-f7c3c8212b56c0d6dae63787efc99093f4383415/Simple_cartesian<Gmpq>/GMP-6.3.0/MPFR-4.2.1";
    std::string inputDigest;
    std::string orientation;
    std::string intersectionError;
    bool boundaryValid = false;
    bool closed = false;
    bool edgeManifold = false;
    bool vertexManifold = false;
    bool consistentlyOriented = false;
    bool componentOrientationsValid = false;
    bool intersectionCompleted = false;
    bool intersectionFree = false;
    std::uint64_t exactPairTests = 0;
    std::vector<int> unreferencedVertices, degenerateFaces, nonManifoldVertices;
    std::vector<std::array<int,2>> boundaryEdges, nonManifoldEdges, inconsistentEdges;
    std::vector<std::array<int,2>> intersectingPairs;
    std::vector<int> faceComponents, componentVolumeSigns;
    std::vector<std::string> componentSignedVolumesExact;
};

// An immutable result computed from represented binary64 coordinates and indexed
// triangles. No report is accepted as input and no artifact/seal status is granted.
class ExactBoundaryAudit {
public:
    const ExactBoundaryFacts& facts() const { return facts_; }
private:
    ExactBoundaryFacts facts_;
    friend ExactBoundaryAudit auditExactTriangleBoundary(
        const Eigen::MatrixXd&, const Eigen::MatrixXi&, ExactBoundaryOrientation);
};

// Invalid shape, nonfinite/empty input and invalid/repeated face indices throw
// invalid_argument before geometric predicates. Unreferenced vertices and exact
// zero-area faces produce invalid facts. Open valid-index meshes retain diagnostics.
// AnyNonzero does not establish nesting, containment, material union or thickness.
ExactBoundaryAudit auditExactTriangleBoundary(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces,
    ExactBoundaryOrientation orientation);

} // namespace verification
