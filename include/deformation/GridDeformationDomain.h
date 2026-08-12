#pragma once

#include <Eigen/Core>

#include <array>


class VoxelGrid;


namespace deformation
{


/** Trilinear coordinates of one point in an active voxel-grid cell. */
struct TrilinearEmbedding
{
    /** Local deformation-node indices in x-fastest binary corner order. */
    std::array<int, 8> nodes{};

    /** Non-negative trilinear coefficients corresponding to nodes. */
    std::array<double, 8> coeffs{};

    int cellLinearIndex = -1;
};


/**
 * Fixed-connectivity graph and interpolation domain for Section 5.
 *
 * Nodes are the union of corners belonging to cells that are inside the
 * input or can intersect its surface according to the signed-distance
 * half-diagonal test.  Local node numbering follows ascending global corner
 * numbering, making bindings deterministic across runs.
 */
class GridDeformationDomain
{
public:
    static GridDeformationDomain fromVoxelGrid(
        const VoxelGrid& voxelGrid
    );

    Eigen::Index nodeCount() const noexcept;
    Eigen::Index edgeCount() const noexcept;
    Eigen::Index cellCount() const noexcept;
    Eigen::Index activeCellCount() const noexcept;

    const Eigen::MatrixXd& restNodes() const noexcept;
    const Eigen::MatrixXi& edges() const noexcept;
    const Eigen::VectorXi& localToGlobalCorner() const noexcept;
    const Eigen::VectorXi& localToGlobal() const noexcept;
    const Eigen::VectorXi& globalToLocalCorners() const noexcept;
    const Eigen::VectorXi& globalToLocal() const noexcept;
    const Eigen::VectorXi& globalToLocalCorner() const noexcept;
    const Eigen::MatrixXi& cellNodes() const noexcept;
    const Eigen::VectorXi& activeCellMask() const noexcept;
    const Eigen::Vector3i& dimensions() const noexcept;
    const Eigen::Vector3d& voxelSize() const noexcept;
    const Eigen::MatrixXd& fullRestCorners() const noexcept;

    int componentCount() const noexcept;
    const Eigen::VectorXi& nodeComponents() const noexcept;
    bool connected() const noexcept;

    int globalToLocalCorner(int globalCorner) const;
    int localToGlobalCorner(int localNode) const;
    int localNodeForGlobalCorner(int globalCorner) const;
    int globalCornerForLocalNode(int localNode) const;

    int cellLinearIndex(int x, int y, int z) const;
    Eigen::Vector3i cellCoordinates(int linearIndex) const;
    bool isActiveCell(int linearIndex) const;
    std::array<int, 8> cellNodes(int linearIndex) const;

    /** Full-grid corner indices in x/y/z binary-bit order. */
    std::array<int, 8> cellGlobalCorners(int linearIndex) const;

    /** Embed a finite point in a deterministically selected active cell. */
    TrilinearEmbedding embed(const Eigen::Vector3d& point) const;

    /** Return the nearest active rest node; ties select the lowest index. */
    int nearestActiveNode(const Eigen::Vector3d& point) const;

    /**
     * Trilinearly transfer partition-of-unity node weights to arbitrary
     * points.  Inputs outside the active domain are rejected.
     */
    Eigen::MatrixXd interpolateWeights(
        const Eigen::MatrixXd& points,
        const Eigen::MatrixXd& nodeWeights
    ) const;

    /** Average the eight deformed corners of every full-grid cell. */
    Eigen::MatrixXd deformedCellCenters(
        const Eigen::MatrixXd& fullDeformedCorners
    ) const;

private:
    Eigen::MatrixXd restNodes_;
    Eigen::MatrixXi edges_;
    Eigen::VectorXi localToGlobalCorner_;
    Eigen::VectorXi globalToLocalCorner_;
    Eigen::MatrixXi cellNodes_;
    Eigen::MatrixXi cellGlobalCorners_;
    Eigen::VectorXi activeCellMask_;
    Eigen::VectorXi nodeComponents_;
    Eigen::Vector3i dimensions_ = Eigen::Vector3i::Zero();
    Eigen::Vector3d voxelSize_ = Eigen::Vector3d::Zero();
    Eigen::MatrixXd fullRestCorners_;
    Eigen::Vector3d gridMinimum_ = Eigen::Vector3d::Zero();
    Eigen::Vector3d gridMaximum_ = Eigen::Vector3d::Zero();
    double geometricTolerance_ = 0.0;
    int activeCellCount_ = 0;
    int componentCount_ = 0;
};


} // namespace deformation
