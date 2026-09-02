#pragma once

#include <Eigen/Core>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>


namespace voxel
{


/**
 * Explicit axis-aligned lattice contract used by material-occupancy imports.
 *
 * Cell payloads use canonical x-fast order:
 * ``x + nx * (y + ny * z)``.  This differs from several historical cache
 * formats on purpose; legacy adapters must perform an explicit coordinate
 * remap instead of silently reinterpreting identifiers.
 */
struct GridSpec
{
    Eigen::Vector3i dimensions = Eigen::Vector3i::Zero();
    Eigen::Vector3d minimum = Eigen::Vector3d::Zero();
    Eigen::Vector3d maximum = Eigen::Vector3d::Ones();
};


struct VoxelResourceEstimate
{
    std::uint64_t cellCount = 0U;
    std::uint64_t cornerCount = 0U;
    std::uint64_t materialCellCount = 0U;
    std::uint64_t activeCellCount = 0U;
    std::uint64_t fillStateBytes = 0U;
    std::uint64_t domainPayloadBytes = 0U;
};


/**
 * Immutable, shareable topology and masks for a regular material lattice.
 *
 * The object is intentionally independent from the legacy VoxelGrid facade.
 * New production code can share one domain across restarts while each
 * candidate owns only a compact VoxelFillState bitset.
 */
class VoxelDomain
{
public:
    static std::shared_ptr<const VoxelDomain> create(
        const GridSpec& grid,
        const std::vector<unsigned char>& materialMaskXFast,
        const std::vector<unsigned char>& activeMaskXFast = {},
        const std::vector<unsigned char>& carveableMaskXFast = {}
    );

    const GridSpec& grid() const noexcept;
    const Eigen::Vector3d& spacing() const noexcept;
    std::size_t cellCount() const noexcept;
    std::size_t cornerCount() const noexcept;

    std::size_t cellIndex(int x, int y, int z) const;
    Eigen::Vector3i cellCoordinates(std::size_t index) const;
    std::size_t cornerIndex(int x, int y, int z) const;
    Eigen::Vector3d cellCenter(int x, int y, int z) const;
    Eigen::Vector3d cornerPosition(int x, int y, int z) const;

    bool material(std::size_t index) const;
    bool active(std::size_t index) const;
    bool carveable(std::size_t index) const;
    const std::vector<unsigned char>& materialMask() const noexcept;
    const std::vector<unsigned char>& activeMask() const noexcept;
    const std::vector<unsigned char>& carveableMask() const noexcept;
    const std::vector<int>& activeCells() const noexcept;

    const std::string& artifactId() const noexcept;
    const VoxelResourceEstimate& resourceEstimate() const noexcept;

private:
    VoxelDomain() = default;

    GridSpec grid_;
    Eigen::Vector3d spacing_ = Eigen::Vector3d::Zero();
    std::vector<unsigned char> materialMask_;
    std::vector<unsigned char> activeMask_;
    std::vector<unsigned char> carveableMask_;
    std::vector<int> activeCells_;
    std::string artifactId_;
    VoxelResourceEstimate estimate_;
};


/** Mutable occupancy snapshot backed by a bitset and one shared domain. */
class VoxelFillState
{
public:
    explicit VoxelFillState(std::shared_ptr<const VoxelDomain> domain);

    const std::shared_ptr<const VoxelDomain>& domain() const noexcept;
    bool filled(std::size_t index) const;
    void setFilled(std::size_t index, bool value);
    void refillMaterial();
    std::size_t filledCount() const noexcept;
    std::uint64_t contentHash() const noexcept;

    /** Exact compact storage required for a bitset containing ``bits``. */
    static std::size_t wordCount(std::size_t bits) noexcept;

private:
    void setUnchecked(std::size_t index, bool value) noexcept;

    std::shared_ptr<const VoxelDomain> domain_;
    std::vector<std::uint64_t> words_;
    std::size_t filledCount_ = 0U;
};


} // namespace voxel
