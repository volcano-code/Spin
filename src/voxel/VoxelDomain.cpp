#include "voxel/VoxelDomain.h"

#include "manufacturing/SlicerValidation.h"

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>


namespace voxel
{
namespace
{


std::size_t checkedProduct(
    const Eigen::Vector3i& dimensions,
    int increment,
    const char* label
)
{
    std::size_t product = 1U;
    for (int axis = 0; axis < 3; ++axis)
    {
        const long long side =
            static_cast<long long>(dimensions(axis)) + increment;
        if (side <= 0)
        {
            throw std::invalid_argument(
                std::string(label) + " dimensions must be positive"
            );
        }
        const std::size_t unsignedSide = static_cast<std::size_t>(side);
        if (product > std::numeric_limits<std::size_t>::max() / unsignedSide)
        {
            throw std::overflow_error(
                std::string(label) + " dimensions overflow size_t"
            );
        }
        product *= unsignedSide;
    }
    return product;
}


void appendU64(std::string& bytes, std::uint64_t value)
{
    for (int shift = 0; shift < 64; shift += 8)
    {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffU));
    }
}


void appendDouble(std::string& bytes, double value)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t));
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    appendU64(bytes, bits);
}


void validateMask(
    const std::vector<unsigned char>& mask,
    std::size_t expected,
    const char* label,
    bool allowEmpty
)
{
    if ((!allowEmpty || !mask.empty()) && mask.size() != expected)
    {
        throw std::invalid_argument(
            std::string(label) + " must contain exactly nx*ny*nz entries"
        );
    }
    if (std::any_of(mask.begin(), mask.end(), [](unsigned char value)
        {
            return value > 1U;
        }))
    {
        throw std::invalid_argument(
            std::string(label) + " entries must be zero or one"
        );
    }
}


} // namespace


std::shared_ptr<const VoxelDomain> VoxelDomain::create(
    const GridSpec& grid,
    const std::vector<unsigned char>& materialMaskXFast,
    const std::vector<unsigned char>& activeMaskXFast,
    const std::vector<unsigned char>& carveableMaskXFast
)
{
    const std::size_t cells = checkedProduct(grid.dimensions, 0, "grid");
    const std::size_t corners = checkedProduct(grid.dimensions, 1, "corner");
    if (
        !grid.minimum.allFinite() ||
        !grid.maximum.allFinite() ||
        (grid.maximum.array() <= grid.minimum.array()).any()
    )
    {
        throw std::invalid_argument(
            "grid bounds must be finite and strictly increasing"
        );
    }
    if (cells > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error(
            "voxel domain exceeds supported integer cell indexing"
        );
    }
    validateMask(materialMaskXFast, cells, "material mask", false);
    validateMask(activeMaskXFast, cells, "active mask", true);
    validateMask(carveableMaskXFast, cells, "carveable mask", true);
    if (std::none_of(
            materialMaskXFast.begin(),
            materialMaskXFast.end(),
            [](unsigned char value) { return value != 0U; }
        ))
    {
        throw std::invalid_argument("material mask contains no occupied cell");
    }

    std::shared_ptr<VoxelDomain> result(new VoxelDomain());
    result->grid_ = grid;
    result->spacing_ = (grid.maximum - grid.minimum).cwiseQuotient(
        grid.dimensions.cast<double>()
    );
    result->materialMask_ = materialMaskXFast;
    result->activeMask_ = activeMaskXFast.empty()
        ? materialMaskXFast
        : activeMaskXFast;
    result->carveableMask_ = carveableMaskXFast.empty()
        ? materialMaskXFast
        : carveableMaskXFast;

    for (std::size_t index = 0; index < cells; ++index)
    {
        if (
            result->carveableMask_[index] != 0U &&
            result->materialMask_[index] == 0U
        )
        {
            throw std::invalid_argument(
                "carveable mask must be a subset of the material mask"
            );
        }
        if (
            result->materialMask_[index] != 0U &&
            result->activeMask_[index] == 0U
        )
        {
            throw std::invalid_argument(
                "active mask must contain every material cell"
            );
        }
        if (result->activeMask_[index] != 0U)
        {
            result->activeCells_.push_back(static_cast<int>(index));
        }
    }

    static constexpr char schemaTag[] = "make-it-stand-voxel-domain-v1";
    std::string canonical(schemaTag, sizeof(schemaTag));
    for (int axis = 0; axis < 3; ++axis)
    {
        appendU64(canonical, static_cast<std::uint64_t>(grid.dimensions(axis)));
        appendDouble(canonical, grid.minimum(axis));
        appendDouble(canonical, grid.maximum(axis));
    }
    canonical.append(
        reinterpret_cast<const char*>(result->materialMask_.data()),
        result->materialMask_.size()
    );
    canonical.append(
        reinterpret_cast<const char*>(result->activeMask_.data()),
        result->activeMask_.size()
    );
    canonical.append(
        reinterpret_cast<const char*>(result->carveableMask_.data()),
        result->carveableMask_.size()
    );
    result->artifactId_ = manufacturing::computeContentSha256(canonical);

    result->estimate_.cellCount = static_cast<std::uint64_t>(cells);
    result->estimate_.cornerCount = static_cast<std::uint64_t>(corners);
    result->estimate_.materialCellCount = static_cast<std::uint64_t>(
        std::count(result->materialMask_.begin(), result->materialMask_.end(), 1U)
    );
    result->estimate_.activeCellCount = static_cast<std::uint64_t>(
        result->activeCells_.size()
    );
    result->estimate_.fillStateBytes = static_cast<std::uint64_t>(
        VoxelFillState::wordCount(cells) * sizeof(std::uint64_t)
    );
    result->estimate_.domainPayloadBytes = static_cast<std::uint64_t>(
        result->materialMask_.size() + result->activeMask_.size() +
        result->carveableMask_.size() +
        result->activeCells_.size() * sizeof(int)
    );
    return result;
}


const GridSpec& VoxelDomain::grid() const noexcept { return grid_; }
const Eigen::Vector3d& VoxelDomain::spacing() const noexcept { return spacing_; }
std::size_t VoxelDomain::cellCount() const noexcept { return materialMask_.size(); }


std::size_t VoxelDomain::cornerCount() const noexcept
{
    return
        static_cast<std::size_t>(grid_.dimensions(0) + 1) *
        static_cast<std::size_t>(grid_.dimensions(1) + 1) *
        static_cast<std::size_t>(grid_.dimensions(2) + 1);
}


std::size_t VoxelDomain::cellIndex(int x, int y, int z) const
{
    if (
        x < 0 || x >= grid_.dimensions(0) ||
        y < 0 || y >= grid_.dimensions(1) ||
        z < 0 || z >= grid_.dimensions(2)
    )
    {
        throw std::out_of_range("voxel cell coordinates are outside the grid");
    }
    return static_cast<std::size_t>(x) +
        static_cast<std::size_t>(grid_.dimensions(0)) *
        (static_cast<std::size_t>(y) +
         static_cast<std::size_t>(grid_.dimensions(1)) *
         static_cast<std::size_t>(z));
}


Eigen::Vector3i VoxelDomain::cellCoordinates(std::size_t index) const
{
    if (index >= cellCount())
    {
        throw std::out_of_range("voxel cell index is outside the grid");
    }
    const std::size_t nx = static_cast<std::size_t>(grid_.dimensions(0));
    const std::size_t ny = static_cast<std::size_t>(grid_.dimensions(1));
    const int x = static_cast<int>(index % nx);
    index /= nx;
    const int y = static_cast<int>(index % ny);
    const int z = static_cast<int>(index / ny);
    return Eigen::Vector3i(x, y, z);
}


std::size_t VoxelDomain::cornerIndex(int x, int y, int z) const
{
    const Eigen::Vector3i cornerDimensions =
        grid_.dimensions + Eigen::Vector3i::Ones();
    if (
        x < 0 || x >= cornerDimensions(0) ||
        y < 0 || y >= cornerDimensions(1) ||
        z < 0 || z >= cornerDimensions(2)
    )
    {
        throw std::out_of_range("voxel corner coordinates are outside the grid");
    }
    return static_cast<std::size_t>(x) +
        static_cast<std::size_t>(cornerDimensions(0)) *
        (static_cast<std::size_t>(y) +
         static_cast<std::size_t>(cornerDimensions(1)) *
         static_cast<std::size_t>(z));
}


Eigen::Vector3d VoxelDomain::cellCenter(int x, int y, int z) const
{
    static_cast<void>(cellIndex(x, y, z));
    return grid_.minimum + spacing_.cwiseProduct(Eigen::Vector3d(
        static_cast<double>(x) + 0.5,
        static_cast<double>(y) + 0.5,
        static_cast<double>(z) + 0.5
    ));
}


Eigen::Vector3d VoxelDomain::cornerPosition(int x, int y, int z) const
{
    static_cast<void>(cornerIndex(x, y, z));
    return grid_.minimum + spacing_.cwiseProduct(Eigen::Vector3d(
        static_cast<double>(x),
        static_cast<double>(y),
        static_cast<double>(z)
    ));
}


bool VoxelDomain::material(std::size_t index) const
{
    return materialMask_.at(index) != 0U;
}


bool VoxelDomain::active(std::size_t index) const
{
    return activeMask_.at(index) != 0U;
}


bool VoxelDomain::carveable(std::size_t index) const
{
    return carveableMask_.at(index) != 0U;
}


const std::vector<unsigned char>& VoxelDomain::materialMask() const noexcept
{
    return materialMask_;
}


const std::vector<unsigned char>& VoxelDomain::activeMask() const noexcept
{
    return activeMask_;
}


const std::vector<unsigned char>& VoxelDomain::carveableMask() const noexcept
{
    return carveableMask_;
}


const std::vector<int>& VoxelDomain::activeCells() const noexcept
{
    return activeCells_;
}


const std::string& VoxelDomain::artifactId() const noexcept
{
    return artifactId_;
}


const VoxelResourceEstimate& VoxelDomain::resourceEstimate() const noexcept
{
    return estimate_;
}


std::size_t VoxelFillState::wordCount(std::size_t bits) noexcept
{
    return (bits + 63U) / 64U;
}


VoxelFillState::VoxelFillState(std::shared_ptr<const VoxelDomain> domain)
    : domain_(std::move(domain))
{
    if (!domain_)
    {
        throw std::invalid_argument("voxel fill state requires a domain");
    }
    words_.assign(wordCount(domain_->cellCount()), 0U);
    refillMaterial();
}


const std::shared_ptr<const VoxelDomain>& VoxelFillState::domain() const noexcept
{
    return domain_;
}


bool VoxelFillState::filled(std::size_t index) const
{
    if (index >= domain_->cellCount())
    {
        throw std::out_of_range("voxel fill index is outside the domain");
    }
    return (words_[index / 64U] & (UINT64_C(1) << (index % 64U))) != 0U;
}


void VoxelFillState::setUnchecked(std::size_t index, bool value) noexcept
{
    const std::uint64_t bit = UINT64_C(1) << (index % 64U);
    std::uint64_t& word = words_[index / 64U];
    const bool previous = (word & bit) != 0U;
    if (previous == value)
    {
        return;
    }
    if (value)
    {
        word |= bit;
        ++filledCount_;
    }
    else
    {
        word &= ~bit;
        --filledCount_;
    }
}


void VoxelFillState::setFilled(std::size_t index, bool value)
{
    if (index >= domain_->cellCount())
    {
        throw std::out_of_range("voxel fill index is outside the domain");
    }
    if (value && !domain_->material(index))
    {
        throw std::invalid_argument(
            "fill state cannot add material outside its immutable domain"
        );
    }
    setUnchecked(index, value);
}


void VoxelFillState::refillMaterial()
{
    std::fill(words_.begin(), words_.end(), UINT64_C(0));
    filledCount_ = 0U;
    for (std::size_t index = 0; index < domain_->cellCount(); ++index)
    {
        if (domain_->material(index))
        {
            setUnchecked(index, true);
        }
    }
}


std::size_t VoxelFillState::filledCount() const noexcept
{
    return filledCount_;
}


std::uint64_t VoxelFillState::contentHash() const noexcept
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    for (const unsigned char byte : domain_->artifactId())
    {
        hash ^= static_cast<std::uint64_t>(byte);
        hash *= UINT64_C(1099511628211);
    }
    // Delimit the domain identity from the packed fill payload.
    hash ^= UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
    for (const std::uint64_t word : words_)
    {
        hash ^= word;
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}


} // namespace voxel
