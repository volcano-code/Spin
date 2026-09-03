#include "spin2014/AdaptiveVoxelDomain.h"

#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace
{


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-it-2014-adaptive-domain "
              << message << '\n';
    return EXIT_FAILURE;
}


std::shared_ptr<const voxel::VoxelDomain> fixture()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(4, 4, 4);
    grid.minimum = Eigen::Vector3d(-1.0, -2.0, -3.0);
    grid.maximum = Eigen::Vector3d(1.0, 2.0, 3.0);
    std::vector<unsigned char> material(64U, 1U);
    std::vector<unsigned char> carveable(64U, 1U);
    // Protect one finest cell as fixed shell material.
    carveable[0] = 0U;
    return voxel::VoxelDomain::create(grid, material, {}, carveable);
}


} // namespace


int main()
{
    const auto source = fixture();
    spin2014::AdaptiveVoxelDomain domain(source, 1, 2);
    if (domain.leaves().size() != 8U)
    {
        return fail("mid-level initialization did not create eight leaves");
    }
    if (
        std::abs(domain.solidMoments().mass() - 48.0) > 1e-12 ||
        std::abs(domain.variableMoments().mass() - 47.25) > 1e-12
    )
    {
        return fail("analytic cuboid moment sums are not conserved");
    }

    physics::RawMassMoments10 sum;
    for (const spin2014::AdaptiveVoxelCell& cell : domain.leaves())
    {
        sum += cell.moments;
    }
    if (!sum.values().isApprox(domain.variableMoments().values(), 1e-12))
    {
        return fail("leaf moments do not partition the carveable domain");
    }

    const std::string initialHash = domain.partitionHash();
    Eigen::VectorXd fill = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(domain.leaves().size())
    );
    fill(0) = 0.5;
    if (
        domain.splitFractionalLeaves(fill, 0.05) != 1U ||
        domain.leaves().size() != 14U ||
        initialHash == domain.partitionHash()
    )
    {
        return fail("deterministic split did not replace one leaf by children");
    }
    sum = physics::RawMassMoments10();
    for (const spin2014::AdaptiveVoxelCell& cell : domain.leaves())
    {
        sum += cell.moments;
    }
    if (!sum.values().isApprox(domain.variableMoments().values(), 1e-12))
    {
        return fail("split failed to conserve raw moments");
    }

    Eigen::VectorXd nearBinary = fill;
    nearBinary.setConstant(1e-3);
    const std::string refinedHash = domain.partitionHash();
    if (
        domain.mergeBinarySiblings(nearBinary, 0.05) != 0U ||
        domain.partitionHash() != refinedHash
    )
    {
        return fail("near-binary siblings were rounded during merge");
    }

    fill.setZero();
    const std::size_t merged = domain.mergeBinarySiblings(fill, 1e-9);
    if (merged != 1U || domain.leaves().size() != 8U)
    {
        return fail("binary-equivalent siblings did not merge");
    }

    const Eigen::MatrixXd laplacian = domain.uniformLaplacian();
    if (
        laplacian.rows() != 8 ||
        !laplacian.isApprox(laplacian.transpose(), 1e-14) ||
        !(laplacian * Eigen::VectorXd::Ones(8)).isZero(1e-14)
    )
    {
        return fail("uniform face-neighbor Laplacian is not symmetric");
    }
    if (!domain.validatePartition())
    {
        return fail("final adaptive partition failed its audit");
    }

    spin2014::AdaptiveVoxelDomain repeated(source, 1, 2);
    if (repeated.partitionHash() != domain.partitionHash())
    {
        return fail("adaptive initialization is not deterministic");
    }

    std::cout << "[PASS] scenario=spin-it-2014-adaptive-domain"
              << " leaves=" << domain.leaves().size()
              << " hash=" << domain.partitionHash() << '\n';
    return EXIT_SUCCESS;
}
