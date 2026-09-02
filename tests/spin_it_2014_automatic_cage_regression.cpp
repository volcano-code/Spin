#include "spin2014/AutomaticCageBinding.h"

#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace
{


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-it-2014-automatic-cage "
              << message << '\n';
    return EXIT_FAILURE;
}


std::shared_ptr<const voxel::VoxelDomain> cubeDomain()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(4, 4, 4);
    grid.minimum = Eigen::Vector3d::Zero();
    grid.maximum = Eigen::Vector3d::Constant(2.0);
    return voxel::VoxelDomain::create(
        grid,
        std::vector<unsigned char>(64U, 1U)
    );
}


Eigen::MatrixXd cubeVertices(double minimum, double maximum)
{
    Eigen::MatrixXd vertices(8, 3);
    for (int bit = 0; bit < 8; ++bit)
    {
        vertices(bit, 0) = (bit & 1) != 0 ? maximum : minimum;
        vertices(bit, 1) = (bit & 2) != 0 ? maximum : minimum;
        vertices(bit, 2) = (bit & 4) != 0 ? maximum : minimum;
    }
    return vertices;
}


} // namespace


int main()
{
    spin2014::AutomaticCageBindingRequest request;
    request.domain = cubeDomain();
    request.outerVertices = cubeVertices(0.0, 2.0);
    request.voidCellVertices = {cubeVertices(0.75, 1.25)};
    request.maximumHandleCount = 12;
    request.weightOptions.maxIterations = 250;
    request.weightOptions.tolerance = 1e-9;

    const spin2014::AutomaticCageBinding first =
        spin2014::buildAutomaticCageBinding(request);
    const spin2014::AutomaticCageBinding repeated =
        spin2014::buildAutomaticCageBinding(request);
    if (
        !first.weightDiagnostics.converged ||
        first.handles.size() < 8U || first.handles.size() > 12U ||
        first.outerWeights.rows() != request.outerVertices.rows() ||
        first.outerWeights.cols() !=
            static_cast<Eigen::Index>(first.handles.size()) ||
        first.voidCellWeights.size() != 1U ||
        first.voidCellWeights.front().rows() != 8 ||
        first.bindingHash.size() != 64U ||
        first.bindingHash != repeated.bindingHash ||
        !first.outerWeights.isApprox(repeated.outerWeights, 0.0)
    )
    {
        return fail("automatic BBW binding is incomplete or nondeterministic");
    }
    const double outerRowError =
        (first.outerWeights.rowwise().sum().array() - 1.0).abs().maxCoeff();
    const double cellRowError =
        (first.voidCellWeights.front().rowwise().sum().array() - 1.0)
            .abs().maxCoeff();
    if (
        outerRowError > 1e-8 || cellRowError > 1e-8 ||
        first.weightDiagnostics.maximumBoundViolation > 1e-8
    )
    {
        return fail("interpolated cage weights violate their partition");
    }
    std::cout << "[PASS] scenario=spin-it-2014-automatic-cage handles="
              << first.handles.size()
              << " hash=" << first.bindingHash << '\n';
    return EXIT_SUCCESS;
}
