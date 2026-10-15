#include "spin2014/AdaptiveVoxelSpinSolver.h"
#include "spin2014/DualMaterialSpinSolver.h"

#include "voxel/VoxelDomain.h"
#include "source_surface_test_utils.h"
#include "../experiments/spin_it_reproduction/strict/dual_material_geometry_export.h"

#include <Eigen/Core>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace
{


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-it-2014-dual-material "
              << message << '\n';
    return EXIT_FAILURE;
}


std::shared_ptr<const voxel::VoxelDomain> fixture()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(4, 4, 4);
    grid.minimum = Eigen::Vector3d(-1.0, -1.0, 0.0);
    grid.maximum = Eigen::Vector3d(1.0, 1.0, 2.0);
    std::vector<unsigned char> material(64U, 1U);
    std::vector<unsigned char> carveable(64U, 0U);
    for (int z = 1; z < 3; ++z)
    {
        for (int y = 1; y < 3; ++y)
        {
            for (int x = 1; x < 3; ++x)
            {
                carveable[static_cast<std::size_t>(x + 4 * (y + 4 * z))] = 1U;
            }
        }
    }
    return voxel::VoxelDomain::create(grid, material, {}, carveable);
}


std::shared_ptr<const voxel::VoxelDomain> asymmetricFixture()
{
    // Frozen before solving: a 4x4x6 body with an offset upper-side extension.
    // The carveable core remains one cell away from the exterior. This creates
    // nonzero cross moments without altering or repairing the output geometry.
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(7, 6, 6);
    grid.minimum = Eigen::Vector3d(-1.75, -1.5, 0.0);
    grid.maximum = Eigen::Vector3d(1.75, 1.5, 3.0);
    std::vector<unsigned char> material(7U * 6U * 6U, 0U);
    std::vector<unsigned char> carveable(material.size(), 0U);
    const auto index = [&grid](int x, int y, int z) {
        return static_cast<std::size_t>(x + grid.dimensions.x() *
            (y + grid.dimensions.y() * z));
    };
    for (int z = 0; z < 6; ++z)
        for (int y = 1; y <= 4; ++y)
            for (int x = 1; x <= 4; ++x)
                material[index(x, y, z)] = 1U;
    for (int z = 1; z <= 5; ++z)
        for (int y = 3; y <= 5; ++y)
            for (int x = 4; x <= 6; ++x)
                material[index(x, y, z)] = 1U;
    for (int z = 1; z <= 4; ++z)
        for (int y = 2; y <= 3; ++y)
            for (int x = 2; x <= 3; ++x)
                carveable[index(x, y, z)] = 1U;
    return voxel::VoxelDomain::create(grid, material, {}, carveable);
}


std::string arrayMember(const std::string& json, const std::string& name)
{
    const std::size_t key = json.find("\"" + name + "\"");
    if (key == std::string::npos) throw std::runtime_error("missing Sumone asset field: " + name);
    const std::size_t begin = json.find('[', key);
    if (begin == std::string::npos) throw std::runtime_error("missing Sumone array: " + name);
    int depth = 0;
    for (std::size_t end = begin; end < json.size(); ++end)
    {
        if (json[end] == '[') ++depth;
        else if (json[end] == ']' && --depth == 0)
            return json.substr(begin, end - begin + 1);
    }
    throw std::runtime_error("unterminated Sumone array: " + name);
}


std::vector<double> parseNumbers(const std::string& text)
{
    std::vector<double> result;
    const char* cursor = text.c_str();
    while (*cursor)
    {
        if ((*cursor >= '0' && *cursor <= '9') || *cursor == '-' || *cursor == '+')
        {
            char* end = nullptr;
            const double value = std::strtod(cursor, &end);
            if (end == cursor) throw std::runtime_error("invalid numeric Sumone asset field");
            result.push_back(value);
            cursor = end;
        }
        else ++cursor;
    }
    return result;
}


std::vector<unsigned char> decodeRle(const std::string& json, const std::string& name,
                                     std::size_t expected)
{
    const std::vector<double> values = parseNumbers(arrayMember(json, name));
    if (values.size() % 2U != 0U) throw std::runtime_error("odd Sumone RLE pair count");
    std::vector<unsigned char> result;
    result.reserve(expected);
    for (std::size_t i = 0; i < values.size(); i += 2U)
    {
        const int value = static_cast<int>(values[i]);
        const int count = static_cast<int>(values[i + 1U]);
        if ((value != 0 && value != 1) || count <= 0 ||
            static_cast<double>(value) != values[i] || static_cast<double>(count) != values[i + 1U] ||
            static_cast<std::size_t>(count) > expected - result.size())
            throw std::runtime_error("invalid Sumone binary RLE pair");
        result.insert(result.end(), static_cast<std::size_t>(count), static_cast<unsigned char>(value));
    }
    if (result.size() != expected) throw std::runtime_error("Sumone RLE cell count mismatch");
    return result;
}


std::shared_ptr<const voxel::VoxelDomain> sumoneFigureFiveFixture()
{
    std::ifstream input("experiments/spin_it_reproduction/spin2014/paper-assets/sumone-dual-density.vox.json");
    if (!input) throw std::runtime_error("cannot open project-generated Sumone Figure 5 asset");
    const std::string json((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    const auto dimensions = parseNumbers(arrayMember(json, "dimensions"));
    const auto minimum = parseNumbers(arrayMember(json, "minimum"));
    const auto maximum = parseNumbers(arrayMember(json, "maximum"));
    if (dimensions.size() != 3U || minimum.size() != 3U || maximum.size() != 3U)
        throw std::runtime_error("invalid Sumone grid metadata");
    voxel::GridSpec grid;
    for (int axis = 0; axis < 3; ++axis)
    {
        grid.dimensions[axis] = static_cast<int>(dimensions[static_cast<std::size_t>(axis)]);
        grid.minimum[axis] = minimum[static_cast<std::size_t>(axis)];
        grid.maximum[axis] = maximum[static_cast<std::size_t>(axis)];
    }
    const std::size_t count = static_cast<std::size_t>(grid.dimensions.prod());
    const auto material = decodeRle(json, "material_rle", count);
    const auto carveable = decodeRle(json, "carveable_rle", count);
    return voxel::VoxelDomain::create(grid, material, {}, carveable);
}


spin2014::AdaptiveVoxelSpinRequest request()
{
    spin2014::AdaptiveVoxelSpinRequest result;
    result.domain = fixture();
    result.mode = spin::SpinMode::StaticBalance;
    result.initialLevel = 1;
    result.maximumLevel = 2;
    result.maximumAdaptivePasses = 3;
    result.initialFill = 0.5;
    result.laplacianWeight = 1e-3;
    result.slqp.maximumIterations = 160;
    return result;
}


} // namespace


int main(int argc,char** argv) try
{
    std::filesystem::path exportDirectory;
    if(argc==3&&std::string(argv[1])=="--export-directory") {
        exportDirectory=argv[2];
        if(!std::filesystem::create_directory(exportDirectory))return fail("export directory must be new");
        auto registration=dual_geometry_export::file(exportDirectory/"protocol.json");
        registration<<"{\"classification\":\"SYNTHETIC_REGRESSION\",\"protocol\":\"existing frozen dual-material regression inputs; no parameter changes\",\"original_paper_case\":false}\n";
    } else if(argc!=1)return fail("expected --export-directory DIR");
    const spin2014::AdaptiveVoxelSpinResult oneMaterial =
        spin2014::solveAdaptiveVoxelSpin(request());
    spin2014::DualMaterialSpinRequest equalDensity;
    equalDensity.hollowing = request();
    equalDensity.boundaryDensity = 1.0;
    equalDensity.interiorDensity = 1.0;
    const spin2014::DualMaterialSpinResult equal =
        spin2014::solveDualMaterialSpin(equalDensity);
    if (
        oneMaterial.status != spin2014::AdaptiveSpinStatus::Converged ||
        equal.hollowing.status != spin2014::AdaptiveSpinStatus::Converged ||
        equal.hollowing.outputHash != oneMaterial.outputHash ||
        equal.hollowing.binaryVoidMask != oneMaterial.binaryVoidMask ||
        !equal.hollowing.finalMoments.values().isApprox(
            oneMaterial.finalMoments.values(), 0.0
        )
    )
    {
        return fail("equal-density limit does not reduce exactly to hollowing");
    }

    spin2014::DualMaterialSpinRequest denseInterior = equalDensity;
    denseInterior.boundaryDensity = 1.17;
    denseInterior.interiorDensity = 8.1;
    const spin2014::DualMaterialSpinResult dense =
        spin2014::solveDualMaterialSpin(denseInterior);

    // Independent material-label integration: integrate each emitted base
    // grid cell as an analytic cuboid instead of subtracting one solver moment
    // aggregate from another. This catches swapped densities, wrong labels,
    // and a non-geometric interiorContribution.
    surface_test::Ten independentlyIntegratedBoundary =
        surface_test::Ten::Zero();
    surface_test::Ten independentlyIntegratedInterior =
        surface_test::Ten::Zero();
    const auto& domain = *denseInterior.hollowing.domain;
    const auto& grid = domain.grid();
    const Eigen::Vector3d spacing = domain.spacing();
    for (std::size_t index = 0; index < dense.materialLabels.size(); ++index)
    {
        const unsigned char label = dense.materialLabels[index];
        if (label == 0U)
        {
            continue;
        }
        const Eigen::Vector3i coordinates = domain.cellCoordinates(index);
        const Eigen::Vector3d lower = grid.minimum +
            coordinates.cast<double>().cwiseProduct(spacing);
        const Eigen::Vector3d upper = lower + spacing;
        const surface_test::Ten cell = surface_test::boxMoments(lower, upper);
        if (label == 1U)
        {
            independentlyIntegratedBoundary +=
                denseInterior.boundaryDensity * cell;
        }
        else if (label == 2U)
        {
            independentlyIntegratedInterior +=
                denseInterior.interiorDensity * cell;
        }
        else
        {
            return fail("dual-material label is outside {0,1,2}");
        }
    }
    if (
        (dense.boundaryContribution.values() -
         independentlyIntegratedBoundary).cwiseAbs().maxCoeff() > 1e-11 ||
        (dense.interiorContribution.values() -
         independentlyIntegratedInterior).cwiseAbs().maxCoeff() > 1e-11 ||
        (dense.hollowing.finalMoments.values() -
         (independentlyIntegratedBoundary + independentlyIntegratedInterior))
                .cwiseAbs().maxCoeff() > 1e-11
    )
    {
        return fail("independent label-based density-weighted moments disagree");
    }
    if (
        dense.hollowing.status != spin2014::AdaptiveSpinStatus::Converged ||
        !dense.hollowing.finalMoments.values().isApprox(
            (dense.boundaryContribution + dense.interiorContribution).values(),
            1e-11
        ) ||
        dense.materialLabels.size() != denseInterior.hollowing.domain->cellCount() ||
        dense.outputHash.size() != 64U
    )
    {
        return fail("density-weighted material contributions are inconsistent");
    }
    bool sawBoundary = false;
    bool sawInterior = false;
    for (unsigned char label : dense.materialLabels)
    {
        sawBoundary = sawBoundary || label == 1U;
        sawInterior = sawInterior || label == 2U;
        if (label > 2U)
        {
            return fail("dual-material label is outside {0,1,2}");
        }
    }
    if (!sawBoundary || !sawInterior)
    {
        return fail("dual-material output omitted a retained material class");
    }

    // The static paper formulation has no free transverse angle, but it still
    // uses the published Givens objective/constraint evaluation path.  This
    // deliberately unequal-density fixture verifies that the dual-material
    // wrapper no longer rejects PAPER_GIVENS or silently changes formulation.
    spin2014::DualMaterialSpinRequest paperDense = denseInterior;
    paperDense.hollowing.formulation = spin2014::Spin2014Formulation::PaperGivens;
    const spin2014::DualMaterialSpinResult paper =
        spin2014::solveDualMaterialSpin(paperDense);
    if (
        paper.hollowing.status != spin2014::AdaptiveSpinStatus::Converged ||
        paper.hollowing.formulation != spin2014::Spin2014Formulation::PaperGivens ||
        !paper.hollowing.finalMoments.values().isApprox(
            (paper.boundaryContribution + paper.interiorContribution).values(),
            1e-11
        ) ||
        paper.materialLabels != dense.materialLabels
    )
    {
        return fail("paper-Givens dual-material product is inconsistent");
    }

    spin2014::AdaptiveVoxelSpinRequest asymmetricHollowing = request();
    asymmetricHollowing.domain = asymmetricFixture();
    asymmetricHollowing.maximumLevel = 3;
    asymmetricHollowing.maximumAdaptivePasses = 4;
    spin2014::DualMaterialSpinRequest asymmetricRequest;
    asymmetricRequest.hollowing = asymmetricHollowing;
    asymmetricRequest.boundaryDensity = 1.25;
    asymmetricRequest.interiorDensity = 6.75;
    const spin2014::DualMaterialSpinResult asymmetric =
        spin2014::solveDualMaterialSpin(asymmetricRequest);
    if (asymmetric.hollowing.status != spin2014::AdaptiveSpinStatus::RelaxationDidNotConverge)
    {
        return fail("frozen asymmetric negative-control status changed: " +
            asymmetric.hollowing.diagnostic);
    }
    std::cout << "[EXPECTED_NEGATIVE] scenario=asymmetric-dual-material"
              << " status=" << spin2014::adaptiveSpinStatusToken(asymmetric.hollowing.status)
              << " input_hash=" << asymmetric.hollowing.inputHash
              << " diagnostic=" << asymmetric.hollowing.diagnostic << '\n';

    // This is the project's reconstructed Figure 5 Sumone asset. Its topology,
    // volume and transform are frozen by the checked-in occupancy hash; the
    // published material densities are used, but the asset is not author-exact.
    spin2014::AdaptiveVoxelSpinRequest sumoneHollowing;
    sumoneHollowing.domain = sumoneFigureFiveFixture();
    sumoneHollowing.mode = spin::SpinMode::SpinningTop;
    sumoneHollowing.initialLevel = 2;
    sumoneHollowing.maximumLevel = 4;
    sumoneHollowing.maximumAdaptivePasses = 8;
    sumoneHollowing.initialFill = 0.5;
    sumoneHollowing.laplacianWeight = 1e-3;
    sumoneHollowing.binaryConstraintTolerance = 0.005;
    sumoneHollowing.slqp.maximumIterations = 240;
    const spin2014::AdaptiveVoxelSpinResult sumoneOneMaterial =
        spin2014::solveAdaptiveVoxelSpin(sumoneHollowing);
    spin2014::DualMaterialSpinRequest sumoneEqualRequest;
    sumoneEqualRequest.hollowing = sumoneHollowing;
    sumoneEqualRequest.boundaryDensity = 1.0;
    sumoneEqualRequest.interiorDensity = 1.0;
    const spin2014::DualMaterialSpinResult sumoneEqual =
        spin2014::solveDualMaterialSpin(sumoneEqualRequest);
    if (sumoneOneMaterial.status != spin2014::AdaptiveSpinStatus::Converged ||
        sumoneEqual.hollowing.status != spin2014::AdaptiveSpinStatus::Converged ||
        sumoneEqual.hollowing.outputHash != sumoneOneMaterial.outputHash ||
        sumoneEqual.hollowing.binaryVoidMask != sumoneOneMaterial.binaryVoidMask)
    {
        return fail("Sumone equal-density limit does not reproduce one-material optimization");
    }
    spin2014::DualMaterialSpinRequest sumoneDualRequest;
    sumoneDualRequest.hollowing = sumoneHollowing;
    sumoneDualRequest.boundaryDensity = 1.17;
    sumoneDualRequest.interiorDensity = 8.1;
    const spin2014::DualMaterialSpinResult sumoneDual =
        spin2014::solveDualMaterialSpin(sumoneDualRequest);
    if (sumoneDual.hollowing.status != spin2014::AdaptiveSpinStatus::Converged ||
        sumoneDual.materialLabels.size() != sumoneHollowing.domain->cellCount())
    {
        return fail("reconstructed Sumone dual-density optimization did not produce labels: " +
            sumoneDual.hollowing.diagnostic);
    }
    surface_test::Ten sumoneBoundaryMoments = surface_test::Ten::Zero();
    surface_test::Ten sumoneInteriorMoments = surface_test::Ten::Zero();
    bool sumoneSawBoundary = false;
    bool sumoneSawInterior = false;
    std::size_t sumoneBoundaryCellCount = 0;
    std::size_t sumoneInteriorCellCount = 0;
    const Eigen::Vector3d sumoneSpacing = sumoneHollowing.domain->spacing();
    for (std::size_t cellIndex = 0; cellIndex < sumoneDual.materialLabels.size(); ++cellIndex)
    {
        const unsigned char label = sumoneDual.materialLabels[cellIndex];
        if (label == 0U) continue;
        const Eigen::Vector3i cell = sumoneHollowing.domain->cellCoordinates(cellIndex);
        const Eigen::Vector3d lower = sumoneHollowing.domain->grid().minimum +
            cell.cast<double>().cwiseProduct(sumoneSpacing);
        const surface_test::Ten moments = surface_test::boxMoments(lower, lower + sumoneSpacing);
        if (label == 1U)
        {
            sumoneBoundaryMoments += 1.17 * moments;
            sumoneSawBoundary = true;
            ++sumoneBoundaryCellCount;
        }
        else if (label == 2U)
        {
            sumoneInteriorMoments += 8.1 * moments;
            sumoneSawInterior = true;
            ++sumoneInteriorCellCount;
        }
        else return fail("Sumone material label is outside {0,1,2}");
    }
    const surface_test::Ten sumoneExpected = sumoneBoundaryMoments + sumoneInteriorMoments;
    const double sumoneMomentError =
        (sumoneDual.hollowing.finalMoments.values() - sumoneExpected).cwiseAbs().maxCoeff();
    const double sumoneBoundaryError =
        (sumoneDual.boundaryContribution.values() - sumoneBoundaryMoments).cwiseAbs().maxCoeff();
    const double sumoneInteriorError =
        (sumoneDual.interiorContribution.values() - sumoneInteriorMoments).cwiseAbs().maxCoeff();
    const double sumoneCompositionError =
        (sumoneDual.hollowing.finalMoments.values() -
         (sumoneDual.boundaryContribution + sumoneDual.interiorContribution).values())
            .cwiseAbs().maxCoeff();
    // On this frozen structural counterpart and locked level-4 profile, the
    // optimizer removes every carveable cell. Preserve that outcome as a
    // negative paper-case result: the independent material accounting is
    // consistent, but this run does not demonstrate a retained dense insert.
    if (!sumoneSawBoundary || sumoneSawInterior || sumoneInteriorCellCount != 0U ||
        sumoneBoundaryCellCount != 192U || sumoneMomentError > 1e-11 ||
        !sumoneDual.hollowing.finalMoments.values().isApprox(
            (sumoneDual.boundaryContribution + sumoneDual.interiorContribution).values(), 1e-11))
    {
        return fail("independent Sumone label integration disagrees: boundary=" +
            std::to_string(sumoneBoundaryError) + " interior=" + std::to_string(sumoneInteriorError) +
            " total=" + std::to_string(sumoneMomentError) +
            " composition=" + std::to_string(sumoneCompositionError) +
            " boundary_cells=" + std::to_string(sumoneBoundaryCellCount) +
            " interior_cells=" + std::to_string(sumoneInteriorCellCount) +
            " saw_boundary=" + std::to_string(sumoneSawBoundary) +
            " saw_interior=" + std::to_string(sumoneSawInterior));
    }

    if(!exportDirectory.empty()) {
        dual_geometry_export::exportCase(exportDirectory/"dense_regression",denseInterior,dense);
        dual_geometry_export::exportCase(exportDirectory/"sumone_reconstructed",sumoneDualRequest,sumoneDual);
    }
    std::cout << "[PASS] scenario=spin-it-2014-dual-material"
              << " mass=" << dense.hollowing.finalMoments.mass()
              << " paper_givens=validated"
              << " independent_label_moments=validated"
              << " sumone_density_ratio=" << (8.1 / 1.17)
              << " sumone_independent_moment_error=" << sumoneMomentError
              << " sumone_fixed_shell_cells=" << sumoneBoundaryCellCount
              << " sumone_retained_dense_cells=" << sumoneInteriorCellCount
              << " sumone_scientific_hash=" << sumoneDual.outputHash
              << " sumone_geometry_hash=" << sumoneDual.hollowing.geometryHash
              << " hash=" << dense.outputHash
              << '\n';
    return EXIT_SUCCESS;
}
catch(const std::exception& error){return fail(error.what());}
