#include "api/BalanceModel.h"
#include "stability/SupportRegion.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>


namespace
{


static_assert(
    std::is_copy_constructible<api::BalanceResult>::value,
    "BalanceResult snapshots must remain cheaply copyable"
);
static_assert(
    !std::is_move_constructible<api::BalanceResult>::value,
    "BalanceResult move must not expose an empty moved-from snapshot"
);
static_assert(
    !std::is_move_assignable<api::BalanceResult>::value,
    "BalanceResult move assignment must not expose empty snapshots"
);
static_assert(
    !std::is_move_constructible<api::BalanceModel>::value,
    "BalanceModel move must not expose a null moved-from facade"
);
static_assert(
    !std::is_move_assignable<api::BalanceModel>::value,
    "BalanceModel move assignment must not expose a null facade"
);


void makeCube(Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces)
{
    vertices.resize(8, 3);
    vertices <<
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0,
        0, 0, 1,
        1, 0, 1,
        1, 1, 1,
        0, 1, 1;
    faces.resize(12, 3);
    faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
}


template <typename Function>
bool throws(Function function)
{
    try
    {
        function();
    }
    catch (const std::exception&)
    {
        return true;
    }
    return false;
}


int fail(const std::string& reason)
{
    std::cerr << "[FAIL] scenario=r9-balance-model reason="
              << reason << '\n';
    return EXIT_FAILURE;
}


} // namespace


int main()
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeCube(vertices, faces);

    if (!throws([&]
    {
        Eigen::MatrixXd invalid = vertices;
        invalid(0, 0) = std::numeric_limits<double>::quiet_NaN();
        api::BalanceModel model(invalid, faces);
    }))
    {
        return fail("non-finite source vertices were accepted");
    }

    api::BalanceModelConfig modelConfig;
    modelConfig.voxelResolution = 8;
    api::BalanceModel model(vertices, faces, modelConfig);

    Eigen::MatrixXd unorderedSupport(7, 3);
    unorderedSupport <<
        1, 0, 1,
        0, 0, 0,
        0.5, 0, 0.5,
        1, 0, 0,
        0, 0, 1,
        1, 0, 1,
        0, 0, 0;
    const stability::SupportRegion canonical =
        stability::SupportRegion::fromPolygon(
            unorderedSupport,
            Eigen::Vector3d(0.0, -5.0, 0.0)
        );
    if (
        canonical.supportPolygon2D().rows() != 4 ||
        canonical.basisU().cross(canonical.basisV()).dot(
            canonical.gravity()
        ) < 1.0 - 1.0e-12
    )
    {
        return fail("explicit support hull was not canonical CCW");
    }

    Eigen::MatrixXi invalidFaces = faces;
    invalidFaces(0, 0) = static_cast<int>(vertices.rows());
    if (!throws([&]
    {
        api::BalanceModel invalidModel(vertices, invalidFaces, modelConfig);
    }))
    {
        return fail("out-of-range face index was accepted");
    }

    Eigen::MatrixXd badSupport(3, 3);
    badSupport << 0, 0, 0, 1, 0, 0, 0, 0.1, 1;
    if (!throws([&] { model.setSupportPolygon(badSupport, 0.0); }))
    {
        return fail("non-coplanar support was accepted");
    }

    Eigen::MatrixXd support(4, 3);
    support << 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1;
    Eigen::MatrixXd fictionalSupport = support;
    fictionalSupport.row(0) << 0.25, 0.0, 0.25;
    if (!throws([&]
    {
        model.setSupportPolygon(fictionalSupport, 0.0);
    }))
    {
        return fail("support point not present in true mesh contacts was accepted");
    }
    model.setSupportPolygon(support, 0.0);
    if (!throws([&]
    {
        model.setGravity(Eigen::Vector3d(0.0, 0.0, -1.0));
    }))
    {
        return fail("gravity changed after a target was established");
    }
    if (!throws([&]
    {
        model.addHandle(
            Eigen::Vector3d(0.5, 0.5, 0.5),
            2.0 * Eigen::Matrix3d::Identity()
        );
    }))
    {
        return fail("improper handle rotation was accepted");
    }

    api::BalanceModel badHandleModel(vertices, faces, modelConfig);
    badHandleModel.setSupportPolygon(support, 0.0);
    badHandleModel.addHandle(Eigen::Vector3d(-10.0, -10.0, -10.0));
    if (!throws([&] { badHandleModel.optimize(); }))
    {
        return fail("non-embeddable user handle was accepted");
    }

    api::BalanceModel conflictingHandleModel(vertices, faces, modelConfig);
    conflictingHandleModel.setSupportPolygon(support, 0.0);
    conflictingHandleModel.addHandle(Eigen::Vector3d(0.0, 0.0, 0.0));
    if (!throws([&] { conflictingHandleModel.optimize(); }))
    {
        return fail("user handle overlapping the fixed support was accepted");
    }

    api::BalanceRunConfig runConfig;
    runConfig.optimizer.maximumCarvingBlocks = 2;
    runConfig.optimizer.maximumDeformationIterationsPerBlock = 4;
    const api::BalanceResult first = model.optimize(runConfig);
    const api::BalanceResult second = model.optimize(runConfig);
    if (
        !first.succeeded() ||
        first.statusToken() != "stable_after_carving" ||
        first.statusToken() != second.statusToken() ||
        first.result().geometryHash != second.result().geometryHash ||
        first.diagnosticsJson() != second.diagnosticsJson() ||
        first.outerVertices().rows() != vertices.rows() ||
        first.outerFaces().rows() != faces.rows() ||
        first.diagnosticsJson().find(
            "make-it-stand-r9-balance-result-v1"
        ) == std::string::npos ||
        first.diagnosticsJson().find("\"manufacturing\"") ==
            std::string::npos
    )
    {
        return fail("stable cube run was not deterministic and successful");
    }

    const std::filesystem::path outputRoot =
        std::filesystem::temp_directory_path() /
        "make-it-stand-r9-balance-model-regression";
    std::error_code cleanupError;
    std::filesystem::remove_all(outputRoot, cleanupError);
    std::filesystem::create_directories(outputRoot);
    const std::filesystem::path printable = outputRoot / "printable.stl";
    const std::filesystem::path outer = outputRoot / "outer.obj";
    const std::filesystem::path inner = outputRoot / "inner.obj";
    const std::filesystem::path diagnostics = outputRoot / "diagnostics.json";
    first.savePrintableMesh(printable.string());
    first.saveOuterMesh(outer.string());
    first.saveInnerMesh(inner.string());
    first.saveDiagnostics(diagnostics.string());
    std::ifstream diagnosticInput(diagnostics, std::ios::binary);
    const std::string savedDiagnostics{
        std::istreambuf_iterator<char>(diagnosticInput),
        std::istreambuf_iterator<char>()
    };
    if (
        !std::filesystem::is_regular_file(printable) ||
        std::filesystem::file_size(printable) == 0U ||
        !std::filesystem::is_regular_file(outer) ||
        !std::filesystem::is_regular_file(inner) ||
        savedDiagnostics != first.diagnosticsJson()
    )
    {
        return fail("value-owned result artifacts were not exported atomically");
    }

    const std::filesystem::path symlinkSentinel =
        outputRoot / "symlink-sentinel.txt";
    const std::filesystem::path symlinkOutput =
        outputRoot / "diagnostics-link.json";
    {
        std::ofstream sentinel(
            symlinkSentinel,
            std::ios::binary | std::ios::trunc
        );
        sentinel << "do-not-follow";
    }
    std::filesystem::create_symlink(
        symlinkSentinel.filename(),
        symlinkOutput
    );
    first.saveDiagnostics(symlinkOutput.string());
    std::ifstream symlinkSentinelInput(symlinkSentinel, std::ios::binary);
    const std::string symlinkSentinelContents{
        std::istreambuf_iterator<char>(symlinkSentinelInput),
        std::istreambuf_iterator<char>()
    };
    std::ifstream symlinkOutputInput(symlinkOutput, std::ios::binary);
    const std::string symlinkOutputContents{
        std::istreambuf_iterator<char>(symlinkOutputInput),
        std::istreambuf_iterator<char>()
    };
    if (
        symlinkSentinelContents != "do-not-follow" ||
        std::filesystem::is_symlink(symlinkOutput) ||
        symlinkOutputContents != first.diagnosticsJson()
    )
    {
        return fail("atomic diagnostics export followed a destination symlink");
    }

    const std::filesystem::path concurrent =
        outputRoot / "concurrent-diagnostics.json";
    std::vector<std::exception_ptr> writeFailures(12U);
    std::vector<std::thread> writers;
    writers.reserve(writeFailures.size());
    for (std::size_t index = 0; index < writeFailures.size(); ++index)
    {
        writers.emplace_back([&, index]
        {
            try
            {
                first.saveDiagnostics(concurrent.string());
            }
            catch (...)
            {
                writeFailures[index] = std::current_exception();
            }
        });
    }
    for (std::thread& writer : writers)
    {
        writer.join();
    }
    for (const std::exception_ptr& failure : writeFailures)
    {
        if (failure)
        {
            return fail("concurrent atomic diagnostics export collided");
        }
    }
    std::ifstream concurrentInput(concurrent, std::ios::binary);
    const std::string concurrentContents{
        std::istreambuf_iterator<char>(concurrentInput),
        std::istreambuf_iterator<char>()
    };
    if (concurrentContents != first.diagnosticsJson())
    {
        return fail("concurrent atomic diagnostics output was corrupted");
    }
    for (const auto& entry : std::filesystem::directory_iterator(outputRoot))
    {
        if (
            entry.path().filename().string().find(".make_it_stand_r9.") !=
            std::string::npos
        )
        {
            return fail("atomic diagnostics export leaked a temporary file");
        }
    }

    const std::filesystem::path blockedTarget = outputRoot / "blocked-target";
    std::filesystem::create_directory(blockedTarget);
    if (!throws([&] { first.saveDiagnostics(blockedTarget.string()); }))
    {
        return fail("atomic diagnostics export replaced a directory target");
    }
    if (!std::filesystem::is_directory(blockedTarget))
    {
        return fail("failed atomic diagnostics export damaged its target");
    }
    for (const auto& entry : std::filesystem::directory_iterator(outputRoot))
    {
        if (
            entry.path().filename().string().find(".make_it_stand_r9.") !=
            std::string::npos
        )
        {
            return fail("failed atomic diagnostics export leaked a temporary");
        }
    }
    std::filesystem::remove_all(outputRoot, cleanupError);

    api::BalanceModelConfig unsafeModelConfig = modelConfig;
    unsafeModelConfig.voxelResolution = 6;
    api::BalanceModel unsafe(vertices, faces, unsafeModelConfig);
    unsafe.setSupportPolygon(support, 1.4);
    api::BalanceRunConfig unsafeRun;
    unsafeRun.optimizer.maximumCarvingBlocks = 1;
    unsafeRun.optimizer.maximumDeformationIterationsPerBlock = 1;
    const api::BalanceResult rejected = unsafe.optimize(unsafeRun);
    const std::filesystem::path guarded = outputRoot / "guarded.stl";
    std::filesystem::create_directories(outputRoot);
    {
        std::ofstream sentinel(guarded, std::ios::binary | std::ios::trunc);
        sentinel << "preserve-me";
    }
    if (
        rejected.manufacturingReady() ||
        !throws([&] { rejected.savePrintableMesh(guarded.string()); })
    )
    {
        return fail("non-ready result bypassed the sealed R8 exporter");
    }
    std::ifstream sentinelInput(guarded, std::ios::binary);
    const std::string sentinelContents{
        std::istreambuf_iterator<char>(sentinelInput),
        std::istreambuf_iterator<char>()
    };
    if (sentinelContents != "preserve-me")
    {
        return fail("failed guarded export damaged its existing target");
    }
    std::filesystem::remove_all(outputRoot, cleanupError);

    std::cout << "[PASS] scenario=r9-balance-model status="
              << first.statusToken()
              << " manufacturing_ready="
              << first.manufacturingReady()
              << " geometry_hash=" << first.result().geometryHash
              << '\n';
    return EXIT_SUCCESS;
}
