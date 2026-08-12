#include "manufacturing/SlicerValidation.h"

#include <algorithm>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <system_error>
#include <vector>

#include <unistd.h>


namespace
{


using manufacturing::AutomaticRepairStatus;
using manufacturing::SlicerValidationConfig;
using manufacturing::SlicerValidationReport;
using manufacturing::SlicerValidationStatus;


constexpr const char* EXPECTED_VERSION = "PrusaSlicer 2.9.0-fake";


int fail(const std::string& scenario, const std::string& reason)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << " reason=" << reason << '\n';
    return EXIT_FAILURE;
}


int pass(const std::string& scenario)
{
    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


class ScopedTemporaryDirectory
{
public:
    explicit ScopedTemporaryDirectory(const std::string& scenario)
    {
        const std::filesystem::path base =
            std::filesystem::temp_directory_path();
        for (int attempt = 0; attempt < 100; ++attempt)
        {
            path_ = base / (
                "make_it_stand_slicer_" + scenario + "_" +
                std::to_string(static_cast<long long>(::getpid())) + "_" +
                std::to_string(attempt)
            );
            std::error_code error;
            if (std::filesystem::create_directory(path_, error))
            {
                return;
            }
            if (error)
            {
                throw std::runtime_error(
                    "failed to create slicer test directory"
                );
            }
        }
        throw std::runtime_error("could not allocate slicer test directory");
    }

    ~ScopedTemporaryDirectory()
    {
        std::error_code ignored;
        std::filesystem::remove_all(path_, ignored);
    }

    const std::filesystem::path& path() const noexcept
    {
        return path_;
    }

private:
    std::filesystem::path path_;
};


class ScopedEnvironment
{
public:
    ScopedEnvironment(const char* name, const std::string& value)
        : name_(name)
    {
        const char* old = std::getenv(name_.c_str());
        if (old != nullptr)
        {
            hadOld_ = true;
            old_ = old;
        }
        if (::setenv(name_.c_str(), value.c_str(), 1) != 0)
        {
            throw std::runtime_error("failed to set fake slicer environment");
        }
    }

    ~ScopedEnvironment()
    {
        if (hadOld_)
        {
            ::setenv(name_.c_str(), old_.c_str(), 1);
        }
        else
        {
            ::unsetenv(name_.c_str());
        }
    }

private:
    std::string name_;
    std::string old_;
    bool hadOld_ = false;
};


void writeText(const std::filesystem::path& path, const std::string& text)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output << text;
    if (!output)
    {
        throw std::runtime_error("failed to write slicer test fixture");
    }
}


std::string readText(const std::filesystem::path& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input)
    {
        throw std::runtime_error("failed to read slicer test artifact");
    }
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    );
}


struct Fixture
{
    ScopedTemporaryDirectory temporary;
    std::filesystem::path executable;
    std::filesystem::path profile;
    std::filesystem::path input;
    std::filesystem::path output;
    SlicerValidationConfig config;

    Fixture(
        const std::string& scenario,
        const std::string& fakeExecutable,
        const std::string& mode = "success",
        bool shellSafeNames = false
    )
        :
        temporary(scenario),
        executable(std::filesystem::absolute(fakeExecutable)),
        profile(
            temporary.path() /
            (shellSafeNames
                ? "profile with spaces ; dollar$ paren() quote'.ini"
                : "profile.ini")
        ),
        input(
            temporary.path() /
            (shellSafeNames
                ? "input with spaces ; dollar$ paren() quote'.stl"
                : "input.stl")
        ),
        output(
            temporary.path() /
            (shellSafeNames
                ? "output with spaces ; dollar$ paren() quote'.gcode"
                : "output.gcode")
        )
    {
        writeText(profile, "# fake exported profile\nmode=" + mode + "\n");
        writeText(
            input,
            "solid fixture\n"
            "facet normal 0 0 1\n"
            "outer loop\n"
            "vertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\n"
            "endloop\nendfacet\nendsolid fixture\n"
        );

        config.executableAbsolutePath = executable.string();
        config.expectedExecutableSha256 =
            manufacturing::computeFileSha256(executable.string());
        config.profileAbsolutePath = profile.string();
        config.expectedProfileSha256 =
            manufacturing::computeFileSha256(profile.string());
        config.expectedVersion = EXPECTED_VERSION;
        config.timeoutMilliseconds = 1500;
        config.maximumDiagnosticBytes = 512;
    }

    SlicerValidationReport validate() const
    {
        return manufacturing::validateWithPinnedSlicer(
            input.string(),
            output.string(),
            config
        );
    }
};


bool exactPinsVerified(const SlicerValidationReport& report)
{
    return
        report.configurationComplete &&
        report.executableAvailable &&
        report.profileAvailable &&
        report.executableHashVerified &&
        report.profileHashVerified;
}


int runNotConfigured()
{
    const std::string scenario = "not-configured";
    const SlicerValidationReport report =
        manufacturing::validateWithPinnedSlicer("", "", {});
    if (
        report.status != SlicerValidationStatus::NotConfigured ||
        report.configurationComplete ||
        report.versionProbeAttempted ||
        report.attempted ||
        report.completed ||
        report.passed ||
        report.automaticRepair != AutomaticRepairStatus::Unknown
    )
    {
        return fail(scenario, "opt-in configuration gate failed");
    }
    return pass(scenario);
}


int runPinMismatch(const std::string& fakeExecutable, bool executableMismatch)
{
    const std::string scenario = executableMismatch
        ? "executable-hash-mismatch"
        : "profile-hash-mismatch";
    Fixture fixture(scenario, fakeExecutable);
    if (executableMismatch)
    {
        fixture.config.expectedExecutableSha256 = std::string(64, '0');
    }
    else
    {
        fixture.config.expectedProfileSha256 = std::string(64, 'f');
    }

    const SlicerValidationReport report = fixture.validate();
    const SlicerValidationStatus expected = executableMismatch
        ? SlicerValidationStatus::ExecutableHashMismatch
        : SlicerValidationStatus::ProfileHashMismatch;
    if (
        report.status != expected ||
        report.versionProbeAttempted ||
        report.versionVerified ||
        report.attempted ||
        report.completed ||
        report.passed ||
        std::filesystem::exists(fixture.output)
    )
    {
        return fail(scenario, "pin mismatch did not fail before execution");
    }
    return pass(scenario);
}


int runVersionMismatch(const std::string& fakeExecutable)
{
    const std::string scenario = "version-mismatch";
    Fixture fixture(scenario, fakeExecutable);
    const ScopedEnvironment mode(
        "MAKE_IT_STAND_FAKE_SLICER_VERSION_MODE",
        "mismatch"
    );
    const SlicerValidationReport report = fixture.validate();
    if (
        report.status != SlicerValidationStatus::VersionMismatch ||
        !exactPinsVerified(report) ||
        !report.versionProbeAttempted ||
        report.versionVerified ||
        report.observedVersion == EXPECTED_VERSION ||
        report.attempted ||
        report.completed ||
        report.passed ||
        std::filesystem::exists(fixture.output)
    )
    {
        return fail(scenario, "version mismatch did not stop before slicing");
    }
    return pass(scenario);
}


int runSuccess(const std::string& fakeExecutable, bool shellSafe)
{
    const std::string scenario = shellSafe ? "shell-safe" : "success";
    Fixture fixture(scenario, fakeExecutable, "success", shellSafe);
    const SlicerValidationReport report = fixture.validate();
    if (
        report.status != SlicerValidationStatus::Passed ||
        !exactPinsVerified(report) ||
        !report.versionProbeAttempted ||
        !report.versionVerified ||
        !report.attempted ||
        !report.completed ||
        !report.exitedNormally ||
        !report.exitSuccess ||
        report.exitCode != 0 ||
        report.timedOut ||
        !report.passed ||
        report.automaticRepair != AutomaticRepairStatus::NoneReported ||
        report.automaticRepairMessageCount != 0 ||
        !report.toolpath.completed ||
        !report.toolpath.nonempty ||
        !report.toolpath.text ||
        !report.toolpath.hasLayers ||
        !report.toolpath.hasDeposition ||
        report.toolpath.layerCount < 2 ||
        report.toolpath.depositionCommandCount < 2 ||
        report.inputMeshSha256.size() != 64 ||
        report.outputToolpathSha256.size() != 64 ||
        report.invocationSha256.size() != 64 ||
        report.outputToolpathSha256 != report.toolpath.sha256 ||
        !std::filesystem::is_regular_file(fixture.output) ||
        readText(fixture.output).find("G1 X10.0") == std::string::npos
    )
    {
        return fail(scenario, "pinned fake slicer success contract failed");
    }
    return pass(scenario);
}


int runRepair(const std::string& fakeExecutable)
{
    const std::string scenario = "repair";
    Fixture fixture(scenario, fakeExecutable, "repair");
    const SlicerValidationReport report = fixture.validate();
    if (
        report.status != SlicerValidationStatus::AutoRepairDetected ||
        !report.attempted ||
        !report.completed ||
        !report.exitSuccess ||
        report.passed ||
        report.automaticRepair != AutomaticRepairStatus::Reported ||
        report.automaticRepairMessageCount <= 0 ||
        report.stderrExcerpt.find("repair") == std::string::npos
    )
    {
        return fail(scenario, "known repair diagnostic was not fail-closed");
    }
    return pass(scenario);
}


int runInvalidToolpath(const std::string& fakeExecutable)
{
    const std::string scenario = "invalid-toolpath";
    Fixture empty(scenario + "-empty", fakeExecutable, "empty");
    const SlicerValidationReport emptyReport = empty.validate();
    Fixture noDeposition(
        scenario + "-no-deposition",
        fakeExecutable,
        "no-deposition"
    );
    const SlicerValidationReport depositionReport = noDeposition.validate();
    Fixture noLayers(
        scenario + "-no-layers",
        fakeExecutable,
        "no-layers"
    );
    const SlicerValidationReport layerReport = noLayers.validate();

    if (
        emptyReport.status != SlicerValidationStatus::InvalidToolpath ||
        emptyReport.passed ||
        !emptyReport.toolpath.completed ||
        emptyReport.toolpath.nonempty ||
        depositionReport.status != SlicerValidationStatus::InvalidToolpath ||
        depositionReport.passed ||
        !depositionReport.toolpath.hasLayers ||
        depositionReport.toolpath.hasDeposition ||
        layerReport.status != SlicerValidationStatus::InvalidToolpath ||
        layerReport.passed ||
        layerReport.toolpath.hasLayers ||
        !layerReport.toolpath.hasDeposition
    )
    {
        return fail(scenario, "invalid G-code semantic gate failed");
    }
    return pass(scenario);
}


int runNonzero(const std::string& fakeExecutable, bool signal)
{
    const std::string scenario = signal ? "signal" : "nonzero";
    Fixture fixture(
        scenario,
        fakeExecutable,
        signal ? "signal" : "nonzero"
    );
    const SlicerValidationReport report = fixture.validate();
    if (
        report.status != SlicerValidationStatus::NonzeroExit ||
        !report.attempted ||
        !report.completed ||
        report.exitSuccess ||
        report.passed ||
        report.automaticRepair != AutomaticRepairStatus::Unknown ||
        (!signal && (!report.exitedNormally || report.exitCode != 23)) ||
        (signal &&
            (report.exitedNormally || report.terminationSignal != SIGTERM))
    )
    {
        return fail(scenario, "abnormal slicer termination was misclassified");
    }
    return pass(scenario);
}


int runTimeout(const std::string& fakeExecutable)
{
    const std::string scenario = "timeout";
    Fixture fixture(scenario, fakeExecutable, "sleep");
    fixture.config.timeoutMilliseconds = 100;
    const auto start = std::chrono::steady_clock::now();
    const SlicerValidationReport report = fixture.validate();
    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start
    ).count();
    if (
        report.status != SlicerValidationStatus::TimedOut ||
        !report.attempted ||
        !report.timedOut ||
        report.exitSuccess ||
        report.passed ||
        report.automaticRepair != AutomaticRepairStatus::Unknown ||
        elapsed > 3000
    )
    {
        return fail(scenario, "slicer timeout was not bounded and fail-closed");
    }
    return pass(scenario);
}


int runReportJson(const std::string& fakeExecutable)
{
    const std::string scenario = "report-json";
    Fixture fixture(scenario, fakeExecutable);
    const SlicerValidationReport report = fixture.validate();
    if (!report.passed)
    {
        return fail(scenario, "report fixture did not pass slicing");
    }

    const std::string first = report.toJson();
    const std::string second = report.toJson();
    const std::filesystem::path destination =
        fixture.temporary.path() / "slicer report.json";
    writeText(destination, "preexisting sentinel\n");
    manufacturing::writeSlicerValidationReport(report, destination.string());

    bool temporaryLeak = false;
    for (const auto& entry :
         std::filesystem::directory_iterator(fixture.temporary.path()))
    {
        if (entry.path().filename().string().find(".r8tmp.") !=
            std::string::npos)
        {
            temporaryLeak = true;
        }
    }

    if (
        first != second ||
        readText(destination) != first ||
        first.find("\"status\": \"passed\"") == std::string::npos ||
        first.find("\"automatic_repair\": \"none_reported\"") ==
            std::string::npos ||
        first.find("\"passed\": true") == std::string::npos ||
        first.find(report.outputToolpathSha256) == std::string::npos ||
        temporaryLeak
    )
    {
        return fail(scenario, "stable atomic slicer JSON contract failed");
    }
    return pass(scenario);
}


} // namespace


int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3)
    {
        std::cerr
            << "usage: slicer-validation-regression <scenario> "
            << "[absolute-fake-slicer]\n";
        return EXIT_FAILURE;
    }

    const std::string scenario = argv[1];
    if (scenario == "not-configured")
    {
        return runNotConfigured();
    }
    if (argc != 3)
    {
        return fail(scenario, "scenario requires the fake slicer path");
    }

    const std::string fakeExecutable = argv[2];
    if (scenario == "executable-hash-mismatch")
        return runPinMismatch(fakeExecutable, true);
    if (scenario == "profile-hash-mismatch")
        return runPinMismatch(fakeExecutable, false);
    if (scenario == "version-mismatch")
        return runVersionMismatch(fakeExecutable);
    if (scenario == "success") return runSuccess(fakeExecutable, false);
    if (scenario == "repair") return runRepair(fakeExecutable);
    if (scenario == "invalid-toolpath")
        return runInvalidToolpath(fakeExecutable);
    if (scenario == "nonzero") return runNonzero(fakeExecutable, false);
    if (scenario == "signal") return runNonzero(fakeExecutable, true);
    if (scenario == "timeout") return runTimeout(fakeExecutable);
    if (scenario == "shell-safe") return runSuccess(fakeExecutable, true);
    if (scenario == "report-json") return runReportJson(fakeExecutable);

    return fail(scenario, "unknown scenario");
}
