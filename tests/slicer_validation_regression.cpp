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
#include <thread>
#include <vector>

#include <fcntl.h>
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


class ScopedClosedStandardDescriptors
{
public:
    ScopedClosedStandardDescriptors()
    {
        savedInput_ = ::fcntl(STDIN_FILENO, F_DUPFD_CLOEXEC, 3);
        savedOutput_ = ::fcntl(STDOUT_FILENO, F_DUPFD_CLOEXEC, 3);
        if (savedInput_ < 0 || savedOutput_ < 0)
        {
            if (savedInput_ >= 0) ::close(savedInput_);
            if (savedOutput_ >= 0) ::close(savedOutput_);
            throw std::runtime_error(
                "failed to save standard descriptors for regression"
            );
        }
        ::close(STDIN_FILENO);
        ::close(STDOUT_FILENO);
    }

    ~ScopedClosedStandardDescriptors()
    {
        if (savedInput_ >= 0)
        {
            ::dup2(savedInput_, STDIN_FILENO);
            ::close(savedInput_);
        }
        if (savedOutput_ >= 0)
        {
            ::dup2(savedOutput_, STDOUT_FILENO);
            ::close(savedOutput_);
        }
    }

    ScopedClosedStandardDescriptors(
        const ScopedClosedStandardDescriptors&
    ) = delete;
    ScopedClosedStandardDescriptors& operator=(
        const ScopedClosedStandardDescriptors&
    ) = delete;

private:
    int savedInput_ = -1;
    int savedOutput_ = -1;
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


int runSigchldOwnership(const std::string& fakeExecutable)
{
    const std::string scenario = "sigchld-ownership";
    Fixture fixture(scenario, fakeExecutable);
    struct sigaction previous{};
    struct sigaction ignored{};
    ignored.sa_handler = SIG_IGN;
    ::sigemptyset(&ignored.sa_mask);
    ignored.sa_flags = 0;
    if (::sigaction(SIGCHLD, &ignored, &previous) != 0)
    {
        return fail(scenario, "could not install SIGCHLD test disposition");
    }
    SlicerValidationReport report;
    try
    {
        report = fixture.validate();
    }
    catch (...)
    {
        ::sigaction(SIGCHLD, &previous, nullptr);
        throw;
    }
    ::sigaction(SIGCHLD, &previous, nullptr);
    if (report.passed || report.completed || report.exitSuccess)
    {
        return fail(scenario, "unsupported SIGCHLD ownership was accepted");
    }
    return pass(scenario);
}


int runClosedStdioLaunchFailure(const std::string& fakeExecutable)
{
    const std::string scenario = "closed-stdio-launch-failure";
    ScopedTemporaryDirectory executableDirectory(scenario + "_executable");
    const std::filesystem::path privateExecutable =
        executableDirectory.path() / "fake-slicer";
    std::filesystem::copy_file(fakeExecutable, privateExecutable);
    Fixture fixture(scenario, privateExecutable.string());
    ScopedEnvironment versionMode(
        "MAKE_IT_STAND_FAKE_SLICER_VERSION_MODE",
        "unlink-self"
    );

    SlicerValidationReport report;
    {
        ScopedClosedStandardDescriptors closedDescriptors;
        report = fixture.validate();
    }
    if (
        !report.versionVerified ||
        !report.attempted ||
        report.status != SlicerValidationStatus::LaunchFailed ||
        report.passed ||
        report.exitSuccess ||
        report.diagnostic.find("pre-exec stage 4") == std::string::npos
    )
    {
        return fail(
            scenario,
            "closed stdio corrupted pre-exec failure classification"
        );
    }
    return pass(scenario);
}


int runCaptureLimit(const std::string& fakeExecutable)
{
    const std::string scenario = "capture-limit";
    Fixture fixture(scenario, fakeExecutable, "noisy");
    fixture.config.maximumCaptureBytes = 128;

    SlicerValidationReport report;
    try
    {
        report = fixture.validate();
    }
    catch (const std::exception& error)
    {
        return fail(
            scenario,
            std::string("capture failure escaped the report: ") + error.what()
        );
    }
    bool temporaryArtifact = false;
    for (const std::filesystem::directory_entry& entry :
         std::filesystem::directory_iterator(fixture.temporary.path()))
    {
        if (entry.path().filename().string().find(".r8.") != std::string::npos)
        {
            temporaryArtifact = true;
        }
    }
    if (
        report.status != SlicerValidationStatus::CaptureFailed ||
        !report.versionVerified ||
        !report.attempted ||
        !report.completed ||
        !report.exitedNormally ||
        report.exitCode != 0 ||
        report.exitSuccess ||
        report.passed ||
        report.diagnostic.find("exceeds limit") == std::string::npos ||
        std::filesystem::exists(fixture.output) ||
        temporaryArtifact
    )
    {
        return fail(scenario, "capture limit was not transactional");
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


int runVersionTruncation(const std::string& fakeExecutable)
{
    const std::string scenario = "version-truncation";
    Fixture fixture(scenario, fakeExecutable);
    fixture.config.maximumDiagnosticBytes =
        std::char_traits<char>::length(EXPECTED_VERSION);
    const ScopedEnvironment mode(
        "MAKE_IT_STAND_FAKE_SLICER_VERSION_MODE",
        "long-suffix"
    );
    const SlicerValidationReport report = fixture.validate();
    if (
        report.status != SlicerValidationStatus::VersionMismatch ||
        report.versionVerified || report.attempted || report.passed ||
        report.observedVersion.size() <= fixture.config.maximumDiagnosticBytes
    )
    {
        return fail(scenario, "truncated version excerpt was used as identity");
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


int runStableInvocation(const std::string& fakeExecutable)
{
    const std::string scenario = "stable-invocation";
    Fixture fixture(scenario, fakeExecutable);
    const SlicerValidationReport first = fixture.validate();
    const SlicerValidationReport second = fixture.validate();
    if (!first.passed || !second.passed ||
        first.invocationSha256 != second.invocationSha256)
    {
        return fail(scenario, "equivalent argv fingerprint was not stable");
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
    Fixture garbage(
        scenario + "-garbage",
        fakeExecutable,
        "garbage-gcode"
    );
    const SlicerValidationReport garbageReport = garbage.validate();
    Fixture extrusionOnly(
        scenario + "-extrusion-only",
        fakeExecutable,
        "extrusion-only"
    );
    const SlicerValidationReport extrusionOnlyReport =
        extrusionOnly.validate();

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
        !layerReport.toolpath.hasDeposition ||
        garbageReport.status != SlicerValidationStatus::InvalidToolpath ||
        garbageReport.passed || garbageReport.toolpath.hasLayers ||
        garbageReport.toolpath.hasDeposition ||
        extrusionOnlyReport.status != SlicerValidationStatus::InvalidToolpath ||
        extrusionOnlyReport.passed ||
        !extrusionOnlyReport.toolpath.hasLayers ||
        extrusionOnlyReport.toolpath.hasDeposition
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

    Fixture descendantFixture(
        scenario + "-descendant",
        fakeExecutable,
        "fork-grandchild"
    );
    descendantFixture.config.timeoutMilliseconds = 100;
    const std::filesystem::path lateMarker =
        descendantFixture.profile.string() + ".late-marker";
    const SlicerValidationReport descendantReport =
        descendantFixture.validate();
    std::this_thread::sleep_for(std::chrono::milliseconds(800));
    if (descendantReport.status != SlicerValidationStatus::TimedOut ||
        !descendantReport.timedOut || descendantReport.passed ||
        std::filesystem::exists(lateMarker))
    {
        return fail(scenario, "timed-out slicer descendant survived");
    }

    Fixture orphanFixture(
        scenario + "-orphan",
        fakeExecutable,
        "orphan-grandchild"
    );
    orphanFixture.config.timeoutMilliseconds = 1000;
    const std::filesystem::path orphanMarker =
        orphanFixture.profile.string() + ".late-marker";
    const SlicerValidationReport orphanReport = orphanFixture.validate();
    std::this_thread::sleep_for(std::chrono::milliseconds(800));
    if (orphanReport.passed || std::filesystem::exists(orphanMarker))
    {
        return fail(scenario, "slicer descendant outlived leader");
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
        first.find(
            "\"output_toolpath_sha256\": \"" +
            report.outputToolpathSha256 + "\""
        ) == std::string::npos ||
        first.find(
            "\"expected_executable_sha256\": \"" +
            report.expectedExecutableSha256 + "\""
        ) == std::string::npos ||
        first.find(
            "\"version_stdout_sha256\": \"" +
            report.versionStdoutSha256 + "\""
        ) == std::string::npos ||
        temporaryLeak
    )
    {
        return fail(scenario, "stable atomic slicer JSON contract failed");
    }
    return pass(scenario);
}


int runRelativePathsAndCollisions(const std::string& fakeExecutable)
{
    const std::string scenario = "paths";
    Fixture fixture(scenario, fakeExecutable);
    const std::filesystem::path original = std::filesystem::current_path();
    std::filesystem::current_path(fixture.temporary.path());
    std::filesystem::create_directory(
        fixture.temporary.path() / "nested"
    );
    SlicerValidationReport relative;
    try
    {
        relative = manufacturing::validateWithPinnedSlicer(
            "input.stl",
            "nested/output.gcode",
            fixture.config
        );
    }
    catch (...)
    {
        std::filesystem::current_path(original);
        throw;
    }
    std::filesystem::current_path(original);
    if (
        relative.status != SlicerValidationStatus::Passed ||
        !std::filesystem::is_regular_file(
            fixture.temporary.path() / "nested/output.gcode"
        )
    )
    {
        return fail(scenario, "valid relative paths were not canonicalized");
    }

    const SlicerValidationReport inputCollision =
        manufacturing::validateWithPinnedSlicer(
            fixture.input.string(),
            fixture.input.string(),
            fixture.config
        );
    const SlicerValidationReport profileCollision =
        manufacturing::validateWithPinnedSlicer(
            fixture.input.string(),
            fixture.profile.string(),
            fixture.config
        );
    const SlicerValidationReport executableCollision =
        manufacturing::validateWithPinnedSlicer(
            fixture.input.string(),
            fixture.executable.string(),
            fixture.config
        );
    const std::filesystem::path hardLink =
        fixture.temporary.path() / "input-hardlink.gcode";
    std::filesystem::create_hard_link(fixture.input, hardLink);
    const SlicerValidationReport hardLinkCollision =
        manufacturing::validateWithPinnedSlicer(
            fixture.input.string(),
            hardLink.string(),
            fixture.config
        );
    if (
        inputCollision.status != SlicerValidationStatus::InvalidConfiguration ||
        profileCollision.status != SlicerValidationStatus::InvalidConfiguration ||
        executableCollision.status != SlicerValidationStatus::InvalidConfiguration ||
        hardLinkCollision.status != SlicerValidationStatus::InvalidConfiguration ||
        readText(fixture.input).find("solid fixture") == std::string::npos ||
        readText(fixture.profile).find("mode=success") == std::string::npos
    )
    {
        return fail(scenario, "destructive output collision was accepted");
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
    if (scenario == "version-truncation")
        return runVersionTruncation(fakeExecutable);
    if (scenario == "success") return runSuccess(fakeExecutable, false);
    if (scenario == "stable-invocation")
        return runStableInvocation(fakeExecutable);
    if (scenario == "repair") return runRepair(fakeExecutable);
    if (scenario == "invalid-toolpath")
        return runInvalidToolpath(fakeExecutable);
    if (scenario == "nonzero") return runNonzero(fakeExecutable, false);
    if (scenario == "signal") return runNonzero(fakeExecutable, true);
    if (scenario == "timeout") return runTimeout(fakeExecutable);
    if (scenario == "shell-safe") return runSuccess(fakeExecutable, true);
    if (scenario == "report-json") return runReportJson(fakeExecutable);
    if (scenario == "paths") return runRelativePathsAndCollisions(fakeExecutable);
    if (scenario == "sigchld-ownership")
        return runSigchldOwnership(fakeExecutable);
    if (scenario == "closed-stdio-launch-failure")
        return runClosedStdioLaunchFailure(fakeExecutable);
    if (scenario == "capture-limit")
        return runCaptureLimit(fakeExecutable);

    return fail(scenario, "unknown scenario");
}
