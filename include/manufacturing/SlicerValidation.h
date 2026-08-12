#pragma once

#include <cstdint>
#include <string>
#include <vector>


namespace manufacturing
{


/** Terminal state of one explicitly configured external-slicer validation. */
enum class SlicerValidationStatus
{
    NotConfigured,
    InvalidConfiguration,
    InputUnavailable,
    ExecutableUnavailable,
    ProfileUnavailable,
    ExecutableHashMismatch,
    ProfileHashMismatch,
    VersionProbeFailed,
    VersionMismatch,
    LaunchFailed,
    TimedOut,
    NonzeroExit,
    CaptureFailed,
    AutoRepairUnknown,
    AutoRepairDetected,
    InvalidToolpath,
    InputChanged,
    OutputCommitFailed,
    Passed
};


/**
 * Auto-repair evidence is deliberately four-state.
 *
 * `NoneReported` means that the pinned adapter found no known PrusaSlicer
 * repair diagnostic in the captured streams.  It is not physical evidence and
 * it must never be synthesized from an absent or failed process.
 */
enum class AutomaticRepairStatus
{
    Unsupported,
    Unknown,
    NoneReported,
    Reported
};


/** Exact opt-in pins for the one supported PrusaSlicer command profile. */
struct SlicerValidationConfig
{
    /** Must be an absolute path; PATH is never searched. */
    std::string executableAbsolutePath;

    /** Lower/upper-case hexadecimal is accepted; the observed digest is lower-case. */
    std::string expectedExecutableSha256;

    /** Must be an absolute path to one self-contained exported PrusaSlicer config. */
    std::string profileAbsolutePath;
    std::string expectedProfileSha256;

    /** Exact trimmed output of `<executable> --version`. */
    std::string expectedVersion;

    /** Applied independently to the version probe and slicing process. */
    int timeoutMilliseconds = 300000;

    /** Maximum bytes retained in each human-readable stdout/stderr excerpt. */
    std::uintmax_t maximumDiagnosticBytes = 4096;

    /** Maximum bytes accepted from each complete captured process stream. */
    std::uintmax_t maximumCaptureBytes = 16U * 1024U * 1024U;
};


/** Minimal semantic inspection of a generated FFF G-code file. */
struct GcodeToolpathAudit
{
    bool completed = false;
    bool nonempty = false;
    bool text = false;
    bool hasLayers = false;
    bool hasDeposition = false;

    std::uintmax_t bytes = 0;
    std::uintmax_t lineCount = 0;
    std::uintmax_t layerCount = 0;
    std::uintmax_t depositionCommandCount = 0;
    std::string sha256;
};


/** Reproducible evidence from one pinned PrusaSlicer smoke test. */
struct SlicerValidationReport
{
    SlicerValidationStatus status =
        SlicerValidationStatus::NotConfigured;
    AutomaticRepairStatus automaticRepair =
        AutomaticRepairStatus::Unknown;

    bool configurationComplete = false;
    bool executableAvailable = false;
    bool profileAvailable = false;
    bool executableHashVerified = false;
    bool profileHashVerified = false;
    bool versionProbeAttempted = false;
    bool versionVerified = false;
    bool attempted = false;
    bool completed = false;
    bool exitedNormally = false;
    bool exitSuccess = false;
    bool timedOut = false;
    bool passed = false;

    int exitCode = -1;
    int terminationSignal = 0;

    std::string executableAbsolutePath;
    std::string expectedExecutableSha256;
    std::string observedExecutableSha256;
    std::string profileAbsolutePath;
    std::string expectedProfileSha256;
    std::string observedProfileSha256;
    std::string expectedVersion;
    std::string observedVersion;

    std::string inputMeshSha256;
    std::string outputToolpathSha256;
    std::string invocationSha256;

    std::uintmax_t versionStdoutBytes = 0;
    std::uintmax_t versionStderrBytes = 0;
    std::string versionStdoutSha256;
    std::string versionStderrSha256;
    std::string versionStdoutExcerpt;
    std::string versionStderrExcerpt;

    std::uintmax_t stdoutBytes = 0;
    std::uintmax_t stderrBytes = 0;
    std::string stdoutSha256;
    std::string stderrSha256;
    std::string stdoutExcerpt;
    std::string stderrExcerpt;

    int automaticRepairMessageCount = 0;
    GcodeToolpathAudit toolpath;
    std::string diagnostic;
    std::vector<std::string> failures;

    /** Stable JSON; no timestamp, ambient PATH, or unpinned host setting. */
    std::string toJson() const;
};


const char* slicerValidationStatusToken(
    SlicerValidationStatus status
) noexcept;


const char* automaticRepairStatusToken(
    AutomaticRepairStatus status
) noexcept;


/** Compute a lower-case SHA-256 digest; throws when the file cannot be read. */
std::string computeFileSha256(const std::string& filename);


/** Compute a lower-case SHA-256 digest of an in-memory byte string. */
std::string computeContentSha256(const std::string& content);


/**
 * Run one argv-based, pinned-binary/profile PrusaSlicer FFF smoke test.
 *
 * The adapter executes exactly two argv-based commands, never through a shell:
 * `<exe> --version`, followed by an isolated `--load ... --export-gcode` run.
 * It never searches PATH, installs software, uploads a toolpath, or talks to a
 * printer. The child still inherits the host process environment and dynamic
 * library closure. This fork/exec adapter must be called while the process is
 * single-threaded, with normal exclusive SIGCHLD ownership: no
 * SIG_IGN/SA_NOCLDWAIT and no concurrent reaper for this child. `Passed` means
 * only that the exact binary/profile/input produced a
 * parseable nonempty G-code artifact without a known reported auto-repair.
 */
SlicerValidationReport validateWithPinnedSlicer(
    const std::string& inputMesh,
    const std::string& outputToolpath,
    const SlicerValidationConfig& config
);


/** Atomically replace one report using a same-directory temporary file. */
void writeSlicerValidationReport(
    const SlicerValidationReport& report,
    const std::string& filename
);


} // namespace manufacturing
