#include "manufacturing/SlicerValidation.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>


namespace manufacturing
{
namespace
{


class Sha256
{
public:
    Sha256()
        : state_{{
            0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
            0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U
        }}, block_{}, blockSize_(0), byteCount_(0)
    {
    }

    void update(const unsigned char* data, std::size_t size)
    {
        byteCount_ += size;
        while (size > 0)
        {
            const std::size_t count = std::min(
                size,
                block_.size() - blockSize_
            );
            std::memcpy(block_.data() + blockSize_, data, count);
            data += count;
            size -= count;
            blockSize_ += count;
            if (blockSize_ == block_.size())
            {
                transform();
                blockSize_ = 0;
            }
        }
    }

    std::string finish()
    {
        const std::uint64_t bits =
            static_cast<std::uint64_t>(byteCount_) * UINT64_C(8);
        block_[blockSize_++] = 0x80U;
        if (blockSize_ > 56U)
        {
            std::fill(block_.begin() + blockSize_, block_.end(), 0U);
            transform();
            blockSize_ = 0;
        }
        std::fill(block_.begin() + blockSize_, block_.begin() + 56, 0U);
        for (int byte = 0; byte < 8; ++byte)
        {
            block_[63 - byte] = static_cast<unsigned char>(
                bits >> (8U * static_cast<unsigned int>(byte))
            );
        }
        transform();

        std::ostringstream output;
        output.imbue(std::locale::classic());
        output << std::hex << std::setfill('0');
        for (const std::uint32_t value : state_)
        {
            output << std::setw(8) << value;
        }
        return output.str();
    }

private:
    static std::uint32_t rotate(std::uint32_t value, unsigned int count)
    {
        return (value >> count) | (value << (32U - count));
    }

    void transform()
    {
        static constexpr std::uint32_t k[64] = {
            0x428a2f98U,0x71374491U,0xb5c0fbcfU,0xe9b5dba5U,
            0x3956c25bU,0x59f111f1U,0x923f82a4U,0xab1c5ed5U,
            0xd807aa98U,0x12835b01U,0x243185beU,0x550c7dc3U,
            0x72be5d74U,0x80deb1feU,0x9bdc06a7U,0xc19bf174U,
            0xe49b69c1U,0xefbe4786U,0x0fc19dc6U,0x240ca1ccU,
            0x2de92c6fU,0x4a7484aaU,0x5cb0a9dcU,0x76f988daU,
            0x983e5152U,0xa831c66dU,0xb00327c8U,0xbf597fc7U,
            0xc6e00bf3U,0xd5a79147U,0x06ca6351U,0x14292967U,
            0x27b70a85U,0x2e1b2138U,0x4d2c6dfcU,0x53380d13U,
            0x650a7354U,0x766a0abbU,0x81c2c92eU,0x92722c85U,
            0xa2bfe8a1U,0xa81a664bU,0xc24b8b70U,0xc76c51a3U,
            0xd192e819U,0xd6990624U,0xf40e3585U,0x106aa070U,
            0x19a4c116U,0x1e376c08U,0x2748774cU,0x34b0bcb5U,
            0x391c0cb3U,0x4ed8aa4aU,0x5b9cca4fU,0x682e6ff3U,
            0x748f82eeU,0x78a5636fU,0x84c87814U,0x8cc70208U,
            0x90befffaU,0xa4506cebU,0xbef9a3f7U,0xc67178f2U
        };
        std::uint32_t w[64]{};
        for (int index = 0; index < 16; ++index)
        {
            w[index] =
                (static_cast<std::uint32_t>(block_[4 * index]) << 24U) |
                (static_cast<std::uint32_t>(block_[4 * index + 1]) << 16U) |
                (static_cast<std::uint32_t>(block_[4 * index + 2]) << 8U) |
                static_cast<std::uint32_t>(block_[4 * index + 3]);
        }
        for (int index = 16; index < 64; ++index)
        {
            const std::uint32_t a = w[index - 15];
            const std::uint32_t b = w[index - 2];
            const std::uint32_t s0 = rotate(a, 7) ^ rotate(a, 18) ^ (a >> 3U);
            const std::uint32_t s1 = rotate(b, 17) ^ rotate(b, 19) ^ (b >> 10U);
            w[index] = w[index - 16] + s0 + w[index - 7] + s1;
        }
        std::uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        std::uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
        for (int index = 0; index < 64; ++index)
        {
            const std::uint32_t s1 = rotate(e, 6) ^ rotate(e, 11) ^ rotate(e, 25);
            const std::uint32_t choice = (e & f) ^ (~e & g);
            const std::uint32_t first = h + s1 + choice + k[index] + w[index];
            const std::uint32_t s0 = rotate(a, 2) ^ rotate(a, 13) ^ rotate(a, 22);
            const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
            const std::uint32_t second = s0 + majority;
            h = g; g = f; f = e; e = d + first;
            d = c; c = b; b = a; a = first + second;
        }
        state_[0] += a; state_[1] += b; state_[2] += c; state_[3] += d;
        state_[4] += e; state_[5] += f; state_[6] += g; state_[7] += h;
    }

    std::array<std::uint32_t, 8> state_;
    std::array<unsigned char, 64> block_;
    std::size_t blockSize_;
    std::size_t byteCount_;
};


std::string escapeJson(const std::string& value)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    for (const unsigned char character : value)
    {
        switch (character)
        {
            case '"': output << "\\\""; break;
            case '\\': output << "\\\\"; break;
            case '\b': output << "\\b"; break;
            case '\f': output << "\\f"; break;
            case '\n': output << "\\n"; break;
            case '\r': output << "\\r"; break;
            case '\t': output << "\\t"; break;
            default:
                if (character < 0x20U)
                {
                    output << "\\u" << std::hex << std::setw(4)
                        << std::setfill('0') << static_cast<int>(character)
                        << std::dec;
                }
                else
                {
                    output << character;
                }
        }
    }
    return output.str();
}


std::string normalizedHash(const std::string& value)
{
    std::string result = value;
    std::transform(
        result.begin(),
        result.end(),
        result.begin(),
        [](unsigned char character)
        {
            return static_cast<char>(std::tolower(character));
        }
    );
    return result;
}


bool validHash(const std::string& value)
{
    return value.size() == 64U && std::all_of(
        value.begin(), value.end(),
        [](unsigned char character) { return std::isxdigit(character) != 0; }
    );
}


std::string trim(const std::string& value)
{
    const std::size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
    {
        return std::string();
    }
    const std::size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1U);
}


bool regularFile(const std::string& filename)
{
    std::error_code error;
    return std::filesystem::is_regular_file(filename, error);
}


std::filesystem::path normalizedAbsolutePath(
    const std::filesystem::path& path
)
{
    std::error_code error;
    const std::filesystem::path absolute = std::filesystem::absolute(
        path,
        error
    );
    if (error)
    {
        return std::filesystem::path();
    }
    return absolute.lexically_normal();
}


bool aliasesExistingArtifact(
    const std::filesystem::path& output,
    const std::filesystem::path& artifact
)
{
    if (output == artifact)
    {
        return true;
    }
    std::error_code outputError;
    const bool outputExists = std::filesystem::exists(output, outputError);
    if (outputError || !outputExists)
    {
        return false;
    }
    std::error_code equivalentError;
    return std::filesystem::equivalent(output, artifact, equivalentError) &&
        !equivalentError;
}


std::vector<unsigned char> readBytes(const std::string& filename)
{
    std::ifstream input(filename, std::ios::binary);
    if (!input)
    {
        throw std::runtime_error("cannot read pinned slicer artifact");
    }
    return std::vector<unsigned char>(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    );
}


std::string bytesHash(const std::vector<unsigned char>& bytes)
{
    Sha256 hash;
    if (!bytes.empty())
    {
        hash.update(bytes.data(), bytes.size());
    }
    return hash.finish();
}


struct ProcessCapture
{
    bool launchSucceeded = false;
    bool completed = false;
    bool exitedNormally = false;
    bool timedOut = false;
    int exitCode = -1;
    int signal = 0;
    std::vector<unsigned char> standardOutput;
    std::vector<unsigned char> standardError;
};


std::filesystem::path uniqueTemporary(
    const std::filesystem::path& directory,
    const std::string& stem
)
{
    static std::atomic<unsigned long long> counter{0};
    for (int attempt = 0; attempt < 128; ++attempt)
    {
        const std::filesystem::path candidate = directory /
            ("." + stem + ".r8." +
             std::to_string(static_cast<long long>(::getpid())) + "." +
             std::to_string(counter.fetch_add(1)));
        std::error_code error;
        if (!std::filesystem::exists(candidate, error))
        {
            return candidate;
        }
    }
    throw std::runtime_error("cannot allocate unique slicer temporary path");
}


ProcessCapture runProcess(
    const std::vector<std::string>& arguments,
    int timeoutMilliseconds,
    const std::filesystem::path& directory,
    const std::string& stem
)
{
    ProcessCapture result;
    const std::filesystem::path stdoutPath = uniqueTemporary(directory, stem + ".out");
    const std::filesystem::path stderrPath = uniqueTemporary(directory, stem + ".err");
    const int stdoutDescriptor = ::open(
        stdoutPath.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0600
    );
    const int stderrDescriptor = ::open(
        stderrPath.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0600
    );
    if (stdoutDescriptor < 0 || stderrDescriptor < 0)
    {
        if (stdoutDescriptor >= 0) ::close(stdoutDescriptor);
        if (stderrDescriptor >= 0) ::close(stderrDescriptor);
        throw std::runtime_error("cannot create slicer capture files");
    }

    const pid_t child = ::fork();
    if (child == 0)
    {
        ::dup2(stdoutDescriptor, STDOUT_FILENO);
        ::dup2(stderrDescriptor, STDERR_FILENO);
        ::close(stdoutDescriptor);
        ::close(stderrDescriptor);
        if (::chdir(directory.c_str()) != 0)
        {
            _exit(126);
        }

        std::vector<char*> argv;
        argv.reserve(arguments.size() + 1U);
        for (const std::string& argument : arguments)
        {
            argv.push_back(const_cast<char*>(argument.c_str()));
        }
        argv.push_back(nullptr);
        ::execv(arguments.front().c_str(), argv.data());
        _exit(127);
    }
    ::close(stdoutDescriptor);
    ::close(stderrDescriptor);

    if (child < 0)
    {
        std::error_code ignored;
        std::filesystem::remove(stdoutPath, ignored);
        std::filesystem::remove(stderrPath, ignored);
        return result;
    }
    result.launchSucceeded = true;

    int status = 0;
    const auto deadline = std::chrono::steady_clock::now() +
        std::chrono::milliseconds(timeoutMilliseconds);
    while (true)
    {
        const pid_t state = ::waitpid(child, &status, WNOHANG);
        if (state == child)
        {
            result.completed = true;
            break;
        }
        if (state < 0)
        {
            break;
        }
        if (std::chrono::steady_clock::now() >= deadline)
        {
            result.timedOut = true;
            ::kill(child, SIGKILL);
            ::waitpid(child, &status, 0);
            break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    if (WIFEXITED(status))
    {
        result.exitedNormally = true;
        result.exitCode = WEXITSTATUS(status);
    }
    else if (WIFSIGNALED(status))
    {
        result.signal = WTERMSIG(status);
    }
    result.standardOutput = readBytes(stdoutPath.string());
    result.standardError = readBytes(stderrPath.string());
    std::error_code ignored;
    std::filesystem::remove(stdoutPath, ignored);
    std::filesystem::remove(stderrPath, ignored);
    return result;
}


void assignCapture(
    const ProcessCapture& capture,
    std::uintmax_t maximumExcerpt,
    std::uintmax_t& stdoutBytes,
    std::string& stdoutHash,
    std::string& stdoutExcerpt,
    std::uintmax_t& stderrBytes,
    std::string& stderrHash,
    std::string& stderrExcerpt
)
{
    stdoutBytes = capture.standardOutput.size();
    stderrBytes = capture.standardError.size();
    stdoutHash = bytesHash(capture.standardOutput);
    stderrHash = bytesHash(capture.standardError);
    stdoutExcerpt.assign(
        capture.standardOutput.begin(),
        capture.standardOutput.begin() + static_cast<std::ptrdiff_t>(std::min(
            maximumExcerpt,
            static_cast<std::uintmax_t>(capture.standardOutput.size())
        ))
    );
    stderrExcerpt.assign(
        capture.standardError.begin(),
        capture.standardError.begin() + static_cast<std::ptrdiff_t>(std::min(
            maximumExcerpt,
            static_cast<std::uintmax_t>(capture.standardError.size())
        ))
    );
}


GcodeToolpathAudit auditGcode(const std::string& filename)
{
    GcodeToolpathAudit audit;
    if (!regularFile(filename))
    {
        return audit;
    }
    const std::vector<unsigned char> bytes = readBytes(filename);
    audit.completed = true;
    audit.bytes = bytes.size();
    audit.nonempty = !bytes.empty();
    audit.sha256 = bytesHash(bytes);
    audit.text = std::find(bytes.begin(), bytes.end(), 0U) == bytes.end();
    if (!audit.text)
    {
        return audit;
    }

    std::istringstream input(std::string(bytes.begin(), bytes.end()));
    input.imbue(std::locale::classic());
    std::string line;
    while (std::getline(input, line))
    {
        ++audit.lineCount;
        const std::string stripped = trim(line);
        bool layerMarker = stripped == ";LAYER_CHANGE";
        if (stripped.rfind(";LAYER:", 0) == 0)
        {
            const std::string value = trim(stripped.substr(7));
            layerMarker = !value.empty() && std::all_of(
                value.begin(),
                value.end(),
                [](unsigned char character)
                {
                    return std::isdigit(character) != 0;
                }
            );
        }
        if (layerMarker)
        {
            ++audit.layerCount;
        }

        const std::size_t comment = line.find(';');
        const std::string command = trim(line.substr(0, comment));
        std::istringstream fields(command);
        fields.imbue(std::locale::classic());
        std::string opcode;
        fields >> opcode;
        bool finitePositiveExtrusion = false;
        std::string field;
        while (fields >> field)
        {
            if (field.size() <= 1U || field.front() != 'E')
            {
                continue;
            }
            std::istringstream valueInput(field.substr(1));
            valueInput.imbue(std::locale::classic());
            double value = 0.0;
            char trailing = '\0';
            if ((valueInput >> value) && !(valueInput >> trailing) &&
                std::isfinite(value) && value > 0.0)
            {
                finitePositiveExtrusion = true;
            }
        }
        if ((opcode == "G0" || opcode == "G00" || opcode == "G1" ||
             opcode == "G01") && finitePositiveExtrusion)
        {
            ++audit.depositionCommandCount;
        }
    }
    audit.hasLayers = audit.layerCount > 0;
    audit.hasDeposition = audit.depositionCommandCount > 0;
    return audit;
}


int repairMessageCount(const std::vector<unsigned char>& bytes)
{
    std::string text(bytes.begin(), bytes.end());
    std::transform(
        text.begin(), text.end(), text.begin(),
        [](unsigned char character)
        {
            return static_cast<char>(std::tolower(character));
        }
    );
    const std::array<const char*, 4> tokens{{
        "auto-repair", "auto repair", "repaired", "fixing model"
    }};
    int count = 0;
    for (const char* token : tokens)
    {
        std::size_t position = 0;
        while ((position = text.find(token, position)) != std::string::npos)
        {
            ++count;
            position += std::strlen(token);
        }
    }
    return count;
}


std::vector<unsigned char> stringBytes(const std::string& value)
{
    return std::vector<unsigned char>(value.begin(), value.end());
}


void writeAtomic(const std::string& filename, const std::vector<unsigned char>& bytes)
{
    if (filename.empty())
    {
        throw std::invalid_argument("atomic slicer output filename is empty");
    }
    const std::filesystem::path destination(filename);
    const std::filesystem::path directory = destination.has_parent_path()
        ? destination.parent_path() : std::filesystem::path(".");
    const std::filesystem::path temporary = uniqueTemporary(
        directory,
        destination.filename().string()
    );
    int descriptor = ::open(
        temporary.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0666
    );
    if (descriptor < 0)
    {
        throw std::runtime_error("cannot create slicer output temporary");
    }
    try
    {
        std::size_t offset = 0;
        while (offset < bytes.size())
        {
            const ssize_t written = ::write(
                descriptor,
                bytes.data() + offset,
                bytes.size() - offset
            );
            if (written < 0 && errno == EINTR) continue;
            if (written <= 0) throw std::runtime_error("slicer output write failed");
            offset += static_cast<std::size_t>(written);
        }
        if (::fsync(descriptor) != 0)
        {
            throw std::runtime_error("slicer output sync failed");
        }
        if (::close(descriptor) != 0)
        {
            descriptor = -1;
            throw std::runtime_error("slicer output close failed");
        }
        descriptor = -1;
        if (::rename(temporary.c_str(), destination.c_str()) != 0)
        {
            throw std::runtime_error("slicer output commit failed");
        }
    }
    catch (...)
    {
        if (descriptor >= 0)
        {
            ::close(descriptor);
        }
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}


std::string invocationFingerprint(const std::vector<std::string>& arguments)
{
    Sha256 hash;
    for (const std::string& argument : arguments)
    {
        hash.update(
            reinterpret_cast<const unsigned char*>(argument.data()),
            argument.size()
        );
        const unsigned char separator = 0U;
        hash.update(&separator, 1U);
    }
    return hash.finish();
}


} // namespace


const char* slicerValidationStatusToken(SlicerValidationStatus status) noexcept
{
    switch (status)
    {
        case SlicerValidationStatus::NotConfigured: return "not_configured";
        case SlicerValidationStatus::InvalidConfiguration: return "invalid_configuration";
        case SlicerValidationStatus::InputUnavailable: return "input_unavailable";
        case SlicerValidationStatus::ExecutableUnavailable: return "executable_unavailable";
        case SlicerValidationStatus::ProfileUnavailable: return "profile_unavailable";
        case SlicerValidationStatus::ExecutableHashMismatch: return "executable_hash_mismatch";
        case SlicerValidationStatus::ProfileHashMismatch: return "profile_hash_mismatch";
        case SlicerValidationStatus::VersionProbeFailed: return "version_probe_failed";
        case SlicerValidationStatus::VersionMismatch: return "version_mismatch";
        case SlicerValidationStatus::LaunchFailed: return "launch_failed";
        case SlicerValidationStatus::TimedOut: return "timed_out";
        case SlicerValidationStatus::NonzeroExit: return "nonzero_exit";
        case SlicerValidationStatus::AutoRepairUnknown: return "auto_repair_unknown";
        case SlicerValidationStatus::AutoRepairDetected: return "auto_repair_detected";
        case SlicerValidationStatus::InvalidToolpath: return "invalid_toolpath";
        case SlicerValidationStatus::InputChanged: return "input_changed";
        case SlicerValidationStatus::OutputCommitFailed: return "output_commit_failed";
        case SlicerValidationStatus::Passed: return "passed";
    }
    return "unknown";
}


const char* automaticRepairStatusToken(AutomaticRepairStatus status) noexcept
{
    switch (status)
    {
        case AutomaticRepairStatus::Unsupported: return "unsupported";
        case AutomaticRepairStatus::Unknown: return "unknown";
        case AutomaticRepairStatus::NoneReported: return "none_reported";
        case AutomaticRepairStatus::Reported: return "reported";
    }
    return "unknown";
}


std::string computeFileSha256(const std::string& filename)
{
    return bytesHash(readBytes(filename));
}


std::string SlicerValidationReport::toJson() const
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << std::boolalpha
        << "{\n"
        << "  \"schema\": \"make-it-stand-r8-slicer-v1\",\n"
        << "  \"claim_ceiling\": \"pinned local slicer smoke test only\",\n"
        << "  \"status\": \"" << slicerValidationStatusToken(status) << "\",\n"
        << "  \"automatic_repair\": \""
        << automaticRepairStatusToken(automaticRepair) << "\",\n"
        << "  \"configuration_complete\": " << configurationComplete << ",\n"
        << "  \"attempted\": " << attempted << ",\n"
        << "  \"completed\": " << completed << ",\n"
        << "  \"passed\": " << passed << ",\n"
        << "  \"exit_code\": " << exitCode << ",\n"
        << "  \"termination_signal\": " << terminationSignal << ",\n"
        << "  \"timed_out\": " << timedOut << ",\n"
        << "  \"executable_path\": \"" << escapeJson(executableAbsolutePath) << "\",\n"
        << "  \"executable_sha256\": \"" << observedExecutableSha256 << "\",\n"
        << "  \"profile_path\": \"" << escapeJson(profileAbsolutePath) << "\",\n"
        << "  \"profile_sha256\": \"" << observedProfileSha256 << "\",\n"
        << "  \"expected_version\": \"" << escapeJson(expectedVersion) << "\",\n"
        << "  \"observed_version\": \"" << escapeJson(observedVersion) << "\",\n"
        << "  \"input_mesh_sha256\": \"" << inputMeshSha256 << "\",\n"
        << "  \"invocation_sha256\": \"" << invocationSha256 << "\",\n"
        << "  \"stdout_sha256\": \"" << stdoutSha256 << "\",\n"
        << "  \"stderr_sha256\": \"" << stderrSha256 << "\",\n"
        << "  \"stdout_excerpt\": \"" << escapeJson(stdoutExcerpt) << "\",\n"
        << "  \"stderr_excerpt\": \"" << escapeJson(stderrExcerpt) << "\",\n"
        << "  \"toolpath\": {\"bytes\": " << toolpath.bytes
        << ", \"layers\": " << toolpath.layerCount
        << ", \"deposition_commands\": " << toolpath.depositionCommandCount
        << ", \"sha256\": \"" << toolpath.sha256 << "\"},\n"
        << "  \"diagnostic\": \"" << escapeJson(diagnostic) << "\",\n"
        << "  \"failures\": [";
    for (std::size_t index = 0; index < failures.size(); ++index)
    {
        if (index > 0) output << ", ";
        output << '"' << escapeJson(failures[index]) << '"';
    }
    output << "]\n}\n";
    return output.str();
}


SlicerValidationReport validateWithPinnedSlicer(
    const std::string& inputMesh,
    const std::string& outputToolpath,
    const SlicerValidationConfig& config
)
{
    SlicerValidationReport report;
    report.executableAbsolutePath = config.executableAbsolutePath;
    report.expectedExecutableSha256 = normalizedHash(
        config.expectedExecutableSha256
    );
    report.profileAbsolutePath = config.profileAbsolutePath;
    report.expectedProfileSha256 = normalizedHash(config.expectedProfileSha256);
    report.expectedVersion = config.expectedVersion;

    const bool anyConfigured =
        !config.executableAbsolutePath.empty() ||
        !config.expectedExecutableSha256.empty() ||
        !config.profileAbsolutePath.empty() ||
        !config.expectedProfileSha256.empty() ||
        !config.expectedVersion.empty();
    if (!anyConfigured)
    {
        report.status = SlicerValidationStatus::NotConfigured;
        report.failures.push_back("slicer_not_configured");
        return report;
    }

    report.configurationComplete =
        std::filesystem::path(config.executableAbsolutePath).is_absolute() &&
        std::filesystem::path(config.profileAbsolutePath).is_absolute() &&
        validHash(config.expectedExecutableSha256) &&
        validHash(config.expectedProfileSha256) &&
        !config.expectedVersion.empty() &&
        config.timeoutMilliseconds > 0;
    if (!report.configurationComplete || outputToolpath.empty())
    {
        report.status = SlicerValidationStatus::InvalidConfiguration;
        report.failures.push_back("invalid_slicer_configuration");
        return report;
    }
    const std::filesystem::path inputPath = normalizedAbsolutePath(inputMesh);
    const std::filesystem::path executablePath = normalizedAbsolutePath(
        config.executableAbsolutePath
    );
    const std::filesystem::path profilePath = normalizedAbsolutePath(
        config.profileAbsolutePath
    );
    const std::filesystem::path outputPath = normalizedAbsolutePath(
        outputToolpath
    );
    if (
        inputPath.empty() || executablePath.empty() || profilePath.empty() ||
        outputPath.empty() || aliasesExistingArtifact(outputPath, inputPath) ||
        aliasesExistingArtifact(outputPath, executablePath) ||
        aliasesExistingArtifact(outputPath, profilePath)
    )
    {
        report.status = SlicerValidationStatus::InvalidConfiguration;
        report.failures.push_back("slicer_artifact_path_collision");
        return report;
    }
    if (!regularFile(inputPath.string()))
    {
        report.status = SlicerValidationStatus::InputUnavailable;
        report.failures.push_back("input_mesh_unavailable");
        return report;
    }
    report.inputMeshSha256 = computeFileSha256(inputPath.string());
    if (!regularFile(executablePath.string()) ||
        ::access(executablePath.c_str(), X_OK) != 0)
    {
        report.status = SlicerValidationStatus::ExecutableUnavailable;
        report.failures.push_back("slicer_executable_unavailable");
        return report;
    }
    report.executableAvailable = true;
    report.observedExecutableSha256 = computeFileSha256(
        executablePath.string()
    );
    if (report.observedExecutableSha256 != report.expectedExecutableSha256)
    {
        report.status = SlicerValidationStatus::ExecutableHashMismatch;
        report.failures.push_back("slicer_executable_hash_mismatch");
        return report;
    }
    report.executableHashVerified = true;

    if (!regularFile(profilePath.string()))
    {
        report.status = SlicerValidationStatus::ProfileUnavailable;
        report.failures.push_back("slicer_profile_unavailable");
        return report;
    }
    report.profileAvailable = true;
    report.observedProfileSha256 = computeFileSha256(profilePath.string());
    if (report.observedProfileSha256 != report.expectedProfileSha256)
    {
        report.status = SlicerValidationStatus::ProfileHashMismatch;
        report.failures.push_back("slicer_profile_hash_mismatch");
        return report;
    }
    report.profileHashVerified = true;

    const std::filesystem::path directory = outputPath.has_parent_path()
        ? outputPath.parent_path() : std::filesystem::path(".");
    if (!std::filesystem::is_directory(directory))
    {
        report.status = SlicerValidationStatus::InvalidConfiguration;
        report.failures.push_back("slicer_output_directory_unavailable");
        return report;
    }

    report.versionProbeAttempted = true;
    const ProcessCapture version = runProcess(
        {executablePath.string(), "--version"},
        config.timeoutMilliseconds,
        directory,
        "slicer-version"
    );
    assignCapture(
        version,
        config.maximumDiagnosticBytes,
        report.versionStdoutBytes,
        report.versionStdoutSha256,
        report.versionStdoutExcerpt,
        report.versionStderrBytes,
        report.versionStderrSha256,
        report.versionStderrExcerpt
    );
    if (!version.launchSucceeded || version.timedOut ||
        !version.exitedNormally || version.exitCode != 0)
    {
        report.status = SlicerValidationStatus::VersionProbeFailed;
        report.failures.push_back("slicer_version_probe_failed");
        return report;
    }
    report.observedVersion = trim(report.versionStdoutExcerpt);
    if (report.observedVersion != config.expectedVersion)
    {
        report.status = SlicerValidationStatus::VersionMismatch;
        report.failures.push_back("slicer_version_mismatch");
        return report;
    }
    report.versionVerified = true;

    const std::filesystem::path stagedOutput = uniqueTemporary(
        directory,
        outputPath.filename().string() + ".gcode"
    );
    const std::vector<std::string> arguments{
        executablePath.string(),
        "--load", profilePath.string(),
        "--export-gcode",
        "--output", stagedOutput.string(),
        inputPath.string()
    };
    report.invocationSha256 = invocationFingerprint(arguments);
    report.attempted = true;
    const ProcessCapture slice = runProcess(
        arguments,
        config.timeoutMilliseconds,
        directory,
        "slicer-run"
    );
    assignCapture(
        slice,
        config.maximumDiagnosticBytes,
        report.stdoutBytes,
        report.stdoutSha256,
        report.stdoutExcerpt,
        report.stderrBytes,
        report.stderrSha256,
        report.stderrExcerpt
    );
    report.completed = slice.completed;
    report.exitedNormally = slice.exitedNormally;
    report.exitCode = slice.exitCode;
    report.terminationSignal = slice.signal;
    report.timedOut = slice.timedOut;
    report.exitSuccess = slice.exitedNormally && slice.exitCode == 0;
    if (!slice.launchSucceeded)
    {
        report.status = SlicerValidationStatus::LaunchFailed;
        report.failures.push_back("slicer_launch_failed");
    }
    else if (slice.timedOut)
    {
        report.status = SlicerValidationStatus::TimedOut;
        report.failures.push_back("slicer_timed_out");
    }
    else if (!report.exitSuccess)
    {
        report.status = SlicerValidationStatus::NonzeroExit;
        report.failures.push_back("slicer_nonzero_exit");
    }
    else
    {
        report.automaticRepairMessageCount =
            repairMessageCount(slice.standardOutput) +
            repairMessageCount(slice.standardError);
        report.automaticRepair = report.automaticRepairMessageCount > 0
            ? AutomaticRepairStatus::Reported
            : AutomaticRepairStatus::NoneReported;
        if (report.automaticRepair == AutomaticRepairStatus::Reported)
        {
            report.status = SlicerValidationStatus::AutoRepairDetected;
            report.failures.push_back("slicer_auto_repair_detected");
        }
        else
        {
            report.toolpath = auditGcode(stagedOutput.string());
            if (!report.toolpath.completed || !report.toolpath.nonempty ||
                !report.toolpath.text || !report.toolpath.hasLayers ||
                !report.toolpath.hasDeposition)
            {
                report.status = SlicerValidationStatus::InvalidToolpath;
                report.failures.push_back("slicer_toolpath_invalid");
            }
            else if (computeFileSha256(inputPath.string()) != report.inputMeshSha256)
            {
                report.status = SlicerValidationStatus::InputChanged;
                report.failures.push_back("input_mesh_changed_during_slicing");
            }
            else
            {
                try
                {
                    writeAtomic(outputToolpath, readBytes(stagedOutput.string()));
                    report.outputToolpathSha256 = computeFileSha256(outputToolpath);
                    report.status = SlicerValidationStatus::Passed;
                    report.passed = true;
                }
                catch (const std::exception& error)
                {
                    report.status = SlicerValidationStatus::OutputCommitFailed;
                    report.diagnostic = error.what();
                    report.failures.push_back("slicer_output_commit_failed");
                }
            }
        }
    }
    std::error_code ignored;
    std::filesystem::remove(stagedOutput, ignored);
    return report;
}


void writeSlicerValidationReport(
    const SlicerValidationReport& report,
    const std::string& filename
)
{
    writeAtomic(filename, stringBytes(report.toJson()));
}


} // namespace manufacturing
