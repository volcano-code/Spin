#include <chrono>
#include <csignal>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>

#include <sys/types.h>
#include <unistd.h>


namespace
{


constexpr const char* FAKE_VERSION = "PrusaSlicer 2.9.0-fake";


std::string environment(const char* name)
{
    const char* value = std::getenv(name);
    return value == nullptr ? std::string() : std::string(value);
}


std::string profileMode(const std::string& filename)
{
    std::ifstream input(filename, std::ios::binary);
    if (!input)
    {
        return "profile-unavailable";
    }

    std::string line;
    while (std::getline(input, line))
    {
        constexpr const char* prefix = "mode=";
        if (line.rfind(prefix, 0) == 0)
        {
            return line.substr(std::char_traits<char>::length(prefix));
        }
    }
    return "success";
}


bool writeToolpath(const std::string& filename, const std::string& contents)
{
    std::ofstream output(filename, std::ios::binary | std::ios::trunc);
    output << contents;
    return static_cast<bool>(output);
}


int runVersionProbe(const char* executablePath)
{
    const std::string mode = environment(
        "MAKE_IT_STAND_FAKE_SLICER_VERSION_MODE"
    );
    if (mode == "nonzero")
    {
        std::cerr << "fake version probe failed\n";
        return 41;
    }
    if (mode == "sleep")
    {
        std::this_thread::sleep_for(std::chrono::seconds(5));
    }
    if (mode == "signal")
    {
        std::raise(SIGTERM);
        return 128 + SIGTERM;
    }
    if (mode == "mismatch")
    {
        std::cout << "PrusaSlicer 0.0.0-wrong\n";
        return EXIT_SUCCESS;
    }
    if (mode == "long-suffix")
    {
        std::cout << FAKE_VERSION << std::string(4096, 'x') << '\n';
        return EXIT_SUCCESS;
    }
    if (mode == "unlink-self")
    {
        std::cout << FAKE_VERSION << '\n';
        std::cout.flush();
        if (::unlink(executablePath) != 0)
        {
            std::cerr << "fake version probe could not unlink itself\n";
            return 42;
        }
        return EXIT_SUCCESS;
    }

    std::cout << FAKE_VERSION << '\n';
    return EXIT_SUCCESS;
}


int runSlice(int argc, char** argv)
{
    std::string profile;
    std::string output;
    std::string input;
    bool exportRequested = false;

    for (int index = 1; index < argc; ++index)
    {
        const std::string argument = argv[index];
        if (argument == "--load" && index + 1 < argc)
        {
            profile = argv[++index];
        }
        else if (
            (argument == "--output" || argument == "-o") &&
            index + 1 < argc
        )
        {
            output = argv[++index];
        }
        else if (argument.rfind("--output=", 0) == 0)
        {
            output = argument.substr(9);
        }
        else if (argument == "--export-gcode")
        {
            exportRequested = true;
        }
        else if (!argument.empty() && argument[0] != '-')
        {
            input = argument;
        }
    }

    if (profile.empty() || output.empty() || input.empty() || !exportRequested)
    {
        std::cerr << "fake slicer received incomplete Prusa-style argv\n";
        return 64;
    }

    std::ifstream mesh(input, std::ios::binary);
    if (!mesh)
    {
        std::cerr << "fake slicer cannot read input mesh\n";
        return 66;
    }

    const std::string mode = profileMode(profile);
    if (mode == "profile-unavailable")
    {
        std::cerr << "fake slicer cannot read profile\n";
        return 66;
    }
    if (mode == "nonzero")
    {
        std::cout << "fake slicing started\n";
        std::cerr << "fake slicing failed with requested nonzero exit\n";
        return 23;
    }
    if (mode == "sleep")
    {
        std::this_thread::sleep_for(std::chrono::seconds(5));
    }
    if (mode == "fork-grandchild")
    {
        const pid_t descendant = ::fork();
        if (descendant < 0)
        {
            return 74;
        }
        if (descendant == 0)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(600));
            writeToolpath(profile + ".late-marker", "survived timeout\n");
            _exit(EXIT_SUCCESS);
        }
        std::this_thread::sleep_for(std::chrono::seconds(5));
    }
    if (mode == "orphan-grandchild")
    {
        const pid_t descendant = ::fork();
        if (descendant < 0)
        {
            return 74;
        }
        if (descendant == 0)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(600));
            writeToolpath(profile + ".late-marker", "survived leader exit\n");
            _exit(EXIT_SUCCESS);
        }
        return EXIT_SUCCESS;
    }
    if (mode == "signal")
    {
        std::raise(SIGTERM);
        return 128 + SIGTERM;
    }
    if (mode == "empty")
    {
        return writeToolpath(output, "") ? EXIT_SUCCESS : 73;
    }
    if (mode == "no-deposition")
    {
        return writeToolpath(
            output,
            "; generated by fake slicer\n"
            ";LAYER_CHANGE\n"
            ";Z:0.20\n"
            "G1 X10.0 Y10.0 F1200\n"
        ) ? EXIT_SUCCESS : 73;
    }
    if (mode == "no-layers")
    {
        return writeToolpath(
            output,
            "; generated by fake slicer\n"
            "M82\n"
            "G1 X10.0 Y10.0 E0.40 F1200\n"
        ) ? EXIT_SUCCESS : 73;
    }
    if (mode == "garbage-gcode")
    {
        return writeToolpath(
            output,
            ";LAYER:this is not a numeric layer\n"
            "G0ARBAGE EVIDENCE\n"
        ) ? EXIT_SUCCESS : 73;
    }
    if (mode == "extrusion-only")
    {
        return writeToolpath(
            output,
            ";LAYER_CHANGE\n"
            "G1 E0.40 F1200\n"
        ) ? EXIT_SUCCESS : 73;
    }
    if (mode == "relative-absolute-retraction")
    {
        return writeToolpath(
            output,
            ";LAYER_CHANGE\n"
            "M83\n"
            "G1 E10\n"
            "M82\n"
            "G1 X1 E1\n"
        ) ? EXIT_SUCCESS : 73;
    }
    if (mode == "stationary-planar")
    {
        return writeToolpath(
            output,
            ";LAYER_CHANGE\n"
            "M82\n"
            "G92 X1 Y2 E0\n"
            "G1 X1 Y2 E1\n"
        ) ? EXIT_SUCCESS : 73;
    }
    if (mode == "g90-ambiguous-extrusion")
    {
        return writeToolpath(
            output,
            ";LAYER_CHANGE\n"
            "M83\n"
            "G1 E-10\n"
            "G90\n"
            "G1 X1 E-1\n"
        ) ? EXIT_SUCCESS : 73;
    }
    if (mode == "g91-m82-override")
    {
        return writeToolpath(
            output,
            ";LAYER_CHANGE\n"
            "G91\n"
            "M82\n"
            "G1 X1 E1\n"
        ) ? EXIT_SUCCESS : 73;
    }
    if (mode == "overflow-delta")
    {
        return writeToolpath(
            output,
            ";LAYER_CHANGE\n"
            "M82\n"
            "G92 X-1e308 E-1e308\n"
            "G1 X1e308 E1e308\n"
        ) ? EXIT_SUCCESS : 73;
    }

    if (mode == "repair")
    {
        /* Include stable diagnostics from the supported Prusa-style family. */
        std::cerr
            << "Warning: The model was repaired.\n"
            << "1 degenerate facet was removed during auto-repair.\n";
    }
    if (mode == "noisy")
    {
        std::cout << std::string(4096, 'x') << '\n';
    }

    const bool written = writeToolpath(
        output,
        "; generated by fake slicer\n"
        ";LAYER_CHANGE\n"
        ";Z:0.20\n"
        "M82\n"
        "G1 X10.0 Y10.0 E0.40 F1200\n"
        ";LAYER_CHANGE\n"
        ";Z:0.40\n"
        "G1 X20.0 Y10.0 E0.75 F1200\n"
    );
    if (!written)
    {
        std::cerr << "fake slicer cannot write output toolpath\n";
        return 73;
    }
    std::cout << "Slicing completed using the pinned fake profile.\n";
    return EXIT_SUCCESS;
}


} // namespace


int main(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        return runVersionProbe(argv[0]);
    }
    return runSlice(argc, argv);
}
