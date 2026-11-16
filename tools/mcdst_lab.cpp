#include "McDstLab.h"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct CommandLine {
    mcdst_lab::Options options;
    unsigned threads = 0;
    std::size_t max_files = 0;
    std::vector<std::string> paths;
    bool help = false;
};

void Usage() {
    std::cout << "Usage: mcdst_lab [options] INPUT_DIRECTORY OUTPUT_DIRECTORY\n"
                 "  --threads N|all      Workers across files (default all hardware threads)\n"
                 "  --max-files N       Maximum number of input files (default all)\n"
                 "  --hist-only          Write QA-only .qa.root files without McDst tree\n"
                 "  --sqrtsnn GeV        Symmetric CM beam energy (default 3)\n"
                 "  --nucleon-mass GeV   SMASH beam nucleon mass (default 0.938)\n"
                 "  --eta-max value      Plot range [-value,value] (default 6)\n"
                 "  --pt-max GeV         Plot range [0,value] (default 3)\n"
                 "Processes .root files directly in INPUT_DIRECTORY (not recursive).\n"
                 "Creates OUTPUT_DIRECTORY; writes <name>.McDst.root with McDst/Event/Particle.\n"
                 "Input must be in SMASH CM axes (Projectile +z, Target -z).\n"
                 "Output: Target rest, +z opposite Projectile; no particle/event cuts.\n"
                 "FXT multiplicity: known charged particles, -2 < eta_lab < 0, pT > 0.15 GeV/c.\n"
                 "Centrality: 0-10, 10-30, 30-50, 50-80% of the shared right-tail event integral.\n"
                 "Equal FXT multiplicities stay together; NFxt=0 events are included.\n"
                 "An existing output is replaced only after success.\n";
}

std::string NextValue(int& index, int argc, char** argv, const std::string& option) {
    if (++index == argc) {
        throw std::runtime_error("missing value for " + option);
    }
    return argv[index];
}

bool IsUnsignedInteger(const std::string& text) {
    return !text.empty() && text.find_first_not_of("0123456789") == std::string::npos;
}

unsigned ParseThreads(const std::string& text) {
    if (text == "all") {
        return 0;
    }
    const char* message = "--threads requires a positive integer or all";
    if (!IsUnsignedInteger(text)) {
        throw std::runtime_error(message);
    }
    const auto value = std::stoull(text);
    if (value == 0 || value > std::numeric_limits<unsigned>::max()) {
        throw std::runtime_error(message);
    }
    return static_cast<unsigned>(value);
}

std::size_t ParseMaxFiles(const std::string& text) {
    if (!IsUnsignedInteger(text)) {
        throw std::runtime_error("--max-files requires a non-negative integer");
    }
    return std::stoull(text);
}

double& NumericOption(mcdst_lab::Options& options, const std::string& name) {
    if (name == "--sqrtsnn") {
        return options.sqrtsnn;
    }
    if (name == "--nucleon-mass") {
        return options.nucleon_mass;
    }
    if (name == "--eta-max") {
        return options.eta_max;
    }
    if (name == "--pt-max") {
        return options.pt_max;
    }
    throw std::runtime_error("unknown option: " + name);
}

void ParseOption(CommandLine& command, const std::string& name, int& index, int argc, char** argv) {
    if (name == "--hist-only") {
        command.options.hist_only = true;
        return;
    }
    if (name == "--threads") {
        command.threads = ParseThreads(NextValue(index, argc, argv, name));
        return;
    }
    if (name == "--max-files") {
        command.max_files = ParseMaxFiles(NextValue(index, argc, argv, name));
        return;
    }

    double& destination = NumericOption(command.options, name);
    const auto text = NextValue(index, argc, argv, name);
    std::size_t end = 0;
    destination = std::stod(text, &end);
    if (end != text.size()) {
        throw std::runtime_error("invalid value for " + name);
    }
}

CommandLine ParseArguments(int argc, char** argv) {
    CommandLine command;
    bool positional_only = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (!positional_only && (argument == "--help" || argument == "-h")) {
            command.help = true;
            return command;
        }
        if (!positional_only && argument == "--") {
            positional_only = true;
            continue;
        }
        if (!positional_only && !argument.empty() && argument.front() == '-') {
            ParseOption(command, argument, index, argc, argv);
        } else {
            command.paths.push_back(argument);
        }
    }
    return command;
}

} // namespace

int main(int argc, char** argv) {
    CommandLine command;
    try {
        command = ParseArguments(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "Argument error: " << error.what() << '\n';
        return 1;
    }
    if (command.help || command.paths.size() != 2) {
        Usage();
        return command.help ? 0 : 1;
    }

    try {
        const auto result =
            mcdst_lab::ProcessDirectory(command.paths[0], command.paths[1], command.options,
                                        command.threads, command.max_files);
        return result.failed == 0 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "Lab processing failed: " << error.what() << '\n';
        return 2;
    }
}
