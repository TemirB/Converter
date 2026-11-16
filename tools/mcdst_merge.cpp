#include "McDstLab.h"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

struct CommandLine {
    mcdst_lab::MergeOptions options;
    std::vector<std::string> paths;
    bool help = false;
};

void Usage() {
    std::cout
        << "Usage: mcdst_merge [--x-max N] OUTPUT.root INPUT.root|DIRECTORY [...]\n"
           "  --x-max N  Inclusive maximum integer multiplicity (default: input ranges).\n"
           "             Higher counts go to overflow; event trees and calibration stay complete.\n"
           "Merge mcdst_lab results with identical FXT calibration.\n"
           "Input directories are scanned for .root files, not recursively.\n"
           "Events and local histograms are combined; common calibration is kept once.\n"
           "Existing outputs, duplicate sources and incompatible inputs are rejected.\n";
}

int ParseMaximum(const std::string& text) {
    const char* message = "--x-max requires a non-negative integer within the histogram bin limit";
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos) {
        throw std::runtime_error(message);
    }
    const auto value = std::stoull(text);
    if (value > static_cast<unsigned long long>(std::numeric_limits<int>::max() - 3)) {
        throw std::runtime_error(message);
    }
    return static_cast<int>(value);
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
        if (!positional_only && argument == "--x-max") {
            if (++index == argc) {
                throw std::runtime_error("missing value for --x-max");
            }
            command.options.x_max = ParseMaximum(argv[index]);
        } else if (!positional_only && !argument.empty() && argument.front() == '-') {
            throw std::runtime_error("unknown option: " + argument);
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
    if (command.help || command.paths.size() < 2) {
        Usage();
        return command.help ? 0 : 1;
    }
    try {
        const auto& output = command.paths.front();
        const std::vector<std::string> inputs(command.paths.begin() + 1, command.paths.end());
        const auto result = mcdst_lab::Merge(inputs, output, command.options);
        std::cout << "Merged: " << result.events << " events, " << result.particles
                  << " particles -> " << output << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "McDst merge failed: " << error.what() << '\n';
        return 2;
    }
}
