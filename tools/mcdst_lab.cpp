#include "McDstLab.h"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void Usage() {
    std::cout << "Usage: mcdst_lab [options] INPUT_DIRECTORY OUTPUT_DIRECTORY\n"
                 "  --threads N|all      Workers across files (default all hardware threads)\n"
                 "  --hist-only          Write QA-only .qa.root files without McDst tree\n"
                 "  --sqrtsnn GeV        Symmetric CM beam energy (default 3)\n"
                 "  --nucleon-mass GeV   SMASH beam nucleon mass (default 0.938)\n"
                 "  --eta-max value      Plot range [-value,value] (default 6)\n"
                 "  --pt-max GeV         Plot range [0,value] (default 3)\n"
                 "Processes .root files directly in INPUT_DIRECTORY (not recursive).\n"
                 "Creates OUTPUT_DIRECTORY; writes <name>.McDst.root with McDst/Event/Particle.\n"
                 "Input must be in SMASH CM axes (Projectile +z, Target -z).\n"
                 "Output: Target rest, +z opposite Projectile; no particle/event cuts.\n"
                 "An existing output is replaced only after success.\n";
}
} // namespace

int main(int argc, char** argv) {
    mcdst_lab::Options options;
    unsigned threads = 0;
    std::vector<std::string> paths;
    try {
        bool positional = false;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (!positional && (arg == "--help" || arg == "-h")) {
                Usage();
                return 0;
            }
            if (!positional && arg == "--") {
                positional = true;
                continue;
            }
            if (!positional && arg == "--hist-only") {
                options.hist_only = true;
                continue;
            }
            if (!positional && arg == "--threads") {
                if (++i == argc)
                    throw std::runtime_error("missing value for --threads");
                const std::string count = argv[i];
                if (count == "all") {
                    threads = 0;
                } else {
                    if (count.empty() || count.find_first_not_of("0123456789") != std::string::npos)
                        throw std::runtime_error("--threads requires a positive integer or all");
                    const auto value = std::stoull(count);
                    if (value == 0 || value > std::numeric_limits<unsigned>::max())
                        throw std::runtime_error("--threads requires a positive integer or all");
                    threads = static_cast<unsigned>(value);
                }
                continue;
            }
            if (!positional && !arg.empty() && arg.front() == '-') {
                double* value = nullptr;
                if (arg == "--sqrtsnn")
                    value = &options.sqrtsnn;
                else if (arg == "--nucleon-mass")
                    value = &options.nucleon_mass;
                else if (arg == "--eta-max")
                    value = &options.eta_max;
                else if (arg == "--pt-max")
                    value = &options.pt_max;
                else
                    throw std::runtime_error("unknown option: " + arg);
                if (++i == argc)
                    throw std::runtime_error("missing value for " + arg);
                std::size_t end = 0;
                *value = std::stod(argv[i], &end);
                if (end != std::string(argv[i]).size())
                    throw std::runtime_error("invalid value for " + arg);
            } else {
                paths.push_back(arg);
            }
        }
        if (paths.size() != 2) {
            Usage();
            return 1;
        }
    } catch (const std::exception& error) {
        std::cerr << "Argument error: " << error.what() << '\n';
        return 1;
    }
    try {
        const auto result = mcdst_lab::ProcessDirectory(paths[0], paths[1], options, threads);
        return result.failed == 0 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "Lab processing failed: " << error.what() << '\n';
        return 2;
    }
}
