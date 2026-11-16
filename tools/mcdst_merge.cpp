#include "McDstLab.h"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void Usage() {
    std::cout << "Usage: mcdst_merge OUTPUT.root INPUT.root|DIRECTORY [...]\n"
                 "Merge mcdst_lab results with identical FXT calibration.\n"
                 "Input directories are scanned for .root files, not recursively.\n"
                 "Events and local histograms are combined; common calibration is kept once.\n"
                 "Existing outputs, duplicate sources and incompatible inputs are rejected.\n";
}
} // namespace

int main(int argc, char** argv) {
    if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        Usage();
        return 0;
    }
    if (argc < 3) {
        Usage();
        return 1;
    }
    try {
        const std::vector<std::string> inputs(argv + 2, argv + argc);
        const auto result = mcdst_lab::Merge(inputs, argv[1]);
        std::cout << "Merged: " << result.events << " events, " << result.particles
                  << " particles -> " << argv[1] << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "McDst merge failed: " << error.what() << '\n';
        return 2;
    }
}
