#include "OscarConverter.h"
#include "RootConverter.h"

#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    ConversionOptions options;
    bool root_input = false;
    std::vector<std::string> filenames;
    for (int i = 1; i < argc; ++i) {
        const std::string argument(argv[i]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "Usage: " << argv[0]
                      << " [--keep-empty] [--root] <input> <output.McDst.root>\n";
            return 0;
        }
        if (argument == "--keep-empty")
            options.keep_empty_events = true;
        else if (argument == "--root")
            root_input = true;
        else if (!argument.empty() && argument.front() == '-') {
            std::cerr << "Unknown option: " << argument << '\n';
            return 1;
        } else
            filenames.push_back(argument);
    }
    if (filenames.size() != 2 || (root_input && options.keep_empty_events)) {
        std::cerr << "Usage: " << argv[0]
                  << " [--keep-empty] [--root] <input> <output.McDst.root>\n";
        return 1;
    }
    const bool success = root_input ? ConvertRootToMcDst(filenames[0], filenames[1])
                                    : OscarConverter().Convert(filenames[0], filenames[1], options);
    return success ? 0 : 2;
}
