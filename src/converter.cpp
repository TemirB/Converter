#include "ConverterValidation.h"
#include "OscarConverter.h"
#include "OscarConverterV1.h"
#include "RootConverter.h"

#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    ConversionOptions options;
    bool root_input = false;
    bool legacy = false, validate = false;
    const auto usage = [&](std::ostream& output) {
        output << "Usage: " << argv[0] << " [--keep-empty] [--root] <input> <output.McDst.root>\n"
               << "       " << argv[0] << " --v1 <input.oscar> <output.McDst.root>\n"
               << "       " << argv[0] << " --validate <input.oscar> [output-directory]\n";
    };
    std::vector<std::string> filenames;
    for (int i = 1; i < argc; ++i) {
        const std::string argument(argv[i]);
        if (argument == "--help" || argument == "-h") {
            usage(std::cout);
            return 0;
        }
        if (argument == "--keep-empty")
            options.keep_empty_events = true;
        else if (argument == "--root")
            root_input = true;
        else if (argument == "--v1")
            legacy = true;
        else if (argument == "--validate")
            validate = true;
        else if (!argument.empty() && argument.front() == '-') {
            std::cerr << "Unknown option: " << argument << '\n';
            return 1;
        } else
            filenames.push_back(argument);
    }
    const bool invalidModes = (root_input && options.keep_empty_events) ||
                              ((legacy || validate) && (root_input || options.keep_empty_events)) ||
                              (legacy && validate);
    if (invalidModes ||
        (validate ? filenames.empty() || filenames.size() > 2 : filenames.size() != 2)) {
        usage(std::cerr);
        return 1;
    }
    if (validate)
        return static_cast<int>(
            ValidateOscarVersions(filenames[0], filenames.size() == 2 ? filenames[1] : ""));
    if (legacy)
        return OscarConverterV1().Convert(filenames[0], filenames[1]) ? 0 : 2;
    const bool success = root_input ? ConvertRootToMcDst(filenames[0], filenames[1])
                                    : OscarConverter().Convert(filenames[0], filenames[1], options);
    return success ? 0 : 2;
}
