#include "ConverterValidation.h"
#include "OscarConverter.h"
#include "RootConverter.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct CommandLine {
    ConversionOptions conversion;
    bool root_input = false;
    bool version_selected = false;
    bool validate = false;
    bool validate_all = false;
    bool help = false;
    std::vector<std::string> paths;
};

void Usage(std::ostream& output) {
    output << "Usage: converter [--version 0|1|2] [--keep-empty] INPUT.oscar OUTPUT.McDst.root\n"
              "       converter --root INPUT.root OUTPUT.McDst.root\n"
              "       converter --v1 INPUT.oscar OUTPUT.McDst.root\n"
              "       converter --validate INPUT.oscar [OUTPUT_DIRECTORY]\n"
              "       converter --validate-all INPUT.oscar [OUTPUT_DIRECTORY]\n"
              "Versions: 0 = original, 1 = optimized (default), 2 = inelastic interaction types.\n"
              "--v1 retains its historical meaning: the original algorithm (version 0).\n";
}

CommandLine ParseArguments(int argc, char** argv) {
    CommandLine command;
    bool paths_only = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (!paths_only && (argument == "--help" || argument == "-h")) {
            command.help = true;
            return command;
        }
        if (!paths_only && argument == "--") {
            paths_only = true;
            continue;
        }
        if (!paths_only && argument == "--keep-empty")
            command.conversion.keep_empty_events = true;
        else if (!paths_only && argument == "--root")
            command.root_input = true;
        else if (!paths_only && (argument == "--version" || argument == "--v1")) {
            if (command.version_selected)
                throw std::runtime_error("select a converter version only once");
            command.version_selected = true;
            std::string value = "0";
            if (argument == "--version") {
                if (++index == argc)
                    throw std::runtime_error("missing --version value");
                value = argv[index];
            }
            if (value != "0" && value != "1" && value != "2")
                throw std::runtime_error("--version must be 0, 1 or 2");
            command.conversion.version = static_cast<OscarVersion>(value.front() - '0');
        } else if (!paths_only && argument == "--validate")
            command.validate = true;
        else if (!paths_only && argument == "--validate-all")
            command.validate_all = true;
        else if (!paths_only && !argument.empty() && argument.front() == '-')
            throw std::runtime_error("unknown option: " + argument);
        else
            command.paths.push_back(argument);
    }
    if ((command.validate && command.validate_all) ||
        ((command.validate || command.validate_all) &&
         (command.root_input || command.version_selected ||
          command.conversion.keep_empty_events)) ||
        (command.root_input &&
         (command.version_selected || command.conversion.keep_empty_events)) ||
        (command.conversion.version == OscarVersion::Legacy &&
         command.conversion.keep_empty_events))
        throw std::runtime_error("incompatible conversion modes");
    const bool validation = command.validate || command.validate_all;
    if (validation ? command.paths.empty() || command.paths.size() > 2 : command.paths.size() != 2)
        throw std::runtime_error("incorrect number of paths");
    return command;
}

} // namespace

int main(int argc, char** argv) {
    CommandLine command;
    try {
        command = ParseArguments(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "Argument error: " << error.what() << '\n';
        Usage(std::cerr);
        return 1;
    }
    if (command.help) {
        Usage(std::cout);
        return 0;
    }
    if (command.validate || command.validate_all) {
        const auto directory = command.paths.size() == 2 ? command.paths[1] : "";
        return static_cast<int>(command.validate_all
                                    ? ValidateAllOscarVersions(command.paths[0], directory)
                                    : ValidateOscarVersions(command.paths[0], directory));
    }
    const bool success =
        command.root_input
            ? ConvertRootToMcDst(command.paths[0], command.paths[1])
            : OscarConverter().Convert(command.paths[0], command.paths[1], command.conversion);
    return success ? 0 : 2;
}
