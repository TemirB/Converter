#include "OscarConverter.h"
#include "McDstWriter.h"
#include "OscarReader.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

bool OscarConverter::Convert(const std::string& inputFilename, const std::string& outputFilename,
                             const ConversionOptions& options) {
    try {
        if (std::filesystem::weakly_canonical(inputFilename) ==
            std::filesystem::weakly_canonical(outputFilename)) {
            throw std::runtime_error("input and output must be different files");
        }
        std::ifstream input(inputFilename);
        if (!input)
            throw std::runtime_error("cannot open input: " + inputFilename);
        const std::string definitions =
            "Npart=initial nucleons with any hadronic interaction; "
            "Ncoll=binary hadronic scatterings of all species; "
            "NN counters=transport history, not MC Glauber; "
            "McParticle.status: 1=spectator, 0=other; "
            "spectator coordinates=final output; other coordinates=last hadronic interaction; "
            "momenta=final output; genealogy=unknown; event selection=" +
            std::string(options.keep_empty_events ? "all" : "projectile_target_scattering");
        McDstWriter writer(outputFilename, definitions);
        if (!writer.Error().empty())
            throw std::runtime_error(writer.Error());
        OscarReader reader(input);
        OscarEvent event;
        long long input_events = 0, skipped_events = 0;
        while (reader.ReadNext(event)) {
            ++input_events;
            if (!options.keep_empty_events && !event.projectile_target_scattering) {
                ++skipped_events;
                continue;
            }
            if (!writer.WriteEvent(event))
                throw std::runtime_error(writer.Error());
        }
        if (!reader.Error().empty())
            throw std::runtime_error(reader.Error());
        if (!writer.Finish(input_events, skipped_events))
            throw std::runtime_error(writer.Error());
        std::cout << "Converted " << input_events - skipped_events << " of " << input_events
                  << " events; skipped " << skipped_events << " empty events.\n";
        return true;
    } catch (const std::exception& error) {
        std::cerr << "Conversion failed: " << error.what() << '\n';
        return false;
    }
}
