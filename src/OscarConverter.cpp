#include "OscarConverter.h"
#include "Particle.h"

#include "McEvent.h"
#include "McParticle.h"
#include "McRun.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"

#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <locale>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

bool OscarConverter::Convert(const std::string& inputFilename, const std::string& outputFilename,
                             const ConversionOptions& options) {
    constexpr int kInitialParticles = 394;
    constexpr int kMaxParticles = std::numeric_limits<unsigned short>::max();
    constexpr std::string_view kWhitespace = " \t\r\n\f\v";
    enum class Mode { Idle, Initial, History, Final, Skip };
    struct BufferedParticle {
        Particle particle;
        bool written = false;
    };

    std::string temporaryFilename;
    std::size_t lineNumber = 0;
    try {
        std::ifstream input(inputFilename);
        if (!input)
            throw std::runtime_error("cannot open input: " + inputFilename);
        if (outputFilename.empty() || std::filesystem::is_directory(outputFilename))
            throw std::runtime_error("output must be a file path");
        if (std::filesystem::exists(outputFilename) &&
            std::filesystem::equivalent(inputFilename, outputFilename))
            throw std::runtime_error("input and output must be different files");

        std::string temporaryName = outputFilename + ".tmp.XXXXXX";
        const int descriptor = mkstemp(temporaryName.data());
        if (descriptor < 0)
            throw std::runtime_error("cannot create output beside " + outputFilename + ": " +
                                     std::strerror(errno));
        temporaryFilename = std::move(temporaryName);
        if (close(descriptor) != 0)
            throw std::runtime_error("cannot close temporary output: " +
                                     std::string(std::strerror(errno)));

        std::unique_ptr<TFile> outputFile(
            TFile::Open(temporaryFilename.c_str(), "RECREATE", "Oscar to McDst"));
        if (!outputFile || outputFile->IsZombie())
            throw std::runtime_error("cannot create output file: " + outputFilename);
        outputFile->SetCompressionLevel(1);
        TClonesArray events("McEvent", 1), particles("McParticle", 1000);
        auto tree = std::make_unique<TTree>("McDst", "Converted Oscar Data");
        tree->SetDirectory(outputFile.get());
        TClonesArray* eventArray = &events;
        TClonesArray* particleArray = &particles;
        auto* eventBranch = tree->Branch("Event", &eventArray, 65536, 99);
        auto* particleBranch = tree->Branch("Particle", &particleArray, 65536, 99);
        if (!eventBranch || !particleBranch)
            throw std::runtime_error("cannot create McDst branches");
        eventBranch->SetAutoDelete(kFALSE);
        particleBranch->SetAutoDelete(kFALSE);
        std::vector<int> sourceIds;
        if (options.record_source_ids && !tree->Branch("SourceId", &sourceIds))
            throw std::runtime_error("cannot create SourceId branch");

        // Reuse the event object; only its number and impact parameter change.
        auto* event = new (events[0]) McEvent();
        event->setNes(1);
        event->setStepNr(1);
        event->setStepT(200.);
        int child[2] = {-1, -1};
        std::unordered_map<int, BufferedParticle> buffer;
        buffer.reserve(kInitialParticles);
        Mode mode = Mode::Idle;
        int eventNumber = -1, eventEnsemble = 0, remainingParticles = 0;
        std::string line, marker, kind, phase, key;
        std::istringstream fields;
        fields.imbue(std::locale::classic());

        while (std::getline(input, line)) {
            ++lineNumber;
            const auto start = line.find_first_not_of(kWhitespace);
            if (start == std::string::npos)
                continue;

            if (line[start] != '#') {
                if (mode == Mode::Idle)
                    throw std::runtime_error("particle outside an active event");
                if (mode != Mode::History) {
                    if (remainingParticles == 0)
                        throw std::runtime_error("particle block exceeds its declared count");
                    --remainingParticles;
                }
                if (mode == Mode::Skip)
                    continue;
                if (mode == Mode::Final) {
                    // Verify 12 columns, converting only the ID from the final list.
                    std::string_view remainingFields(line), idToken;
                    for (int column = 0; column < 12; ++column) {
                        const auto begin = remainingFields.find_first_not_of(kWhitespace);
                        if (begin == std::string_view::npos)
                            throw std::runtime_error("invalid final particle record");
                        remainingFields.remove_prefix(begin);
                        const auto end = remainingFields.find_first_of(kWhitespace);
                        if (column == 10)
                            idToken = remainingFields.substr(0, end);
                        remainingFields.remove_prefix(
                            end == std::string_view::npos ? remainingFields.size() : end);
                    }
                    if (idToken.front() == '+')
                        idToken.remove_prefix(1);
                    int id = -1;
                    const auto result =
                        std::from_chars(idToken.data(), idToken.data() + idToken.size(), id);
                    if (result.ec != std::errc{} || result.ptr != idToken.data() + idToken.size())
                        throw std::runtime_error("invalid final particle record");

                    const auto found = buffer.find(id);
                    if (found == buffer.end() || found->second.written)
                        continue;
                    Particle& p = found->second.particle;
                    if (id < kInitialParticles)
                        p.t = p.x = p.y = p.z = -999.;
                    const int idx = particles.GetEntriesFast();
                    new (particles[idx]) McParticle(idx, p.pdg, 0, 0, 0, -1, 0, child, p.px, p.py,
                                                    p.pz, p.p0, p.x, p.y, p.z, p.t);
                    if (options.record_source_ids)
                        sourceIds.push_back(id);
                    found->second.written = true;
                    continue;
                }

                fields.clear();
                fields.str(line);
                Particle p;
                if (!(fields >> p.t >> p.x >> p.y >> p.z >> p.mass >> p.p0 >> p.px >> p.py >>
                      p.pz >> p.pdg >> p.id >> p.charge))
                    throw std::runtime_error("invalid OSCAR particle record");
                for (double value : {p.t, p.x, p.y, p.z, p.mass, p.p0, p.px, p.py, p.pz}) {
                    if (!std::isfinite(value))
                        throw std::runtime_error("non-finite particle value");
                }
                buffer.insert_or_assign(p.id, BufferedParticle{p, false});
                continue;
            }

            fields.clear();
            fields.str(line);
            if (!(fields >> marker >> kind) || marker != "#" ||
                (kind != "event" && kind != "interaction"))
                continue;
            if (remainingParticles != 0)
                throw std::runtime_error("particle block ended before its declared count");
            const bool finalSeen = mode == Mode::Final || mode == Mode::Skip;
            if (kind == "interaction") {
                if (mode == Mode::Idle || finalSeen)
                    throw std::runtime_error("interaction outside an active event");
                mode = Mode::History;
                continue;
            }

            int number = -1, ensemble = 0, count = -1;
            if (!(fields >> number >> phase) || number < 0)
                throw std::runtime_error("invalid event header");
            if (phase == "ensemble" && (!(fields >> ensemble >> phase) || ensemble < 0))
                throw std::runtime_error("invalid ensemble header");
            if (!(fields >> count) || count < 0 || count > kMaxParticles)
                throw std::runtime_error("invalid event particle count");

            if (phase == "in") {
                if (mode != Mode::Idle)
                    throw std::runtime_error("new event before previous end");
                eventNumber = number;
                eventEnsemble = ensemble;
                buffer.clear();
                particles.Delete();
                sourceIds.clear();
                if (static_cast<std::size_t>(count) > buffer.bucket_count())
                    buffer.reserve(count);
                remainingParticles = count;
                mode = Mode::Initial;
                continue;
            }
            if (mode == Mode::Idle || number != eventNumber || ensemble != eventEnsemble)
                throw std::runtime_error("event/ensemble header does not match active event");
            if (phase == "out") {
                if (finalSeen)
                    throw std::runtime_error("duplicate final particle block");
                remainingParticles = count;
                mode = count == kInitialParticles && !options.keep_empty_events ? Mode::Skip
                                                                                : Mode::Final;
                continue;
            }
            if (phase != "end")
                throw std::runtime_error("unknown event phase");
            if (!finalSeen)
                throw std::runtime_error("missing final particle block");
            double impactParameter = -1.;
            if (!(fields >> key >> impactParameter) || key != "impact" ||
                !std::isfinite(impactParameter) || impactParameter < 0.)
                throw std::runtime_error("invalid impact parameter");

            const bool skip = mode == Mode::Skip;
            mode = Mode::Idle;
            if (skip)
                continue;
            event->setEventNr(eventNumber);
            event->setB(impactParameter);
            if (tree->Fill() < 0 || outputFile->TestBit(TFile::kWriteError))
                throw std::runtime_error("failed to write event");
        }
        if (input.bad())
            throw std::runtime_error("input read failed");
        if (mode != Mode::Idle)
            throw std::runtime_error("EOF before event end");

        outputFile->cd();
        const auto entries = tree->GetEntries();
        const int runCount =
            entries <= std::numeric_limits<int>::max() ? static_cast<int>(entries) : 0;
        McRun run("SMASH", "Converted from Oscar file", 0, 0, 0., 0, 0, 0., 0., 0., -1, 0., 0., 0.,
                  runCount);
        if (run.Write() <= 0 || outputFile->Write() <= 0 || outputFile->TestBit(TFile::kWriteError))
            throw std::runtime_error("failed to finalize ROOT output");
        tree.reset();
        outputFile->Close();
        if (outputFile->TestBit(TFile::kWriteError))
            throw std::runtime_error("failed to close ROOT output");
        outputFile.reset();
        std::error_code error;
        std::filesystem::rename(temporaryFilename, outputFilename, error);
        if (error)
            throw std::runtime_error("cannot replace output: " + error.message());
        temporaryFilename.clear();
        std::cout << "Conversion completed. Total events: " << entries << '\n';
        return true;
    } catch (const std::exception& error) {
        if (!temporaryFilename.empty()) {
            std::error_code ignored;
            std::filesystem::remove(temporaryFilename, ignored);
        }
        std::cerr << "Conversion failed: ";
        if (lineNumber != 0)
            std::cerr << "line " << lineNumber << ": ";
        std::cerr << error.what() << '\n';
        return false;
    }
}
