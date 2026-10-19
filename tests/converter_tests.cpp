#include "McEvent.h"
#include "McParticle.h"
#include "OscarConverter.h"
#include "OscarReader.h"
#include "RootConverter.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

const std::string kHeader =
    "#!OSCAR2013 full_event_history t x y z mass p0 px py pz pdg ID charge\n";

void Check(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

std::string ParticleLine(int id, int pdg, double time = 0., double x = 0., double px = 0.) {
    const double mass = pdg == 211 ? .13957 : (pdg == 2214 ? 1.232 : .938);
    std::ostringstream line;
    line << time << ' ' << x << " 0 0 " << mass << ' ' << std::sqrt(mass * mass + px * px) << ' '
         << px << " 0 0 " << pdg << ' ' << id << " 1\n";
    return line.str();
}

std::string Initial() {
    // Sparse IDs deliberately break the old ID < 394 assumption.
    return "# event 0 in 3\n" + ParticleLine(10, 2212) + ParticleLine(1001, 2112) +
           ParticleLine(2002, 2212);
}

std::string End(const std::string& flag = "yes") {
    return "# event 0 end 0 impact 7.5 scattering_projectile_target " + flag + "\n";
}

std::string EmptyEvent() {
    return kHeader + Initial() + "# event 0 out 3\n" + ParticleLine(10, 2212, 250.) +
           ParticleLine(1001, 2112, 250.) + ParticleLine(2002, 2212, 250.) + End("no");
}

std::string ElasticEvent() {
    return kHeader + Initial() + "# interaction in 2 out 2 rho 0 weight 1 partial 1 type 1\n" +
           ParticleLine(10, 2212, 10.) + ParticleLine(1001, 2112, 10.) +
           ParticleLine(10, 2212, 10., 1., .2) + ParticleLine(1001, 2112, 10., 2., -.2) +
           "# event 0 out 3\n" + ParticleLine(10, 2212, 250.) + ParticleLine(1001, 2112, 250.) +
           ParticleLine(2002, 2212, 250., 100.) + End();
}

std::string InelasticInverseEvent() {
    return kHeader + Initial() + "# interaction in 2 out 2 rho 0 weight 1 partial 1 type 3\n" +
           ParticleLine(10, 2212, 10.) + ParticleLine(1001, 2112, 10.) +
           ParticleLine(3000, 2212, 10.) + ParticleLine(3001, 2214, 10.) +
           "# interaction in 2 out 2 rho 0 weight 1 partial 1 type 3\n" +
           ParticleLine(3000, 2212, 20.) + ParticleLine(3001, 2214, 20.) +
           ParticleLine(3002, 2212, 20.) + ParticleLine(3003, 2112, 20.) + "# event 0 out 3\n" +
           ParticleLine(2002, 2212, 250.) + ParticleLine(3002, 2212, 250.) +
           ParticleLine(3003, 2112, 250.) + End();
}

std::string PionEvent() {
    return kHeader + Initial() + "# interaction in 2 out 3 rho 0 weight 1 partial 1 type 4\n" +
           ParticleLine(10, 2212, 10.) + ParticleLine(1001, 2112, 10.) +
           ParticleLine(3000, 2212, 10.) + ParticleLine(3001, 2112, 10.) +
           ParticleLine(3002, 211, 10., 2., .1) +
           "# interaction in 2 out 2 rho 0 weight 1 partial 1 type 1\n" +
           ParticleLine(3000, 2212, 20.) + ParticleLine(3002, 211, 20., 3., .1) +
           ParticleLine(3000, 2212, 20.) + ParticleLine(3002, 211, 20., 4., .2) +
           "# event 0 out 4\n" + ParticleLine(3002, 211, 250., 40., .3) +
           ParticleLine(3001, 2112, 250.) + ParticleLine(2002, 2212, 250., 100.) +
           ParticleLine(3000, 2212, 250.) + End();
}

std::string Replace(std::string input, const std::string& from, const std::string& to) {
    std::size_t at = 0;
    bool replaced = false;
    while ((at = input.find(from, at)) != std::string::npos) {
        input.replace(at, from.size(), to);
        at += to.size();
        replaced = true;
    }
    Check(replaced, "test replacement was not found");
    return input;
}

OscarEvent Read(const std::string& contents) {
    std::istringstream input(contents);
    OscarReader reader(input);
    OscarEvent event;
    Check(reader.ReadNext(event), reader.Error());
    OscarEvent next;
    Check(!reader.ReadNext(next) && reader.Error().empty(),
          "unexpected trailing event or error: " + reader.Error());
    return event;
}

void Reject(const std::string& contents, const std::string& expected) {
    std::istringstream input(contents);
    OscarReader reader(input);
    OscarEvent event;
    while (reader.ReadNext(event)) {
    }
    Check(reader.Error().find(expected) != std::string::npos,
          "expected '" + expected + "', received '" + reader.Error() + "'");
}

struct TempDirectory {
    fs::path path =
        fs::temp_directory_path() / ("converter-tests-" + std::to_string(std::random_device{}()));
    TempDirectory() {
        Check(fs::create_directory(path), "cannot create test directory");
    }
    ~TempDirectory() {
        std::error_code ignored;
        fs::remove_all(path, ignored);
    }
};

void Write(const fs::path& path, const std::string& contents) {
    std::ofstream file(path);
    file << contents;
    Check(bool(file), "cannot write fixture");
}

void CheckRoot(const fs::path& path, int expected_entries, int expected_participants = 2,
               int expected_collisions = 2, int expected_spectators = 1, bool pion = true) {
    TFile file(path.c_str(), "READ");
    Check(!file.IsZombie(), "ROOT output cannot be reopened");
    auto* tree = dynamic_cast<TTree*>(file.Get("McDst"));
    Check(tree && tree->GetEntries() == expected_entries, "incorrect output event count");
    auto* definitions = dynamic_cast<TNamed*>(file.Get("ConverterDefinitions"));
    Check(definitions &&
              std::string(definitions->GetTitle()).find("not MC Glauber") != std::string::npos,
          "missing participant definitions");
    if (!expected_entries)
        return;
    TClonesArray events("McEvent", 1), particles("McParticle", 10);
    TClonesArray* event_pointer = &events;
    TClonesArray* particle_pointer = &particles;
    std::vector<int> ids;
    auto* id_pointer = &ids;
    int spectators = -1;
    tree->SetBranchAddress("Event", &event_pointer);
    tree->SetBranchAddress("Particle", &particle_pointer);
    tree->SetBranchAddress("OscarParticleId", &id_pointer);
    tree->SetBranchAddress("SpectatorNucleons", &spectators);
    for (int entry = 0; entry < expected_entries; ++entry) {
        Check(tree->GetEntry(entry) > 0, "cannot read output entry");
        auto* event = static_cast<McEvent*>(events.At(0));
        Check(event->npart() == expected_participants && event->ncoll() == expected_collisions,
              "incorrect ROOT participant/collision count");
        Check(event->stepT() == 250. && event->b() == 7.5, "incorrect ROOT event time/impact");
        Check(spectators == expected_spectators, "incorrect ROOT spectator count");
        Check(ids.size() == static_cast<std::size_t>(particles.GetEntries()),
              "missing original IDs");
        bool found_pion = false, found_spectator = false;
        for (int i = 0; i < particles.GetEntries(); ++i) {
            auto* particle = static_cast<McParticle*>(particles.At(i));
            Check(particle->firstChild() == 65535 && particle->lastChild() == 65535 &&
                      particle->parent() == 65535,
                  "fake genealogy remained");
            if (particle->status() == 1) {
                found_spectator = true;
                Check(particle->t() == 250., "spectator coordinates were replaced by a sentinel");
            }
            if (particle->pdg() == 211) {
                found_pion = true;
                Check(particle->t() == 20. && particle->x() == 4. &&
                          std::abs(particle->px() - .3) < 1e-6,
                      "incorrect pion emission/final momentum");
                Check(ids.at(i) == 3002, "OSCAR pion ID was lost");
            }
        }
        Check(found_spectator == (expected_spectators > 0) && found_pion == pion,
              "incorrect output particle classification");
    }
    tree->ResetBranchAddresses();
}

void TestRootInput(const fs::path& path, int count = 2) {
    TFile input(path.c_str(), "RECREATE");
    TTree tree("Particles", "legacy ROOT fixture");
    double impact = 7.5;
    std::vector<double> t{20., 20.}, x{1., 2.}, y{0., 0.}, z{0., 0.};
    std::vector<double> mass{.13957, .13957}, energy{.3, .3}, px{.2, -.2}, py{0., 0.}, pz{0., 0.};
    std::vector<int> pdg{211, -211}, ids{5, 6}, charge{1, -1};
    tree.Branch("npart", &count);
    tree.Branch("impactParameter", &impact);
    tree.Branch("t", &t);
    tree.Branch("x", &x);
    tree.Branch("y", &y);
    tree.Branch("z", &z);
    tree.Branch("mass", &mass);
    tree.Branch("p0", &energy);
    tree.Branch("px", &px);
    tree.Branch("py", &py);
    tree.Branch("pz", &pz);
    tree.Branch("pdg", &pdg);
    tree.Branch("ID", &ids);
    tree.Branch("charge", &charge);
    tree.Fill();
    tree.Write();
}

} // namespace

int main() {
    try {
        const auto event = Read(PionEvent());
        Check(event.initial_nucleons == 3 && event.participants == 2 && event.spectators == 1,
              "sparse initial IDs were misclassified");
        Check(event.binary_collisions == 2 && event.nn_inelastic_collisions == 1 &&
                  event.nn_elastic_collisions == 0,
              "incorrect transport collision types");
        Check(event.particles.back().t == 20. && event.particles.back().x == 4. &&
                  event.particles.back().px == .3,
              "emission coordinates replaced final momentum");
        const auto elastic = Read(ElasticEvent());
        Check(elastic.participants == 2 && elastic.nn_elastic_collisions == 1 &&
                  elastic.projectile_target_scattering,
              "elastic participants were lost");
        const auto inverse = Read(InelasticInverseEvent());
        Check(inverse.particles.size() == 3 && inverse.participants == 2 &&
                  inverse.nn_inelastic_collisions == 1,
              "inverse reaction history was lost");
        const auto empty = Read(EmptyEvent());
        Check(empty.participants == 0 && empty.spectators == 3 &&
                  !empty.projectile_target_scattering,
              "empty event was misclassified");
        const auto ensemble = Read(Replace(PionEvent(), "# event 0 ", "# event 0 ensemble 7 "));
        Check(ensemble.ensemble == 7 && ensemble.participants == 2,
              "ensemble headers were not parsed");
        Check(Read(Replace(ElasticEvent(), "scattering_projectile_target yes", "empty no"))
                  .projectile_target_scattering,
              "legacy empty flag is inverted");
        Check(Read(Replace(ElasticEvent(), " scattering_projectile_target yes", ""))
                  .projectile_target_scattering,
              "legacy scattering fallback failed");
        Check(Read(Replace(PionEvent(), "# event", "  # event")).participants == 2,
              "leading whitespace broke headers");

        // Extended records must have exactly the declared number of columns.
        std::istringstream plain(PionEvent());
        std::ostringstream extended;
        std::string extended_line;
        while (std::getline(plain, extended_line)) {
            if (extended_line.rfind("#!", 0) == 0)
                extended_line += " Ncoll";
            else if (!extended_line.empty() && extended_line.front() != '#')
                extended_line += " 0";
            extended << extended_line << '\n';
        }
        Check(Read(extended.str()).participants == 2, "extended OSCAR records failed");

        const auto decay = Read(
            kHeader + Initial() + "# interaction in 2 out 2 rho 0 weight 1 partial 1 type 3\n" +
            ParticleLine(10, 2212, 10.) + ParticleLine(1001, 2112, 10.) +
            ParticleLine(3000, 2212, 10.) + ParticleLine(3001, 2214, 10.) +
            "# interaction in 1 out 2 rho 0 weight 1 partial 1 type 5\n" +
            ParticleLine(3001, 2214, 250.) + ParticleLine(3002, 2112, 250., 2.) +
            ParticleLine(3003, 211, 250., 2., .2) + "# event 0 out 4\n" +
            ParticleLine(2002, 2212, 250.) + ParticleLine(3000, 2212, 250.) +
            ParticleLine(3002, 2112, 250.) + ParticleLine(3003, 211, 250., 2., .2) + End());
        Check(decay.binary_collisions == 1 && decay.particles.back().t == 250.,
              "final decays affected binary counts or emission times");

        const auto wall = Read(kHeader + "# event 0 in 1\n" + ParticleLine(10, 2212) +
                               "# interaction in 1 out 1 rho 0 weight 0 partial 0 type 6\n" +
                               ParticleLine(10, 2212, 10.) + ParticleLine(10, 2212, 10., 1.) +
                               "# event 0 out 1\n" + ParticleLine(10, 2212, 250.) + End("no"));
        Check(wall.participants == 0 && wall.binary_collisions == 0 && wall.spectators == 1,
              "wall crossing created a participant");

        Reject(PionEvent().substr(0, PionEvent().find("# event 0 end")), "EOF before event end");
        Reject(Replace(PionEvent(), "# event 0 out 4", "# event 0 out 5"), "declared count");
        Reject(Replace(PionEvent(), "# event 0 out 4", "# event 0 out -1"), "final particle count");
        Reject(Replace(PionEvent(), ParticleLine(3002, 211, 250., 40., .3),
                       ParticleLine(9999, 211, 250., 40., .3)),
               "live history");
        Reject(Replace(PionEvent(), "# event 0 end", "# event 1 end"), "does not match");
        Reject(Replace(PionEvent(), "impact 7.5", "impact -1"), "impact parameter");
        Reject(Replace(PionEvent(), "scattering_projectile_target yes",
                       "scattering_projectile_target maybe"),
               "scattering flag");
        Reject(Replace(PionEvent(), "#!OSCAR2013", "#!OSCAR1999"), "expected OSCAR2013");
        Reject(Replace(PionEvent(), "t x y z mass", "x t y z mass"), "particle columns");
        Reject(Replace(PionEvent(), "type 4", "type 3"), "interaction type");
        Reject(kHeader + Initial() + End(), "missing final particle block");
        Reject(Replace(EmptyEvent(), "scattering_projectile_target no",
                       "scattering_projectile_target yes"),
               "without interaction history");
        Reject(kHeader + Initial() + Initial(), "interleaved ensembles");
        Reject(kHeader + "# event 0 in 1\nnot a particle\n", "particle record");
        Reject(kHeader + "# event 0 in 2\n" + ParticleLine(4, 2212) + ParticleLine(4, 2212),
               "duplicate particle ID");
        Reject(kHeader + "# event 0 in 1\n" + ParticleLine(4, 0), "PDG");

        TempDirectory temp;
        const auto source = temp.path / "input.oscar", output = temp.path / "output.root";
        Write(source, PionEvent());
        Check(OscarConverter().Convert(source.string(), output.string()),
              "valid conversion failed");
        CheckRoot(output, 1);
        Write(source, ElasticEvent());
        Check(OscarConverter().Convert(source.string(), output.string()),
              "elastic event was rejected");
        CheckRoot(output, 1, 2, 1, 1, false);
        Write(source, InelasticInverseEvent());
        Check(OscarConverter().Convert(source.string(), output.string()),
              "inelastic event with original multiplicity was rejected");
        CheckRoot(output, 1, 2, 2, 1, false);
        Write(source, EmptyEvent());
        Check(OscarConverter().Convert(source.string(), output.string()),
              "empty conversion failed");
        CheckRoot(output, 0);
        Check(OscarConverter().Convert(source.string(), output.string(), {true}),
              "keep-empty conversion failed");
        CheckRoot(output, 1, 0, 0, 3, false);
        Write(source, Replace(PionEvent(), "# event 0 ", "# event 0 ensemble 7 "));
        Check(OscarConverter().Convert(source.string(), output.string()),
              "ensemble conversion failed");
        CheckRoot(output, 1);
        Write(source, PionEvent() +
                          Replace(PionEvent().substr(kHeader.size()), "# event 0 ", "# event 9 "));
        Check(OscarConverter().Convert(source.string(), output.string()),
              "multiple event conversion failed");
        CheckRoot(output, 2);
        const auto bytes_before = fs::file_size(output);
        Write(source, PionEvent() + Initial());
        Check(!OscarConverter().Convert(source.string(), output.string()),
              "corrupt input returned success");
        Check(fs::file_size(output) == bytes_before, "failed conversion overwrote existing output");
        CheckRoot(output, 2);
        Check(!OscarConverter().Convert(source.string(), source.string()),
              "input/output identity was accepted");
        Check(!OscarConverter().Convert((temp.path / "missing.oscar").string(), output.string()),
              "missing input was accepted");
        for (const auto& entry : fs::directory_iterator(temp.path))
            Check(entry.path().filename().string().find(".tmp.") == std::string::npos,
                  "staging file leaked");
        TestRootInput(temp.path / "legacy.root");
        Check(ConvertRootToMcDst((temp.path / "legacy.root").string(),
                                 (temp.path / "legacy-out.root").string()),
              "legacy ROOT conversion failed");
        TFile legacy((temp.path / "legacy-out.root").c_str(), "READ");
        auto* legacy_tree = dynamic_cast<TTree*>(legacy.Get("McDst"));
        Check(legacy_tree && legacy_tree->GetEntries() == 1, "legacy ROOT event was lost");
        TClonesArray legacy_events("McEvent", 1);
        auto* legacy_events_pointer = &legacy_events;
        legacy_tree->SetBranchAddress("Event", &legacy_events_pointer);
        legacy_tree->GetEntry(0);
        const auto* legacy_event = static_cast<McEvent*>(legacy_events.At(0));
        Check(legacy_event->npart() == -1 && legacy_event->ncoll() == -1 &&
                  legacy_event->stepT() == 20.,
              "legacy input invented participant counts or output time");
        legacy_tree->ResetBranchAddresses();
        const auto legacy_bytes = fs::file_size(temp.path / "legacy-out.root");
        TestRootInput(temp.path / "invalid-legacy.root", 3);
        Check(!ConvertRootToMcDst((temp.path / "invalid-legacy.root").string(),
                                  (temp.path / "legacy-out.root").string()),
              "mismatched ROOT vectors were accepted");
        Check(fs::file_size(temp.path / "legacy-out.root") == legacy_bytes,
              "invalid ROOT conversion overwrote existing output");
        std::cout << "Parser, physics classification, ROOT round trips and atomic failure checks "
                     "passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "TEST FAILURE: " << error.what() << '\n';
        return 1;
    }
}
