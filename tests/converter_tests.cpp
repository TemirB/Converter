#include "McEvent.h"
#include "McParticle.h"
#include "OscarConverter.h"
#include "OscarReader.h"
#include "RootConverter.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
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
    // Sparse IDs exercise the fixed ID boundary independently of array positions.
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

std::string AuSpectatorEvent(bool inelastic) {
    auto pdg = [](int id) { return id % 197 < 79 ? 2212 : 2112; };
    std::string contents = kHeader + "# event 0 in 394\n";
    for (int id = 0; id < 394; ++id)
        contents += ParticleLine(id, pdg(id));
    // Both nuclei's protons keep their IDs through any number of elastic scatterings.
    for (double time : {10., 20.}) {
        contents += "# interaction in 2 out 2 rho 0 weight 1 partial 1 type 1\n";
        contents += ParticleLine(0, 2212, time) + ParticleLine(197, 2212, time);
        contents += ParticleLine(0, 2212, time, 1., .2) + ParticleLine(197, 2212, time, 2., -.2);
    }
    if (inelastic) {
        contents += "# interaction in 2 out 2 rho 0 weight 1 partial 1 type 3\n";
        contents += ParticleLine(0, 2212, 30.) + ParticleLine(197, 2212, 30.);
        contents += ParticleLine(394, 2212, 30.) + ParticleLine(395, 2214, 30.);
    }
    contents += "# event 0 out 394\n";
    for (int id = 0; id < 394; ++id) {
        if (inelastic && (id == 0 || id == 197))
            continue;
        const double px = id == 0 ? .2 : (id == 197 ? -.2 : 0.);
        contents += ParticleLine(id, pdg(id), 250., 100., px);
    }
    if (inelastic)
        contents += ParticleLine(394, 2212, 250.) + ParticleLine(395, 2214, 250.);
    return contents + End();
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

struct ExpectedParticle {
    int pdg;
    double t, x = 0., px = 0.;
};

struct ExpectedEvent {
    int number;
    std::vector<ExpectedParticle> particles;
};

void CheckRoot(const fs::path& path, const std::vector<ExpectedEvent>& expected) {
    TFile file(path.c_str(), "READ");
    Check(!file.IsZombie(), "ROOT output cannot be reopened");
    auto* tree = dynamic_cast<TTree*>(file.Get("McDst"));
    Check(tree && tree->GetEntries() == static_cast<Long64_t>(expected.size()),
          "incorrect output event count");
    Check(tree->GetListOfBranches()->GetEntries() == 2 && tree->GetBranch("Event") &&
              tree->GetBranch("Particle"),
          "legacy output must contain only Event and Particle branches");
    TClonesArray events("McEvent", 1), particles("McParticle", 10);
    TClonesArray* event_pointer = &events;
    TClonesArray* particle_pointer = &particles;
    tree->SetBranchAddress("Event", &event_pointer);
    tree->SetBranchAddress("Particle", &particle_pointer);
    for (std::size_t entry = 0; entry < expected.size(); ++entry) {
        Check(tree->GetEntry(static_cast<Long64_t>(entry)) > 0, "cannot read output entry");
        Check(events.GetEntries() == 1, "event array was not cleared between events");
        auto* event = static_cast<McEvent*>(events.At(0));
        Check(event->eventNr() == expected[entry].number && event->stepT() == 200. &&
                  event->b() == 7.5 && event->phi() == 0. && event->numberOfSteps() == 1 &&
                  event->stepNumber() == 1,
              "legacy event metadata changed");
        Check(event->npart() == -1 && event->ncoll() == -1,
              "legacy output invented participant/collision counts");
        TString comment;
        event->comment(comment);
        Check(comment.IsNull(), "legacy event comment must be empty");
        const auto& wanted = expected[entry].particles;
        Check(particles.GetEntries() == static_cast<int>(wanted.size()),
              "incorrect particle intersection or particle array was not cleared");
        std::vector<bool> matched(wanted.size(), false);
        for (int i = 0; i < particles.GetEntries(); ++i) {
            auto* particle = static_cast<McParticle*>(particles.At(i));
            Check(particle->index() == i && particle->status() == 0 && particle->parent() == 0 &&
                      particle->parentDecay() == 0 && particle->mate() == 65535 &&
                      particle->decay() == 0 && particle->firstChild() == 65535 &&
                      particle->lastChild() == 65535,
                  "legacy particle indices, status or genealogy changed");
            bool found = false;
            for (std::size_t j = 0; j < wanted.size(); ++j) {
                const auto& p = wanted[j];
                if (matched[j] || particle->pdg() != p.pdg || particle->t() != p.t ||
                    particle->x() != p.x || std::abs(particle->px() - p.px) > 1e-6)
                    continue;
                const double yz = p.t == -999. ? -999. : 0.;
                Check(particle->y() == yz && particle->z() == yz && particle->py() == 0. &&
                          particle->pz() == 0.,
                      "legacy particle coordinates or momenta changed");
                const double mass = p.pdg == 211 ? .13957 : (p.pdg == 2214 ? 1.232 : .938);
                Check(std::abs(particle->e() - std::sqrt(mass * mass + p.px * p.px)) < 1e-5,
                      "energy did not come from the buffered history particle");
                matched[j] = found = true;
                break;
            }
            Check(found, "unexpected particle or history coordinates/momenta were overwritten");
        }
    }
    tree->ResetBranchAddresses();
}

std::string Multiplicity394Event() {
    std::string contents = kHeader + "# event 0 in 394\n";
    for (int id = 0; id < 394; ++id)
        contents += ParticleLine(id, 2212);
    contents += "# event 0 out 394\n";
    for (int id = 0; id < 394; ++id)
        contents += ParticleLine(id, 2212, 250., 100., .5);
    return contents + End();
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
        Check(event.initial_nucleons == 3 && event.participants == 3 && event.spectators == 0,
              "reader did not apply the fixed spectator ID boundary");
        Check(event.binary_collisions == 2 && event.nn_inelastic_collisions == 1 &&
                  event.nn_elastic_collisions == 0,
              "incorrect transport collision types");
        Check(event.particles.back().t == 20. && event.particles.back().x == 4. &&
                  event.particles.back().px == .3,
              "emission coordinates replaced final momentum");
        const auto elastic = Read(ElasticEvent());
        Check(elastic.participants == 2 && elastic.spectators == 1 &&
                  elastic.particles.front().spectator && elastic.nn_elastic_collisions == 1 &&
                  elastic.projectile_target_scattering,
              "elastic scattering changed the ID-based spectator selection");
        const auto au_elastic = Read(AuSpectatorEvent(false));
        Check(au_elastic.initial_nucleons == 394 && au_elastic.participants == 0 &&
                  au_elastic.spectators == 394 && au_elastic.binary_collisions == 2 &&
                  au_elastic.nn_elastic_collisions == 2 &&
                  std::all_of(au_elastic.particles.begin(), au_elastic.particles.end(),
                              [](const Particle& p) { return p.spectator; }),
              "repeated elastic scattering removed Au spectators");
        const auto au_inelastic = Read(AuSpectatorEvent(true));
        Check(au_inelastic.participants == 2 && au_inelastic.spectators == 392 &&
                  au_inelastic.binary_collisions == 3 &&
                  au_inelastic.nn_inelastic_collisions == 1 &&
                  std::all_of(au_inelastic.particles.begin(), au_inelastic.particles.end(),
                              [](const Particle& p) { return p.spectator == (p.id < 394); }),
              "reader did not distinguish retained initial IDs from inelastic products");
        const auto boundary_ids =
            Read(kHeader + "# event 0 in 2\n" + ParticleLine(393, 2212) + ParticleLine(394, 2212) +
                 "# event 0 out 2\n" + ParticleLine(393, 2212, 250.) +
                 ParticleLine(394, 2212, 250.) + End("no"));
        Check(boundary_ids.spectators == 1 && boundary_ids.participants == 1 &&
                  boundary_ids.particles[0].spectator && !boundary_ids.particles[1].spectator,
              "spectator ID boundary must include 393 and exclude 394");
        const auto inverse = Read(InelasticInverseEvent());
        Check(inverse.particles.size() == 3 && inverse.participants == 3 &&
                  inverse.nn_inelastic_collisions == 1,
              "inverse reaction history was lost");
        const auto empty = Read(EmptyEvent());
        Check(empty.participants == 2 && empty.spectators == 1 &&
                  !empty.projectile_target_scattering,
              "spectator selection unexpectedly depends on interaction history");
        const auto ensemble = Read(Replace(PionEvent(), "# event 0 ", "# event 0 ensemble 7 "));
        Check(ensemble.ensemble == 7 && ensemble.participants == 3,
              "ensemble headers were not parsed");
        Check(Read(Replace(ElasticEvent(), "scattering_projectile_target yes", "empty no"))
                  .projectile_target_scattering,
              "legacy empty flag is inverted");
        Check(Read(Replace(ElasticEvent(), " scattering_projectile_target yes", ""))
                  .projectile_target_scattering,
              "legacy scattering fallback failed");
        Check(Read(Replace(PionEvent(), "# event", "  # event")).participants == 3,
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
        Check(Read(extended.str()).participants == 3, "extended OSCAR records failed");

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
        const std::vector<ExpectedParticle> pion_particles = {
            {211, 20., 4., .2}, {2112, 10.}, {2212, 0.}, {2212, 20.}};
        const std::vector<ExpectedParticle> empty_particles = {
            {2212, -999., -999.}, {2112, 0.}, {2212, 0.}};
        Write(source, PionEvent());
        Check(OscarConverter().Convert(source.string(), output.string()),
              "valid conversion failed");
        CheckRoot(output, {{0, pion_particles}});
        Write(source, ElasticEvent());
        Check(OscarConverter().Convert(source.string(), output.string()),
              "elastic event was rejected");
        CheckRoot(output, {{0, {{2212, -999., -999., .2}, {2112, 10., 2., -.2}, {2212, 0.}}}});
        Write(source, InelasticInverseEvent());
        Check(OscarConverter().Convert(source.string(), output.string()),
              "inelastic event with original multiplicity was rejected");
        CheckRoot(output, {{0, {{2212, 0.}, {2212, 20.}, {2112, 20.}}}});
        Write(source, EmptyEvent());
        Check(OscarConverter().Convert(source.string(), output.string()),
              "scattering flag incorrectly rejected a non-394 event");
        CheckRoot(output, {{0, empty_particles}});
        Check(OscarConverter().Convert(source.string(), output.string(), {true}),
              "keep-empty conversion failed");
        CheckRoot(output, {{0, empty_particles}});

        const std::string count394 = Multiplicity394Event();
        Write(source, count394);
        Check(OscarConverter().Convert(source.string(), output.string()),
              "394-particle conversion failed");
        CheckRoot(output, {});
        Check(OscarConverter().Convert(source.string(), output.string(), {true}),
              "keep-empty did not retain a 394-particle event");
        CheckRoot(output, {{0, std::vector<ExpectedParticle>(394, {2212, -999., -999.})}});

        // No format header, low-ID non-nucleon, repeated/unknown final IDs, and
        // interaction counts that deliberately differ from the records below.
        const std::string boundary =
            "# event 0 ensemble 4 in 3\n" + ParticleLine(393, 211, 1., 2., .1) +
            ParticleLine(394, 2212, 1., 3., .2) + ParticleLine(1000, 2112, 1., 4., .3) +
            "#\n# a service comment\n# interaction in 2 out 3 rho 0 type 4\n" +
            ParticleLine(393, 211, 2., 5., .4) + ParticleLine(394, 2212, 3., 6., .5) +
            "#\n# another service comment\n# event 0 ensemble 4 out 4\n" +
            ParticleLine(394, 2212, 250., 60., .9) + ParticleLine(393, 211, 250., 50., .8) +
            ParticleLine(394, 2212, 250., 70., .7) + ParticleLine(9999, 211, 250.) +
            "# event 0 ensemble 4 end 0 impact 7.5 scattering_projectile_target no\n";
        Write(source, boundary);
        Check(OscarConverter().Convert(source.string(), output.string()),
              "headerless history/intersection conversion failed");
        CheckRoot(output, {{0, {{211, -999., -999., .4}, {2212, 3., 6., .5}}}});

        Write(source, extended.str());
        Check(OscarConverter().Convert(source.string(), output.string()),
              "extended record conversion failed");
        CheckRoot(output, {{0, pion_particles}});
        Write(source, Replace(PionEvent(), "# event 0 ", "# event 0 ensemble 7 "));
        Check(OscarConverter().Convert(source.string(), output.string()),
              "ensemble conversion failed");
        CheckRoot(output, {{0, pion_particles}});

        // The second out block includes an ID that existed only in the first
        // event. Neither that ID nor the first event's clone-array slots may leak.
        const std::string second = "# event 9 in 1\n" + ParticleLine(394, 2112, 2., 6., .4) +
                                   "# event 9 out 2\n" + ParticleLine(3002, 211, 250.) +
                                   ParticleLine(394, 2112, 250., 60., .9) +
                                   "# event 9 end 0 impact 7.5 scattering_projectile_target no\n";
        Write(source, PionEvent() + second);
        Check(OscarConverter().Convert(source.string(), output.string()),
              "multiple event conversion failed");
        const std::vector<ExpectedEvent> multiple = {{0, pion_particles},
                                                     {9, {{2112, 2., 6., .4}}}};
        CheckRoot(output, multiple);
        const auto bytes_before = fs::file_size(output);
        auto reject_conversion = [&](const std::string& contents, const std::string& message) {
            Write(source, contents);
            Check(!OscarConverter().Convert(source.string(), output.string()), message);
            Check(fs::file_size(output) == bytes_before,
                  "failed conversion overwrote existing output");
            CheckRoot(output, multiple);
        };
        reject_conversion(PionEvent() + Initial(), "truncated event returned success");
        reject_conversion(Replace(PionEvent(), "# event 0 in 3", "# event 0 in 4"),
                          "incomplete initial block returned success");
        reject_conversion(Replace(PionEvent(), "# event 0 out 4", "# event 0 out 5"),
                          "incomplete final block returned success");
        reject_conversion(Replace(PionEvent(), "# event 0 out 4", "# event 0 out 3"),
                          "extra final particle returned success");
        reject_conversion(Replace(PionEvent(), "# event 0 end", "# event 1 end"),
                          "mismatched event number returned success");
        reject_conversion(Replace(PionEvent(), "# event 0 out", "# event 0 ensemble 2 out"),
                          "mismatched ensemble returned success");
        reject_conversion(
            Replace(PionEvent(), ParticleLine(3002, 211, 20., 4., .2), "not a particle\n"),
            "invalid interaction particle returned success");
        reject_conversion(Replace(PionEvent(), ParticleLine(3002, 211, 250., 40., .3),
                                  "250 40 0 0 .13957 .3 .3 0 0 211 3002\n"),
                          "short final particle record returned success");
        reject_conversion(Replace(count394, ParticleLine(393, 2212, 250., 100., .5), ""),
                          "incomplete skipped final block returned success");
        reject_conversion(Replace(count394, End(), ParticleLine(394, 2212, 250.) + End()),
                          "extra skipped final particle returned success");
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
        std::cout << "Parser, legacy conversion, ROOT round trips and atomic failure checks "
                     "passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "TEST FAILURE: " << error.what() << '\n';
        return 1;
    }
}
