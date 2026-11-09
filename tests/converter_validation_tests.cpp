#include "ConverterValidation.h"
#include "McEvent.h"
#include "McParticle.h"
#include "OscarConverter.h"
#include "OscarConverterV1.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TNamed.h"
#include "TTree.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

void Check(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

struct TempDirectory {
    fs::path path = fs::temp_directory_path() /
                    ("converter validation " + std::to_string(std::random_device{}()));
    TempDirectory() {
        fs::create_directories(path);
    }
    ~TempDirectory() {
        std::error_code ignored;
        fs::remove_all(path, ignored);
    }
};

struct FixtureParticle {
    int id, pdg = 211;
    std::array<double, 6> values = {1., 2., 3., 4., 5., 6.};
};
struct FixtureEvent {
    int number;
    double impact;
    std::vector<FixtureParticle> particles;
};

void WriteMcDst(const fs::path& path, const std::vector<FixtureEvent>& input) {
    TFile file(path.c_str(), "RECREATE");
    TClonesArray events("McEvent", 1), particles("McParticle", 10);
    TClonesArray* eventPointer = &events;
    TClonesArray* particlePointer = &particles;
    std::vector<int> ids;
    TTree tree("McDst", "validation fixture");
    tree.Branch("Event", &eventPointer);
    tree.Branch("Particle", &particlePointer);
    tree.Branch("SourceId", &ids);
    int children[2] = {-1, -1};
    for (const auto& source : input) {
        events.Delete();
        particles.Delete();
        ids.clear();
        auto* event = new (events[0]) McEvent();
        event->setEventNr(source.number);
        event->setB(source.impact);
        for (const auto& particle : source.particles) {
            const auto& v = particle.values;
            const int i = particles.GetEntriesFast();
            new (particles[i]) McParticle(i, particle.pdg, 0, 0, 0, -1, 0, children, v[0], v[1],
                                          v[2], 10., v[3], v[4], v[5], 20.);
            ids.push_back(particle.id);
        }
        tree.Fill();
    }
    tree.Write();
}

using Row = std::map<std::string, std::string>;

std::vector<std::string> Split(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream input(line);
    std::string field;
    while (std::getline(input, field, ','))
        fields.push_back(field);
    if (!line.empty() && line.back() == ',')
        fields.emplace_back();
    return fields;
}

std::vector<Row> ReadCsv(const fs::path& path) {
    std::ifstream input(path);
    Check(static_cast<bool>(input), "CSV was not created");
    std::string line;
    std::getline(input, line);
    const auto header = Split(line);
    std::vector<Row> result;
    while (std::getline(input, line)) {
        const auto fields = Split(line);
        Check(fields.size() == header.size(), "CSV row has the wrong number of columns");
        Row row;
        for (std::size_t i = 0; i < fields.size(); ++i)
            row[header[i]] = fields[i];
        result.push_back(std::move(row));
    }
    return result;
}

TH1& Histogram(TFile& file, const char* name) {
    auto* histogram = dynamic_cast<TH1*>(file.Get(name));
    Check(histogram != nullptr, "ROOT histogram is missing");
    return *histogram;
}

void CheckHistogram(TFile& file, const char* name, double entries, double mean) {
    auto& histogram = Histogram(file, name);
    Check(histogram.GetEntries() == entries && histogram.Integral() == entries,
          "wrong histogram population or entries outside the axis range");
    Check(std::abs(histogram.GetMean() - mean) < 1e-10, "wrong histogram delta or sign");
}

void TestComparison(const fs::path& directory) {
    const auto current = directory / "current.root", v1 = directory / "v1.root";
    const auto report = directory / "diff.csv";
    WriteMcDst(current, {{7, 6., {{401}, {400}}}, {7, 8., {}}});
    WriteMcDst(v1, {{7, 6., {{400}, {401}}}, {7, 8., {}}});
    Check(CompareMcDstVersions(current.string(), v1.string(), report.string()) ==
              ValidationResult::Equal,
          "particle reordering must not produce differences");
    auto rows = ReadCsv(report);
    Check(rows.size() == 4 && rows[1].at("current_index") != rows[1].at("v1_index"),
          "fixture did not exercise reordered particles");
    Check(rows.back().at("event_occurrence") == "1" &&
              rows.back().at("current_multiplicity") == "0",
          "repeated event numbers or empty events were lost");
    for (const auto& row : rows) {
        Check(row.at("status") == "equal", "identical data marked different");
        Check(row.at("diff_multiplicity") == "0" && row.at("diff_impact_parameter") == "0",
              "identical event metadata differs");
    }
    {
        TFile histograms((directory / "diff.root").c_str(), "READ");
        Check(!histograms.IsZombie(), "ROOT histogram file was not created");
        for (const char* name : {"diff_px", "diff_py", "diff_pz", "diff_x", "diff_y", "diff_z",
                                 "diff_multiplicity", "diff_impact_parameter"})
            CheckHistogram(histograms, name, 2., 0.);
        auto& events = Histogram(histograms, "event_status");
        auto& particles = Histogram(histograms, "particle_status");
        Check(events.GetEntries() == 2 && events.GetBinContent(1) == 2 &&
                  particles.GetEntries() == 2 && particles.GetBinContent(1) == 2,
              "reordered or empty events changed histogram counts");
        Check(std::string(events.GetXaxis()->GetBinLabel(3)) == "current_only" &&
                  dynamic_cast<TNamed*>(histograms.Get("comparison")) &&
                  dynamic_cast<TNamed*>(histograms.Get("current_input")),
              "histogram labels or comparison metadata missing");
    }

    WriteMcDst(current, {{7, 6., {{401}, {400}, {402}, {404}}}, {9, 8., {{900}}}});
    WriteMcDst(v1, {{7, 2., {{400, -211, {0., 4., 1., 8., 2., 9.}}, {403}}}, {10, 9., {{901}}}});
    Check(CompareMcDstVersions(current.string(), v1.string(), report.string()) ==
              ValidationResult::Different,
          "changed values must produce differences");
    rows = ReadCsv(report);
    bool foundParticle = false, foundCurrentEvent = false, foundV1Event = false;
    for (const auto& row : rows) {
        if (row.at("particle_id") == "400") {
            const std::array<std::string, 6> names = {"px", "py", "pz", "x", "y", "z"};
            const std::array<double, 6> expected = {1., -2., 2., -4., 3., -3.};
            for (std::size_t i = 0; i < names.size(); ++i)
                Check(std::stod(row.at("diff_" + names[i])) == expected[i],
                      "wrong particle delta or sign");
            Check(row.at("diff_multiplicity") == "2" && row.at("diff_impact_parameter") == "4" &&
                      row.at("status") == "different" && row.at("v1_pdg") == "-211",
                  "wrong event deltas or PDG comparison");
            foundParticle = true;
        }
        if (row.at("particle_id") == "403")
            Check(row.at("status") == "v1_only" && row.at("diff_px").empty() &&
                      row.at("current_px").empty(),
                  "missing particle was replaced by numeric zero");
        if (row.at("row_type") == "event" && row.at("event_nr") == "9") {
            Check(row.at("status") == "current_only" && row.at("diff_impact_parameter").empty(),
                  "missing v1 event not recorded");
            foundCurrentEvent = true;
        }
        if (row.at("row_type") == "event" && row.at("event_nr") == "10") {
            Check(row.at("status") == "v1_only", "missing current event not recorded");
            foundV1Event = true;
        }
    }
    Check(foundParticle && foundCurrentEvent && foundV1Event, "comparison omitted records");
    {
        TFile histograms((directory / "diff.root").c_str(), "READ");
        const std::array<const char*, 8> names = {
            "diff_px", "diff_py", "diff_pz",           "diff_x",
            "diff_y",  "diff_z",  "diff_multiplicity", "diff_impact_parameter"};
        const std::array<double, 8> expected = {1., -2., 2., -4., 3., -3., 2., 4.};
        for (std::size_t i = 0; i < names.size(); ++i)
            CheckHistogram(histograms, names[i], 1., expected[i]);
        auto& events = Histogram(histograms, "event_status");
        auto& particles = Histogram(histograms, "particle_status");
        Check(events.GetEntries() == 3 && events.GetBinContent(2) == 1 &&
                  events.GetBinContent(3) == 1 && events.GetBinContent(4) == 1 &&
                  particles.GetEntries() == 7 && particles.GetBinContent(2) == 1 &&
                  particles.GetBinContent(3) == 4 && particles.GetBinContent(4) == 2,
              "missing records were lost or inserted as zero differences");
    }
    Check(CompareMcDstVersions(current.string(), v1.string(), current.string()) ==
              ValidationResult::Error,
          "comparison must not overwrite its own input");
    Check(CompareMcDstVersions(current.string(), v1.string(),
                               (directory / "current.csv").string()) == ValidationResult::Error,
          "histogram output must not overwrite a McDst input");
}

void TestConverters(const fs::path& directory) {
    const auto source = directory / "input.oscar";
    const auto output = directory / "validation output";
    const std::string event =
        "# event 17 ensemble 3 in 2\n"
        "0 1 2 3 .938 4 .1 .2 .3 2212 393 1\n"
        "0 4 5 6 .139 4 .4 .5 .6 211 394 1\n"
        "# interaction in 1 out 1 rho 0 weight 1 partial 1 type 1\n"
        "1 7 8 9 .139 4 .7 .8 .9 211 394 1\n"
        "# event 17 ensemble 3 out 2\n"
        "200 40 50 60 .139 4 1 2 3 211 394 1\n"
        "200 10 20 30 .938 4 4 5 6 2212 393 1\n"
        "# event 17 ensemble 3 end 0 impact 7.5 scattering_projectile_target yes\n";
    {
        std::ofstream file(source);
        file << event;
    }
    Check(ValidateOscarVersions(source.string(), output.string()) == ValidationResult::Equal,
          "legacy and current conversion disagree on ensemble fixture");
    const auto rows = ReadCsv(output / "diff.csv");
    Check(rows.size() == 3 && rows[1].at("particle_id") == "393" &&
              rows[1].at("current_x") == "-999" &&
              std::abs(std::stod(rows[2].at("current_px")) - .7) < 1e-6,
          "history or spectator values changed");
    Check(OscarConverterV1().Convert(source.string(), (directory / "legacy.root").string()),
          "standalone v1 conversion failed");
    Check(OscarConverter().Convert(source.string(), (directory / "normal.root").string()),
          "standalone current conversion failed");
    for (const char* name : {"legacy.root", "normal.root"}) {
        TFile file((directory / name).c_str(), "READ");
        auto* tree = dynamic_cast<TTree*>(file.Get("McDst"));
        Check(tree && tree->GetEntries() == 1 && !tree->GetBranch("SourceId"),
              "validation instrumentation leaked into normal output");
    }

    // A skipped event must not leak IDs or clone-array entries into the next event.
    std::ostringstream skipped;
    skipped << "# event 16 ensemble 3 in 394\n";
    for (int id = 0; id < 394; ++id)
        skipped << "0 1 2 3 .938 4 .1 .2 .3 2212 " << id << " 1\n";
    skipped << "# event 16 ensemble 3 out 394\n";
    for (int id = 0; id < 394; ++id)
        skipped << "200 10 20 30 .938 4 1 2 3 2212 " << id << " 1\n";
    skipped << "# event 16 ensemble 3 end 0 impact 9 scattering_projectile_target no\n";
    {
        std::ofstream file(source);
        file << skipped.str() << event
             << "# event 18 ensemble 3 in 0\n# event 18 ensemble 3 out 0\n"
                "# event 18 ensemble 3 end 0 impact 4 scattering_projectile_target no\n";
    }
    Check(ValidateOscarVersions(source.string(), output.string()) == ValidationResult::Equal,
          "skipped or empty event handling changed");
    const auto multiple = ReadCsv(output / "diff.csv");
    Check(multiple.size() == 4 && multiple.front().at("event_nr") == "17" &&
              multiple.back().at("current_multiplicity") == "0",
          "394-particle selection or empty-event reporting changed");

    // Keep v1's original parser: headers without ensemble are not repaired.
    std::string plain = event;
    std::size_t pos;
    while ((pos = plain.find(" ensemble 3")) != std::string::npos)
        plain.erase(pos, 11);
    {
        std::ofstream file(source);
        file << plain;
    }
    Check(ValidateOscarVersions(source.string()) == ValidationResult::Different,
          "v1's original header handling was changed");
    Check(fs::exists(source.string() + ".validation/diff.csv"), "default report path is wrong");
    {
        TFile histograms((source.string() + ".validation/diff.root").c_str(), "READ");
        CheckHistogram(histograms, "diff_px", 0., 0.);
        CheckHistogram(histograms, "diff_multiplicity", 0., 0.);
        CheckHistogram(histograms, "diff_impact_parameter", 0., 0.);
        Check(Histogram(histograms, "event_status").GetBinContent(3) == 1 &&
                  Histogram(histograms, "particle_status").GetBinContent(3) == 2,
              "unmatched data was not recorded in ROOT status histograms");
    }
    Check(ValidateOscarVersions((directory / "missing.oscar").string()) == ValidationResult::Error,
          "missing input was accepted");
    {
        std::ofstream file(source);
        file << "# event 1 in 1\ninvalid particle\n";
    }
    Check(ValidateOscarVersions(source.string(), (directory / "invalid").string()) ==
                  ValidationResult::Error &&
              !fs::exists(directory / "invalid/diff.csv") &&
              !fs::exists(directory / "invalid/diff.root"),
          "failed conversion produced a comparison report");
}

} // namespace

int main() {
    try {
        TempDirectory temp;
        TestComparison(temp.path);
        TestConverters(temp.path);
        std::cout
            << "Legacy conversion, ID matching, signed deltas and missing-record checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "TEST FAILURE: " << error.what() << '\n';
        return 1;
    }
}
