#include "McDstLabInternal.h"

#include <TBranch.h>
#include <TBufferJSON.h>
#include <TFileMerger.h>
#include <TKey.h>
#include <TNamed.h>

#include <algorithm>
#include <filesystem>
#include <set>
#include <tuple>
#include <vector>

namespace mcdst_lab {
namespace {

using namespace detail;
namespace fs = std::filesystem;

// These objects describe the shared calibration, not a contribution from each input file.
constexpr std::array<const char*, 15> kSharedMetadata = {
    "LabFrame",
    "LabSqrtsNN",
    "LabNucleonMass",
    "LabBeta",
    "LabGamma",
    "LabHistOnly",
    "FxtDefinition",
    "FxtEtaMin",
    "FxtEtaMax",
    "FxtPtMin",
    "FxtCentralityDefinition",
    "FxtCentralityScope",
    "FxtCentralityEvents",
    "FxtCentralityFiles",
    "FxtCentralityFailedFiles",
};

constexpr std::array<const char*, 10> kHistograms = {
    "hEtaPtBefore",   "hEtaPtAfter", "hNAll",       "hNCh",        "hNFxt",
    "hFxtCentrality", "hNFxt_0_10",  "hNFxt_10_30", "hNFxt_30_50", "hNFxt_50_80",
};

constexpr std::array<const char*, 10> kCounters = {
    "LabEvents",
    "LabParticles",
    "LabUnknownPDG",
    "LabMaskedPositions",
    "hEtaPtBefore_undefined_eta",
    "hEtaPtBefore_outside_plot",
    "hEtaPtBefore_axial",
    "hEtaPtAfter_undefined_eta",
    "hEtaPtAfter_outside_plot",
    "hEtaPtAfter_axial",
};

template <class T> T& Get(TDirectory& directory, const char* name) {
    auto* object = dynamic_cast<T*>(directory.Get(name));
    Require(object != nullptr, "missing or invalid object: " + std::string(name));
    return *object;
}

std::unique_ptr<TFile> OpenInput(const fs::path& path) {
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    Require(file && !file->IsZombie() && !file->TestBit(TFile::kRecovered),
            "cannot read complete input: " + path.string());
    Get<TNamed>(*file, "LabFrame");
    Require(!file->Get("MergeInputs"),
            "input already merged; use original per-file results: " + path.string());
    return file;
}

std::vector<fs::path> CollectInputs(const std::vector<std::string>& inputs) {
    std::vector<fs::path> files;
    std::set<fs::path> unique;
    auto add = [&](const fs::path& path) {
        const auto resolved = fs::canonical(path);
        Require(fs::is_regular_file(resolved), "input is not a file: " + path.string());
        Require(unique.insert(resolved).second, "duplicate input: " + path.string());
        files.push_back(resolved);
    };
    for (const auto& input : inputs) {
        if (fs::is_directory(input)) {
            for (const auto& entry : fs::directory_iterator(input)) {
                if (entry.is_regular_file() && entry.path().extension() == ".root") {
                    add(entry.path());
                }
            }
        } else {
            add(input);
        }
    }
    Require(!files.empty(), "no input .root files");
    std::sort(files.begin(), files.end());
    return files;
}

using CalibrationRows = std::vector<std::tuple<int, Long64_t, double, double, double, int>>;

CalibrationRows ReadCalibration(TFile& file) {
    auto& tree = Get<TTree>(file, "FxtCentralityMap");
    int multiplicity = 0;
    Long64_t events = 0;
    double low = 0.0;
    double high = 0.0;
    double percentile = 0.0;
    int category = -1;
    ResetTreeAddresses reset(tree);
    auto bind = [&](const char* name, auto* value) {
        auto* branch = tree.GetBranch(name);
        Require(branch && branch->GetEntries() == tree.GetEntries(),
                "incomplete calibration branch: " + std::string(name));
        Require(tree.SetBranchAddress(name, value) >= 0,
                "invalid calibration branch: " + std::string(name));
    };
    bind("NFxt", &multiplicity);
    bind("Events", &events);
    bind("PercentileLow", &low);
    bind("PercentileHigh", &high);
    bind("Percentile", &percentile);
    bind("Class", &category);
    CalibrationRows rows;
    int previous = -1;
    Long64_t total = 0;
    for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry) {
        Require(tree.GetEntry(entry) > 0, "cannot read FXT calibration");
        Require(multiplicity > previous && events > 0,
                "invalid or duplicated FXT calibration rows; use original per-file results");
        rows.emplace_back(multiplicity, events, low, high, percentile, category);
        previous = multiplicity;
        total += events;
    }
    Require(total == Get<TParameter<Long64_t>>(file, "FxtCentralityEvents").GetVal(),
            "FXT calibration event count mismatch");
    return rows;
}

using Schema = std::map<std::string, std::string>;

void AddBranchSchema(TBranch& branch, Long64_t entries, Schema& schema) {
    Require(branch.GetEntries() == entries,
            "incomplete event branch: " + std::string(branch.GetFullName().Data()));
    schema[branch.GetFullName().Data()] =
        std::string(branch.GetClassName()) + ":" + branch.GetTitle();
    TIter next(branch.GetListOfBranches());
    while (auto* child = static_cast<TBranch*>(next())) {
        AddBranchSchema(*child, entries, schema);
    }
}

Schema TreeSchema(TTree& tree) {
    Schema schema;
    TIter next(tree.GetListOfBranches());
    while (auto* branch = static_cast<TBranch*>(next())) {
        AddBranchSchema(*branch, tree.GetEntries(), schema);
    }
    return schema;
}

bool SameAxes(const TH1& first, const TH1& next) {
    if (first.GetDimension() != next.GetDimension()) {
        return false;
    }
    const std::array<const TAxis*, 3> a{first.GetXaxis(), first.GetYaxis(), first.GetZaxis()};
    const std::array<const TAxis*, 3> b{next.GetXaxis(), next.GetYaxis(), next.GetZaxis()};
    for (int axis = 0; axis < first.GetDimension(); ++axis) {
        if (a[axis]->GetNbins() != b[axis]->GetNbins() ||
            a[axis]->GetXmin() != b[axis]->GetXmin() || a[axis]->GetXmax() != b[axis]->GetXmax()) {
            return false;
        }
    }
    return true;
}

void CheckCompatible(TFile& reference, TFile& file, const CalibrationRows& calibration) {
    for (const char* name : kSharedMetadata) {
        const auto& expected = Get<TObject>(reference, name);
        const auto& actual = Get<TObject>(file, name);
        Require(TBufferJSON::ConvertToJSON(&expected) == TBufferJSON::ConvertToJSON(&actual),
                "incompatible metadata: " + std::string(name));
    }
    Require(ReadCalibration(file) == calibration,
            "incompatible FXT calibration; rerun mcdst_lab on the combined input sample");
    for (const char* name : kHistograms) {
        Require(SameAxes(Get<TH1>(reference, name), Get<TH1>(file, name)),
                "incompatible histogram axes: " + std::string(name));
    }
    auto& expected = Get<TH1>(reference, "hNFxtReference");
    auto& actual = Get<TH1>(file, "hNFxtReference");
    Require(SameAxes(expected, actual), "incompatible hNFxtReference axes");
    for (int bin = 0; bin < expected.GetNcells(); ++bin) {
        Require(expected.GetBinContent(bin) == actual.GetBinContent(bin),
                "incompatible hNFxtReference contents");
    }
    const bool hist_only = Get<TParameter<int>>(reference, "LabHistOnly").GetVal() != 0;
    const Long64_t events = Get<TParameter<Long64_t>>(file, "LabEvents").GetVal();
    for (const char* name : {"McDst", "LabMultiplicity"}) {
        if (hist_only && std::string(name) == "McDst") {
            continue;
        }
        auto& tree = Get<TTree>(file, name);
        Require(tree.GetEntries() == events &&
                    TreeSchema(tree) == TreeSchema(Get<TTree>(reference, name)),
                "incompatible event tree: " + std::string(name));
    }
}

Result ValidateInputs(const std::vector<fs::path>& paths, TFile& reference) {
    const auto calibration = ReadCalibration(reference);
    std::set<std::string> sources;
    Result result;
    for (const auto& path : paths) {
        auto file = OpenInput(path);
        CheckCompatible(reference, *file, calibration);
        const std::string source = Get<TNamed>(*file, "LabSource").GetTitle();
        Require(sources.insert(source).second, "duplicate LabSource: " + source);
        result.events += Get<TParameter<Long64_t>>(*file, "LabEvents").GetVal();
        result.particles += Get<TParameter<Long64_t>>(*file, "LabParticles").GetVal();
        result.unknown_pdg += Get<TParameter<Long64_t>>(*file, "LabUnknownPDG").GetVal();
    }
    return result;
}

void MergeLocalObjects(const std::vector<fs::path>& paths, const std::string& output) {
    TFileMerger merger(false, false);
    merger.SetMaxOpenedFiles(16);
    Require(merger.OutputFile(output.c_str(), "RECREATE", 1), "cannot create merge output");
    for (const char* name : {"McDst", "LabMultiplicity"}) {
        merger.AddObjectNames(name);
    }
    for (const char* name : kHistograms) {
        merger.AddObjectNames(name);
    }
    for (const char* name : kCounters) {
        merger.AddObjectNames(name);
    }
    for (const auto& path : paths) {
        Require(merger.AddFile(path.c_str()), "cannot add merge input: " + path.string());
    }
    Require(merger.PartialMerge(TFileMerger::kAll | TFileMerger::kOnlyListed |
                                TFileMerger::kKeepCompression),
            "ROOT failed to merge input files");
    // Some ROOT versions close the output inside PartialMerge(kRegular).
    if (auto* file = merger.GetOutputFile()) {
        Require(!file->TestBit(TFile::kWriteError), "ROOT merge write error");
        merger.CloseOutputFile();
    }
}

void CopyDirectory(TDirectory& source, TDirectory& destination) {
    std::set<std::string> copied;
    TIter next(source.GetListOfKeys());
    while (auto* key = static_cast<TKey*>(next())) {
        if (!copied.insert(key->GetName()).second) {
            continue;
        }
        std::unique_ptr<TObject> object(key->ReadObj());
        Require(object != nullptr, "cannot read source metadata");
        destination.cd();
        if (auto* directory = dynamic_cast<TDirectory*>(object.get())) {
            auto* child = destination.mkdir(key->GetName());
            Require(child != nullptr, "cannot create source metadata directory");
            CopyDirectory(*directory, *child);
        } else {
            Require(!dynamic_cast<TTree*>(object.get()), "unexpected tree in source metadata");
            Write(*object, key->GetName());
        }
    }
}

void CopyFileMetadata(TFile& input, TDirectory& destination) {
    std::set<std::string> copied;
    TIter next(input.GetListOfKeys());
    while (auto* key = static_cast<TKey*>(next())) {
        const std::string type = key->GetClassName();
        if ((type != "McRun" && type != "TNamed" && type.rfind("TParameter<", 0) != 0) ||
            !copied.insert(key->GetName()).second) {
            continue;
        }
        std::unique_ptr<TObject> object(key->ReadObj());
        Require(object != nullptr, "cannot read file metadata");
        destination.cd();
        Write(*object, key->GetName());
    }
    if (auto* source = dynamic_cast<TDirectory*>(input.Get("SourceMetadata"))) {
        auto* target = destination.mkdir("SourceMetadata");
        Require(target != nullptr, "cannot archive SourceMetadata");
        CopyDirectory(*source, *target);
    }
}

void WriteProvenance(const std::vector<fs::path>& paths, TFile& output) {
    auto* archive = output.mkdir("MergedSources");
    Require(archive != nullptr, "cannot create MergedSources");
    output.cd();
    TTree inputs("MergeInputs", "Input files and their ranges in the merged event trees");
    std::string input_path;
    std::string lab_source;
    Long64_t first_entry = 0;
    Long64_t events = 0;
    inputs.Branch("InputFile", &input_path);
    inputs.Branch("LabSource", &lab_source);
    inputs.Branch("FirstEntry", &first_entry);
    inputs.Branch("Events", &events);
    for (std::size_t index = 0; index < paths.size(); ++index) {
        auto file = OpenInput(paths[index]);
        input_path = paths[index].string();
        lab_source = Get<TNamed>(*file, "LabSource").GetTitle();
        events = Get<TParameter<Long64_t>>(*file, "LabEvents").GetVal();
        auto* directory = archive->mkdir(std::to_string(index).c_str());
        Require(directory != nullptr, "cannot create input metadata archive");
        CopyFileMetadata(*file, *directory);
        Require(inputs.Fill() >= 0, "cannot write merge input manifest");
        first_entry += events;
    }
    output.cd();
    Write(inputs);
    inputs.SetDirectory(nullptr);
}

void RestoreSharedObjects(TFile& reference, TFile& output) {
    output.cd();
    for (const char* name : kSharedMetadata) {
        Write(Get<TObject>(reference, name), name);
    }
    Write(Get<TH1>(reference, "hNFxtReference"));
    auto& source = Get<TTree>(reference, "FxtCentralityMap");
    std::unique_ptr<TTree> mapping(source.CloneTree(-1, "fast"));
    Require(mapping != nullptr, "cannot copy shared FXT calibration");
    mapping->SetDirectory(&output);
    Write(*mapping);

    // Keep the first run header at its traditional location; every source is archived separately.
    TIter next(reference.GetListOfKeys());
    std::set<std::string> copied;
    while (auto* key = static_cast<TKey*>(next())) {
        if (std::string(key->GetClassName()) == "McRun" && copied.insert(key->GetName()).second) {
            Write(Get<TObject>(reference, key->GetName()), key->GetName());
        }
    }
    if (auto* original = dynamic_cast<TDirectory*>(reference.Get("SourceMetadata"))) {
        // Partial merges can leave an empty directory even with kOnlyListed.
        auto* directory = output.GetDirectory("SourceMetadata");
        if (directory) {
            directory->Delete("*;*");
        } else {
            directory = output.mkdir("SourceMetadata");
        }
        Require(directory != nullptr, "cannot create SourceMetadata");
        CopyDirectory(*original, *directory);
    }
    output.cd();
    TNamed source_name("LabSource", "merged outputs; see MergeInputs and MergedSources");
    Write(source_name);
}

} // namespace

Result Merge(const std::vector<std::string>& inputs, const std::string& output) {
    InitializeRoot();
    Require(!output.empty() && !fs::exists(output),
            "output already exists or is empty; choose a new output path");
    const auto paths = CollectInputs(inputs);
    auto reference = OpenInput(paths.front());
    const auto result = ValidateInputs(paths, *reference);
    TemporaryFile temporary(output);
    MergeLocalObjects(paths, temporary.Path());
    {
        TFile merged(temporary.Path().c_str(), "UPDATE");
        Require(!merged.IsZombie() && !merged.TestBit(TFile::kRecovered),
                "cannot reopen complete merge output");
        RestoreSharedObjects(*reference, merged);
        WriteProvenance(paths, merged);
        Require(Get<TTree>(merged, "LabMultiplicity").GetEntries() == result.events,
                "merged event count mismatch");
        Require(Get<TParameter<Long64_t>>(merged, "LabEvents").GetVal() == result.events,
                "merged LabEvents mismatch");
        CloseOutput(merged, "cannot finalize merged output");
    }
    temporary.Commit(output);
    return result;
}

} // namespace mcdst_lab
