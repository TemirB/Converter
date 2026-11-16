#include "McDstLabInternal.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <thread>
#include <unordered_set>
#include <vector>

namespace mcdst_lab {
namespace {

using namespace detail;
namespace fs = std::filesystem;

struct FileJob {
    fs::path input;
    fs::path output;
    std::unique_ptr<TemporaryFile> staged;
    FxtCounts fxt_counts;
};

std::string OutputFilename(const fs::path& input, bool hist_only) {
    std::string stem = input.stem().string();
    constexpr std::string_view suffix = ".mcdst";
    if (stem.size() >= suffix.size()) {
        std::string ending = stem.substr(stem.size() - suffix.size());
        std::transform(ending.begin(), ending.end(), ending.begin(),
                       [](unsigned char character) { return std::tolower(character); });
        if (ending == suffix) {
            stem.resize(stem.size() - suffix.size());
        }
    }
    return stem + (hist_only ? ".qa.root" : ".mcDst.root");
}

std::vector<fs::path> FindInputFiles(const fs::path& directory, std::size_t max_files) {
    std::vector<fs::path> files;
    for (const auto& file : fs::directory_iterator(directory)) {
        if (file.is_regular_file() && file.path().extension() == ".root") {
            files.push_back(file.path());
        }
    }
    Require(!files.empty(), "no .root files in input directory (scan is not recursive)");
    std::sort(files.begin(), files.end());
    if (max_files != 0 && files.size() > max_files) {
        files.resize(max_files);
    }
    return files;
}

std::vector<FileJob> PlanBatch(const std::string& input_name, const std::string& output_name,
                               bool hist_only, std::size_t max_files) {
    Require(fs::is_directory(input_name), "input must be a directory");
    const auto input = fs::canonical(input_name);
    const auto output = fs::weakly_canonical(output_name);
    Require(!output_name.empty() && input != output, "input and output directories must differ");
    if (fs::exists(output)) {
        Require(fs::is_directory(output), "output must be a directory");
    }

    std::vector<FileJob> jobs;
    std::unordered_set<std::string> output_names;
    for (const auto& file : FindInputFiles(input, max_files)) {
        const auto filename = OutputFilename(file, hist_only);
        Require(output_names.insert(filename).second,
                "multiple input files map to the same output filename: " + filename);
        jobs.push_back({file, output / filename, nullptr, {}});
    }
    fs::create_directories(output);
    return jobs;
}

unsigned WorkerCount(unsigned requested, std::size_t files) {
    const unsigned available =
        requested ? requested : std::max(1u, std::thread::hardware_concurrency());
    return static_cast<unsigned>(std::min<std::size_t>(available, files));
}

Result StageFile(FileJob& job, const Options& options) {
    CheckDestination(job.input.string(), job.output.string());
    auto temporary = std::make_unique<TemporaryFile>(job.output.string());
    auto processed = WriteTransformedFile(job.input.string(), temporary->Path(), options);
    job.fxt_counts = std::move(processed.fxt_counts);
    job.staged = std::move(temporary);
    return processed.summary;
}

void ReportFile(const FileJob& job, const Result& result) {
    std::cout << job.input.filename().string() << " -> " << job.output.filename().string() << ": "
              << result.events << " events, " << result.particles << " particles";
    if (result.unknown_pdg != 0) {
        std::cout << "; unknown PDG records=" << result.unknown_pdg
                  << " (NCh/NFxt may be incomplete; see NUnknownPDG)";
    }
    std::cout << '\n';
}

template <class Action> void RunWorkers(unsigned count, const Action& action) {
    std::vector<std::thread> threads;
    try {
        for (unsigned index = 1; index < count; ++index) {
            threads.emplace_back(action);
        }
    } catch (...) {
        // Finish already started workers if the OS refuses another thread.
        for (auto& thread : threads) {
            thread.join();
        }
        throw;
    }
    action();
    for (auto& thread : threads) {
        thread.join();
    }
}

BatchResult StageFiles(std::vector<FileJob>& jobs, const Options& options, unsigned workers) {
    BatchResult result;
    result.workers = workers;
    std::atomic<std::size_t> next{0};
    std::mutex reporting;

    RunWorkers(workers, [&] {
        while (true) {
            const std::size_t index = next.fetch_add(1);
            if (index >= jobs.size()) {
                return;
            }
            auto& job = jobs[index];
            try {
                const auto summary = StageFile(job, options);
                std::lock_guard<std::mutex> lock(reporting);
                ++result.succeeded;
                ReportFile(job, summary);
            } catch (const std::exception& error) {
                std::lock_guard<std::mutex> lock(reporting);
                ++result.failed;
                std::cerr << job.input.filename().string() << ": " << error.what() << '\n';
            }
        }
    });
    return result;
}

CentralityCalibration CalibrateBatch(const std::vector<FileJob>& jobs) {
    FxtCounts combined;
    for (const auto& job : jobs) {
        if (!job.staged) {
            continue;
        }
        for (const auto& [multiplicity, count] : job.fxt_counts) {
            combined[multiplicity] += count;
        }
    }
    return CalibrateCentrality(combined);
}

void AnnotateAndPublish(std::vector<FileJob>& jobs, const CentralityCalibration& calibration,
                        const BatchResult& result) {
    const CalibrationScope scope{
        result.succeeded,
        result.failed,
        "selected directory inputs; successfully staged files",
    };

    // Keep every output staged until all centrality annotations succeed.
    for (const auto& job : jobs) {
        if (job.staged) {
            AddCentrality(job.staged->Path(), calibration, scope);
        }
    }
    for (auto& job : jobs) {
        if (job.staged) {
            job.staged->Commit(job.output.string());
        }
    }
}

void ReportBatch(const CentralityCalibration& calibration, const BatchResult& result) {
    std::cout << "FXT centrality: " << calibration.events << " events across " << result.succeeded
              << " files; classes 0-10, 10-30, 30-50, 50-80%; "
                 "equal multiplicities stay together\n";
    if (result.failed != 0) {
        std::cerr << "FXT calibration excludes " << result.failed << " failed input files\n";
    }
    std::cout << "Completed: " << result.succeeded << "; failed: " << result.failed << '\n';
}

} // namespace

BatchResult ProcessDirectory(const std::string& input, const std::string& output,
                             const Options& options, unsigned threads, std::size_t max_files) {
    const Frame frame(options);
    auto jobs = PlanBatch(input, output, options.hist_only, max_files);
    InitializeRoot();
    const unsigned workers = WorkerCount(threads, jobs.size());
    std::cout << "Files: " << jobs.size() << "; workers: " << workers << "; beta=" << frame.beta
              << "; gamma=" << frame.gamma << '\n';

    const auto result = StageFiles(jobs, options, workers);
    const auto calibration = CalibrateBatch(jobs);
    AnnotateAndPublish(jobs, calibration, result);
    ReportBatch(calibration, result);
    return result;
}

} // namespace mcdst_lab
