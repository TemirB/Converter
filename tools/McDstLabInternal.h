#pragma once

#include "McDstLab.h"

#include <TFile.h>
#include <TH1D.h>
#include <TObject.h>
#include <TParameter.h>
#include <TTree.h>

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string_view>

namespace mcdst_lab::detail {

constexpr double kFxtEtaMin = -2.0;
constexpr double kFxtEtaMax = 0.0;
constexpr double kFxtPtMin = 0.15;
constexpr int kMultiplicityBins = 65536;
constexpr double kMultiplicityMin = -0.5;
constexpr double kMultiplicityMax = 65535.5;

struct CentralityClass {
    double upper_percentile;
    const char* label;
    const char* histogram_name;
};

inline constexpr std::array<CentralityClass, 4> kCentralityClasses = {{
    {10.0, "0-10", "hNFxt_0_10"},
    {30.0, "10-30", "hNFxt_10_30"},
    {50.0, "30-50", "hNFxt_30_50"},
    {80.0, "50-80", "hNFxt_50_80"},
}};

using FxtCounts = std::map<int, Long64_t>;

struct StagedResult {
    Result summary;
    FxtCounts fxt_counts;
};

struct CentralityBin {
    Long64_t events = 0;
    double percentile_low = 0.0;
    double percentile_high = 0.0;
    double percentile = 0.0;
    int class_index = -1; // -1 denotes the peripheral 80-100% remainder.
};

struct CentralityCalibration {
    std::map<int, CentralityBin> bins;
    Long64_t events = 0;
};

struct CalibrationScope {
    std::size_t files = 0;
    std::size_t failed_files = 0;
    const char* description = "";
};

// The caller owns the staging file until every calculation and write succeeds.
class TemporaryFile {
  public:
    explicit TemporaryFile(const std::string& destination);
    ~TemporaryFile();

    TemporaryFile(const TemporaryFile&) = delete;
    TemporaryFile& operator=(const TemporaryFile&) = delete;

    const std::string& Path() const;
    void Commit(const std::string& destination);

  private:
    std::string path_;
};

class ResetTreeAddresses {
  public:
    explicit ResetTreeAddresses(TTree& tree) : tree_(tree) {}
    ~ResetTreeAddresses() {
        tree_.ResetBranchAddresses();
    }

    ResetTreeAddresses(const ResetTreeAddresses&) = delete;
    ResetTreeAddresses& operator=(const ResetTreeAddresses&) = delete;

  private:
    TTree& tree_;
};

inline void Require(bool condition, std::string_view message) {
    if (!condition) {
        throw std::runtime_error(std::string(message));
    }
}

inline void CheckFloat(double value) {
    Require(std::isfinite(value) && std::abs(value) <= std::numeric_limits<float>::max(),
            "non-finite value or value outside mcDst Float_t range");
}

inline void Write(TObject& object, const char* name = nullptr) {
    Require(object.Write(name, TObject::kOverwrite) > 0,
            "cannot write ROOT object: " + std::string(object.GetName()));
}

template <class T> void WriteParameter(const char* name, T value) {
    TParameter<T> parameter(name, value);
    Write(parameter);
}

inline void CloseOutput(TFile& file, const char* message) {
    file.Close();
    Require(!file.TestBit(TFile::kWriteError), message);
}

inline std::unique_ptr<TH1D> MakeMultiplicityHistogram(const char* name, const char* title) {
    auto histogram =
        std::make_unique<TH1D>(name, title, kMultiplicityBins, kMultiplicityMin, kMultiplicityMax);
    histogram->SetDirectory(nullptr);
    return histogram;
}

void InitializeRoot();
void CheckDestination(const std::string& input, const std::string& output);
StagedResult WriteTransformedFile(const std::string& input, const std::string& output,
                                  const Options& options);

CentralityCalibration CalibrateCentrality(const FxtCounts& counts);
void AddCentrality(const std::string& path, const CentralityCalibration& calibration,
                   const CalibrationScope& scope);

} // namespace mcdst_lab::detail
