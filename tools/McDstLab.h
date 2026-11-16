#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace mcdst_lab {

struct Options {
    double sqrtsnn = 3.0;
    double nucleon_mass = 0.938;
    double eta_max = 6.0;
    double pt_max = 3.0;
    bool hist_only = false;
};

// Symmetric beams in the SMASH CM frame: Projectile +z, Target -z.
struct Frame {
    explicit Frame(const Options& options);
    double gamma;
    double beta;

    double Energy(double energy, double pz) const;

    // Includes the rotation about y: output +z is opposite to Projectile.
    double Longitudinal(double energy, double pz) const;
};

struct Result {
    long long events = 0;
    long long particles = 0;
    long long unknown_pdg = 0;
};

struct BatchResult {
    unsigned workers = 0;
    std::size_t succeeded = 0;
    std::size_t failed = 0;
};

// Includes NFxt: all known charged particles with -2 < eta_lab < 0, pT > 0.15 GeV/c.
// Single-file centrality is calibrated on this file, including NFxt=0 events.
// Throws on failure. The destination is replaced only after a complete write.
Result Process(const std::string& input, const std::string& output, const Options& options = {});

// Processes *.root in this directory (not recursively), writing *.mcDst.root.
// The QA-only hist_only mode writes *.qa.root without an McDst tree.
// threads == 0 selects hardware_concurrency; failures do not stop other files.
// Centrality uses one shared right-tail integral across successfully processed selected files.
// Equal NFxt values share a percentile midpoint and stay in the same class.
// max_files == 0 selects all input files.
BatchResult ProcessDirectory(const std::string& input, const std::string& output,
                             const Options& options = {}, unsigned threads = 0,
                             std::size_t max_files = 0);

// Merge per-file lab/QA results with identical FXT calibration. Input directories are expanded
// non-recursively. Shared calibration is copied once; existing outputs are never overwritten.
Result Merge(const std::vector<std::string>& inputs, const std::string& output);

} // namespace mcdst_lab
