#pragma once

#include <string>

namespace mcdst_lab {

struct Options {
    double sqrtsnn = 3.;
    double nucleon_mass = .938;
    double eta_max = 6.;
    double pt_max = 3.;
    bool hist_only = false;
};

// Symmetric beams in the SMASH CM frame: Projectile +z, Target -z.
struct Frame {
    explicit Frame(const Options& options);
    double gamma, beta;
    double Energy(double energy, double pz) const;
    // Includes the rotation about y: output +z is opposite to Projectile.
    double Longitudinal(double energy, double pz) const;
};

struct Result {
    long long events = 0, particles = 0, unknown_pdg = 0;
};

struct BatchResult {
    unsigned workers = 0;
    std::size_t succeeded = 0, failed = 0;
};

// Throws on failure. The destination is replaced only after a complete write.
Result Process(const std::string& input, const std::string& output, const Options& options = {});

// Processes *.root in this directory (not recursively), preserving filenames.
// threads == 0 selects hardware_concurrency; failures do not stop other files.
BatchResult ProcessDirectory(const std::string& input, const std::string& output,
                             const Options& options = {}, unsigned threads = 0);

} // namespace mcdst_lab
