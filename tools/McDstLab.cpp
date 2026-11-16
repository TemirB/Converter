#include "McDstLabInternal.h"

#include <cerrno>
#include <cstring>
#include <filesystem>
#include <unistd.h>

namespace mcdst_lab::detail {

TemporaryFile::TemporaryFile(const std::string& destination) : path_(destination + ".tmp.XXXXXX") {
    const int descriptor = mkstemp(path_.data());
    Require(descriptor >= 0,
            "cannot create temporary output: " + std::string(std::strerror(errno)));
    if (close(descriptor) != 0) {
        std::filesystem::remove(path_);
        throw std::runtime_error("cannot close temporary output");
    }
}

TemporaryFile::~TemporaryFile() {
    if (!path_.empty()) {
        std::error_code ignored;
        std::filesystem::remove(path_, ignored);
    }
}

const std::string& TemporaryFile::Path() const {
    return path_;
}

void TemporaryFile::Commit(const std::string& destination) {
    std::filesystem::rename(path_, destination);
    path_.clear();
}

void CheckDestination(const std::string& input, const std::string& output) {
    Require(!output.empty() && !std::filesystem::is_directory(output),
            "output must be a file path");
    Require(std::filesystem::weakly_canonical(input) != std::filesystem::weakly_canonical(output) &&
                (!std::filesystem::exists(output) || !std::filesystem::equivalent(input, output)),
            "input and output must be different files");
}

} // namespace mcdst_lab::detail

namespace mcdst_lab {

Frame::Frame(const Options& options) {
    detail::Require(std::isfinite(options.sqrtsnn) && std::isfinite(options.nucleon_mass) &&
                        options.nucleon_mass > 0.0 && options.sqrtsnn > 2.0 * options.nucleon_mass,
                    "sqrtsnn must exceed twice the positive nucleon mass (GeV)");
    gamma = options.sqrtsnn / (2.0 * options.nucleon_mass);
    beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    detail::Require(std::isfinite(gamma) && beta < 1.0, "beam energy is too large for this boost");
}

// The boost stops Target; the subsequent rotation makes +z opposite to Projectile.
double Frame::Energy(double energy, double pz) const {
    return gamma * (energy + beta * pz);
}

double Frame::Longitudinal(double energy, double pz) const {
    return -gamma * (pz + beta * energy);
}

Result Process(const std::string& input, const std::string& output, const Options& options) {
    detail::CheckDestination(input, output);
    detail::TemporaryFile temporary(output);

    const auto staged = detail::WriteTransformedFile(input, temporary.Path(), options);
    const auto calibration = detail::CalibrateCentrality(staged.fxt_counts);
    const detail::CalibrationScope scope{1, 0, "single input file"};
    detail::AddCentrality(temporary.Path(), calibration, scope);

    temporary.Commit(output);
    return staged.summary;
}

} // namespace mcdst_lab
