#pragma once

#include "Particle.h"

#include <cstddef>
#include <istream>
#include <string>
#include <vector>

struct OscarEvent {
    int number = -1;
    int ensemble = 0;
    double impact_parameter = 0.;
    double final_time = 0.;
    bool projectile_target_scattering = false;
    int initial_nucleons = 0;
    int participants = 0;
    int spectators = 0;
    int binary_collisions = 0;
    int nn_elastic_collisions = 0;
    int nn_inelastic_collisions = 0;
    std::vector<Particle> particles;
};

// Streaming reader of complete, serial SMASH collision histories.
// False means EOF or an error; Error() distinguishes these outcomes.
class OscarReader {
  public:
    explicit OscarReader(std::istream& input) : input_(input) {}
    bool ReadNext(OscarEvent& event);
    const std::string& Error() const {
        return error_;
    }

  private:
    bool ReadLine(std::string& line);
    bool ReadParticles(int count, std::vector<Particle>& particles);
    bool Fail(const std::string& message);

    std::istream& input_;
    std::size_t line_number_ = 0;
    std::size_t particle_columns_ = 0;
    bool format_seen_ = false;
    std::string error_;
};
