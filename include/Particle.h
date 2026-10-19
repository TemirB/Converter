#pragma once

// OSCAR2013: positions in fm, masses and momenta in GeV (c = 1).
struct Particle {
    double t = 0., x = 0., y = 0., z = 0.;
    double mass = 0., p0 = 0., px = 0., py = 0., pz = 0.;
    int pdg = 0, id = -1, charge = 0;
    bool spectator = false;
};

inline bool IsNucleon(int pdg) {
    return pdg == 2212 || pdg == 2112;
}
