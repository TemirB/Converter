#pragma once

#include<cstdint>

class alignas(8) Particle {
public:
    double t, x, y, z, mass, p0, px, py, pz; // 8*9 = 72 bytes
    int pdg, ID, charge; // 4*3 = 12 bytes
    // 72 + 12 = 84 -> 88 (padding to allign)

    Particle();
    Particle(double t, double x, double y, double z, double mass, double p0, double px, double py, double pz, int pdg, int ID, int charge);
    void SetSpectator();
};
