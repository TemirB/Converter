#pragma once

class Particle {
public:
    double t, x, y, z, mass, p0, px, py, pz; // 9*8 = 72 byte
    int pdg, ID, charge; // 3*4 = 12 byte
    // 75 bytes -> 

    Particle();
    Particle(double t, double x, double y, double z, double mass, double p0, double px, double py, double pz, int pdg, int ID, int charge);
    void SetSpectator();
};
