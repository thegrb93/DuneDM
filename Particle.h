#pragma once
// class for lorentz transformation

class Particle{
	
public:
    Particle();
    Particle(double);
    double   m, px, py, pz, E;
    void     FourMomentum(double, double, double, double);
    void     Lorentz(Particle);
    double   norm, normpx, normpy, normpz, pt, theta, phi;
    void calcOptionalKinematics();
};

