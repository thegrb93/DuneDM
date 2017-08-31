#pragma once
// class detector
class Particle;
class DUNEDetector {
    double posx, posy, posz;
    double minsx, minsy, minsz;
    double maxsx, maxsy, maxsz;
public:
    DUNEDetector();
	double thetaenter (double);
	double thetaexit (double);
	double Lenter (double);
	double Lexit (double);
	double Ldet (Particle &DM);
	bool intersect (double dx, double dy, double dz, double& tmin, double& tmax);
};

