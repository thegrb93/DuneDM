#pragma once
// class detector
class Particle;
class DUNEDetector {

	double theta;
public:
	double thetaenter (double);
	double thetaexit (double);
	double Lenter (double);
	double Lexit (double);
	double Ldet (double);
	void intersect (int &, int &, Particle&);
};

