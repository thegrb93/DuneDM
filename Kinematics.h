#pragma once
// class Kinematics

class Kinematics {

	double px, py, pz, p0, a, b, c, p;
public:
	double pAbs (double,double,double,double);
	double pt (double,double,double,double);
	double theta (double,double,double,double);
	double phi (double,double,double,double);
	double invariantmass (double,double,double,double);
	double lambda (double,double,double);
};

