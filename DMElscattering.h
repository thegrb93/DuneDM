#pragma once
// class DMscattering

class Particle;
class DMscattering {
    std::vector<double> neutrino_sigma, neutrino_energy;
public:
    DMscattering();
	double EeTheta (double,double,double);
	double EeTMax (double,double);
	double EeTMin (double,double);
	double F1 (double,double,double,double);
	double dsigmadEe (double,double,double,double,double,double);
	double F2 (double,double,double,double);
	double sigma (double,double,double,double,double);
	bool probscatter(double,double,double,double,Particle&,double LXdet);
	void scatterevent(double,double,double,double,Particle&,Particle &);
	bool probscatterNeutrino(Particle &DM, double LXdet);
	void scattereventNeutrino(Particle &DM, Particle &electron);

    double nudSigmadEe(double nE, double theta);

    double nuSigma(double nE);
};

