#pragma once
// class DMscattering

class Particle;
class DMscattering {

	double Ee, Thetael;
	double EDM, MDM;
	double MDP;
	double kappa, alphaD;
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
	void probscatter(int &,int &,double &,double,double,double,double,Particle&);
	void scatterevent(int &,int &,double,double,double,double,Particle&,Particle &);
	void probscatterNeutrino(int &dswitch, int &Nscat, double &pMax, Particle &DM);

	void scattereventNeutrino(int &dswitch, int &Nelec, Particle &DM, Particle &electron);

    double nudSigmadEe(double nE, double theta);

    double nuSigma(double nE);
};

