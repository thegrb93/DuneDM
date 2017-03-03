#pragma once
// class DMscattering

class Particle;
class DMscattering {

	double Ee, Thetael;
	double EDM, MDM;
	double MDP;
	double kappa, alphaD;
public:
	double EeTheta (double,double,double);
	double EeTMax (double,double);
	double EeTMin (double,double);
	double F1 (double,double,double,double);
	double dsigmadEe (double,double,double,double,double,double);
	double F2 (double,double,double,double);
	double sigma (double,double,double,double,double);
	void probscatter(int &,int &,double &,double,double,double,double,Particle&);
	void scatterevent(int &,int &,double,double,double,double,Particle&,Particle &);
};
