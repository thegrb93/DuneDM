#include <iostream>
#include <cmath>
#include <vector>

#include "DMElscattering.h"
#include "Particle.h"
#include "Random.h"
#include "Kinematics.h"
#include "DUNEdet.h"

// Electron energy as a function of 
// electron scattering angle cross and dark matter energy
double DMscattering::EeTheta (double EDM, double Thetael, double MDM) {
	double rEeTheta;
	double Me = 0.000511; 
	double EeThetaN, EeThetaD;
	EeThetaN = Me*((EDM+Me)*(EDM+Me)+(EDM*EDM-MDM*MDM)*cos(Thetael)*cos(Thetael));
	EeThetaD = (EDM+Me)*(EDM+Me)-(EDM*EDM-MDM*MDM)*cos(Thetael)*cos(Thetael);
	rEeTheta = EeThetaN/EeThetaD;
	return(rEeTheta);
}
// Maximum electron energy as a function of 
// dark matter energy and mass
double DMscattering::EeTMax (double EDM, double MDM) {
	double rEeTMax;
	double ThetaelMax = 0.0;
	rEeTMax = EeTheta(EDM,ThetaelMax,MDM);
        //std::cout<<"Eemax"<<rEeTMax<<std::endl;
	return(rEeTMax);
}
// Minimum electron energy as a function of 
// dark matter energy and mass
double DMscattering::EeTMin (double EDM, double MDM) {
	double rEeTMin;
	double Pi = 3.141592653589793;
	double ThetaelMin = Pi/2.0;
	rEeTMin = EeTheta(EDM,ThetaelMin,MDM);
        std::cout<<"Eemin "<<rEeTMin<<std::endl;
	return(rEeTMin);
}
// Function F1
// dsigma/dEe = 4*Pi*kappa*kappa*alpha*alphaD*F1
double DMscattering::F1 (double Ee, double EDM, double MDM, double MDP) {
	double rF1;
	double Me = 0.000511;
	double F1N, F1D;
	F1N = 2.0*Me*EDM*EDM-(2.0*Me*EDM+MDM*MDM)*(Ee-Me);
	F1D = (EDM*EDM-MDM*MDM)*(MDP*MDP+2*Me*Ee-2*Me*Me)*(MDP*MDP+2*Me*Ee-2*Me*Me);
	rF1 = F1N/F1D;
	return(rF1);
}
//  differential DM - electron scattering cross section dsigma/dEe
double DMscattering::dsigmadEe (double Ee, double EDM, double MDM, double MDP, double kappa, double alphaD) {
	double rdsig;
	double Pi = 3.141592653589793;
	double alphaEM = 1/137.035999074;
	double coef;
	coef = 4*Pi*kappa*kappa*alphaEM*alphaD;
	rdsig = coef*F1(Ee,EDM,MDM,MDP);
	return(rdsig);
}
// Function F2(Ee) 
// Total DM - electron scattering cross section equals
// sigma =  4*Pi*kappa*kappa*alpha*alphaD*( F2(EeMax)- F2(EeMin) ) 
double DMscattering::F2 (double Ee, double EDM, double MDM, double MDP) {
	double rF2;
	double Me = 0.000511;
	double F2N1, F2N2, F2D;
	F2N1 = (4*EDM*EDM*Me*Me+2*EDM*Me*MDP*MDP + MDP*MDP*MDM*MDM)/(2*Ee*Me-2*Me*Me+MDP*MDP);
	F2N2 = (2*EDM*Me+MDM*MDM)*log(2*Ee*Me-2*Me*Me+MDP*MDP);
	F2D  = 4*Me*Me*(EDM*EDM-MDM*MDM);
	rF2  = -(F2N1+F2N2)/F2D;
	return(rF2);
}
//
double DMscattering::sigma (double EDM, double MDM, double MDP, double kappa, double alphaD) {
	double rsig;
	double Pi = 3.141592653589793;
	double alphaEM = 1/137.035999074;
	double coef;
	coef = 4*Pi*kappa*kappa*alphaEM*alphaD;
	double EeMaxA, EeMinA;
	EeMaxA = EeTMax(EDM,MDM);
	EeMinA = EeTMin(EDM,MDM);
	rsig = coef*(F2(EeMaxA,EDM,MDM,MDP)-F2(EeMinA,EDM,MDM,MDP));
	return(rsig);
}
//
void DMscattering::probscatter (int &dswitch, int &Nscat, double &pMax, double MDP, double MDM, double kap, double alD, Particle& DM) {
	double pscat, Rscat;
	double LXdet, XS;
	double prob;
	double pMax0 = 1.0e-10;
	double ne = 5.1e+23;
	//int Nscatter;	
	double convmcm, convGeV2cm2;
	convGeV2cm2 = 3.89e-28;
	convmcm = 100.0;
	pscat = Random::Flat(0,1);
	Kinematics kin;
	DUNEDetector det;
	if (dswitch == 1)
	{

		LXdet = det.Ldet(DM);
		LXdet = LXdet*convmcm;
		XS = sigma(DM.E,MDM,MDP,kap,alD);
		XS = XS*convGeV2cm2;
		//std::cout<<DM.E<<"\t"<<(XS*pow(10,39))<<std::endl;

		prob = XS*ne*LXdet;

		if (prob > pMax0)
		{
			pMax = prob;

		}
		Rscat = prob/pMax0;

		if (Rscat > pscat)
		{
			Nscat = Nscat+1;
			dswitch = 2;

		}
	}


}

void DMscattering::probscatterNeutrino (int &dswitch, int &Nscat, double &pMax, std::vector<double>& energies, std::vector<double>& xsections, Particle& DM) {
	double pscat, Rscat;
	double LXdet, XS;
	double prob;
	double pMax0 = 1.0e-10;
	double ne = 5.1e+23;
	//int Nscatter;
	double convmcm, convGeV2cm2;
	convGeV2cm2 = 1e-46;
	convmcm = 100.0;
	pscat = Random::Flat(0,1);
	Kinematics kin;
	DUNEDetector det;
	if (dswitch == 1)
	{

		LXdet = det.Ldet(DM);
		LXdet = LXdet*convmcm;
        XS = -1;
        if(DM.E >= energies[0])
		    for(int i = 1; i<energies.size(); ++i)
                if(DM.E <= energies[i])
                {
                    XS = (DM.E - energies[i-1])/(energies[i] - energies[i-1])*(xsections[i]-xsections[i-1])+xsections[i-1];
                }
        if(XS < 0) return;
		XS = XS*convGeV2cm2;
		//std::cout<<DM.E<<"\t"<<(XS*pow(10,39))<<std::endl;

		prob = XS*ne*LXdet;

		if (prob > pMax0)
		{
			pMax = prob;

		}
		Rscat = prob/pMax0;

		if (Rscat > pscat)
		{
			Nscat = Nscat+1;
			dswitch = 2;

		}
	}


}
//
void DMscattering::scatterevent (int &dswitch, int &Nelec, double MDP, double MDM, double kap, double alD, Particle& DM, Particle &electron) {
    double Pi = 3.141592653589793;
    double Me = 0.000511;
    double EeMin, EeMax;
    double xe, ye, Thetae, Phie, Ee;
    double pe, pex, pey, pez;
    double dsig, sig, psig;
    double dsigMax, psigMax;
    double probe, Re;
    int eswitch;
    //std::cout <<"value is" << pex <<"\t"<< pey <<"\t"<< pez <<"\t"<< Ee <<std::endl;
    if (dswitch == 2)
    {


        eswitch = 0;
        EeMax = EeTMax(DM.E,MDM);
        EeMin = EeTMin(DM.E,MDM);
        sig = sigma(DM.E,MDM,MDP,kap,alD);
        dsigMax = dsigmadEe(EeMin,DM.E,MDM,MDP,kap,alD);
        psigMax =(EeMax-EeMin)*dsigMax/sig;

        while (eswitch == 0)
        {
            probe = Random::Flat(0,1);
            xe = Random::Flat(0,1);
            Thetae = xe*Pi;

            Ee = EeMin + xe*(EeMax-EeMin);
            dsig = dsigmadEe(Ee,DM.E,MDM,MDP,kap,alD);
            psig = (EeMax-EeMin)*dsig/sig;
            Re = psig/psigMax;
            if (Re > probe)
            {
                eswitch = 1;
                ye = Random::Flat(0,1);
                Phie = ye*2.0*Pi;
                pe = sqrt(Ee*Ee-Me*Me);
                pex = pe*sin(Thetae)*cos(Phie);
                pey = pe*sin(Thetae)*sin(Phie);
                pez = pe*cos(Thetae);
                electron.FourMomentum(pex,pey,pez,Ee);
                Nelec = Nelec+1;

            }
        }
    }
}
//
/*void DMscattering::scattereventNeutrino (int &dswitch, int &Nelec, double MDM, Particle& DM, Particle &electron) {
    double Pi = 3.141592653589793;
    double Me = 0.000511;
    double EeMin, EeMax;
    double xe, ye, Thetae, Phie, Ee;
    double pe, pex, pey, pez;
    double dsig, sig, psig;
    double dsigMax, psigMax;
    double probe, Re;
    int eswitch;
    //std::cout <<"value is" << pex <<"\t"<< pey <<"\t"<< pez <<"\t"<< Ee <<std::endl;
    if (dswitch == 2)
    {


        eswitch = 0;
        EeMax = EeTMax(DM.E,0);
        EeMin = EeTMin(DM.E,0);
        sig = sigma(DM.E,MDM,MDP,kap,alD);
        dsigMax = dsigmadEe(EeMin,DM.E,MDM,MDP,kap,alD);
        psigMax =(EeMax-EeMin)*dsigMax/sig;

        while (eswitch == 0)
        {
            probe = Random::Flat(0,1);
            xe = Random::Flat(0,1);
            Thetae = xe*Pi;

            Ee = EeMin + xe*(EeMax-EeMin);
            dsig = dsigmadEe(Ee,DM.E,MDM,MDP,kap,alD);
            psig = (EeMax-EeMin)*dsig/sig;
            Re = psig/psigMax;
            if (Re > probe)
            {
                eswitch = 1;
                ye = Random::Flat(0,1);
                Phie = ye*2.0*Pi;
                pe = sqrt(Ee*Ee-Me*Me);
                pex = pe*sin(Thetae)*cos(Phie);
                pey = pe*sin(Thetae)*sin(Phie);
                pez = pe*cos(Thetae);
                electron.FourMomentum(pex,pey,pez,Ee);
                Nelec = Nelec+1;

            }
        }
    }
}*/
