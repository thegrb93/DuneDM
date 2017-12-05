#pragma once
#include <string>

class TFile;
class TH1D;
class TH2D;

class DMHistograms
{
public:
    DMHistograms(std::string);
    ~DMHistograms();

    void AddProductionDM(TH1D* E, TH1D* Px, TH1D* Py, TH1D* Pz, TH1D* Th, TH1D* Phi);
    void AddProductionNu(double E, double Pz, double Pt, double Th, double W);
    void AddDetectorDM(double E, double x, double y, double Px, double Py, double Pz, double Pt, double Th, double Phi);
    void AddDetectorNu(double E, double Pz, double Pt, double Th, double W);
    void AddScatterDM(double E, double Pz, double Pt, double Th);
    void AddScatterNu(double E, double Pz, double Pt, double Th, double W);
    void AddScatterSigElectron(double E, double Pz, double Pt, double Th);
    void AddScatterBgElectron(double E, double Pz, double Pt, double Th, double W);
    void ScaleDarkmatter(double scale, double detscale);
    void ScaleNeutrinos(double scale, double detscale);
    void NormalizeHistograms();
    void SaveHistograms();

    static std::string run_name;
    TFile* dm_output, *nu_output;
    bool found_nu_output;
    //Misc
    TH1D *dmtime;
    TH2D *dmxy2;
    //Input dark matter distributions
    TH1D *dmpx1, *dmpy1, *dmpz1, *dmpt1, *dme1, *dmthe1, *dmphi1, *dmx2, *dmy2, *dmpx2, *dmpy2, *dmpz2, *dmpt2, *dme2, *dmthe2, *dmphi2, *dmpz3, *dme3, *dmthe3;
    //Input neutrino distributions
    TH1D *nupz1, *nupt1, *nue1, *nuthe1, *nupz2, *nupt2, *nue2, *nuthe2, *nupz3, *nupt3, *nue3, *nuthe3;
    //Signal-electron distributions
    TH1D *dm_epz1, *dm_ept1, *dm_ethe1, *dm_ee1, *dm_ee1smear, *dm_ee1smearr;
    //Background-electron distributions
    TH1D *nu_epz1, *nu_ept1, *nu_ethe1, *nu_ee1, *nu_ee1smear, *nu_ee1smearr;
    //Neutrino Analysis
    TH1D *nue_e, *nuebar_e, *numu_e, *numubar_e, *nue_e_ratio, *nuebar_e_ratio, *numu_e_ratio, *numubar_e_ratio;
};
