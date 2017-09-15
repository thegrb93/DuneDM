#pragma once

class TFile;
class TH1D;

class DMHistograms
{
public:
    DMHistograms();
    ~DMHistograms();

    void AddProductionDM(double E, double Pz, double Pt, double Th);
    void AddProductionNu(double E, double Pz, double Pt, double Th, int t, double W);
    void AddDetectorDM(double E, double Pz, double Pt, double Th);
    void AddDetectorNu(double E, double Pz, double Pt, double Th, double W);
    void AddScatterDM(double E, double Pz, double Pt, double Th);
    void AddScatterNu(double E, double Pz, double Pt, double Th, double W);
    void AddScatterSigElectron(double E, double Pz, double Pt, double Th);
    void AddScatterBgElectron(double E, double Pz, double Pt, double Th, double W);
    void ScaleNeutrinos(double scale);
    void NormalizeHistograms();
    void SaveHistograms();

    TFile* output;
    //Misc
    TH1D *dmtime;
    //Input dark matter distributions
    TH1D *dmpz1, *dmpt1, *dme1, *dmthe1, *dmpz2, *dmpt2, *dme2, *dmthe2, *dmpz3, *dme3, *dme3smear, *dme3smearr;
    //Input neutrino distributions
    TH1D *nupz1, *nupt1, *nue1, *nuthe1, *nupz2, *nupt2, *nue2, *nuthe2, *nupz3, *nue3, *nue3smear, *nue3smearr;
    //Signal-electron distributions
    TH1D *dm_epz1, *dm_ept1, *dm_ethe1, *dm_ee1, *dm_ee1smear, *dm_ee1smearr;
    //Background-electron distributions
    TH1D *nu_epz1, *nu_ept1, *nu_ethe1, *nu_ee1, *nu_ee1smear, *nu_ee1smearr;
    //Neutrino Analysis
    TH1D *nue_e, *nuebar_e, *numu_e, *numubar_e;
};
