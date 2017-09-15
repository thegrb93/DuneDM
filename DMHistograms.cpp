#include "DMHistograms.h"

#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TF1.h>
#include <TCanvas.h>
#include <THStack.h>
#include <sstream>
#include <sys/stat.h>


DMHistograms::DMHistograms()
{
    output = new TFile("output.root", "RECREATE");

    dmpz1 = new TH1D("dmpz1","\\chi P_{z};P_{z} (GeV/c)", 100, 0, 60);
    dmpt1 = new TH1D("dmpt1","\\chi P_{t};P_{t} (GeV/c)", 100, 0, 0);
    dme1 = new TH1D("dme1","\\chi E;E (GeV)", 100, 0, 60);
    dmtime = new TH1D("dmtime","\\chi vs \\nu  arrival time;Time (s)", 100, 0, 1e-9);
    dmthe1 = new TH1D("dmthe1","\\chi \\theta;\\theta", 100, 0, M_PI);
    dmpz2 = new TH1D("dmpz2","\\chi in Detector P_{z};P_{z} (GeV/c)", 100, 0, 60);
    dmpt2 = new TH1D("dmpt2","\\chi in Detector P_{t};P_{t} (GeV/c)", 100, 0, 0);
    dme2 = new TH1D("dme2","\\chi in Detector E;E (GeV)", 100, 0, 60);
    dmthe2 = new TH1D("dmthe2","\\chi in Detector \\theta;\\theta", 100, 0, 0);
    //dmpz3 = new TH1D("dmpz3","\\chi Scatter P_{z};P_{z} (GeV/c)", 100, 0, 0);
    //dme3 = new TH1D("dme3","\\chi Scatter E;E (GeV)", 100, 0, 0);
    //dme3smear = new TH1D("dme3smear","\\chi  Smeared Scatter E;E (GeV)", 100, 0, 0);
    //dme3smearr = new TH1D("dme3smearr","\\chi  Scatter Energy Smearing Ratio;E (GeV)", 100, 0, 0);
    dm_epz1 = new TH1D("dm_epz1","\\chi - e^{-} Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
    dm_ept1 = new TH1D("dm_ept1","\\chi - e^{-} Scatter P_{t};P_{t} (GeV/c)", 100, 0, 6);
    dm_ethe1 = new TH1D("dm_ethe1","\\chi - e^{-} Scatter \\theta;\\theta)", 100, 0, 3.2);
    dm_ee1 = new TH1D("dm_ee1","\\chi - e^{-} Scatter E;E (GeV)", 100, 0, 6);
    dm_ee1smear = new TH1D("dm_ee1smear","\\chi - e^{-} Smeared Scatter E;E (GeV)", 100, 0, 6);
    dm_ee1smearr = new TH1D("dm_ee1smearr","\\chi - e^{-} Scatter Energy Smearing Ratio;E (GeV)", 40, 0, 6);

    nupz1 = new TH1D("nupz1","\\nu P_{z};P_{z} (GeV/c)", 100, 0, 60);
    nupt1 = new TH1D("nupt1","\\nu P_{t};P_{t} (GeV/c)", 100, 0, 0);
    nue1 = new TH1D("nue1","\\nu E;E (GeV)", 100, 0, 60);
    nuthe1 = new TH1D("nuthe1","\\nu \\theta;\\theta", 100, 0, M_PI);
    nupz2 = new TH1D("nupz2","\\nu in Detector P_{z};P_{z} (GeV/c)", 100, 0, 60);
    nupt2 = new TH1D("nupt2","\\nu in Detector P_{t};P_{t} (GeV/c)", 100, 0, 0);
    nue2 = new TH1D("nue2","\\nu in Detector E;E (GeV)", 100, 0, 60);
    nuthe2 = new TH1D("nuthe2","\\nu in Detector \\theta;\\theta", 100, 0, 0);
    //nupz3 = new TH1D("nupz3","\\nu Scatter P_{z};P_{z} (GeV/c)", 100, 0, 0);
    //nue3 = new TH1D("nue3","\\nu Scatter E;E (GeV)", 100, 0, 0);
    //nue3smear = new TH1D("nue3smear","\\nu Smeared Scatter E;E (GeV)", 100, 0, 0);
    //nue3smearr = new TH1D("nue3smearr","\\nu Scatter Energy Smearing Ratio;E (GeV)", 100, 0, 0);
    nu_epz1 = new TH1D("nu_epz1","\\nu - e^{-} Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
    nu_ept1 = new TH1D("nu_ept1","\\nu - e^{-} Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
    nu_ethe1 = new TH1D("nu_ethe1","\\nu - e^{-} Scatter \\theta;\\theta)", 100, 0, 3.2);
    nu_ee1 = new TH1D("nu_ee1","\\nu - e^{-} Scatter E;E (GeV)", 100, 0, 6);
    nu_ee1smear = new TH1D("nu_ee1smear","\\nu - e^{-} Smeared Scatter E;E (GeV)", 100, 0, 6);
    nu_ee1smearr = new TH1D("nu_ee1smearr","\\nu - e^{-} Scatter Energy Smearing Ratio;E (GeV)", 40, 0, 6);

    nue_e = new TH1D("nue_e","\\nu_{e} Energy;E (GeV)", 100, 0, 80);
    nuebar_e = new TH1D("nuebar_e","\\bar{\\nu}_{e} Energy;E (GeV)", 100, 0, 80);
    numu_e = new TH1D("numu_e","\\nu_{\\mu} Energy;E (GeV)", 100, 0, 80);
    numubar_e = new TH1D("numubar_e","\\bar{\\nu}_{\\mu} Energy;E (GeV)", 100, 0, 80);
}

DMHistograms::~DMHistograms()
{
    delete dmpz1;
    delete dmpt1;
    delete dme1;
    delete dmpz2;
    delete dmpt2;
    delete dme2;
    delete dmtime;
    delete dmthe1;
    delete dmthe2;
    //delete dmpz3;
    //delete dme3;
    //delete dme3smear;
    //delete dme3smearr;
    delete dm_epz1;
    delete dm_ept1;
    delete dm_ethe1;
    delete dm_ee1;
    delete dm_ee1smear;

    delete nupz1;
    delete nupt1;
    delete nue1;
    delete nupz2;
    delete nupt2;
    delete nue2;
    delete nuthe1;
    delete nuthe2;
    //delete nupz3;
    //delete nue3;
    //delete nue3smear;
    //delete nue3smearr;
    delete nu_epz1;
    delete nu_ept1;
    delete nu_ethe1;
    delete nu_ee1;
    delete nu_ee1smear;
    delete nu_ee1smearr;

    delete nue_e;
    delete nuebar_e;
    delete numu_e;
    delete numubar_e;

    delete output;
}

void saveComparison(TCanvas* canvas, const char* savename, const char* canvastitle, TH1D* hist1, TH1D* hist2, const char* histn1, const char* histn2) {
    canvas->SetFillColor(33);
    canvas->SetFrameFillColor(17);
    canvas->SetGrid();
    canvas->SetLogy(0);

    THStack *hs = new THStack("hs",canvastitle);

    // draw the legend
    TLegend *legend=new TLegend(0.6,0.75,0.88,0.85);
    legend->SetTextFont(72);
    legend->SetTextSize(0.02);

    hist1->SetLineColor(1);
    //hist1->SetFillColor(1);
    //hist1->SetMarkerStyle(21);
    //hist1->SetMarkerSize(0.8);
    //hist1->SetStats(0);
    hs->Add(hist1);

    hist2->SetLineColor(2);
    //hist2->SetFillColor(2);
    //hist1->SetMarkerStyle(21);
    //hist2->SetMarkerSize(0.8);
    //hist2->SetStats(0);
    hs->Add(hist2);
    hs->Draw("hist nostack");

    legend->AddEntry(hist1,histn1);
    legend->AddEntry(hist2,histn2);
    legend->Draw();

    canvas->Update();
    //canvas->Draw();
    canvas->SaveAs(("histograms/"+std::string(savename)+".png").c_str());

    canvas->SetLogy();
    canvas->Update();
    //canvas->Draw();
    canvas->SaveAs(("histograms/"+std::string(savename)+"_log.png").c_str());

    delete legend;
}

void normalizeHisto(TH1* hist){
    hist->Scale(1/hist->Integral());
};

void saveHistogram(TCanvas* canvas, TH1D* hist) {
    std::stringstream filen;
    filen << "histograms/" << hist->GetName() << ".png";
    hist->Draw();
    canvas->Update();
    canvas->SaveAs(filen.str().c_str());
}

void DMHistograms::NormalizeHistograms() {
    normalizeHisto(dmpz1);
    normalizeHisto(dmpt1);
    normalizeHisto(dme1);
    normalizeHisto(dmpz2);
    normalizeHisto(dmpt2);
    normalizeHisto(dme2);
    normalizeHisto(dmtime);
    normalizeHisto(dmthe1);
    normalizeHisto(dmthe2);
    normalizeHisto(dmpz3);
    normalizeHisto(dme3);
    normalizeHisto(dme3smear);
    normalizeHisto(dme3smearr);
    normalizeHisto(dm_epz1);
    normalizeHisto(dm_ept1);
    normalizeHisto(dm_ethe1);
    normalizeHisto(dm_ee1);
    normalizeHisto(dm_ee1smear);
    normalizeHisto(dm_ee1smearr);

    normalizeHisto(nupz1);
    normalizeHisto(nue1);
    normalizeHisto(nupz2);
    normalizeHisto(nue2);
    normalizeHisto(nuthe1);
    normalizeHisto(nuthe2);
    normalizeHisto(nupz3);
    normalizeHisto(nue3);
    normalizeHisto(nue3smear);
    normalizeHisto(nue3smearr);
    normalizeHisto(nu_epz1);
    normalizeHisto(nu_ept1);
    normalizeHisto(nu_ethe1);
    normalizeHisto(nu_ee1);
    normalizeHisto(nu_ee1smear);
    normalizeHisto(nu_ee1smearr);

    normalizeHisto(nue_e);
    normalizeHisto(nuebar_e);
    normalizeHisto(numu_e);
    normalizeHisto(numubar_e);
}

void DMHistograms::SaveHistograms() {
    dmpz1->BufferEmpty();
    dmpt1->BufferEmpty();
    dmthe1->BufferEmpty();
    dme1->BufferEmpty();
    dmpz2->BufferEmpty();
    dmpt2->BufferEmpty();
    dme2->BufferEmpty();
    dmtime->BufferEmpty();
    dmthe2->BufferEmpty();
    //dmpz3->BufferEmpty();
    //dme3->BufferEmpty();
    //dme3smear->BufferEmpty();
    //dme3smearr->BufferEmpty();
    dm_epz1->BufferEmpty();
    dm_ept1->BufferEmpty();
    dm_ethe1->BufferEmpty();
    dm_ee1->BufferEmpty();
    dm_ee1smear->BufferEmpty();
    dm_ee1smearr->BufferEmpty();

    nupz1->BufferEmpty();
    nue1->BufferEmpty();
    nupz2->BufferEmpty();
    nue2->BufferEmpty();
    nuthe1->BufferEmpty();
    nuthe2->BufferEmpty();
    //nupz3->BufferEmpty();
    //nue3->BufferEmpty();
    //nue3smear->BufferEmpty();
    //nue3smearr->BufferEmpty();
    nu_epz1->BufferEmpty();
    nu_ept1->BufferEmpty();
    nu_ethe1->BufferEmpty();
    nu_ee1->BufferEmpty();
    nu_ee1smear->BufferEmpty();
    nu_ee1smearr->BufferEmpty();

    double pare[3] = {1,1,1};
    //TF1 *ge = new TF1("ge","gaus(0)",0.8,1.2);
    TF1 *ge = new TF1("ge","gaus(0)",dm_ee1smearr->GetXaxis()->GetXmin(),dm_ee1smearr->GetXaxis()->GetXmax());
    ge->SetParameters(pare);
    dm_ee1smearr->Fit(ge,"R");

    output->Write();

    TCanvas* canvas = new TCanvas("c1","canvas",1024,1024);
    mkdir("histograms", S_IRWXU | S_IRWXG | S_IRWXO);

    saveHistogram(canvas, dmpz1);
    saveHistogram(canvas, dmpt1);
    saveHistogram(canvas, dme1);
    saveHistogram(canvas, dmpz2);
    saveHistogram(canvas, dmpt2);
    saveHistogram(canvas, dme2);
    saveHistogram(canvas, dmtime);
    saveHistogram(canvas, dmthe1);
    saveHistogram(canvas, dmthe2);
    //saveHistogram(canvas, dmpz3);
    //saveHistogram(canvas, dme3);
    //saveHistogram(canvas, dme3smear);
    //saveHistogram(canvas, dme3smearr);
    saveHistogram(canvas, nupz1);
    saveHistogram(canvas, nupt1);
    saveHistogram(canvas, nue1);
    saveHistogram(canvas, nupz2);
    saveHistogram(canvas, nupt2);
    saveHistogram(canvas, nue2);
    saveHistogram(canvas, nuthe1);
    saveHistogram(canvas, nuthe2);
    //saveHistogram(canvas, nupz3);
    //saveHistogram(canvas, nue3);
    //saveHistogram(canvas, nue3smear);
    //saveHistogram(canvas, nue3smearr);
    saveHistogram(canvas, dm_epz1);
    saveHistogram(canvas, dm_ept1);
    saveHistogram(canvas, dm_ethe1);
    saveHistogram(canvas, dm_ee1);
    saveHistogram(canvas, dm_ee1smear);
    saveHistogram(canvas, dm_ee1smearr);
    saveHistogram(canvas, nu_epz1);
    saveHistogram(canvas, nu_ept1);
    saveHistogram(canvas, nu_ethe1);
    saveHistogram(canvas, nu_ee1);
    saveHistogram(canvas, nu_ee1smear);
    saveHistogram(canvas, nu_ee1smearr);

    saveHistogram(canvas, nue_e);
    saveHistogram(canvas, nuebar_e);
    saveHistogram(canvas, numu_e);
    saveHistogram(canvas, numubar_e);

    saveComparison(canvas, "nu_dmpz1", "DM vs Neutrino P_{z};P_{z} (GeV/c)", dmpz1, nupz1, "Dark matter P_{z}", "Neutrino P_{z}");
    saveComparison(canvas, "nu_dmpt1", "DM vs Neutrino P_{t};P_{t} (GeV/c)", dmpt1, nupt1, "Dark matter P_{t}", "Neutrino P_{t}");
    saveComparison(canvas, "nu_dme1", "DM vs Neutrino E;E (GeV)", dme1, nue1, "Dark matter E", "Neutrino E");
    saveComparison(canvas, "nu_dmet1", "DM vs Neutrino \\theta;\\theta  (rad)", dmthe1, nuthe1, "Dark matter \\theta", "Neutrino \\theta");
    saveComparison(canvas, "nu_dmet2", "DM vs Neutrino Intersections \\theta;\\theta  (rad)", dmthe2, nuthe2, "Dark matter \\theta", "Neutrino \\theta");
    saveComparison(canvas, "nu_dmpz2", "DM vs Neutrino Intersections P_{z};P_{z} (GeV/c)", dmpz2, nupz2, "Dark matter P_{z}", "Neutrino P_{z}");
    saveComparison(canvas, "nu_dmpt2", "DM vs Neutrino Intersections P_{t};P_{t} (GeV/c)", dmpt2, nupt2, "Dark matter P_{t}", "Neutrino P_{t}");
    saveComparison(canvas, "nu_dme2", "DM vs Neutrino Intersections E;E (GeV)", dme2, nue2, "Dark matter E", "Neutrino E");
    saveComparison(canvas, "nu_dm_epz", "DM vs Neutrino Electron Scatter P_{z};P_{z} (GeV/c)", dm_epz1, nu_epz1, "Electron - Dark matter P_{z}", "Electron - Neutrino P_{z}");
    saveComparison(canvas, "nu_dm_ept", "DM vs Neutrino Electron Scatter P_{t};P_{t} (GeV/c)", dm_ept1, nu_ept1, "Electron - Dark matter P_{t}", "Electron - Neutrino P_{t}");
    saveComparison(canvas, "nu_dm_ee", "DM vs Neutrino Electron Scatter E;E (GeV)", dm_ee1, nu_ee1, "Electron - Dark matter E", "Electron - Neutrino E");
    saveComparison(canvas, "nu_dm_etheta", "DM vs Neutrino Electron Scatter \\theta;\\theta (rad)", dm_ethe1, nu_ethe1, "Electron - Dark matter \\theta", "Electron - Neutrino \\theta");

    delete canvas;
}

void DMHistograms::ScaleNeutrinos(double scale) {
    nupz1->Scale(scale);
    nupt1->Scale(scale);
    nuthe1->Scale(scale);
    nue1->Scale(scale);
    nupz2->Scale(scale);
    nupt2->Scale(scale);
    nue2->Scale(scale);
    nuthe2->Scale(scale);
    nu_epz1->Scale(scale);
    nu_ept1->Scale(scale);
    nu_ee1->Scale(scale);
    nu_ethe1->Scale(scale);
}

void DMHistograms::AddProductionDM(double E, double Pz, double Pt, double Th) {
    dme1->Fill(E);
    dmpz1->Fill(Pz);
    dmpt1->Fill(Pt);
    dmthe1->Fill(Th);
}
void DMHistograms::AddDetectorDM(double E, double Pz, double Pt, double Th) {
    dme2->Fill(E);
    dmpz2->Fill(Pz);
    dmpt2->Fill(Pt);
    dmthe2->Fill(Th);
}
void DMHistograms::AddScatterDM(double E, double Pz, double Pt, double Th) {
    dme3->Fill(E);
    dmpz3->Fill(Pz);
    //dmpt3->Fill(Pt);
    //dmthe3->Fill(Th);

    /*if (smear_sigma > 0) {
        E += Random::Gauss(smear_mean, smear_sigma / sqrt(E));
        dme3smear->Fill(darkmatter1.E);
    }*/
}
void DMHistograms::AddScatterSigElectron(double E, double Pz, double Pt, double Th) {
    dm_ee1->Fill(E);
    dm_epz1->Fill(Pz);
    dm_ept1->Fill(Pt);
    dm_ethe1->Fill(Th);

    /*if (smear_sigma > 0) {
        double Es += Random::Gauss(smear_mean, smear_sigma / sqrt(E));
        dm_ee1smear->Fill(Es);
        dm_ee1smearr->Fill(E / Es);
    }*/
}

void DMHistograms::AddProductionNu(double E, double Pz, double Pt, double Th, int t, double W) {
    nue1->Fill(E, W);
    nupz1->Fill(Pz, W);
    nupt1->Fill(Pt, W);
    nuthe1->Fill(Th, W);

    switch(t)
    {
        case 53:
            nue_e->Fill(E/*, W*/);
            break;
        case 52:
            nuebar_e->Fill(E/*, W*/);
            break;
        case 56:
            numu_e->Fill(E/*, W*/);
            break;
        case 55:
            numubar_e->Fill(E/*, W*/);
            break;
        default:
            break;
    }
}

void DMHistograms::AddDetectorNu(double E, double Pz, double Pt, double Th, double W) {
    nue2->Fill(E, W);
    nupz2->Fill(Pz, W);
    nupt2->Fill(Pt, W);
    nuthe2->Fill(Th, W);
}
void DMHistograms::AddScatterNu(double E, double Pz, double Pt, double Th, double W) {
    nue3->Fill(E, W);
    nupz3->Fill(Pz, W);
    //nupt3->Fill(Pt, W);
    //nuthe3->Fill(Th, W);

    /*if (smear_sigma > 0) {
        E += Random::Gauss(smear_mean, smear_sigma / sqrt(E));

        dme3smearr->Fill(darkmatter1.E / darkmatter1E);
        dm_ee1smearr->Fill(electron1.E / electron1E);
    }*/
}
void DMHistograms::AddScatterBgElectron(double E, double Pz, double Pt, double Th, double W) {
    nu_ee1->Fill(E, W);
    nu_epz1->Fill(Pz, W);
    nu_ept1->Fill(Pt, W);
    nu_ethe1->Fill(Th, W);

    /*if (smear_sigma > 0) {
        E += Random::Gauss(smear_mean, smear_sigma / sqrt(E));

        dme3smearr->Fill(darkmatter1.E / darkmatter1E);
        dm_ee1smearr->Fill(electron1.E / electron1E);
    }*/
}
