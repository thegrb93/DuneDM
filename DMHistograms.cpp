#include "DMHistograms.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TF1.h>
#include <TCanvas.h>
#include <THStack.h>
#include <sstream>
#include <sys/stat.h>

std::string DMHistograms::run_name;

DMHistograms::DMHistograms(std::string rname) {
    run_name = rname;
    mkdir("histograms", S_IRWXU | S_IRWXG | S_IRWXO);
    mkdir("root", S_IRWXU | S_IRWXG | S_IRWXO);
    mkdir(("histograms/"+rname).c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    
    dm_output = new TFile(("root/"+rname+".root").c_str(), "RECREATE");

    dmpx1 = new TH1D("dmpx1","#chi Production P_{x};P_{x} (GeV/c)", 100, -20, 20);
    dmpy1 = new TH1D("dmpy1","#chi Production P_{y};P_{y} (GeV/c)", 100, -20, 20);
    dmpz1 = new TH1D("dmpz1","#chi Production P_{z};P_{z} (GeV/c)", 100, 0, 60);
    //dmpt1 = new TH1D("dmpt1","#chi Production P_{t};P_{t} (GeV/c)", 100, 0, 0);
    dme1 = new TH1D("dme1","#chi Production Energy;E (GeV)", 100, 0, 60);
    //dmtime = new TH1D("dmtime","#chi vs #nu  arrival time;Time (s)", 100, 0, 1e-9);
    dmthe1 = new TH1D("dmthe1","#chi Production #theta;#theta", 100, 0, M_PI);
    dmphi1 = new TH1D("dmphi1","#chi Production #phi;#phi", 100, -M_PI, M_PI);
    dmxy2 = new TH2D("dmxy2", "#chi in Detector x/y; x (m); y (m)", 80, -2, 2, 60, -1.5, 1.5);
    dmx2 = new TH1D("dmx2","#chi in Detector x;x (m)", 100, -1.5, 1.5);
    dmy2 = new TH1D("dmy2","#chi in Detector y;y (m)", 100, -1, 1);
    dmpx2 = new TH1D("dmpx2","#chi in Detector P_{x};P_{x} (GeV/c)", 100, -0.4, 0.4);
    dmpy2 = new TH1D("dmpy2","#chi in Detector P_{y};P_{y} (GeV/c)", 100, -0.4, 0.4);
    dmpz2 = new TH1D("dmpz2","#chi in Detector P_{z};P_{z} (GeV/c)", 100, 0, 60);
    dmpt2 = new TH1D("dmpt2","#chi in Detector P_{t};P_{t} (GeV/c)", 100, 0, 0);
    dme2 = new TH1D("dme2","#chi in Detector E;E (GeV)", 100, 0, 100);
    dmthe2 = new TH1D("dmthe2","#chi in Detector #theta;#theta", 100, 0, 0.006);
    dmphi2 = new TH1D("dmphi2","#chi in Detector #phi;#phi", 100, -M_PI, M_PI);
    dmpz3 = new TH1D("dmpz3","#chi Scatter P_{z};P_{z} (GeV/c)", 100, 0, 0);
    dme3 = new TH1D("dme3","#chi Scatter E;E (GeV)", 100, 0, 0);
    dmthe3 = new TH1D("dmthe3","#chi Scatter #theta;#theta", 100, 0, 0.006);
    dm_epz1 = new TH1D("dm_epz1","#chi + e^{-} Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
    dm_ept1 = new TH1D("dm_ept1","#chi + e^{-} Scatter P_{t};P_{t} (GeV/c)", 100, 0, 6);
    dm_ethe1 = new TH1D("dm_ethe1","#chi + e^{-} Scatter #theta;#theta)", 100, 0, 3.2);
    dm_ee1 = new TH1D("dm_ee1","#chi + e^{-} Scatter E;E (GeV)", 100, 0, 6);
    dm_ee1smear = new TH1D("dm_ee1smear","#chi + e^{-} Smeared Scatter E;E (GeV)", 100, 0, 6);
    dm_ee1smearr = new TH1D("dm_ee1smearr","#chi + e^{-} Scatter Energy Smearing Ratio;E (GeV)", 40, 0, 6);

    const char* nu_filen = "data/nu_output.root";
    struct stat buffer;
    if(stat(nu_filen, &buffer) == 0)
    {
        nu_output = new TFile(nu_filen, "READ");
        found_nu_output = true;
        nu_output->GetObject("nupz1", nupz1);
        nu_output->GetObject("nupt1", nupt1);
        nu_output->GetObject("nue1", nue1);
        nu_output->GetObject("nuthe1", nuthe1);
        nu_output->GetObject("nupz2", nupz2);
        nu_output->GetObject("nupt2", nupt2);
        nu_output->GetObject("nue2", nue2);
        nu_output->GetObject("nuthe2", nuthe2);
        nu_output->GetObject("nu_epz1", nu_epz1);
        nu_output->GetObject("nu_ept1", nu_ept1);
        nu_output->GetObject("nu_ethe1", nu_ethe1);
        nu_output->GetObject("nu_ee1", nu_ee1);
        nu_output->GetObject("nu_ee1smear", nu_ee1smear);
        nu_output->GetObject("nu_ee1smearr", nu_ee1smearr);
    }
    else
    {
        std::cout << "Creating new nu_data!\n";
        found_nu_output = false;
        nu_output = new TFile(nu_filen, "RECREATE");
        
        nupz1 = new TH1D("nupz1","#nu Production P_{z};P_{z} (GeV/c)", 100, 0, 60);
        nupt1 = new TH1D("nupt1","#nu Production P_{t};P_{t} (GeV/c)", 100, 0, 0);
        nue1 = new TH1D("nue1","#nu Production Energy;E (GeV)", 100, 0, 60);
        nuthe1 = new TH1D("nuthe1","#nu Production #theta;#theta", 100, 0, M_PI);
        nupz2 = new TH1D("nupz2","#nu in Detector P_{z};P_{z} (GeV/c)", 100, 0, 60);
        nupt2 = new TH1D("nupt2","#nu in Detector P_{t};P_{t} (GeV/c)", 100, 0, 0);
        nue2 = new TH1D("nue2","#nu in Detector E;E (GeV)", 100, 0, 60);
        nuthe2 = new TH1D("nuthe2","#nu in Detector #theta;#theta", 100, 0, 0);
        nupz3 = new TH1D("nupz3","#nu Scatter P_{z};P_{z} (GeV/c)", 100, 0, 0);
        nupt3 = new TH1D("nupt3","#nu Scatter P_{t};P_{t} (GeV/c)", 100, 0, 0);
        nue3 = new TH1D("nue3","#nu Scatter E;E (GeV)", 100, 0, 0);
        nuthe3 = new TH1D("nuthe3","#nu in Detector #theta;#theta", 100, 0, 0);
        nu_epz1 = new TH1D("nu_epz1","#nu + e^{-} Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
        nu_ept1 = new TH1D("nu_ept1","#nu + e^{-} Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
        nu_ethe1 = new TH1D("nu_ethe1","#nu + e^{-} Scatter #theta;#theta)", 100, 0, 3.2);
        nu_ee1 = new TH1D("nu_ee1","#nu + e^{-} Scatter E;E (GeV)", 100, 0, 6);
        nu_ee1smear = new TH1D("nu_ee1smear","#nu + e^{-} Smeared Scatter E;E (GeV)", 100, 0, 6);
        nu_ee1smearr = new TH1D("nu_ee1smearr","#nu + e^{-} Scatter Energy Smearing Ratio;E (GeV)", 40, 0, 6);
    }

    /*
    // more complex binning for passing to fastMC
    std::vector<double> bins;
    // 0.125 GeV bins up to 8 GeV
    for(int i = 0; i<8/0.125; i++)
        bins.push_back(i*.125);
    // 0.5 GeV bins up to 20 GeV
    for(int i = 0; i<(20-8)/0.5; i++)
        bins.push_back(8.0+i*.5);
    // 2.0 GeV bins up to 120 GeV
    for(int i = 0; i<(120-20)/2.0; i++)
        bins.push_back(20.0+i*2.0);

    nue_e = new TH1D("nue_e","#nu_{e} Energy;E (GeV)", bins.size()-1, bins.data());
    nuebar_e = new TH1D("nuebar_e","#bar{#nu}_{e} Energy;E (GeV)", bins.size()-1, bins.data());
    numu_e = new TH1D("numu_e","#nu_{#mu} Energy;E (GeV)", bins.size()-1, bins.data());
    numubar_e = new TH1D("numubar_e","#bar{#nu}_{#mu} Energy;E (GeV)", bins.size()-1, bins.data());

    nue_e->Sumw2();
    nuebar_e->Sumw2();
    numu_e->Sumw2();
    numubar_e->Sumw2();

    nue_e_ratio = new TH1D("nue_e_r", "", bins.size()-1, bins.data());
    nuebar_e_ratio = new TH1D("nuebar_e_r", "", bins.size()-1, bins.data());
    numu_e_ratio = new TH1D("numu_e_r", "", bins.size()-1, bins.data());
    numubar_e_ratio = new TH1D("numubar_e_r", "", bins.size()-1, bins.data());*/
}

DMHistograms::~DMHistograms() {
    delete dmpx1;
    delete dmpy1;
    delete dmpz1;
    //delete dmpt1;
    delete dme1;
    //delete dmtime;
    delete dmthe1;
    delete dmphi1;
    delete dmxy2;
    delete dmx2;
    delete dmy2;
    delete dmpx2;
    delete dmpy2;
    delete dmpz2;
    delete dmpt2;
    delete dme2;
    delete dmthe2;
    delete dmphi2;
    delete dmpz3;
    delete dme3;
    delete dmthe3;
    delete dm_epz1;
    delete dm_ept1;
    delete dm_ethe1;
    delete dm_ee1;
    delete dm_ee1smear;
    delete dm_ee1smearr;

    delete nupz1;
    delete nupt1;
    delete nue1;
    delete nuthe1;
    delete nupz2;
    delete nupt2;
    delete nue2;
    delete nuthe2;
    delete nupz3;
    delete nupt3;
    delete nue3;
    delete nuthe3;
    delete nu_epz1;
    delete nu_ept1;
    delete nu_ethe1;
    delete nu_ee1;
    delete nu_ee1smear;
    delete nu_ee1smearr;

    /*delete nue_e;
    delete nuebar_e;
    delete numu_e;
    delete numubar_e;
    */

    delete dm_output;
    delete nu_output;
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
    canvas->SaveAs(("histograms/"+DMHistograms::run_name+"/"+std::string(savename)+".png").c_str());

    canvas->SetLogy();
    canvas->Update();
    //canvas->Draw();
    canvas->SaveAs(("histograms/"+DMHistograms::run_name+"/"+std::string(savename)+"_log.png").c_str());

    delete legend;
}

void normalizeHisto(TH1* hist){
    hist->Scale(1/hist->Integral());
};

void saveHistogram(TCanvas* canvas, TH1D* hist) {
    std::stringstream filen;
    filen << "histograms/"+DMHistograms::run_name+"/" << hist->GetName() << ".png";
    hist->Draw();
    canvas->Update();
    canvas->SaveAs(filen.str().c_str());
}
void saveHistogram2D(TCanvas* canvas, TH2D* hist) {
    std::stringstream filen;
    filen << "histograms/"+DMHistograms::run_name+"/" << hist->GetName() << ".png";
    hist->Draw("COLZ");
    canvas->Update();
    canvas->SaveAs(filen.str().c_str());
}

void DMHistograms::NormalizeHistograms() {
    normalizeHisto(dmpx1);
    normalizeHisto(dmpy1);
    normalizeHisto(dmpz1);
    //normalizeHisto(dmpt1);
    normalizeHisto(dme1);
    //normalizeHisto(dmtime);
    normalizeHisto(dmthe1);
    normalizeHisto(dmphi1);
    normalizeHisto(dmxy2);
    normalizeHisto(dmx2);
    normalizeHisto(dmy2);
    normalizeHisto(dmpx2);
    normalizeHisto(dmpy2);
    normalizeHisto(dmpz2);
    normalizeHisto(dmpt2);
    normalizeHisto(dme2);
    normalizeHisto(dmthe2);
    normalizeHisto(dmphi2);
    normalizeHisto(dmpz3);
    normalizeHisto(dme3);
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
    normalizeHisto(nupt3);
    normalizeHisto(nue3);
    normalizeHisto(nu_epz1);
    normalizeHisto(nu_ept1);
    normalizeHisto(nu_ethe1);
    normalizeHisto(nu_ee1);
    normalizeHisto(nu_ee1smear);
    normalizeHisto(nu_ee1smearr);

    /*normalizeHisto(nue_e);
    normalizeHisto(nuebar_e);
    normalizeHisto(numu_e);
    normalizeHisto(numubar_e);*/
}

void DMHistograms::SaveHistograms() {
    dmpx1->BufferEmpty();
    dmpy1->BufferEmpty();
    dmpz1->BufferEmpty();
    //dmpt1->BufferEmpty();
    dme1->BufferEmpty();
    dmthe1->BufferEmpty();
    dmphi1->BufferEmpty();
    //dmtime->BufferEmpty();
    dmx2->BufferEmpty();
    dmy2->BufferEmpty();
    dmxy2->BufferEmpty();
    dmpx2->BufferEmpty();
    dmpy2->BufferEmpty();
    dmpz2->BufferEmpty();
    dmpt2->BufferEmpty();
    dme2->BufferEmpty();
    dmthe2->BufferEmpty();
    dmphi2->BufferEmpty();
    dme3->BufferEmpty();
    dmpz3->BufferEmpty();
    dme3->BufferEmpty();
    dmthe3->BufferEmpty();
    dm_epz1->BufferEmpty();
    dm_ept1->BufferEmpty();
    dm_ethe1->BufferEmpty();
    dm_ee1->BufferEmpty();
    dm_ee1smear->BufferEmpty();
    dm_ee1smearr->BufferEmpty();

    nupz1->BufferEmpty();
    nue1->BufferEmpty();
    nupz2->BufferEmpty();
    nuthe1->BufferEmpty();
    nupz2->BufferEmpty();
    nue2->BufferEmpty();
    nuthe2->BufferEmpty();
    nupz3->BufferEmpty();
    nupt3->BufferEmpty();
    nue3->BufferEmpty();
    nuthe3->BufferEmpty();
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

    std::cout << "Writing DM RootFile.\n";
    dm_output->Write();
    if(!found_nu_output){
        std::cout << "Writing NU RootFile.\n";
        nu_output->Write();
    }

    TCanvas* canvas = new TCanvas("c1","canvas",1024,576);

    saveHistogram(canvas, dmpx1);
    saveHistogram(canvas, dmpy1);
    saveHistogram(canvas, dmpz1);
    //saveHistogram(canvas, dmpt1);
    saveHistogram(canvas, dme1);
    //saveHistogram(canvas, dmtime);
    saveHistogram(canvas, dmthe1);
    saveHistogram(canvas, dmphi1);
    saveHistogram(canvas, dmx2);
    saveHistogram(canvas, dmy2);
    saveHistogram2D(canvas, dmxy2);
    saveHistogram(canvas, dmpx2);
    saveHistogram(canvas, dmpy2);
    saveHistogram(canvas, dmpz2);
    saveHistogram(canvas, dmpt2);
    saveHistogram(canvas, dme2);
    saveHistogram(canvas, dmthe2);
    saveHistogram(canvas, dmphi2);
    saveHistogram(canvas, dmpz3);
    saveHistogram(canvas, dme3);
    saveHistogram(canvas, dmthe3);
    
    saveHistogram(canvas, nupz1);
    saveHistogram(canvas, nupt1);
    saveHistogram(canvas, nue1);
    saveHistogram(canvas, nuthe1);
    
    saveHistogram(canvas, nupz2);
    saveHistogram(canvas, nupt2);
    saveHistogram(canvas, nue2);
    saveHistogram(canvas, nuthe2);
    
    saveHistogram(canvas, nupz3);
    saveHistogram(canvas, nupt3);
    saveHistogram(canvas, nue3);
    saveHistogram(canvas, nuthe3);
    
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

    /*saveHistogram(canvas, nue_e);
    saveHistogram(canvas, nuebar_e);
    saveHistogram(canvas, numu_e);
    saveHistogram(canvas, numubar_e);

    TFile* cdr = new TFile("/home/garrett/Desktop/DUNE2015CDRFluxes/LBNO_Optimized_80GeV_StandardDP/g4lbne_v3r2p4b_FHC_ND.root");
    TH1D* nue_e_cdr, *nuebar_e_cdr, *numu_e_cdr, *numubar_e_cdr;
    cdr->GetObject("nue_flux", nue_e_cdr);
    cdr->GetObject("nuebar_flux", nuebar_e_cdr);
    cdr->GetObject("numu_flux", numu_e_cdr);
    cdr->GetObject("numubar_flux", numubar_e_cdr);

    nue_e_ratio->Divide(nue_e, nue_e_cdr);
    nuebar_e_ratio->Divide(nuebar_e, nue_e_cdr);
    numu_e_ratio->Divide(numu_e, nue_e_cdr);
    numubar_e_ratio->Divide(numubar_e, nue_e_cdr);*/

    //saveComparison(canvas, "nu_dmpx1", "DM vs Neutrino P_{x};P_{x} (GeV/c)", dmpx1, nupx1, "Dark matter P_{x}", "Neutrino P_{x}");
    //saveComparison(canvas, "nu_dmpy1", "DM vs Neutrino P_{y};P_{y} (GeV/c)", dmpy1, nupy1, "Dark matter P_{y}", "Neutrino P_{y}");
    saveComparison(canvas, "nu_dmpz1", "DM vs Neutrino P_{z};P_{z} (GeV/c)", dmpz1, nupz1, "Dark matter P_{z}", "Neutrino P_{z}");
    //saveComparison(canvas, "nu_dmpt1", "DM vs Neutrino P_{t};P_{t} (GeV/c)", dmpt1, nupt1, "Dark matter P_{t}", "Neutrino P_{t}");
    saveComparison(canvas, "nu_dme1", "DM vs Neutrino E;E (GeV)", dme1, nue1, "Dark matter E", "Neutrino E");
    saveComparison(canvas, "nu_dmet1", "DM vs Neutrino #theta;#theta  (rad)", dmthe1, nuthe1, "Dark matter #theta", "Neutrino #theta");
    saveComparison(canvas, "nu_dmet2", "DM vs Neutrino Intersections #theta;#theta  (rad)", dmthe2, nuthe2, "Dark matter #theta", "Neutrino #theta");
    saveComparison(canvas, "nu_dmpz2", "DM vs Neutrino Intersections P_{z};P_{z} (GeV/c)", dmpz2, nupz2, "Dark matter P_{z}", "Neutrino P_{z}");
    saveComparison(canvas, "nu_dmpt2", "DM vs Neutrino Intersections P_{t};P_{t} (GeV/c)", dmpt2, nupt2, "Dark matter P_{t}", "Neutrino P_{t}");
    saveComparison(canvas, "nu_dme2", "DM vs Neutrino Intersections E;E (GeV)", dme2, nue2, "Dark matter E", "Neutrino E");
    saveComparison(canvas, "nu_dm_epz", "DM vs Neutrino Electron Scatter P_{z};P_{z} (GeV/c)", dm_epz1, nu_epz1, "Electron - Dark matter P_{z}", "Electron - Neutrino P_{z}");
    saveComparison(canvas, "nu_dm_ept", "DM vs Neutrino Electron Scatter P_{t};P_{t} (GeV/c)", dm_ept1, nu_ept1, "Electron - Dark matter P_{t}", "Electron - Neutrino P_{t}");
    saveComparison(canvas, "nu_dm_ee", "DM vs Neutrino Electron Scatter E;E (GeV)", dm_ee1, nu_ee1, "Electron - Dark matter E", "Electron - Neutrino E");
    saveComparison(canvas, "nu_dm_etheta", "DM vs Neutrino Electron Scatter #theta;#theta (rad)", dm_ethe1, nu_ethe1, "Electron - Dark matter #theta", "Electron - Neutrino #theta");

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

void DMHistograms::AddProductionDM(TH1D* E, TH1D* Px, TH1D* Py, TH1D* Pz, TH1D* Th, TH1D* Phi) {
    dme1->Add(E);
    dmpx1->Add(Px);
    dmpy1->Add(Py);
    dmpz1->Add(Pz);
    dmthe1->Add(Th);
    dmphi1->Add(Phi);
}
void DMHistograms::AddDetectorDM(double E, double x, double y, double Px, double Py, double Pz, double Pt, double Th, double Phi) {
    dme2->Fill(E);
    dmx2->Fill(x);
    dmy2->Fill(y);
    dmxy2->Fill(x, y);
    dmpx2->Fill(Px);
    dmpy2->Fill(Py);
    dmpz2->Fill(Pz);
    dmpt2->Fill(Pt);
    dmthe2->Fill(Th);
    dmphi2->Fill(Phi);
}
void DMHistograms::AddScatterDM(double E, double Pz, double Pt, double Th) {
    dme3->Fill(E);
    dmpz3->Fill(Pz);
    //dmpt3->Fill(Pt);
    dmthe3->Fill(Th);

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

void DMHistograms::AddProductionNu(double E, double Pz, double Pt, double Th, double W) {
    nue1->Fill(E, W);
    nupz1->Fill(Pz, W);
    nupt1->Fill(Pt, W);
    nuthe1->Fill(Th, W);

    /*switch(t)
    {
        case 53:
            nue_e->Fill(E, W);
            break;
        case 52:
            nuebar_e->Fill(E, W);
            break;
        case 56:
            numu_e->Fill(E, W);
            break;
        case 55:
            numubar_e->Fill(E, W);
            break;
        default:
            break;
    }*/
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

