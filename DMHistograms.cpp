#include "DMHistograms.h"
#include "DMAnalysis.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TColor.h>
#include <THStack.h>
#include <sstream>
#include <sys/stat.h>

std::string DMHistograms::run_name;

DMHistograms::DMHistograms(std::string rname, TFile* dm, TFile* nu, bool nu_exists) {
    run_name = rname;
    dm_output = dm;
    nu_output = nu;
    found_nu_output = nu_exists;

    dm->cd();

    dme1 = new TH1D("dme1","#chi Production Energy;E (GeV)", 100, 0, 60);
    dme2 = new TH1D("dme2","#chi in Detector E;E (GeV)", 100, 0, 100);
    dme3 = new TH1D("dme3","#chi Scatter E;E (GeV)", 100, 0, 100);

    dmpz1 = new TH1D("dmpz1","#chi Production P_{z};P_{z} (GeV/c)", 100, 0, 60);
    dmpz2 = new TH1D("dmpz2","#chi in Detector P_{z};P_{z} (GeV/c)", 100, 0, 100);
    dmpz3 = new TH1D("dmpz3","#chi Scatter P_{z};P_{z} (GeV/c)", 100, 0, 100);

    dmthe1 = new TH1D("dmthe1","#chi Production #theta;#theta", 100, 0, M_PI);
    dmthe2 = new TH1D("dmthe2","#chi in Detector #theta;#theta", 100, 0, 0.009);
    dmthe3 = new TH1D("dmthe3","#chi Scatter #theta;#theta", 100, 0, 0.009);

    dmphi1 = new TH1D("dmphi1","#chi Production #phi;#phi", 100, -M_PI, M_PI);
    dmphi2 = new TH1D("dmphi2","#chi in Detector #phi;#phi", 100, -M_PI, M_PI);

    dmpx1 = new TH1D("dmpx1","#chi Production P_{x};P_{x} (GeV/c)", 100, -20, 20);
    dmpx2 = new TH1D("dmpx2","#chi in Detector P_{x};P_{x} (GeV/c)", 100, -0.4, 0.4);
    dmpy1 = new TH1D("dmpy1","#chi Production P_{y};P_{y} (GeV/c)", 100, -20, 20);
    dmpy2 = new TH1D("dmpy2","#chi in Detector P_{y};P_{y} (GeV/c)", 100, -0.4, 0.4);
    //dmpt1 = new TH1D("dmpt1","#chi Production P_{t};P_{t} (GeV/c)", 100, 0, 0);
    dmpt2 = new TH1D("dmpt2","#chi in Detector P_{t};P_{t} (GeV/c)", 100, 0, 0);


    dmxy2 = new TH2D("dmxy2", "#chi in Detector x/y; x (m); y (m)", 80, -2, 2, 60, -1.5, 1.5);
    dmx2 = new TH1D("dmx2","#chi in Detector x;x (m)", 100, -1.5, 1.5);
    dmy2 = new TH1D("dmy2","#chi in Detector y;y (m)", 100, -1, 1);

    dmtime = new TH1D("dmtime","#chi arrival time;Time (s)", 100, 1.90965e-6, 1.90969e-6);

    dm_ee1 = new TH1D("dm_ee1","#chi + e^{-} Scatter E;E (GeV)", 100, 0, 15);
    dm_epz1 = new TH1D("dm_epz1","#chi + e^{-} Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
    dm_ept1 = new TH1D("dm_ept1","#chi + e^{-} Scatter P_{t};P_{t} (GeV/c)", 100, 0, 6);
    dm_ethe1 = new TH1D("dm_ethe1","#chi + e^{-} Scatter #theta;#theta", 100, 0, 1.6);
    dm_ee1smear = new TH1D("dm_ee1smear","#chi + e^{-} Smeared Scatter E;E (GeV)", 100, 0, 15);
    dm_ee1smearr = new TH1D("dm_ee1smearr","#chi + e^{-} Scatter Energy Smearing Ratio;E (GeV)", 40, 0, 15);

    if(nu_exists)
    {
        nu_output->GetObject("nupz1", nupz1);
        nu_output->GetObject("nupt1", nupt1);
        nu_output->GetObject("nue1", nue1);
        nu_output->GetObject("nuthe1", nuthe1);
        nu_output->GetObject("nupz2", nupz2);
        nu_output->GetObject("nupt2", nupt2);
        nu_output->GetObject("nue2", nue2);
        nu_output->GetObject("nuthe2", nuthe2);
        nu_output->GetObject("nupz3", nupz3);
        nu_output->GetObject("nupt3", nupt3);
        nu_output->GetObject("nue3", nue3);
        nu_output->GetObject("nuthe3", nuthe3);
        nu_output->GetObject("nutime", nutime);
        nu_output->GetObject("nu_epz1", nu_epz1);
        nu_output->GetObject("nu_ept1", nu_ept1);
        nu_output->GetObject("nu_ethe1", nu_ethe1);
        nu_output->GetObject("nu_ee1", nu_ee1);
        nu_output->GetObject("nu_ee1smear", nu_ee1smear);
        nu_output->GetObject("nu_ee1smearr", nu_ee1smearr);
    }
    else
    {
        nu->cd();
        std::cout << "Creating new nu_data!\n";
        nue1 = new TH1D("nue1","#nu Production Energy;E (GeV)", 100, 0, 100);
        nue2 = new TH1D("nue2","#nu in Detector E;E (GeV)", 100, 0, 100);
        nue3 = new TH1D("nue3","#nu Scatter E;E (GeV)", 100, 0, 100);

        nupz1 = new TH1D("nupz1","#nu Production P_{z};P_{z} (GeV/c)", 100, 0, 100);
        nupz2 = new TH1D("nupz2","#nu in Detector P_{z};P_{z} (GeV/c)", 100, 0, 100);
        nupz3 = new TH1D("nupz3","#nu Scatter P_{z};P_{z} (GeV/c)", 100, 0, 100);

        nupt1 = new TH1D("nupt1","#nu Production P_{t};P_{t} (GeV/c)", 100, 0, 0);
        nupt2 = new TH1D("nupt2","#nu in Detector P_{t};P_{t} (GeV/c)", 100, 0, 0);
        nupt3 = new TH1D("nupt3","#nu Scatter P_{t};P_{t} (GeV/c)", 100, 0, 0);

        nuthe1 = new TH1D("nuthe1","#nu Production #theta;#theta", 100, 0, M_PI);
        nuthe2 = new TH1D("nuthe2","#nu in Detector #theta;#theta", 100, 0, 0.009);
        nuthe3 = new TH1D("nuthe3","#nu in Detector #theta;#theta", 100, 0, 0.009);

        nutime = new TH1D("nutime","#nu in Detector time;time", 100, 1.90965e-6, 1.90969e-6);

        nu_ee1 = new TH1D("nu_ee1","#nu + e^{-} Scatter E;E (GeV)", 100, 0, 15);
        nu_epz1 = new TH1D("nu_epz1","#nu + e^{-} Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
        nu_ept1 = new TH1D("nu_ept1","#nu + e^{-} Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
        nu_ethe1 = new TH1D("nu_ethe1","#nu + e^{-} Scatter #theta;#theta", 100, 0, 1.6);

        nu_ee1smear = new TH1D("nu_ee1smear","#nu + e^{-} Smeared Scatter E;E (GeV)", 100, 0, 15);
        nu_ee1smearr = new TH1D("nu_ee1smearr","#nu + e^{-} Scatter Energy Smearing Ratio;E (GeV)", 40, 0, 15);
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
    delete nutime;
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

void saveComparison(TCanvas* canvas, const char* savename, const char* canvastitle, const std::vector<TH1D*>& p1, const std::vector<TH1D*>& p2, const std::vector<std::string>& p1_labels, const std::vector<std::string>& p2_labels) {
    static TColor *c1 = new TColor(9001,1,0.5,0);
    canvas->SetFillColor(33);
    canvas->SetFrameFillColor(17);
    canvas->SetGrid();

    int colors1[] = {1, 2, 9001};
    int colors2[] = {8, 4, 6};

    THStack *hs = new THStack("hs",canvastitle);

    TLegend *legend=new TLegend(0.8, 0.5, 1, 0.8);
    //legend->SetTextFont(72);
    legend->SetTextSize(0.03);

    for(int i = 0; i<std::min<int>(p1.size(), 3); ++i)
    {
        p1[i]->SetLineColor(colors1[i]);
        p1[i]->SetLineWidth(2);
        hs->Add(p1[i]);
        legend->AddEntry(p1[i],p1_labels[i].c_str());
    }
    for(int i = 0; i<std::min<int>(p2.size(), 3); ++i)
    {
        p2[i]->SetLineColor(colors2[i]);
        p2[i]->SetLineWidth(2);
        hs->Add(p2[i]);
        legend->AddEntry(p2[i],p2_labels[i].c_str());
    }

    hs->Draw("hist nostack");
    hs->GetYaxis()->SetTickLength(0);
    legend->Draw();

    canvas->SetLogy(0);
    canvas->Update();
    //canvas->Draw();
    canvas->SaveAs(("histograms/"+DMHistograms::run_name+"/"+std::string(savename)+".png").c_str());

    canvas->SetLogy();
    canvas->Update();
    //canvas->Draw();
    canvas->SaveAs(("histograms/"+DMHistograms::run_name+"/"+std::string(savename)+"_log.png").c_str());

    for(int i = 0; i<std::min<int>(p1.size(), 3); ++i)
        p1[i]->Scale(1.0/p1[i]->Integral());
    for(int i = 0; i<std::min<int>(p2.size(), 3); ++i)
        p2[i]->Scale(1.0/p2[i]->Integral());

    hs->Draw("hist nostack");
    legend->Draw();

    canvas->SetLogy(0);
    canvas->Update();
    //canvas->Draw();
    canvas->SaveAs(("histograms/"+DMHistograms::run_name+"/"+std::string(savename)+"_norm.png").c_str());

    canvas->SetLogy();
    canvas->Update();
    //canvas->Draw();
    canvas->SaveAs(("histograms/"+DMHistograms::run_name+"/"+std::string(savename)+"_norm_log.png").c_str());

    delete legend;
    delete hs;
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
    normalizeHisto(dmtime);
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
    normalizeHisto(nutime);
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

void DMHistograms::SaveNeutrinos() {
    nu_output->Write();
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
    nutime->BufferEmpty();
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

    TCanvas* canvas = new TCanvas("c1","canvas",1920,1080);

    /*saveHistogram(canvas, dmpx1);
    saveHistogram(canvas, dmpy1);
    saveHistogram(canvas, dmpz1);
    //saveHistogram(canvas, dmpt1);
    saveHistogram(canvas, dme1);
    saveHistogram(canvas, dmtime);
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
    saveHistogram(canvas, nutime);
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
    saveHistogram(canvas, nu_ee1smearr);*/

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

    std::vector<std::string> legend_labels_chi{"Production #chi #bar{#chi}", "Fiducial Volume #chi #bar{#chi}", "Scattered #chi #bar{#chi}"};
    std::vector<std::string> legend_labels_nu{"Production #nu_{e}", "Fiducial Volume #nu_{e}", "Scattered #nu_{e}"};
    std::vector<std::string> legend_labels_echi{"e^{-} scattered via #chi #bar{#chi}"};
    std::vector<std::string> legend_labels_enu{"e^{-} scattered via #nu_{e}"};

    //saveComparison(canvas, "nu_dmpx1", "DM vs Neutrino P_{x};P_{x} (GeV/c)", dmpx1, nupx1, "Dark matter P_{x}", "Neutrino P_{x}");
    //saveComparison(canvas, "nu_dmpy1", "DM vs Neutrino P_{y};P_{y} (GeV/c)", dmpy1, nupy1, "Dark matter P_{y}", "Neutrino P_{y}");
    saveComparison(canvas, "nu_dmpz", "DM vs Neutrino P_{z};P_{z} (GeV/c)", std::vector<TH1D*>{dmpz1, dmpz2, dmpz3}, std::vector<TH1D*>{nupz1, nupz2, nupz3}, legend_labels_chi, legend_labels_nu);
    //saveComparison(canvas, "nu_dmpt1", "DM vs Neutrino P_{t};P_{t} (GeV/c)", dmpt1, nupt1, "Dark matter P_{t}", "Neutrino P_{t}");
    saveComparison(canvas, "nu_dme", "DM vs Neutrino Energies;E (GeV)", std::vector<TH1D*>{dme1, dme2, dme3}, std::vector<TH1D*>{nue1, nue2, nue3}, legend_labels_chi, legend_labels_nu);
    saveComparison(canvas, "nu_dmet1", "DM vs Neutrino #theta;#theta (rad)", std::vector<TH1D*>{dmthe1, dmthe2, dmthe3}, std::vector<TH1D*>{nuthe1, nuthe2, nuthe3}, legend_labels_chi, legend_labels_nu);

    saveComparison(canvas, "nu_dm_epz", "DM vs Neutrino Electron Scatter P_{z};P_{z} (GeV/c)", std::vector<TH1D*>{dm_epz1}, std::vector<TH1D*>{nu_epz1}, legend_labels_echi, legend_labels_enu);
    // saveComparison(canvas, "nu_dm_ept", "Electron Scatter P_{t};P_{t} (GeV/c)", dm_ept1, nu_ept1, "Electron - Dark matter P_{t}", "Electron - Neutrino P_{t}");
    saveComparison(canvas, "nu_dm_ee", "DM vs Neutrino Electron Scatter E;E (GeV)", std::vector<TH1D*>{dm_ee1}, std::vector<TH1D*>{nu_ee1}, legend_labels_echi, legend_labels_enu);
    saveComparison(canvas, "nu_dm_etheta", "DM vs Neutrino Electron Scatter #theta;#theta (rad)", std::vector<TH1D*>{dm_ethe1}, std::vector<TH1D*>{nu_ethe1}, legend_labels_echi, legend_labels_enu);
    // saveComparison(canvas, "nu_dm_time", "Electron Time;Time", dmtime, nutime, "Electron - Dark matter Time", "Electron - Neutrino Time");

    delete canvas;
}

void DMHistograms::ScaleDarkmatter(double scale, double detscale) {
    dmpx1->Scale(scale);
    dmpy1->Scale(scale);
    dmpz1->Scale(scale);
    //dmpt1->Scale(scale);
    dme1->Scale(scale);
    dmthe1->Scale(scale);
    dmphi1->Scale(scale);
    dmx2->Scale(scale);
    dmy2->Scale(scale);
    dmxy2->Scale(scale);
    dmpx2->Scale(scale);
    dmpy2->Scale(scale);
    dmpz2->Scale(scale);
    dmpt2->Scale(scale);
    dme2->Scale(scale);
    dmthe2->Scale(scale);
    dmphi2->Scale(scale);
    dme3->Scale(detscale);
    dmpz3->Scale(detscale);
    //dmpt3->Scale(detscale);
    dmthe3->Scale(detscale);
    dm_ee1->Scale(detscale);
    dm_epz1->Scale(detscale);
    dm_ept1->Scale(detscale);
    dm_ethe1->Scale(detscale);
    dmtime->Scale(detscale);
}

void DMHistograms::ScaleNeutrinos(double scale, double detscale) {
    nupz1->Scale(scale);
    nupt1->Scale(scale);
    nuthe1->Scale(scale);
    nue1->Scale(scale);
    nupz2->Scale(scale);
    nupt2->Scale(scale);
    nuthe2->Scale(scale);
    nue2->Scale(scale);
    nupz3->Scale(detscale);
    nupt3->Scale(detscale);
    nuthe3->Scale(detscale);
    nue3->Scale(detscale);
    nu_epz1->Scale(detscale);
    nu_ept1->Scale(detscale);
    nu_ee1->Scale(detscale);
    nu_ethe1->Scale(detscale);
    nutime->Scale(detscale);
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
void DMHistograms::AddScatterSigElectron(double E, double Pz, double Pt, double Th, double time) {
    dm_ee1->Fill(E);
    dm_epz1->Fill(Pz);
    dm_ept1->Fill(Pt);
    dm_ethe1->Fill(Th);
    dmtime->Fill(time);
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
void DMHistograms::AddScatterBgElectron(double E, double Pz, double Pt, double Th, double time, double W) {
    nu_ee1->Fill(E, W);
    nu_epz1->Fill(Pz, W);
    nu_ept1->Fill(Pt, W);
    nu_ethe1->Fill(Th, W);
    nutime->Fill(time, W);
    /*if (smear_sigma > 0) {
        E += Random::Gauss(smear_mean, smear_sigma / sqrt(E));

        dme3smearr->Fill(darkmatter1.E / darkmatter1E);
        dm_ee1smearr->Fill(electron1.E / electron1E);
    }*/
}

