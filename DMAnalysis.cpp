#include "DMAnalysis.h"

#include <iostream>
#include <sstream>
#include <regex>

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <TH2D.h>
#include <TF1.h>
#include <TChain.h>
#include <TGraph2D.h>
#include <TLegend.h>
#include <fstream>
#include <TApplication.h>
#include <sys/stat.h>

// Header files for the classes stored in the TTree if any.
#include "ExRootClasses.h"
#include "TClonesArray.h"

// Detector sensitivity headers
#include "Particle.h"
#include "Random.h"
#include "DUNEdet.h"
#include "DMElscattering.h"

const double emass = 0.000511;

DMAnalysis::DMAnalysis():
branch(0),
nentries(0)
{
}
DMAnalysis::~DMAnalysis()
{
}

int DMAnalysis::Process()
{
	Init();
	for(std::vector<std::string>::iterator fname = files.begin(); fname!=files.end(); ++fname) {
		TFile* file = new TFile(fname->c_str());
		TTree* tree;
		
		if(!file->IsOpen())
		{
			std::cout << "Failed to open file: " << file->GetName() << std::endl;
			return 1;
		}

		file->GetObject("LHEF",tree);
		if(!tree) {
			std::cout << "Root couldn't find LHEF tree in the input root file.\n";
			return 1;
		}
		nentries = tree->GetEntries();	

		branch = tree->GetBranch("Particle");
		if(!branch) {
			std::cout << "Root couldn't find Particle branch in the LHEF tree.\n";
			return 1;
		}

		Analyze(*fname);

		delete file;
	}
	
	UnInit();
	return 0;
}

std::regex param_match(".+?_([\\d\\.]+)_([\\d\\.]+)_([\\d\\.]+)_([\\d\\.]+).*");
int DMParameters(const std::string& filen, double& vpmass, double& chimass, double& kappa, double& alpha/*, double& xsection*/)
{
	std::smatch matches;
	if(!std::regex_match(filen, matches, param_match)) goto error;

	{
		std::stringstream parser;
		parser << matches[1];
		parser >> vpmass;
		if(parser.bad()) goto error;
		parser.clear();
		parser << matches[2];
		parser >> chimass;
		if(parser.bad()) goto error;
		parser.clear();
		parser << matches[3];
		parser >> kappa;
		if(parser.bad()) goto error;
        parser.clear();
        parser << matches[4];
        parser >> alpha;
        if(parser.bad()) goto error;
        /*parser.clear();
        parser << matches[5];
        parser >> xsection;
        if(parser.bad()) goto error;*/
		return 0;
	}
	
	error:;
	std::cout << "Tried to parse invalid filename: " << filen << std::endl;
	return 1;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

StatisticsAnalysis::StatisticsAnalysis()
{
}

DMAnalysis* StatisticsAnalysis::create()
{
	return new StatisticsAnalysis();
}

StatisticsAnalysis::~StatisticsAnalysis()
{
	delete graphchi2;
	delete graphdmee;
    delete graphsig;
    delete graphbg;
	delete neutrino_electron_e;
}

void StatisticsAnalysis::Init()
{
	index = 0;
    graphchi2 = new TGraph2D((int)files.size());
    graphchi2->SetName("chi2");
    graphchi2->SetTitle("Dark Matter #chi^{2};VP mass(GeV);#chi mass(GeV);#chi^{2}");

    graphsig = new TGraph2D((int)files.size());
    graphsig->SetName("signal");
    graphsig->SetTitle("Signal;VP mass(GeV);#chi mass(GeV);Counts");

    graphbg = new TGraph2D((int)files.size());
    graphbg->SetName("bg");
    graphbg->SetTitle("Background;VP mass(GeV);#chi mass(GeV);Counts");
	
	graphdmee = new TGraph2D((int)files.size());
	graphdmee->SetName("dmee");
	graphdmee->SetTitle("Electron Scatter E against dm Mean;VP mass(GeV);#chi mass(GeV);E(GeV)");

	DUNEDetector det;
	DMscattering scatter;

	double probMax = 1e-15;
	
    neutrino_electron_e = new TH1D("nuee3","Nu-Electron Scatter E;E (GeV)", 100, 0, 6);
	TFile* neutrinos = new TFile("data/g4lbne_nudata.root");
    TTree* neutrino_tree = (TTree*)neutrinos->Get("nudata");
    neutrino_tree->SetMakeClass(1);
    float ndxdz, ndydz, ndz;
    double Nimpwt;
    double NenergyN[5];
    double NWtNear[5];
    neutrino_tree->SetBranchAddress("Ndxdz", &ndxdz);
    neutrino_tree->SetBranchAddress("Ndydz", &ndydz);
    neutrino_tree->SetBranchAddress("Npz", &ndz);
    neutrino_tree->SetBranchAddress("NenergyN[5]", &NenergyN);
    neutrino_tree->SetBranchAddress("Nimpwt", &Nimpwt);
    neutrino_tree->SetBranchAddress("NWtNear[5]", &NWtNear);

    long long neutrino_entries = neutrino_tree->GetEntries();
    int neutrino_intersectcount = 0;
    int neutrino_scattercount = 0;
    int neutrino_total = 5000000;
    double total_weight = 0;
    std::cout << "Generating neutrino distribution...\n";
    for(long long i = 0, j = 0; j < neutrino_total; i=(i+1)%neutrino_entries, ++j) {
        neutrino_tree->GetEvent(i);
        Particle neutrino(0);
        Particle electron(emass);
        double weight = Nimpwt*NWtNear[0];
        total_weight += weight;

        neutrino.FourMomentum(ndxdz*ndz, ndydz*ndz, ndz, NenergyN[0]);

        //nupz_1->Fill(neutrino.pz, weight);
        //nue_1->Fill(neutrino.E, weight);

        double pLen, dx, dy, dz, tmin, tmax;
        neutrino.getNorm(pLen, dx, dy, dz);
        if(det.intersect(dx, dy, dz, tmin, tmax)) {

            //nupz_2->Fill(neutrino.pz, weight);
            //nue_2->Fill(neutrino.E, weight);
            if(scatter.probscatterNeutrino(probMax, neutrino, tmax-tmin)) {
                scatter.scattereventNeutrino(neutrino, electron);

                //nupz_3->Fill(neutrino.pz, weight);
                //nue_3->Fill(neutrino.E, weight);

                //nuepz_3->Fill(electron.pz, weight);
                neutrino_electron_e->Fill(electron.E, weight);

                //double pnorm = sqrt(electron.px*electron.px+electron.py*electron.py+electron.pz*electron.pz);
                //nuetheta_3->Fill(acos(std::min<double>(std::max<double>(electron.pz/pnorm, -1), 1));
            }
        }
    }
    double scale = (double)100000/total_weight;
    neutrino_electron_e->Scale(scale);
    
    delete neutrinos;
}

void StatisticsAnalysis::Analyze(const std::string& filen)
{
    //if(index>10) return;
	double vpmass, chimass, kappa, alpha;
	if(DMParameters(filen, vpmass, chimass, kappa, alpha)) return;
	std::cout << filen << std::endl;

	DUNEDetector det;
	DMscattering scatter;
	double probMax = 1e-15;
	
	TH1D* dm_electron_e = new TH1D("dme","DM-Electron Scatter E;E (GeV)", 100, 0, 6);
	TClonesArray* array = new TClonesArray("TRootLHEFParticle", 5);
	branch->SetAddress(&array);
		
	for(Int_t i = 0; i < nentries; ++i)
	{
		branch->GetEntry(i);
		for(Int_t j = 0; j < array->GetEntries(); ++j)
		{
			TRootLHEFParticle* particle = (TRootLHEFParticle*)array->At(j);
			if(particle->PID==33)
			{
				Particle darkmatter1(chimass);
				Particle electron1(emass);
				darkmatter1.FourMomentum(particle->Px,particle->Py,-particle->Pz,particle->E);

                double pLen, dx, dy, dz, tmin, tmax;
                darkmatter1.getNorm(pLen, dx, dy, dz);
				if(det.intersect(dx, dy, dz, tmin, tmax)) {
                    if(scatter.probscatter(probMax, vpmass, chimass, kappa, alpha, darkmatter1, tmax-tmin)) {
                        scatter.scatterevent(vpmass, chimass, kappa, alpha, darkmatter1,
                                             electron1);
                            dm_electron_e->Fill(electron1.E);
                    }
                }
			}
		}
	}
	
	double chi2 = 0;
	for(int i = 1; i<=dm_electron_e->GetNbinsX(); ++i)
	{
		double N_bg = neutrino_electron_e->GetBinContent(i);
		double N_tot = dm_electron_e->GetBinContent(i) + N_bg;
		if(N_bg == 0) continue;
		chi2 += 2*(N_bg*log(N_bg/N_tot) + N_tot - N_bg);
	}

    graphsig->SetPoint(index, vpmass, chimass, dm_electron_e->Integral());
	graphchi2->SetPoint(index, vpmass, chimass, chi2);
	graphdmee->SetPoint(index, vpmass, chimass, dm_electron_e->GetMean(1));
	++index;
	
	delete array;
	delete dm_electron_e;
}

void StatisticsAnalysis::UnInit()
{
	/*mkdir("images", S_IRWXU | S_IRWXG | S_IRWXO);
	std::string params(argv[2]);
	std::string outputfilen(argv[3]);*/

	TCanvas* canvas = new TCanvas("output", "", 1920,1080);
	canvas->SetFillColor(33);
	canvas->SetFrameFillColor(17);
	canvas->SetGrid();
	
	canvas->SetTheta(-90);
	canvas->SetPhi(-0.0001);

    graphchi2->Draw("surf1 Z");
    canvas->Update();
    canvas->SaveAs("chi2.png");

    graphsig->Draw("surf1 Z");
    canvas->Update();
    canvas->SaveAs("sig.png");

	graphdmee->Draw("surf1 Z");
	canvas->Update();
	canvas->SaveAs("dme.png");
	
	neutrino_electron_e->Draw();
	canvas->Update();
	canvas->SaveAs("nue.png");

    /*canvas->Draw();
    gApp->Run();*/

	delete canvas;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DetectorAnalysis::DetectorAnalysis()
{
}

DMAnalysis* DetectorAnalysis::create()
{
	return new DetectorAnalysis();
}

DetectorAnalysis::~DetectorAnalysis()
{
}

void DetectorAnalysis::Init()
{
	//if(gOptions[OPT_DETECTOR])
	//	detectorType = gOptions[OPT_DETECTOR].last()->arg;
	//else
    	detectorType = "DUNE";

	if(detectorType=="DUNE")
	{
		smear_mean = 0;
		smear_sigma = 0.06;
	}

	/*if(gOptions[OPT_DET_SMEAR_SIG])
		smear_sigma = toDouble(gOptions[OPT_DET_SMEAR_SIG].last()->arg);
	else
		smear_sigma = 0;
	
	if(gOptions[OPT_DET_SMEAR_MEAN])
		smear_mean = toDouble(gOptions[OPT_DET_SMEAR_MEAN].last()->arg);
	else
		smear_mean = 0;*/
	
	smear_sigma = 0.06;

    DUNEDetector det;
    DMscattering scatter;

    double probMax = 1e-15;

    output = new TFile("output.root", "RECREATE");
    dmpz_1 = new TH1D("dmpz1","Dark matter P_{z};P_{z} (GeV/c)", 100, 0, 60);
    dmpt_1 = new TH1D("dmpt1","Dark matter P_{t};P_{t} (GeV/c)", 100, 0, 0);
    dme_1 = new TH1D("dme1","Dark matter E;E (GeV)", 100, 0, 60);
    dme_v = new TH1D("dmev","Dark matter velocity;Vel (m/s)", 100, 2.5e8, 3e8);
    dme_time = new TH1D("dmetime","Dark matter arrival time;Time (s)", 100, 1.9e-6, 2.5e-6);
    dme_t1 = new TH1D("dmet1","Dark matter \\theta;\\theta", 100, 0, M_PI);
    dmpz_2 = new TH1D("dmpz2","Dark matter in Detector P_{z};P_{z} (GeV/c)", 100, 0, 60);
    dmpt_2 = new TH1D("dmpt2","Dark matter in Detector P_{t};P_{t} (GeV/c)", 100, 0, 0);
    dme_2 = new TH1D("dme2","Dark matter in Detector E;E (GeV)", 100, 0, 60);
    dme_t2 = new TH1D("dmet2","Dark matter in Detector \\theta;\\theta", 100, 0, 0);
    //dmpz_3 = new TH1D("dmpz3","Dark matter Scatter P_{z};P_{z} (GeV/c)", 100, 0, 0);
    //dme_3 = new TH1D("dme3","Dark matter Scatter E;E (GeV)", 100, 0, 0);
    //dme_4 = new TH1D("dme4","Dark matter Smeared Scatter E;E (GeV)", 100, 0, 0);
    //dme_5 = new TH1D("dme5","Dark matter Scatter Energy Smearing Ratio;E (GeV)", 100, 0, 0);

    nupz_1 = new TH1D("nupz1","Neutrino P_{z};P_{z} (GeV/c)", 100, 0, 60);
    nupt_1 = new TH1D("nupt1","Neutrino P_{t};P_{t} (GeV/c)", 100, 0, 0);
    nue_1 = new TH1D("nue1","Neutrino E;E (GeV)", 100, 0, 60);
    nue_t1 = new TH1D("nuet1","Neutrino \\theta;\\theta", 100, 0, M_PI);
    nupz_2 = new TH1D("nupz2","Neutrino in Detector P_{z};P_{z} (GeV/c)", 100, 0, 60);
    nupt_2 = new TH1D("nupt2","Neutrino in Detector P_{t};P_{t} (GeV/c)", 100, 0, 0);
    nue_2 = new TH1D("nue2","Neutrino in Detector E;E (GeV)", 100, 0, 60);
    nue_t2 = new TH1D("nuet2","Neutrino in Detector \\theta;\\theta", 100, 0, 0);
    //nupz_3 = new TH1D("nupz3","Neutrino Scatter P_{z};P_{z} (GeV/c)", 100, 0, 0);
    //nue_3 = new TH1D("nue3","Neutrino Scatter E;E (GeV)", 100, 0, 0);
    //nue_4 = new TH1D("nue4","Neutrino Smeared Scatter E;E (GeV)", 100, 0, 0);
    //nue_5 = new TH1D("nue5","Neutrino Scatter Energy Smearing Ratio;E (GeV)", 100, 0, 0);

    epz_3 = new TH1D("dmepz3","DM-Electron Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
    ept_3 = new TH1D("dmept3","DM-Electron Scatter P_{t};P_{t} (GeV/c)", 100, 0, 6);
    etheta_3 = new TH1D("dmetheta","DM-Electron Scatter \\theta;\\theta)", 100, 0, 3.2);
    ee_3 = new TH1D("dmee3","DM-Electron Scatter E;E (GeV)", 100, 0, 6);
    ee_4 = new TH1D("dmee4","DM-Electron Smeared Scatter E;E (GeV)", 100, 0, 6);
    ee_5 = new TH1D("dmee5","DM-Electron Scatter Energy Smearing Ratio;E (GeV)", 40, 0, 6);

    nuepz_3 = new TH1D("nuepz3","Nu-Electron Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
    nuept_3 = new TH1D("nuept3","Nu-Electron Scatter P_{z};P_{z} (GeV/c)", 100, 0, 6);
    nuetheta_3 = new TH1D("nuetheta","Nu-Electron Scatter \\theta;\\theta)", 100, 0, 3.2);
    nuee_3 = new TH1D("nuee3","Nu-Electron Scatter E;E (GeV)", 100, 0, 6);
    nuee_4 = new TH1D("nuee4","Nu-Electron Smeared Scatter E;E (GeV)", 100, 0, 6);
    nuee_5 = new TH1D("nuee5","Nu-Electron Scatter Energy Smearing Ratio;E (GeV)", 40, 0, 6);

    /*TFile* neutrinos = new TFile("data/g4lbne_nudata.root");
    TTree* neutrino_tree = (TTree*)neutrinos->Get("nudata");
    neutrino_tree->SetMakeClass(1);
    float ndxdz, ndydz, ndz, ne;
    double Nimpwt;
    double NWtNear[5];
    neutrino_tree->SetBranchAddress("Ndxdz", &ndxdz);
    neutrino_tree->SetBranchAddress("Ndydz", &ndydz);
    neutrino_tree->SetBranchAddress("Npz", &ndz);
    neutrino_tree->SetBranchAddress("Nenergy", &ne);
    neutrino_tree->SetBranchAddress("Nimpwt", &Nimpwt);
    neutrino_tree->SetBranchAddress("NWtNear[5]", &NWtNear);

    long long neutrino_entries = neutrino_tree->GetEntries();
    int neutrino_intersectcount = 0;
    int neutrino_scattercount = 0;
    double total_weight = 0;
    for(long long i = 0, j = 0; j < 10000000; i=(i+1)%neutrino_entries, ++j) {
        neutrino_tree->GetEvent(i);
        Particle neutrino(0);
        Particle electron(emass);
        double weight = Nimpwt*NWtNear[0];
        total_weight += weight;

        neutrino.FourMomentum(ndxdz*ndz, ndydz*ndz, ndz, ne);
        double pLen, dx, dy, dz;
        neutrino.getNorm(pLen, dx, dy, dz);
        double pt = sqrt(neutrino.px*neutrino.px+neutrino.py*neutrino.py);
        double theta = acos(std::min<double>(std::max<double>(dz, -1), 1));

        nupz_1->Fill(neutrino.pz, weight);
        nupt_1->Fill(pt, weight);
        nue_1->Fill(neutrino.E, weight);
        nue_t1->Fill(theta, weight);

        double tmin = 0, tmax = 0;
        if(det.intersect(dx, dy, dz, tmin, tmax)) {
            nupz_2->Fill(neutrino.pz, weight);
            nupt_2->Fill(pt, weight);
            nue_2->Fill(neutrino.E, weight);
            nue_t2->Fill(theta, weight);

            if(scatter.probscatterNeutrino(probMax, neutrino, tmax - tmin)) {
                scatter.scattereventNeutrino(neutrino, electron);

                //nupz_3->Fill(neutrino.pz, weight);
                //nue_3->Fill(neutrino.E, weight);

                nuepz_3->Fill(electron.pz, weight);
                nuept_3->Fill(sqrt(electron.px*electron.px+electron.py*electron.py), weight);
                nuee_3->Fill(electron.E, weight);

                nuetheta_3->Fill(acos(std::min<double>(std::max<double>(electron.pz / sqrt(electron.px * electron.px +
                                                                                           electron.py * electron.py +
                                                                                           electron.pz * electron.pz),
                                                                        -1), 1)));
            }
        }
    }
    double scale = (double)10000000/total_weight;
    nupz_1->Scale(scale);
    nue_1->Scale(scale);
    nupz_2->Scale(scale);
    nue_2->Scale(scale);
    nuepz_3->Scale(scale);
    nuept_3->Scale(scale);
    nuee_3->Scale(scale);
    */
}

void saveComparison(TCanvas* canvas, const char* savename, const char* canvastitle, TH1D* hist1, TH1D* hist2, const char* histn1, const char* histn2)
{
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

void saveHistogram(TCanvas* canvas, TH1D* hist)
{
    std::stringstream filen;
    filen << "histograms/" << hist->GetName() << ".png";
    hist->Draw();
    canvas->Update();
    canvas->SaveAs(filen.str().c_str());
}

void DetectorAnalysis::Analyze(const std::string& filen)
{
	double vpmass, chimass, kappa, alpha;
	if(DMParameters(filen, vpmass, chimass, kappa, alpha)) return;

	DUNEDetector det;
	DMscattering scatter;
	double probMax = 1e-10;

	TClonesArray* array = new TClonesArray("TRootLHEFParticle", 5);
	branch->SetAddress(&array);

    //Particle darkphoton(vpmass);
    Particle darkmatter1(chimass);
    //Particle darkmatter2(chimass);
    Particle electron1(emass);
    //Particle electron2(emass);

    for(Int_t i = 0; i < nentries; ++i)
    {
		branch->GetEntry(i);
		for(Int_t j = 0; j < array->GetEntries(); ++j)
		{
			TRootLHEFParticle* particle = (TRootLHEFParticle*)array->At(j);
			if(particle->PID==33)
			{
                darkmatter1.FourMomentum(particle->Px,particle->Py,-particle->Pz,particle->E);
                double pLen, dx, dy, dz;
                darkmatter1.getNorm(pLen, dx, dy, dz);
                double theta = acos(std::min<double>(std::max<double>(dz, -1), 1));
                double pt = sqrt(darkmatter1.px*darkmatter1.px+darkmatter1.py*darkmatter1.py);

				dmpz_1->Fill(-particle->Pz);
                dmpt_1->Fill(pt);
				dme_1->Fill(particle->E);
                dme_t1->Fill(theta);

                double tmin = 0, tmax = 0;
				if(det.intersect(dx, dy, dz, tmin, tmax))
				{
					dmpz_2->Fill(-particle->Pz);
                    dmpt_2->Fill(pt);
					dme_2->Fill(particle->E);
                    dme_t2->Fill(theta);
                    double vel = pLen/particle->E*2.998e8;
                    dme_v->Fill(vel);
                    dme_time->Fill(tmin/vel);

				    if(scatter.probscatter(probMax,vpmass,chimass,kappa,alpha,darkmatter1,tmax-tmin)) {
                        scatter.scatterevent(vpmass, chimass, kappa, alpha, darkmatter1,
                                             electron1);

                        //dmpz_3->Fill(darkmatter1.pz);
                        //dme_3->Fill(darkmatter1.E);

                        //bool thresh = electron1.E < 0.5;
                        epz_3->Fill(electron1.pz);
                        ept_3->Fill(sqrt(electron1.px*electron1.px+electron1.py*electron1.py));
                        /*if(thresh)*/ ee_3->Fill(electron1.E);
                        double pnorm = sqrt(electron1.px * electron1.px + electron1.py * electron1.py +
                                            electron1.pz * electron1.pz);
                        etheta_3->Fill(acos(std::min<double>(std::max<double>(electron1.pz / pnorm, -1), 1)));

                        double darkmatter1E = darkmatter1.E;
                        double electron1E = electron1.E;

                        if (smear_sigma > 0) {
                            darkmatter1.E += Random::Gauss(smear_mean, smear_sigma / sqrt(darkmatter1.E));
                            electron1.E += Random::Gauss(smear_mean, smear_sigma / sqrt(electron1.E));

                            //dme_4->Fill(darkmatter1.E);
                            /*if(thresh)*/ ee_4->Fill(electron1.E);

                            //dme_5->Fill(darkmatter1.E / darkmatter1E);
                            /*if(thresh)*/ ee_5->Fill(electron1.E / electron1E);
                        }
                    }
				}
			}
		}
	}
	delete array;
}

void DetectorAnalysis::UnInit()
{

    dmpz_1->BufferEmpty();
    dmpt_1->BufferEmpty();
    dme_1->BufferEmpty();
    dmpz_2->BufferEmpty();
    dmpt_2->BufferEmpty();
    dme_2->BufferEmpty();
    dme_v->BufferEmpty();
    dme_time->BufferEmpty();
    dme_t1->BufferEmpty();
    dme_t2->BufferEmpty();
    //dmpz_3->BufferEmpty();
    //dme_3->BufferEmpty();
    //dme_4->BufferEmpty();
    //dme_5->BufferEmpty();
    nupz_1->BufferEmpty();
    nue_1->BufferEmpty();
    nupz_2->BufferEmpty();
    nue_2->BufferEmpty();
    nue_t1->BufferEmpty();
    nue_t2->BufferEmpty();
    //nupz_3->BufferEmpty();
    //nue_3->BufferEmpty();
    //nue_4->BufferEmpty();
    //nue_5->BufferEmpty();

    epz_3->BufferEmpty();
    ept_3->BufferEmpty();
    etheta_3->BufferEmpty();
    ee_3->BufferEmpty();
    ee_4->BufferEmpty();
    ee_5->BufferEmpty();
    nuepz_3->BufferEmpty();
    nuept_3->BufferEmpty();
    nuetheta_3->BufferEmpty();
    nuee_3->BufferEmpty();
    nuee_4->BufferEmpty();
    nuee_5->BufferEmpty();

    double pare[3] = {1,1,1};
    //TF1 *ge = new TF1("ge","gaus(0)",0.8,1.2);
    TF1 *ge = new TF1("ge","gaus(0)",ee_5->GetXaxis()->GetXmin(),ee_5->GetXaxis()->GetXmax());
    ge->SetParameters(pare);
    ee_5->Fit(ge,"R");

    output->Write();

    TCanvas* canvas = new TCanvas("c1","canvas",1024,1024);
    mkdir("histograms", S_IRWXU | S_IRWXG | S_IRWXO);
    saveHistogram(canvas, dmpz_1);
    saveHistogram(canvas, dmpt_1);
    saveHistogram(canvas, dme_1);
    saveHistogram(canvas, dmpz_2);
    saveHistogram(canvas, dmpt_2);
    saveHistogram(canvas, dme_2);
    saveHistogram(canvas, dme_v);
    saveHistogram(canvas, dme_time);
    saveHistogram(canvas, dme_t1);
    saveHistogram(canvas, dme_t2);
    //saveHistogram(canvas, dmpz_3);
    //saveHistogram(canvas, dme_3);
    //saveHistogram(canvas, dme_4);
    //saveHistogram(canvas, dme_5);
    saveHistogram(canvas, nupz_1);
    saveHistogram(canvas, nupt_1);
    saveHistogram(canvas, nue_1);
    saveHistogram(canvas, nupz_2);
    saveHistogram(canvas, nupt_2);
    saveHistogram(canvas, nue_2);
    saveHistogram(canvas, nue_t1);
    saveHistogram(canvas, nue_t2);
    //saveHistogram(canvas, nupz_3);
    //saveHistogram(canvas, nue_3);
    //saveHistogram(canvas, nue_4);
    //saveHistogram(canvas, nue_5);
    saveHistogram(canvas, epz_3);
    saveHistogram(canvas, ept_3);
    saveHistogram(canvas, etheta_3);
    saveHistogram(canvas, ee_3);
    saveHistogram(canvas, ee_4);
    saveHistogram(canvas, ee_5);
    saveHistogram(canvas, nuepz_3);
    saveHistogram(canvas, nuept_3);
    saveHistogram(canvas, nuetheta_3);
    saveHistogram(canvas, nuee_3);
    saveHistogram(canvas, nuee_4);
    saveHistogram(canvas, nuee_5);

    saveComparison(canvas, "nu_dmpz1", "DM vs Neutrino P_{z}", dmpz_1, nupz_1, "Dark matter P_{z}", "Neutrino P_{z}");
    saveComparison(canvas, "nu_dme1", "DM vs Neutrino E", dme_1, nue_1, "Dark matter E", "Neutrino E");
    saveComparison(canvas, "nu_dmet1", "DM vs Neutrino \\theta", dme_t1, nue_t1, "Dark matter \\theta", "Neutrino \\theta");
    saveComparison(canvas, "nu_dmet2", "DM vs Neutrino Intersections \\theta", dme_t2, nue_t2, "Dark matter \\theta", "Neutrino \\theta");
    saveComparison(canvas, "nu_dmet2", "DM vs Neutrino Intersections \\theta", dme_t2, nue_t2, "Dark matter \\theta", "Neutrino \\theta");
    saveComparison(canvas, "nu_dmpz2", "DM vs Neutrino Intersections P_{z}", dmpz_2, nupz_2, "Dark matter P_{z}", "Neutrino P_{z}");
    saveComparison(canvas, "nu_dme2", "DM vs Neutrino Intersections E", dme_2, nue_2, "Dark matter E", "Neutrino E");
    saveComparison(canvas, "nu_dm_epz", "DM vs Neutrino Electron Scatter P_{z}", epz_3, nuepz_3, "Electron - Dark matter P_{z}", "Electron - Neutrino P_{z}");
    saveComparison(canvas, "nu_dm_ept", "DM vs Neutrino Electron Scatter P_{t}", ept_3, nuept_3, "Electron - Dark matter P_{t}", "Electron - Neutrino P_{t}");
    saveComparison(canvas, "nu_dm_etheta", "DM vs Neutrino Electron Scatter \\theta", etheta_3, nuetheta_3, "Electron - Dark matter \\theta", "Electron - Neutrino \\theta");
    saveComparison(canvas, "nu_dm_ee", "DM vs Neutrino Electron Scatter E", ee_3, nuee_3, "Electron - Dark matter E", "Electron - Neutrino E");

    delete dmpz_1;
    delete dmpt_1;
    delete dme_1;
    delete dmpz_2;
    delete dmpt_2;
    delete dme_2;
    delete dme_v;
    delete dme_time;
    delete dme_t1;
    delete dme_t2;
    //delete dmpz_3;
    //delete dme_3;
    //delete dme_4;
    //delete dme_5;
    delete nupz_1;
    delete nupt_1;
    delete nue_1;
    delete nupz_2;
    delete nupt_2;
    delete nue_2;
    delete nue_t1;
    delete nue_t2;
    //delete nupz_3;
    //delete nue_3;
    //delete nue_4;
    //delete nue_5;
    delete epz_3;
    delete ept_3;
    delete etheta_3;
    delete ee_3;
    delete ee_4;
    delete nuepz_3;
    delete nuept_3;
    delete nuetheta_3;
    delete nuee_3;
    delete nuee_4;
    delete nuee_5;
    delete canvas;
    delete output;
}
