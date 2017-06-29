#include "DMAnalysis.h"

#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <regex>

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TF1.h>
#include <TChain.h>
#include <TGraph2D.h>
#include <TLegend.h>
#include <fstream>

// Header files for the classes stored in the TTree if any.
#include "ExRootClasses.h"
#include "TClonesArray.h"
#include "TObject.h"

// Detector sensitivity headers
#include "Particle.h"
#include "Random.h"
#include "Kinematics.h"
#include "DUNEdet.h"
#include "DMElscattering.h"



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

std::regex param_match(".+?_([\\d\\.]+)_([\\d\\.]+)_([\\d\\.]+)_([\\d\\.]+)\\.root");
int DMAnalysis::DMParameters(const std::string& filen, double& vpmass, double& chimass, double& kappa, double& alpha)
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
	delete graph;
}

void StatisticsAnalysis::Init()
{
	index = 0;
	graph = new TGraph2D(files.size());
	graph->SetTitle("#mu_{#chi pz};VP mass(GeV);#chi mass(GeV);#mu_{pz}(GeV/c^{2})");
}

void StatisticsAnalysis::Analyze(const std::string& filen)
{
	double vpmass, chimass, kappa, alpha;
	if(DMParameters(filen, vpmass, chimass, kappa, alpha)) return;
	
	double z = 0;
	
	TClonesArray* array = new TClonesArray("TRootLHEFParticle", 5);
	branch->SetAddress(&array);
		
	for(Int_t i = 0; i < nentries; ++i)
	{
		branch->GetEntry(i);
		for(Int_t j = 0; j < array->GetEntries(); ++j)
		{
			TRootLHEFParticle* particle = (TRootLHEFParticle*)array->At(j);
			switch(particle->PID)
			{
				case 33:
					z += particle->Pz;
				break;
				case -33:
				break;
			}
		}
	}
	
	z /= (double)nentries;
	graph->SetPoint(index, vpmass, chimass, z);
	++index;
	
	delete array;
}

void StatisticsAnalysis::UnInit()
{
	/*mkdir("images", S_IRWXU | S_IRWXG | S_IRWXO);
	std::string params(argv[2]);
	std::string outputfilen(argv[3]);*/

	TCanvas* canvas = new TCanvas("output", "", 1600,900);
	canvas->SetFillColor(33);
	canvas->SetFrameFillColor(17);
	canvas->SetGrid();
	
	//canvas->SetTheta(90);
	//canvas->SetPhi(0);
	graph->Draw("surf1");
	
	canvas->Update();
	
	std::stringstream filename;
	//filename << "images/" << name << "_" << outputfilen << ".png";
	filename << "graph.png";
	
	canvas->SaveAs(filename.str().c_str());

	delete canvas;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

DarkMatterDistribution::DarkMatterDistribution()
{
}

DMAnalysis* DarkMatterDistribution::create()
{
	return new DarkMatterDistribution();
}

DarkMatterDistribution::~DarkMatterDistribution()
{
	delete histo;
	delete output;
}

void DarkMatterDistribution::Init()
{
	std::string particle, attr;
	//if(gOptions[OPT_PARTICLE])
	//	particle = std::string(gOptions[OPT_PARTICLE].last()->arg);
	//else
		particle = "33";
	//if(gOptions[OPT_PARTICLEATTRIBUTE])
    //	attr = std::string(gOptions[OPT_PARTICLEATTRIBUTE].last()->arg);
	//else
		attr = "px";
	
	output = new TFile("out.root","RECREATE");
	
	std::stringstream parser;
	parser << particle;
	parser >> pdgCode;
	if(parser.bad()){std::cout << "Error parsing requested particle pdgCode: " << particle << std::endl;}
	
	std::stringstream title;
	switch(pdgCode)
	{
		case 33:
		title << "#chi ";
	}
	
	if(attr == "px")
	{
		attribute = &TRootLHEFParticle::Px;
		title << "p_{x};p_{x}(GeV/c^{2})";
	}
	else if(attr == "py")
	{
		attribute = &TRootLHEFParticle::Py;
		title << "p_{y};p_{y}(GeV/c^{2})";
	}
	else if(attr == "pz")
	{
		attribute = &TRootLHEFParticle::Pz;
		title << "p_{z};p_{z}(GeV/c^{2})";
	}
	else if(attr == "e")
	{
		attribute = &TRootLHEFParticle::E;
		title << "Energy;E(GeV)";
	}
	else if(attr == "m")
	{
		attribute = &TRootLHEFParticle::M;
		title << "Mass;M(GeV)";
	}
	else
	{
		attribute = 0;
		std::cout << "Invalid attribute: " << attr << std::endl;
	}
	
	histo = new TH1D("Particle",title.str().c_str(),100,0,0);
}

void DarkMatterDistribution::Analyze(const std::string& filen)
{
	if(!attribute){std::cout << "Tried to fill histogram with invalid attribute.\n";}
	TClonesArray* array = new TClonesArray("TRootLHEFParticle", 5);
	branch->SetAddress(&array);
		
	for(Int_t i = 0; i < nentries; ++i)
	{
		branch->GetEntry(i);
		for(Int_t j = 0; j < array->GetEntries(); ++j)
		{
			TRootLHEFParticle* particle = (TRootLHEFParticle*)array->At(j);
			if(particle->PID==pdgCode)
			{
				histo->Fill(particle->*attribute);
			}
		}
	}
	histo->BufferEmpty();
	delete array;
}

void DarkMatterDistribution::UnInit()
{
	output->Write(); 
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
}

void DetectorAnalysis::Analyze(const std::string& filen)
{
	double vpmass, chimass, kappa, alpha;
	if(DMParameters(filen, vpmass, chimass, kappa, alpha)) return;
	
	const double emass = 0.000511;

    Kinematics kin;
	DUNEDetector det;
	DMscattering scatter;
	
	int scatterCount = 0;
	int Nscatter;
	int Nelectron;
	double probMax = 10e-15;
	
	TFile* output = new TFile("output.root", "RECREATE");
    TH1D* dmpz_1 = new TH1D("dmpz1","Dark matter P_{z};P_{z} (GeV/c^{2})", 100, 0, 0);
    TH1D* dme_1 = new TH1D("dme1","Dark matter E;E (GeV)", 100, 0, 0);
    TH1D* dmpz_2 = new TH1D("dmpz2","Intersected Dark matter P_{z};P_{z} (GeV/c^{2})", 100, 0, 0);
    TH1D* dme_2 = new TH1D("dme2","Intersected Dark matter E;E (GeV)", 100, 0, 0);
    TH1D* dmpz_3 = new TH1D("dmpz3","Dark matter Scatter P_{z};P_{z} (GeV/c^{2})", 100, 0, 0);
    TH1D* dme_3 = new TH1D("dme3","Dark matter Scatter E;E (GeV)", 100, 0, 0);
    TH1D* dme_4 = new TH1D("dme4","Dark matter Smeared Scatter E;E (GeV)", 100, 0, 0);
    TH1D* dme_5 = new TH1D("dme5","Dark matter Scatter Energy Smearing Ratio;E (GeV)", 100, 0, 0);

    TH1D* nupz_1 = new TH1D("nupz1","Neutrino P_{z};P_{z} (GeV/c^{2})", 100, 0, 0);
    TH1D* nue_1 = new TH1D("nue1","Neutrino E;E (GeV)", 100, 0, 0);
    TH1D* nupz_2 = new TH1D("nupz2","Intersected Neutrino P_{z};P_{z} (GeV/c^{2})", 100, 0, 0);
    TH1D* nue_2 = new TH1D("nue2","Intersected Neutrino E;E (GeV)", 100, 0, 0);
    TH1D* nupz_3 = new TH1D("nupz3","Neutrino Scatter P_{z};P_{z} (GeV/c^{2})", 100, 0, 0);
    TH1D* nue_3 = new TH1D("nue3","Neutrino Scatter E;E (GeV)", 100, 0, 0);
    TH1D* nue_4 = new TH1D("nue4","Neutrino Smeared Scatter E;E (GeV)", 100, 0, 0);
    TH1D* nue_5 = new TH1D("nue5","Neutrino Scatter Energy Smearing Ratio;E (GeV)", 100, 0, 0);

    TH1D* epz_3 = new TH1D("dmepz3","Electron Scatter P_{z};P_{z} (GeV/c^{2})", 100, 0, 0);
    TH1D* ee_3 = new TH1D("dmee3","Electron Scatter E;E (GeV)", 100, 0, 0);
    TH1D* ee_4 = new TH1D("dmee4","Electron Smeared Scatter E;E (GeV)", 100, 0, 0);
    TH1D* ee_5 = new TH1D("dmee5","Electron Scatter Energy Smearing Ratio;E (GeV)", 40, 0, 2);

    TH1D* nuepz_3 = new TH1D("nuepz3","Electron Scatter P_{z};P_{z} (GeV/c^{2})", 100, 0, 0);
    TH1D* nuee_3 = new TH1D("nuee3","Electron Scatter E;E (GeV)", 100, 0, 0);
    TH1D* nuee_4 = new TH1D("nuee4","Electron Smeared Scatter E;E (GeV)", 100, 0, 0);
    TH1D* nuee_5 = new TH1D("nuee5","Electron Scatter Energy Smearing Ratio;E (GeV)", 40, 0, 2);

    TFile* neutrinos = new TFile("g4lbne_nudata.root");
    TTree* neutrino_tree = (TTree*)neutrinos->Get("nudata");
    neutrino_tree->SetMakeClass(1);
    float ndxdz, ndydz, ndz, ne;
    neutrino_tree->SetBranchAddress("Ndxdz", &ndxdz);
    neutrino_tree->SetBranchAddress("Ndydz", &ndydz);
    neutrino_tree->SetBranchAddress("Npz", &ndz);
    neutrino_tree->SetBranchAddress("Nenergy", &ne);

    long long neutrino_entries = neutrino_tree->GetEntries();
    int neutrino_intersectcount = 0;
    int neutrino_scattercount = 0;
    for(Int_t i = 0; i < neutrino_entries; ++i) {
        neutrino_tree->GetEvent(i);
        Particle neutrino(0);
        Particle electron(emass);

        neutrino.FourMomentum(ndxdz*ndz, ndydz*ndz, ndz, ne);

        nupz_1->Fill(neutrino.pz);
        nue_1->Fill(neutrino.E);

        int Switch = 0;
        det.intersect(Switch,neutrino_intersectcount,neutrino);
        if(Switch == 1) {
            nupz_2->Fill(neutrino.pz);
            nue_2->Fill(neutrino.E);
        }
        //scatter.probscatterNeutrino(Switch,neutrino_scattercount,probMax,neutrino);
        Switch = 2;
        scatter.scattereventNeutrino(Switch,neutrino_scattercount,neutrino,electron);
        if(Switch == 2) {
            nupz_3->Fill(neutrino.pz);
            nue_3->Fill(neutrino.E);

            nuepz_3->Fill(electron.pz);
            nuee_3->Fill(electron.E);
        }
    }
    std::cout << "Neutrino intersections: (" << neutrino_intersectcount << " / " << neutrino_entries << ")\n";
    std::cout << "Neutrino scatters: (" << neutrino_scattercount << " / " << neutrino_intersectcount << ")\n";


	TClonesArray* array = new TClonesArray("TRootLHEFParticle", 5);
	branch->SetAddress(&array);

    Particle darkphoton(vpmass);
    Particle darkmatter1(chimass);
    Particle darkmatter2(chimass);
    Particle electron1(emass);
    Particle electron2(emass);

    for(Int_t i = 0; i < nentries; ++i)
    {
		branch->GetEntry(i);
		for(Int_t j = 0; j < array->GetEntries(); ++j)
		{
			TRootLHEFParticle* particle = (TRootLHEFParticle*)array->At(j);
			if(particle->PID==33)
			{
				dmpz_1->Fill(-particle->Pz);
				dme_1->Fill(particle->E);
				
				darkmatter1.FourMomentum(particle->Px,particle->Py,-particle->Pz,particle->E);
				
				int DMSwitch = 0;
				det.intersect(DMSwitch,scatterCount,darkmatter1);
				if(DMSwitch == 1)
				{
					dmpz_2->Fill(-particle->Pz);
					dme_2->Fill(particle->E);
				}
				scatter.probscatter(DMSwitch,Nscatter,probMax,vpmass,chimass,kappa,alpha,darkmatter1);
				scatter.scatterevent(DMSwitch,Nelectron,vpmass,chimass,kappa,alpha,darkmatter1,electron1);
				if(DMSwitch == 2)
				{
					dmpz_3->Fill(darkmatter1.pz);
					dme_3->Fill(darkmatter1.E);
					
					bool thresh = electron1.E < 0.5;
					epz_3->Fill(electron1.pz);
					if(thresh) ee_3->Fill(electron1.E);
					
					double darkmatter1E = darkmatter1.E;
					double electron1E = electron1.E;
					
					if(smear_sigma>0)
					{
						darkmatter1.E += Random::Gauss(smear_mean, smear_sigma/sqrt(darkmatter1.E));
						electron1.E += Random::Gauss(smear_mean, smear_sigma/sqrt(electron1.E));
						
						dme_4->Fill(darkmatter1.E);
						if(thresh) ee_4->Fill(electron1.E);
						
						dme_5->Fill(darkmatter1.E / darkmatter1E);
						if(thresh) ee_5->Fill(electron1.E / electron1E);
					}
				}
			}
		}
	}
    dmpz_1->BufferEmpty();
    dme_1->BufferEmpty();
    dmpz_2->BufferEmpty();
    dme_2->BufferEmpty();
    dmpz_3->BufferEmpty();
    dme_3->BufferEmpty();
    dme_4->BufferEmpty();
    dme_5->BufferEmpty();
    nupz_1->BufferEmpty();
    nue_1->BufferEmpty();
    nupz_2->BufferEmpty();
    nue_2->BufferEmpty();
    nupz_3->BufferEmpty();
    nue_3->BufferEmpty();
    nue_4->BufferEmpty();
    nue_5->BufferEmpty();

    epz_3->BufferEmpty();
    ee_3->BufferEmpty();
    ee_4->BufferEmpty();
    ee_5->BufferEmpty();
    nuepz_3->BufferEmpty();
    nuee_3->BufferEmpty();
    nuee_4->BufferEmpty();
    nuee_5->BufferEmpty();
	
	double parx[3] = {1,1,1};
	//TF1 *gx = new TF1("gx","gaus(0)",0.999,1.001);
	TF1 *gx = new TF1("gx","gaus(0)",dme_5->GetXaxis()->GetXmin(),dme_5->GetXaxis()->GetXmax());
	gx->SetParameters(parx);
	dme_5->Fit(gx,"R");
	
	double pare[3] = {1,1,1};
	//TF1 *ge = new TF1("ge","gaus(0)",0.8,1.2);
	TF1 *ge = new TF1("ge","gaus(0)",ee_5->GetXaxis()->GetXmin(),ee_5->GetXaxis()->GetXmax());
	ge->SetParameters(pare);
	ee_5->Fit(ge,"R");
	
	output->Write();
	
	delete dmpz_1;
	delete dme_1;
	delete dmpz_2;
	delete dme_2;
	delete dmpz_3;
	delete dme_3;
    delete nupz_1;
    delete nue_1;
    delete nupz_2;
    delete nue_2;
    delete nupz_3;
    delete nue_3;
    delete nue_4;
    delete nue_5;
	delete epz_3;
	delete ee_3;
	delete dme_4;
	delete ee_4;
    delete nuepz_3;
    delete nuee_3;
    delete nuee_4;
    delete nuee_5;
	delete output;
	delete array;
}

void DetectorAnalysis::UnInit()
{

}
