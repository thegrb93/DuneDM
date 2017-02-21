#include "histograms.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TChain.h>
#include <TGraph2D.h>
#include <TLegend.h>

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

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

int getDMParameters(const std::string& filen, double& vpmass, double& chimass, double& kappa, double& alpha)
{
	std::vector<std::string> params = split(filen, '_');
	if(params.size()<5)
	{
		goto error;
	}
	
	{
		std::stringstream parser;
		parser << params[1];
		parser >> vpmass;
		if(parser.bad()) goto error;
		parser.clear();
		parser << params[2];
		parser >> chimass;
		if(parser.bad()) goto error;
		parser.clear();
		parser << params[3];
		parser >> kappa;
		if(parser.bad()) goto error;
		parser.clear();
		parser << params[4];
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

DarkMatterAnalysis::DarkMatterAnalysis(int nfiles) :
	graph(new TGraph2D(nfiles)),
	nfiles(nfiles),
	index(0)
{
	graph->SetTitle("#mu_{#chi pz};VP mass(GeV);#chi mass(GeV);#mu_{pz}(GeV/c^{2})");
}

DarkMatterAnalysis::~DarkMatterAnalysis()
{
	delete graph;
}

void DarkMatterAnalysis::Fill(const std::string& filen, TBranch* branch, int nentries)
{
	double vpmass, chimass, kappa, alpha;
	if(getDMParameters(filen, vpmass, chimass, kappa, alpha)) return;
	
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

void DarkMatterAnalysis::Save()
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

DarkMatterDistribution::DarkMatterDistribution(const std::string& particle, const std::string& attr)
{
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

DarkMatterDistribution::~DarkMatterDistribution()
{
	delete histo;
	delete output;
}

void DarkMatterDistribution::Fill(TBranch* branch, int nentries)
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

void DarkMatterDistribution::Save()
{
	output->Write(); 
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void DetectorAnalysis::Fill(const std::string& filen, TBranch* branch, int nentries)
{
	double vpmass, chimass, kappa, alpha;
	if(getDMParameters(filen, vpmass, chimass, kappa, alpha)) return;
	
	const double emass = 0.000511;
	
	Particle darkphoton(vpmass);
	Particle darkmatter1(chimass);
	Particle darkmatter2(chimass);
	Particle electron1(emass);
	Particle electron2(emass);
	
	Kinematics kin;	
	detector det;
	DMscattering scatter;
	
	int scatterCount = 0;
	int Nscatter;
	int Nelectron;
	double probMax = 10e-15;
	
	TFile* output = new TFile("output.root", "RECREATE");
	TH1D* dmpz = new TH1D("dmpz","Dark matter Scatter P_{z};P_{z} (GeV/c^{2})", 100, 0, 0);
	TH1D* dme = new TH1D("dme","Dark matter Scatter E;E (GeV)", 100, 0, 0);
	TH1D* epz = new TH1D("epz","Electron Scatter P_{z};P_{z} (GeV/c^{2})", 100, 0, 0);
	TH1D* ee = new TH1D("ee","Electron Scatter E;E (GeV)", 100, 0, 0);
	
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
				darkmatter1.FourMomentum(particle->Px,particle->Py,-particle->Pz,particle->E);
				
				int DMSwitch = 0;
				det.intersect(DMSwitch,scatterCount,darkmatter1);
				scatter.probscatter(DMSwitch,Nscatter,probMax,vpmass,chimass,kappa,alpha,darkmatter1);
				scatter.scatterevent(DMSwitch,Nelectron,vpmass,chimass,kappa,alpha,darkmatter1,electron1);
				if(DMSwitch == 2)
				{
					dmpz->Fill(darkmatter1.pz);
					dme->Fill(darkmatter1.E);
					epz->Fill(electron1.pz);
					ee->Fill(electron1.E);
				}	
			}
		}
	}
	dmpz->BufferEmpty();
	dme->BufferEmpty();
	epz->BufferEmpty();
	ee->BufferEmpty();
	output->Write();
	
	delete dmpz;
	delete dme;
	delete epz;
	delete ee;
	delete output;
	delete array;
}

void DetectorAnalysis::Save()
{

}


/*
std::string particle_names[5] = {"Pqd","Pqd1","Chi","Chi bar","V"}; 
for(int i = 2; i<4; ++i)
{
std::string pxname = "X Momentum" + params;
std::string pyname = "Y Momentum " + params;
std::string pzname = "Z momentum " + params;
std::string thetaname = "Theta " + params;
std::string ename = "Energy " + params;

std::stringstream numstr;
numstr << i;

particle_hists[i].px = new TH1D((pxname+numstr.str()).c_str(), pxname.c_str(), 100, -2.5, 2.5);
particle_hists[i].py = new TH1D((pyname+numstr.str()).c_str(), pyname.c_str(), 100, -2.5, 2.5);
particle_hists[i].pz = new TH1D((pzname+numstr.str()).c_str(), pzname.c_str(), 100, -70.0, 0.0);
particle_hists[i].theta = new TH1D((thetaname+numstr.str()).c_str(), thetaname.c_str(), 100, -3.2, 3.2);
particle_hists[i].e = new TH1D((ename+numstr.str()).c_str(), ename.c_str(), 100, 0.0, 60.0);

particle_hists[i].px->SetLineColor(i);
particle_hists[i].px->SetMarkerSize(0.8);
particle_hists[i].px->SetStats(0);
particle_hists[i].py->SetLineColor(i);
particle_hists[i].py->SetMarkerSize(0.8);
particle_hists[i].py->SetStats(0);
particle_hists[i].pz->SetLineColor(i);
particle_hists[i].pz->SetMarkerSize(0.8);
particle_hists[i].pz->SetStats(0);
particle_hists[i].theta->SetLineColor(i);
particle_hists[i].theta->SetMarkerSize(0.8);
particle_hists[i].theta->SetStats(0);
particle_hists[i].e->SetLineColor(i);
particle_hists[i].e->SetMarkerSize(0.8);
particle_hists[i].e->SetStats(0);
}
*/

/*histograms->px->Fill(particle->Px);
histograms->py->Fill(particle->Py);
histograms->pz->Fill(particle->Pz);
histograms->theta->Fill(particle->Phi);
histograms->e->Fill(particle->E);*/
			
