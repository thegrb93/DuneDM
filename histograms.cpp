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

void DarkMatterAnalysis::Fill(const std::string& file, TBranch* branch, int nentries)
{
	std::vector<std::string> params = split(file, '_');
	if(params.size()<5){std::cout << "Tried to parse invalid filename: " << file << std::endl; return;}
	std::stringstream parser;
	
	double x, y;
	double z = 0;
	parser << params[1];
	parser >> x;
	if(parser.bad()){}
	parser.clear();
	parser << params[2];
	parser >> y;
	
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
	graph->SetPoint(index, x, y, z);
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
	
}

void DarkMatterDistribution::Save()
{
	output->Write(); 
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
			
