// main program
#include <iostream>
#include <vector>
#include <sstream>

#include <sys/stat.h>

#include <TApplication.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "histograms.h"

static std::vector<std::string> files;
static std::map<std::string,std::string> options = {
	{"mode", "distributions"},
	{"particle", "33"},
	{"attribute", "px"}
};

static int getParticleBranch(TFile* file, TTree*& tree, TBranch*& branch)
{
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
	Long64_t nentries = tree->GetEntries();	

	branch = tree->GetBranch("Particle");
	if(!branch) {
		std::cout << "Root couldn't find Particle branch in the LHEF tree.\n";
		return 1;
	}
	
	return 0;
}

static void plotDistributions()
{
	TFile* f = new TFile(files[0].c_str());
	TTree* tree;
	TBranch* branch;
	
	if(getParticleBranch(f, tree, branch)) return;
	
	DarkMatterDistribution distr(options["particle"],options["attribute"]);
	int nentries = tree->GetEntries();
	distr.Fill(branch, nentries);
	delete f;
	
	distr.Save();
}

static void plotStatistics()
{
	DarkMatterAnalysis analysis(files.size());
		
	for(std::vector<std::string>::iterator i = files.begin(); i!=files.end(); ++i) {
		TFile* f = new TFile(i->c_str());
		TTree* tree;
		TBranch* branch;
		if(getParticleBranch(f, tree, branch)) return;

		int nentries = tree->GetEntries();
		std::cout<<"File (" << (i-files.begin()+1) <<"/"<< files.size() << "), number of entries: " << nentries << std::endl;
		analysis.Fill(*i, branch, nentries);
		
		delete f;
	}
	std::cout << "Done processing.\n";

	analysis.Save();
}

static void plotSensitivity()
{
	TFile* f = new TFile(files[0].c_str());
	TTree* tree;
	TBranch* branch;
	
	if(getParticleBranch(f, tree, branch)) return;
	
	DetectorAnalysis distr;
	int nentries = tree->GetEntries();
	distr.Fill(files[0], branch, nentries);
	delete f;
	
	distr.Save();
}

int main (int argc, char** argv)
{
	for(int i = 1; i < argc; ++i)
	{
		char* str = argv[i];
		if(str[0]=='-')
		{
			std::string option(str+1);
			if(option=="h")
			{
				std::cout << "Usage: DuneDM [options] <files>\n"
							 "Options:\n"
							 "-mode  --  distributions, statistics, or detector. (def. distributions)\n"
							 "-particle  --   particle pdgcode. (def. 33)\n"
							 "-attribute  --  attribute to plot. (def. px)\n";
				return 0;
			}
			std::map<std::string,std::string>::iterator find = options.find(option);
			if(find!=options.end())
			{
				if(i+1<argc)
				{
					find->second = std::string(argv[i+1]);
					++i;
				}
			}
			else
			{
				std::cout << "Invalid option: " << option << "\n";
				return 0;
			}
		}
		else
			files.push_back(std::string(str));
	}

	if(files.size()>0)
	{
		std::string mode = options["mode"];
		if(mode=="distributions")
			plotDistributions();
		else if(mode=="statistics")
			plotStatistics();
		else if(mode=="detector")
			plotSensitivity();
		else
			std::cout << "Invalid mode.\n";
	}
	else
		std::cout << "Expected input root files.\n";

	return 0;
}

