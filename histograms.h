#pragma once
#include <string>

class TFile;
class TBranch;
class TGraph2D;
class TH1D;
class TRootLHEFParticle;

class DarkMatterAnalysis
{
	TGraph2D* graph;
	int nfiles;
	int index;
	
	public:
	DarkMatterAnalysis(int nfiles);
	~DarkMatterAnalysis();
	
	void Fill(const std::string& file, TBranch* branch, int nentries);
	void Save();
};

class DarkMatterDistribution
{
	TFile* output;
	TH1D* histo;
	int pdgCode;
	double TRootLHEFParticle::*attribute;

	public:
	DarkMatterDistribution(const std::string& particle, const std::string& attr);
	~DarkMatterDistribution();
	
	void Fill(TBranch* branch, int nentries);
	void Save();
};

class DetectorAnalysis
{
	public:
	void Fill(const std::string& filen, TBranch* branch, int nentries);
	void Save();
};

