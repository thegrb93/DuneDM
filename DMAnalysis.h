#pragma once
#include <string>
#include <vector>
#include <TROOT.h>

class TApplication;
extern TApplication* gApp;

class TFile;
class TBranch;
class TGraph2D;
class TH1D;
class TRootLHEFParticle;
class DMHistograms;

class DMAnalysis
{
public:
    std::string folder;
	std::vector<std::string> files;
	DMAnalysis();
	virtual ~DMAnalysis();
	int Process();

protected:
	TBranch* branch;
	Long64_t nentries;
	virtual void Analyze(const std::string& file) = 0;
	virtual void Init() = 0;
	virtual void UnInit() = 0;
};

class StatisticsAnalysis : public DMAnalysis
{
	TGraph2D* graphchi2, *graphdmee, *graphnuee, *graphsig, *graphbg;
	TH1D* neutrino_electron_e;
	int index;
	
public:
	StatisticsAnalysis();
	~StatisticsAnalysis();
	static DMAnalysis* create();

protected:
	void Analyze(const std::string& file);
	void Init();
	void UnInit();
};

class DetectorAnalysis : public DMAnalysis
{
    DMHistograms* hists;
	std::string detectorType;
	double smear_sigma;
	double smear_mean;
    bool normalize_histos;
public:
	DetectorAnalysis();
	~DetectorAnalysis();
	static DMAnalysis* create();
protected:
	void Analyze(const std::string& file);
	void Init();
	void UnInit();
};

