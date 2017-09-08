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

class DMAnalysis
{
public:
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
	std::string detectorType;
	double smear_sigma;
	double smear_mean;
    TFile* output;
    //Input dark matter distributions
    TH1D *dmpz_1, *dmpt_1, *dme_1, *dme_v, *dme_time, *dme_t1, *dmpz_2, *dmpt_2, *dme_2, *dme_t2, *dmpz_3, *dme_3, *dme_4, *dme_5;
    //Input neutrino distributions
    TH1D *nupz_1, *nupt_1, *nue_1, *nue_t1, *nupz_2, *nupt_2, *nue_2, *nue_t2, *nupz_3, *nue_3, *nue_4, *nue_5;
    //Signal-electron distributions
    TH1D *epz_3, *ept_3, *etheta_3, *ee_3, *ee_4, *ee_5;
    //Background-electron distributions
    TH1D *nuepz_3, *nuept_3, *nuetheta_3, *nuee_3, *nuee_4, *nuee_5;

public:
	DetectorAnalysis();
	~DetectorAnalysis();
	static DMAnalysis* create();
protected:
	void Analyze(const std::string& file);
	void Init();
	void UnInit();
};

