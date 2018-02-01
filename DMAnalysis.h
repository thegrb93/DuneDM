#pragma once
#include <string>
#include <vector>
#include <TROOT.h>
#include <TVectorD.h>

class TApplication;
extern TApplication* gApp;

class TFile;
class TBranch;
class TGraph2D;
class TH1D;
class TProfile;
class TRootLHEFParticle;
class DMHistograms;

class DMAnalysis
{
public:
	std::vector<std::string> files;
	DMAnalysis();
	virtual ~DMAnalysis();
	int Process(std::string folder);

protected:
	double vpmass, chimass, kappa, alpha;
	virtual int Analyze(TFile*, TBranch*) = 0;
	virtual void Init(TFile*, TBranch*) = 0;
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

protected:
	int Analyze(TFile*, TBranch*);
	void Init(TFile*, TBranch*);
	void UnInit();
};

class DetectorAnalysis : public DMAnalysis
{
    DMHistograms* hists;
	std::string detectorType;
	double smear_sigma;
	double smear_mean;
	double dm_detector_scale;
	double xsection;
    bool normalize_histos;
public:
	DetectorAnalysis();
	~DetectorAnalysis();
protected:
	int Analyze(TFile*, TBranch*);
	void Init(TFile*, TBranch*);
	void UnInit();
};

class SensitivityAnalysis : public DMAnalysis
{
	double smear_sigma;
	double smear_mean;
	TH1D *dm_energy, *nu_energy;
	TProfile *theta_avg;
	TFile* nu_cache, *dm_cache;
	TVectorD* xsection;
	bool loadedDM;
public:
    double dm_detector_scale, nu_detector_scale;
    double chisqr;
	SensitivityAnalysis();
	~SensitivityAnalysis();
protected:
	int Analyze(TFile*, TBranch*);
	void Init(TFile*, TBranch*);
	void UnInit();
};

class SensitivityScan
{
public:
    SensitivityScan();
	int Process(const std::vector<std::string>& folders);
};

