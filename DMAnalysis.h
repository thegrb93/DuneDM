#pragma once
#include "cxxopts.hpp"
#include <string>
#include <fstream>
#include <vector>
#include <TROOT.h>
#include <TVectorD.h>

class TApplication;
extern TApplication* gApp;
extern cxxopts::Options options;

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
	virtual int Init(TFile*, TBranch*) = 0;
	virtual void UnInit() = 0;
};

class DetectorAnalysis : public DMAnalysis
{
    DMHistograms* hists;
	std::string detectorType;
	double smear_sigma;
	double smear_mean;
	double xsection, totalevents;
    bool normalize_histos;
public:
	DetectorAnalysis();
	~DetectorAnalysis();
protected:
	int Analyze(TFile*, TBranch*);
	int Init(TFile*, TBranch*);
	void UnInit();
};

class SensitivityAnalysis : public DMAnalysis
{
	double smear_sigma;
	double smear_mean;
    double detector_efficiency;
    int nbins; double minbin, maxbin;
	TH1D *dm_energy, *nu_energy;
	TProfile *theta_avg;
	TFile* nu_cache, *dm_cache;
	TVectorD* xsection;
	bool loadedDM;
public:
    double chisqr;
    double binstart_max, binend_max;
	SensitivityAnalysis();
	~SensitivityAnalysis();
protected:
	int Analyze(TFile*, TBranch*);
	int Init(TFile*, TBranch*);
	void UnInit();
};

class SensitivityScan
{
public:
    SensitivityScan();
	int Process(const std::vector<std::string>& folders);
};

