#pragma once
#include <string>
#include <vector>
#include <TROOT.h>
#include "optionparser.h"

enum optionIndex {
OPT_UNKNOWN, OPT_HELP, OPT_MODE, OPT_PARTICLE, OPT_PARTICLEATTRIBUTE,
OPT_DETECTOR, OPT_DET_SMEAR_SIG, OPT_DET_SMEAR_MEAN};

extern option::Option* gOptions;


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
	
	static int DMParameters(const std::string& filen, double& vpmass, double& chimass, double& kappa, double& alpha);
protected:
	TBranch* branch;
	Long64_t nentries;
	virtual void Analyze(const std::string& file) = 0;
	virtual void Init() = 0;
	virtual void UnInit() = 0;
};

class StatisticsAnalysis : public DMAnalysis
{
	TGraph2D* graph;
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

class DarkMatterDistribution : public DMAnalysis
{
	TFile* output;
	TH1D* histo;
	int pdgCode;
	double TRootLHEFParticle::*attribute;

public:
	DarkMatterDistribution();
	~DarkMatterDistribution();
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
public:
	DetectorAnalysis();
	~DetectorAnalysis();
	static DMAnalysis* create();
protected:
	void Analyze(const std::string& file);
	void Init();
	void UnInit();
};
