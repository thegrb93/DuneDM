#include "DMAnalysis.h"
#include "DMHistograms.h"

#include <iostream>
#include <sstream>
#include <regex>

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <TH2D.h>
#include <TF1.h>
#include <TChain.h>
#include <TGraph2D.h>
#include <TLegend.h>
#include <fstream>
#include <TApplication.h>
#include <dirent.h>

// Header files for the classes stored in the TTree if any.
#include "ExRootClasses.h"
#include "TClonesArray.h"

// Detector sensitivity headers
#include "Particle.h"
#include "DUNEdet.h"
#include "DMElscattering.h"

const double emass = 0.000511;

void getFolderFiles(std::string path, std::vector<std::string>& files) {
    DIR *dpdf;
    struct dirent *epdf;
    dpdf = opendir(path.c_str());
    if (dpdf != NULL){
        while (epdf = readdir(dpdf)){
            std::string name(epdf->d_name);
            if(name.size()>=5 && name.substr(name.size()-5,5)==".root")
                files.push_back(path+"/"+std::string(epdf->d_name));
        }
    }
}


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
    std::vector<std::string> files;
    getFolderFiles(folder, files);
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

std::regex param_match(".+?_([\\d\\.]+)_([\\d\\.]+)_([\\d\\.]+)_([\\d\\.]+).*");
int DMParameters(const std::string& filen, double& vpmass, double& chimass, double& kappa, double& alpha/*, double& xsection*/)
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
        /*parser.clear();
        parser << matches[5];
        parser >> xsection;
        if(parser.bad()) goto error;*/
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
	delete graphchi2;
	delete graphdmee;
    delete graphsig;
    delete graphbg;
	delete neutrino_electron_e;
}

void StatisticsAnalysis::Init()
{
	index = 0;
    graphchi2 = new TGraph2D((int)files.size());
    graphchi2->SetName("chi2");
    graphchi2->SetTitle("Dark Matter #chi^{2};VP mass(GeV);#chi mass(GeV);#chi^{2}");

    graphsig = new TGraph2D((int)files.size());
    graphsig->SetName("signal");
    graphsig->SetTitle("Signal;VP mass(GeV);#chi mass(GeV);Counts");

    graphbg = new TGraph2D((int)files.size());
    graphbg->SetName("bg");
    graphbg->SetTitle("Background;VP mass(GeV);#chi mass(GeV);Counts");
	
	graphdmee = new TGraph2D((int)files.size());
	graphdmee->SetName("dmee");
	graphdmee->SetTitle("Electron Scatter E against dm Mean;VP mass(GeV);#chi mass(GeV);E(GeV)");

	DUNEDetector det;
	DMscattering scatter;

	double probMax = 1e-15;

    neutrino_electron_e = new TH1D("nuee3","Nu-Electron Scatter E;E (GeV)", 100, 0, 6);
	TFile* neutrinos = new TFile("data/g4lbne_nudata.root");
    TTree* neutrino_tree = (TTree*)neutrinos->Get("nudata");
    neutrino_tree->SetMakeClass(1);
    float ndxdz, ndydz, ndz;
    double Nimpwt;
    double NenergyN[5];
    double NWtNear[5];
    neutrino_tree->SetBranchAddress("Ndxdz", &ndxdz);
    neutrino_tree->SetBranchAddress("Ndydz", &ndydz);
    neutrino_tree->SetBranchAddress("Npz", &ndz);
    neutrino_tree->SetBranchAddress("NenergyN[5]", &NenergyN);
    neutrino_tree->SetBranchAddress("Nimpwt", &Nimpwt);
    neutrino_tree->SetBranchAddress("NWtNear[5]", &NWtNear);

    long long neutrino_entries = neutrino_tree->GetEntries();
    int neutrino_intersectcount = 0;
    int neutrino_scattercount = 0;
    int neutrino_total = 5000000;
    double total_weight = 0;
    std::cout << "Generating neutrino distribution...\n";
    for(long long i = 0, j = 0; j < neutrino_total; i=(i+1)%neutrino_entries, ++j) {
        neutrino_tree->GetEvent(i);
        Particle neutrino(0);
        Particle electron(emass);
        double weight = Nimpwt*NWtNear[0];
        total_weight += weight;

        neutrino.FourMomentum(ndxdz*ndz, ndydz*ndz, ndz, NenergyN[0]);

        //nupz1->Fill(neutrino.pz, weight);
        //nue1->Fill(neutrino.E, weight);

        double pLen, dx, dy, dz, tmin, tmax;
        neutrino.getNorm(pLen, dx, dy, dz);
        if(det.intersect(dx, dy, dz, tmin, tmax)) {

            //nupz2->Fill(neutrino.pz, weight);
            //nue2->Fill(neutrino.E, weight);
            if(scatter.probscatterNeutrino(neutrino, tmax-tmin)) {
                scatter.scattereventNeutrino(neutrino, electron);

                //nupz3->Fill(neutrino.pz, weight);
                //nue3->Fill(neutrino.E, weight);

                //nu_epz1->Fill(electron.pz, weight);
                neutrino_electron_e->Fill(electron.E, weight);

                //double pnorm = sqrt(electron.px*electron.px+electron.py*electron.py+electron.pz*electron.pz);
                //nu_ethe1->Fill(acos(std::min<double>(std::max<double>(electron.pz/pnorm, -1), 1));
            }
        }
    }
    double scale = (double)100000/total_weight;
    neutrino_electron_e->Scale(scale);
    
    delete neutrinos;
}

void StatisticsAnalysis::Analyze(const std::string& filen)
{
    //if(index>10) return;
	double vpmass, chimass, kappa, alpha;
	if(DMParameters(filen, vpmass, chimass, kappa, alpha)) return;
	std::cout << filen << std::endl;

	DUNEDetector det;
	DMscattering scatter;
	double probMax = 1e-15;
	
	TH1D* dm_electron_e = new TH1D("dme","DM-Electron Scatter E;E (GeV)", 100, 0, 6);
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
				Particle darkmatter1(chimass);
				Particle electron1(emass);
				darkmatter1.FourMomentum(particle->Px,particle->Py,-particle->Pz,particle->E);

                double pLen, dx, dy, dz, tmin, tmax;
                darkmatter1.getNorm(pLen, dx, dy, dz);
				if(det.intersect(dx, dy, dz, tmin, tmax)) {
                    if(scatter.probscatter(vpmass, chimass, kappa, alpha, darkmatter1, tmax-tmin)) {
                        scatter.scatterevent(vpmass, chimass, kappa, alpha, darkmatter1,
                                             electron1);
                            dm_electron_e->Fill(electron1.E);
                    }
                }
			}
		}
	}
	
	double chi2 = 0;
	for(int i = 1; i<=dm_electron_e->GetNbinsX(); ++i)
	{
		double N_bg = neutrino_electron_e->GetBinContent(i);
		double N_tot = dm_electron_e->GetBinContent(i) + N_bg;
		if(N_bg == 0) continue;
		chi2 += 2*(N_bg*log(N_bg/N_tot) + N_tot - N_bg);
	}

    graphsig->SetPoint(index, vpmass, chimass, dm_electron_e->Integral());
	graphchi2->SetPoint(index, vpmass, chimass, chi2);
	graphdmee->SetPoint(index, vpmass, chimass, dm_electron_e->GetMean(1));
	++index;
	
	delete array;
	delete dm_electron_e;
}

void StatisticsAnalysis::UnInit()
{
	/*mkdir("images", S_IRWXU | S_IRWXG | S_IRWXO);
	std::string params(argv[2]);
	std::string outputfilen(argv[3]);*/

	TCanvas* canvas = new TCanvas("output", "", 1920,1080);
	canvas->SetFillColor(33);
	canvas->SetFrameFillColor(17);
	canvas->SetGrid();
	
	canvas->SetTheta(-90);
	canvas->SetPhi(-0.0001);

    graphchi2->Draw("surf1 Z");
    canvas->Update();
    canvas->SaveAs("chi2.png");

    graphsig->Draw("surf1 Z");
    canvas->Update();
    canvas->SaveAs("sig.png");

	graphdmee->Draw("surf1 Z");
	canvas->Update();
	canvas->SaveAs("dme.png");
	
	neutrino_electron_e->Draw();
	canvas->Update();
	canvas->SaveAs("nue.png");

    /*canvas->Draw();
    gApp->Run();*/

	delete canvas;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


DetectorAnalysis::DetectorAnalysis() : hists(new DMHistograms)
{
}

DMAnalysis* DetectorAnalysis::create()
{
	return new DetectorAnalysis();
}

DetectorAnalysis::~DetectorAnalysis()
{
}

inline void calculateNeutrinoNearDetector(int ptype, double Necm, double pdPx, double pdPy, double pdPz, double Vx, double Vy, double Vz, double Nimpwt, double& energy, double& weight) {
    const double detx = 0.0;
    const double dety = 0.0;
    const double detz = 57400.0;
    const double rdet = 100.0;

    double parent_mass;
    switch(ptype)
    {
        case 5:
        case 6:
            parent_mass = 0.105658389; //Muon mass (GeV)
        case 8:
        case 9:
            parent_mass = 0.13957; //Pi mass (GeV)
            break;
        case 10:
            parent_mass = 0.49767; //Keon0 mass (GeV)
        case 11:
        case 12:
            parent_mass = 0.49368; //Keon mass (GeV)
        default:
            return;
    }
    double parent_energy = sqrt(pdPx*pdPx +
                                pdPy*pdPy +
                                pdPz*pdPz +
                                parent_mass*parent_mass);
    double gamma = parent_energy / parent_mass;
    double gamma_sqr = gamma*gamma;
    double beta_mag = sqrt((gamma_sqr-1.)/gamma_sqr);

    double rad = sqrt((detx-Vx)*(detx-Vx) +
                      (dety-Vy)*(dety-Vy) +
                      (detz-Vz)*(detz-Vz));

    double parentp = sqrt((pdPx*pdPx)+
                          (pdPy*pdPy)+
                          (pdPz*pdPz));
    double costh_pardet = (pdPx*(detx-Vx) +
                           pdPy*(dety-Vy) +
                           pdPz*(detz-Vz))/(parentp*rad);

    if (costh_pardet>1.) costh_pardet = 1.;
    else if (costh_pardet<-1.) costh_pardet = -1.;
    double theta_pardet = acos(costh_pardet);

    double emrat = 1./(gamma*(1. - beta_mag * cos(theta_pardet)));
    double sangdet = (rdet*rdet /(rad*rad)/ 4.);

    energy = emrat*Necm;
    weight = sangdet * emrat * emrat * Nimpwt / M_PI;
}

void getNeutrinoDistribution(DMHistograms* hists, double normalize){
    DUNEDetector det;
    DMscattering scatter;

    int neutrino_intersectcount = 0;
    int neutrino_scattercount = 0;
    double total_weight = 0;

    Particle neutrino(0);
    Particle electron(emass);

    std::vector<std::string> neutrino_files;
    getFolderFiles("data/neutrinos", neutrino_files);
    for (std::vector<std::string>::iterator filen = neutrino_files.begin();
         filen != neutrino_files.end(); ++filen) {
        std::cout << *filen << std::endl;
        TFile *neutrinos = new TFile(filen->c_str());
        TTree *neutrino_tree = (TTree *) neutrinos->Get("nudata");
        neutrino_tree->SetMakeClass(1);
        int ptype, Ntype;
        float Necm, pdPx, pdPy, pdPz, Vx, Vy, Vz;
        double Nimpwt;
        neutrino_tree->SetBranchAddress("ptype", &ptype);
        neutrino_tree->SetBranchAddress("Ntype", &Ntype);
        neutrino_tree->SetBranchAddress("Necm", &Necm);
        neutrino_tree->SetBranchAddress("pdPx", &pdPx);
        neutrino_tree->SetBranchAddress("pdPy", &pdPy);
        neutrino_tree->SetBranchAddress("pdPz", &pdPz);
        neutrino_tree->SetBranchAddress("Vx", &Vx);
        neutrino_tree->SetBranchAddress("Vy", &Vy);
        neutrino_tree->SetBranchAddress("Vz", &Vz);
        neutrino_tree->SetBranchAddress("Nimpwt", &Nimpwt);

        long long neutrino_entries = neutrino_tree->GetEntries();

        for (long long i = 0, j = 0; j < 1e5; i = (i + 1) % neutrino_entries, ++j) {
            //for (long long i = 0; i < neutrino_entries; ++i) {
            neutrino_tree->GetEvent(i);
            double energy, weight;//Nimpwt * NWtNear[0];
            calculateNeutrinoNearDetector(ptype, Necm, pdPx, pdPy, pdPz, Vx, Vy, Vz, Nimpwt, energy, weight);

            //total_weight += weight;
            //neutrino.FourMomentum(ndxdz * ndz, ndydz * ndz, ndz, ne);
            //double pLen, dx, dy, dz;
            //neutrino.getNorm(pLen, dx, dy, dz);
            //double pt = sqrt(neutrino.px * neutrino.px + neutrino.py * neutrino.py);
            //double theta = acos(std::min<double>(std::max<double>(dz, -1), 1));

            //hists->AddProductionNu(neutrino.E, neutrino.pz, pt, theta, ntype, weight);
            hists->AddProductionNu(energy, 0, 0, 0, Ntype, weight);

            /*double tmin = 0, tmax = 0;
            if (det.intersect(dx, dy, dz, tmin, tmax)) {
                hists->AddDetectorNu(neutrino.E, neutrino.pz, pt, theta, weight);

                if (scatter.probscatterNeutrino(neutrino, tmax - tmin)) {
                    scatter.scattereventNeutrino(neutrino, electron);

                    /*hists->AddScatterNu(
                        darkmatter.E,
                        darkmatter.pz,
                        sqrt(darkmatter.px*darkmatter.px+darkmatter.py*darkmatter.py),
                        acos(std::min<double>(std::max<double>(electron.pz / sqrt(darkmatter.px * darkmatter.px + darkmatter.py * darkmatter.py + darkmatter.pz * darkmatter.pz), -1), 1)),
                        weight
                    );

                    hists->AddScatterBgElectron(
                            electron.E,
                            electron.pz,
                            sqrt(electron.px * electron.px + electron.py * electron.py),
                            acos(std::min<double>(std::max<double>(electron.pz / sqrt(electron.px * electron.px +
                                                                                      electron.py * electron.py +
                                                                                      electron.pz * electron.pz),
                                                                   -1), 1)),
                            weight
                    );
                }
            }*/
        }
        delete neutrino_tree;
        delete neutrinos;
    }
    if (!normalize) {
        //hists->ScaleNeutrinos(3e16 / total_weight);
    }
}

void DetectorAnalysis::Init() {
	//if(gOptions[OPT_DETECTOR])
	//	detectorType = gOptions[OPT_DETECTOR].last()->arg;
	//else
    	detectorType = "DUNE";

	if(detectorType=="DUNE")
	{
		smear_mean = 0;
		smear_sigma = 0.06;
	}
    normalize_histos = false;

	/*if(gOptions[OPT_DET_SMEAR_SIG])
		smear_sigma = toDouble(gOptions[OPT_DET_SMEAR_SIG].last()->arg);
	else
		smear_sigma = 0;
	
	if(gOptions[OPT_DET_SMEAR_MEAN])
		smear_mean = toDouble(gOptions[OPT_DET_SMEAR_MEAN].last()->arg);
	else
		smear_mean = 0;*/
	
	smear_sigma = 0.06;

    getNeutrinoDistribution(hists, normalize_histos);
}

void DetectorAnalysis::Analyze(const std::string& filen) {
    return;
	double vpmass, chimass, kappa, alpha;
	if(DMParameters(filen, vpmass, chimass, kappa, alpha)) return;

	DUNEDetector det;
	DMscattering scatter;

	TClonesArray* array = new TClonesArray("TRootLHEFParticle", 5);
	branch->SetAddress(&array);

    Particle darkmatter(chimass);
    Particle electron(emass);

    for(Int_t i = 0; i < nentries; ++i)
    {
		branch->GetEntry(i);
        TRootLHEFParticle* particle = (TRootLHEFParticle*)array->At(3);

        double pLen, dx, dy, dz, theta, pt;
        darkmatter.FourMomentum(particle->Px,particle->Py,-particle->Pz,particle->E);
        darkmatter.getNorm(pLen, dx, dy, dz);
        theta = acos(std::min<double>(std::max<double>(dz, -1), 1));
        pt = sqrt(darkmatter.px*darkmatter.px+darkmatter.py*darkmatter.py);

        hists->AddProductionDM(particle->E, -particle->Pz, pt, theta);

        double tmin = 0, tmax = 0;
        if(det.intersect(dx, dy, dz, tmin, tmax))
        {
            hists->AddDetectorDM(darkmatter.E, darkmatter.pz, pt, theta);
            //double vel = pLen/particle->E*299792458;
            //dmv->Fill(vel);
            //dmtime->Fill(tmin/vel - (570-3.2)/299792458);

            if(scatter.probscatter(vpmass,chimass,kappa,alpha,darkmatter,tmax-tmin)) {
                scatter.scatterevent(vpmass, chimass, kappa, alpha, darkmatter, electron);

                /*hists->AddScatterDM(
                        darkmatter.E,
                        darkmatter.pz,
                        sqrt(darkmatter.px*darkmatter.px+darkmatter.py*darkmatter.py),
                        acos(std::min<double>(std::max<double>(electron.pz / sqrt(darkmatter.px * darkmatter.px + darkmatter.py * darkmatter.py + darkmatter.pz * darkmatter.pz), -1), 1)
                        ));*/

                hists->AddScatterSigElectron(
                        electron.E,
                        electron.pz,
                        sqrt(electron.px*electron.px+electron.py*electron.py),
                        acos(std::min<double>(std::max<double>(electron.pz / sqrt(electron.px * electron.px + electron.py * electron.py + electron.pz * electron.pz), -1), 1)
                        ));
            }
        }
	}
	delete array;
}

void DetectorAnalysis::UnInit() {
    if(normalize_histos)
        hists->NormalizeHistograms();
    hists->SaveHistograms();

    delete hists;
}
