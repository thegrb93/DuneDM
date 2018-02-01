#include "DMAnalysis.h"
#include "DMHistograms.h"
#include "chisq_pull.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <regex>
#include <functional>

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TF1.h>
#include <TChain.h>
#include <TGraph2D.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TStyle.h>
#include <fstream>
#include <TApplication.h>
#include <dirent.h>
#include <sys/stat.h>

// Header files for the classes stored in the TTree if any.
#include "ExRootClasses.h"
#include "TClonesArray.h"

// Detector sensitivity headers
#include "Particle.h"
#include "DUNEdet.h"
#include "DMElscattering.h"
#include "Random.h"

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

std::regex param_match(".*?DM_([\\d\\.]+)_([\\d\\.]+)_([\\d\\.]+)_([\\d\\.]+).*");
int DMParameters(const std::string& filen, double& vpmass, double& chimass, double& kappa, double& alpha/*, double& xsection*/) {
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

inline void calculateNeutrinoNearDetector(int ptype, int Ntype, double ppenergy, double ppdxdz, double ppdydz, double pppz, double Necm, double pdPx, double pdPy, double pdPz, double Vx, double Vy, double Vz, double mupare, double muparpx, double muparpy, double muparpz, double Nimpwt, double& nu_energy, double& weight) {
    const double detx = 0.0;
    const double dety = 0.0;
    const double detz = 57400.0;
    const double rdet = 100.0;

    const double mumass  =    0.105658389;
    const double taumass =    1.77682;

    double parent_mass;
    switch(ptype)
    {
        case 5:
        case 6:
            parent_mass = mumass; //Muon mass (GeV)
            break;
        case 8:
        case 9:
            parent_mass = 0.13957; //Pi mass (GeV)
            break;
        case 10:
            parent_mass = 0.49767; //Keon0 mass (GeV)
            break;
        case 11:
        case 12:
            parent_mass = 0.49368; //Keon mass (GeV)
            break;
        default:
            std::cout << "INVALID NEUTRINO PARENT! " << ptype << std::endl;
            weight = 0;
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

    nu_energy = emrat*Necm;
    weight = sangdet * emrat * emrat * Nimpwt / M_PI;

    //done for all except polarized muon
    // in which case need to modify weight
    if (ptype==5 || ptype==6)
    {
        //boost new neutrino to mu decay cm
        double beta[3];
        double p_nu[3]; //nu momentum
        beta[0]=pdPx / parent_energy;
        beta[1]=pdPy / parent_energy;
        beta[2]=pdPz / parent_energy;

        p_nu[0] = (detx - Vx) * nu_energy / rad;
        p_nu[1] = (dety - Vy) * nu_energy / rad;
        p_nu[2] = (detz - Vz) * nu_energy / rad;

        double partial = gamma*(beta[0]*p_nu[0]+
                                beta[1]*p_nu[1]+
                                beta[2]*p_nu[2]);
        partial = nu_energy-partial / (gamma+1.);
        double p_dcm_nu[4];
        for (int i=0;i<3;i++) p_dcm_nu[i]=p_nu[i]-beta[i]*gamma*partial;
        p_dcm_nu[3]=0.;
        for (int i=0;i<3;i++) p_dcm_nu[3]+=p_dcm_nu[i]*p_dcm_nu[i];
        p_dcm_nu[3]=sqrt(p_dcm_nu[3]);

        //boost parent of mu to mu production cm
        gamma= ppenergy / parent_mass;
        beta[0] = ppdxdz * pppz / ppenergy;
        beta[1] = ppdydz * pppz / ppenergy;
        beta[2] =                  pppz / ppenergy;
        partial = gamma*(beta[0]*muparpx+
                         beta[1]*muparpy+
                         beta[2]*muparpz);
        partial = mupare - partial / (gamma+1.);
        double p_pcm_mp[4];
        p_pcm_mp[0]=muparpx-beta[0]*gamma*partial;
        p_pcm_mp[1]=muparpy-beta[1]*gamma*partial;
        p_pcm_mp[2]=muparpz-beta[2]*gamma*partial;
        p_pcm_mp[3]=0.;
        for (int i=0;i<3;i++) p_pcm_mp[3]+=p_pcm_mp[i]*p_pcm_mp[i];
        p_pcm_mp[3]=sqrt(p_pcm_mp[3]);

        double wt_ratio = 1.;
        //have to check p_pcm_mp
        //it can be 0 if mupar..=0. (I guess muons created in target??)
        if (p_pcm_mp[3] != 0. ) {
            //calc new decay angle w.r.t. (anti)spin direction
            double costh = (p_dcm_nu[0]*p_pcm_mp[0]+
                            p_dcm_nu[1]*p_pcm_mp[1]+
                            p_dcm_nu[2]*p_pcm_mp[2])/(p_dcm_nu[3]*p_pcm_mp[3]);

            if (costh>1.) costh = 1.;
            else if (costh<-1.) costh = -1.;

            //calc relative weight due to angle difference
            if (Ntype == 53 || Ntype == 52)
            {
                wt_ratio = 1.-costh;
            }
            else if (Ntype == 56 || Ntype == 55)
            {
                double xnu = 2.* Necm / mumass;
                wt_ratio = ( (3.-2.*xnu) - (1.-2.*xnu)*costh ) / (3.-2.*xnu);
            }
            else if (Ntype == 58 || Ntype == 59)
            {
                double xnu = 2.* Necm / taumass;
                wt_ratio = ( (3.-2.*xnu) - (1.-2.*xnu)*costh ) / (3.-2.*xnu);
                std::cout << "calculating weight for tau neutrino; this may not be correct" << std::endl;
            }
            else
            {
                std::cout << "eventRates:: Bad neutrino type = " << Ntype << std::endl;
            }
        }
        weight *= wt_ratio;
    }
}

void iterateNeutrino(int detector_scale,
    std::function<void(Particle&, double)> fproduction,
    std::function<void(Particle&, double, double, double)> fdetector,
    std::function<void(Particle&, Particle&, Particle&, double)> fscatter
){
    DUNEDetector det;
    DMscattering scatter;

    Particle neutrino(0);
    Particle electron(emass);

    std::vector<std::string> neutrino_files;
    getFolderFiles("data/neutrinos", neutrino_files);
    for (auto filen = neutrino_files.begin(); filen != neutrino_files.end(); ++filen) {
        TFile *neutrinos = new TFile(filen->c_str());
        std::cout << neutrinos->GetPath() << std::endl;
        TTree *neutrino_tree = (TTree *) neutrinos->Get("nudata");
        neutrino_tree->SetMakeClass(1);
        int ptype, Ntype;
        float Necm, ndxdz, ndydz, ndz, pdPx, pdPy, pdPz, Vx, Vy, Vz, ppenergy, ppdxdz, ppdydz, pppz, mupare, muparpx, muparpy, muparpz;
        double Nimpwt;
        double NenergyN[5];
        double NWtNear[5];
        neutrino_tree->SetBranchAddress("ptype", &ptype);
        neutrino_tree->SetBranchAddress("Ntype", &Ntype);
        neutrino_tree->SetBranchAddress("Necm", &Necm);
        neutrino_tree->SetBranchAddress("Ndxdz", &ndxdz);
        neutrino_tree->SetBranchAddress("Ndydz", &ndydz);
        neutrino_tree->SetBranchAddress("Npz", &ndz);
        neutrino_tree->SetBranchAddress("NenergyN[5]", &NenergyN);
        neutrino_tree->SetBranchAddress("NWtNear[5]", &NWtNear);
        neutrino_tree->SetBranchAddress("pdPx", &pdPx);
        neutrino_tree->SetBranchAddress("pdPy", &pdPy);
        neutrino_tree->SetBranchAddress("pdPz", &pdPz);
        neutrino_tree->SetBranchAddress("Vx", &Vx);
        neutrino_tree->SetBranchAddress("Vy", &Vy);
        neutrino_tree->SetBranchAddress("Vz", &Vz);
        neutrino_tree->SetBranchAddress("ppenergy", &ppenergy);
        neutrino_tree->SetBranchAddress("ppdxdz", &ppdxdz);
        neutrino_tree->SetBranchAddress("ppdydz", &ppdydz);
        neutrino_tree->SetBranchAddress("pppz", &pppz);
        neutrino_tree->SetBranchAddress("mupare", &mupare);
        neutrino_tree->SetBranchAddress("muparpx", &muparpx);
        neutrino_tree->SetBranchAddress("muparpy", &muparpy);
        neutrino_tree->SetBranchAddress("muparpz", &muparpz);
        neutrino_tree->SetBranchAddress("Nimpwt", &Nimpwt);

        long long neutrino_entries = neutrino_tree->GetEntries();
        for (long long i = 0; i < neutrino_entries; ++i) {
        //for (long long i = 0; i < neutrino_entries; ++i) {
        //for(long long i = 0, j = 0; j < neutrino_total; i=(i+1)%neutrino_entries, ++j) {
            neutrino_tree->GetEvent(i);
            //double energy, weight;
            //calculateNeutrinoNearDetector(ptype, Ntype, ppenergy, ppdxdz, ppdydz, pppz, Necm, pdPx, pdPy, pdPz, Vx, Vy, Vz, mupare, muparpx, muparpy, muparpz, Nimpwt, energy, weight);
            double weight = Nimpwt * NWtNear[0];
            neutrino.FourMomentum(ndxdz*ndz, ndydz*ndz, ndz, NenergyN[0]);
            neutrino.calcOptionalKinematics();

            fproduction(neutrino, weight);
            double tmin = 0, tmax = 0;
            if (det.intersect(neutrino.normpx, neutrino.normpy, neutrino.normpz, tmin, tmax)) {
                fdetector(neutrino, tmin, tmax, weight);
                double pathlength = tmax-tmin;
                for(int j = 0; j<detector_scale; ++j){
                    if (scatter.probscatterNeutrino(neutrino, pathlength)) {
                        Particle neutrino2 = neutrino;
                        scatter.scattereventNeutrino(neutrino2, electron);
                        neutrino2.calcOptionalKinematics();
                        electron.calcOptionalKinematics();
                        fscatter(neutrino, neutrino2, electron, weight);
                    }
                }
            }
        }
        delete neutrino_tree;
        delete neutrinos;
    }
}

void iterateDarkmatter(int detitr, TBranch* branch, double vpmass, double chimass, double kappa, double alpha,
    std::function<void(Particle&, double, double)> fdetector,
    std::function<void(Particle&, Particle&, Particle&)> fscatter
){
    DUNEDetector det;
	DMscattering scatter;

	TClonesArray* array = new TClonesArray("TRootLHEFParticle", 5);
	branch->SetAddress(&array);

    Particle darkmatter(chimass);
    Particle electron(emass);

    long long nentries = branch->GetEntries();
    if(nentries==0){std::cout << "The root file contains no entries!" << std::endl; return;}
    branch->GetEntry(0);
    int darkmatter_index = -1, antidarkmatter_index = -1;
    for(Int_t i = 0; i < array->GetSize(); ++i)
    {
        TRootLHEFParticle* particle = (TRootLHEFParticle*)array->At(i);
        if(particle){
            int PID = particle->PID;
            if(PID == 33) darkmatter_index = i;
            if(PID == -33) antidarkmatter_index = i;
        }
    }
    if(darkmatter_index==-1){std::cout << "The array does not contain darkmatter!" << std::endl; return;}
    if(antidarkmatter_index==-1){std::cout << "The array does not contain antidarkmatter!" << std::endl; return;}
    for(Int_t i = 0; i < nentries; ++i)
    {
		branch->GetEntry(i);
        TRootLHEFParticle* lhefdarkmatter = (TRootLHEFParticle*)array->At(darkmatter_index);
        TRootLHEFParticle* lhefantidarkmatter = (TRootLHEFParticle*)array->At(antidarkmatter_index);
        /*std::cout << ((TRootLHEFParticle*)array->At(0))->PID << std::endl;
        std::cout << ((TRootLHEFParticle*)array->At(1))->PID << std::endl;
        std::cout << ((TRootLHEFParticle*)array->At(2))->PID << std::endl;
        std::cout << ((TRootLHEFParticle*)array->At(3))->PID << std::endl;
        std::cout << ((TRootLHEFParticle*)array->At(4))->PID << std::endl;*/
        
        double tmin, tmax;
        darkmatter.FourMomentum(lhefdarkmatter->Px,lhefdarkmatter->Py,-lhefdarkmatter->Pz,lhefdarkmatter->E);
        darkmatter.calcOptionalKinematics();
        if(det.intersect(darkmatter.normpx, darkmatter.normpy, darkmatter.normpz, tmin, tmax))
        {
            fdetector(darkmatter, tmin, tmax);
            double pathlength = tmax-tmin;
            for(int j = 0; j<detitr; ++j){
                if(scatter.probscatter(vpmass,chimass,kappa,alpha,darkmatter,pathlength)) {
                    Particle darkmatter2 = darkmatter;
                    scatter.scatterevent(vpmass, chimass, kappa, alpha, darkmatter2, electron);
                    darkmatter2.calcOptionalKinematics();
                    electron.calcOptionalKinematics();
                    fscatter(darkmatter, darkmatter2, electron);
                }
            }
        }
	}
	delete array;
}

DMAnalysis::DMAnalysis() {
}
DMAnalysis::~DMAnalysis() {
}

int DMAnalysis::Process(std::string folder) {
    if(DMParameters(folder, vpmass, chimass, kappa, alpha)){
		std::cout << "Couldn't parse DM parameters.\n";
		UnInit();
		return 1;
    }
    
    getFolderFiles(folder, files);
    if(files.size()==0) {
        std::cout << "Didn't find any root files.\n";
        UnInit();
        return 1;
    }
    
    bool first = true;
	for(auto fname = files.begin(); fname!=files.end(); ++fname) {
		TFile file(fname->c_str());
		TTree* tree;
		TBranch* branch;
		
		if(!file.IsOpen()) {
			std::cout << "Failed to open file: " << file.GetName() << std::endl;
			continue;
		}

		file.GetObject("LHEF",tree);
		if(!tree) {
			std::cout << "Root couldn't find LHEF tree in the input root file.\n";
			continue;
		}
		branch = tree->GetBranch("Particle");
		if(!branch) {
			std::cout << "Root couldn't find Particle branch in the LHEF tree.\n";
			continue;
		}
		if(first){
		    first = false;
		    Init(&file, branch);
		}
		if(Analyze(&file, branch))
		    continue;
	}
	
	UnInit();
	return 0;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

StatisticsAnalysis::StatisticsAnalysis(){
}

StatisticsAnalysis::~StatisticsAnalysis(){
	delete graphchi2;
	delete graphdmee;
    delete graphsig;
    delete graphbg;
	delete neutrino_electron_e;
}

void StatisticsAnalysis::Init(TFile*, TBranch*){
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

    neutrino_electron_e = new TH1D("nuee3","Nu-Electron Scatter E;E (GeV)", 100, 0, 6);
	
	iterateNeutrino(1,
	[=](Particle& neutrino, double weight){
	//nupz1->Fill(neutrino.pz, weight);
        //nue1->Fill(neutrino.E, weight);
    },
    [=](Particle& neutrino, double tmin, double tmax, double weight){
        //nupz2->Fill(neutrino.pz, weight);
        //nue2->Fill(neutrino.E, weight);
    },
    [=](Particle&, Particle& neutrino, Particle& electron, double weight){
        //nupz3->Fill(neutrino.pz, weight);
        //nue3->Fill(neutrino.E, weight);

        //nu_epz1->Fill(electron.pz, weight);
        //neutrino_electron_e->Fill(electron.E, weight);

        //double pnorm = sqrt(electron.px*electron.px+electron.py*electron.py+electron.pz*electron.pz);
        //nu_ethe1->Fill(acos(std::min<double>(std::max<double>(electron.pz/pnorm, -1), 1));
    });
}

int StatisticsAnalysis::Analyze(TFile* file, TBranch* branch){
	TH1D* dm_electron_e = new TH1D("dme","DM-Electron Scatter E;E (GeV)", 100, 0, 6);
	
	iterateDarkmatter(1, branch, vpmass, chimass, kappa, alpha,
	[=](Particle& darkmatter, double tmin, double tmax){
	},
	[=](Particle&, Particle&, Particle& electron){
	    dm_electron_e->Fill(electron.E);
	});
	
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
	
	delete dm_electron_e;
	return 0;
}

void StatisticsAnalysis::UnInit(){
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
DetectorAnalysis::DetectorAnalysis(){
    dm_detector_scale = 100;
    xsection = 0;
}

DetectorAnalysis::~DetectorAnalysis(){
}

void DetectorAnalysis::Init(TFile* file, TBranch* branch) {
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
	smear_mean = 0;
	smear_sigma = 0.06;
	
	std::stringstream runname;
	runname << std::fixed << std::setprecision(6) << "DM_" << vpmass << "_" << chimass << "_" << kappa << "_" << alpha;
	hists = new DMHistograms(runname.str());
	
	if(!hists->found_nu_output){
	    double detector_scale = 1000;
	    iterateNeutrino((int)detector_scale,
	    [=](Particle& neutrino, double weight){
            hists->AddProductionNu(neutrino.E, neutrino.pz, neutrino.pt, neutrino.theta, /*ntype,*/ weight);
	    },
	    [=](Particle& neutrino, double tmin, double tmax, double weight){
	        hists->AddDetectorNu(neutrino.E, neutrino.pz, neutrino.pt, neutrino.theta, weight);
	    },
	    [=](Particle&, Particle& neutrino, Particle& electron, double weight){
	        hists->AddScatterNu(neutrino.E, neutrino.pz, neutrino.pt, neutrino.theta, weight);
            hists->AddScatterBgElectron(electron.E, electron.pz, electron.pt, electron.theta, weight);
	    });
	    
        double desired_pot = 6e20;
        double nu_pot = 250000000;
        double scale = desired_pot/nu_pot;
	    hists->ScaleNeutrinos(scale, scale/detector_scale);
    }
}

int DetectorAnalysis::Analyze(TFile* file, TBranch* branch) {
    TVectorD *v;
	file->GetObject("xsection", v);
	if(v)
	    xsection += v->operator()(0);
	else
		std::cout << "Root couldn't find crosssection object in file.\n";
	
    {
        TH1D* E, *Px, *Py, *Pz, *Th, *Phi;
        file->GetObject("energies", E);
        file->GetObject("xmomentum", Px);
        file->GetObject("ymomentum", Py);
        file->GetObject("zmomentum", Pz);
        file->GetObject("theta", Th);
        file->GetObject("phi", Phi);
        if(E) hists->AddProductionDM(E, Px, Py, Pz, Th, Phi);
    }

    iterateDarkmatter((int)dm_detector_scale, branch, vpmass, chimass, kappa, alpha,
    [=](Particle& darkmatter, double tmin, double tmax){
        hists->AddDetectorDM(darkmatter.E, darkmatter.normpx*tmin, darkmatter.normpy*tmin, darkmatter.px, darkmatter.py, darkmatter.pz, darkmatter.pt, darkmatter.theta, darkmatter.phi);
    },
    [=](Particle&, Particle& darkmatter, Particle& electron){
        hists->AddScatterDM(darkmatter.E, darkmatter.pz, darkmatter.pt, darkmatter.theta);
        hists->AddScatterSigElectron(electron.E, electron.pz, electron.pt, electron.theta);
    });
    return 0;
}

void DetectorAnalysis::UnInit() {
    double desired_pot = 6e20;
    if(xsection!=0){
        double dm_pot = std::pow(files.size(),2)/(12*xsection*1e-40*100*6.022140857e23);
        double dm_scale = desired_pot/dm_pot;
        hists->ScaleDarkmatter(dm_scale, dm_scale/dm_detector_scale);
    }
    hists->SaveHistograms();

    delete hists;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

SensitivityAnalysis::SensitivityAnalysis(){	
	smear_mean = 0;
	smear_sigma = 0.06;
	chisqr = 0;
	xsection = 0;
	dm_detector_scale = 100;
	nu_detector_scale = 1000;
	loadedDM = false;
}

SensitivityAnalysis::~SensitivityAnalysis(){
}

void SensitivityAnalysis::Init(TFile* file, TBranch* branch) {
    mkdir("detector", S_IRWXU | S_IRWXG | S_IRWXO);
	const char* nufile = "detector/nu_energies.root";
	struct stat buffer;
    if(stat(nufile, &buffer) == 0)
    {
        nu_cache = new TFile(nufile);
        nu_cache->GetObject("nuenergy", nu_energy);
    }
    else
    {
        nu_cache = new TFile(nufile, "NEW");
        nu_energy = new TH1D("nuenergy","nuenergy",20,0.3,2);
        iterateNeutrino((int)nu_detector_scale,
        [=](Particle&, double){},
        [=](Particle&, double, double, double){},
        [=](Particle&, Particle& neutrino, Particle& electron, double weight){
            double E = electron.E;
            E += Random::Gauss(smear_mean, smear_sigma / std::sqrt(E));
            if(E>0){
                nu_energy->Fill(E, weight);
            }
        });
        nu_energy->Scale(1/nu_detector_scale);
        nu_cache->Write();
    }
    
	std::stringstream dmfile;
	dmfile << "detector/DM_" << std::fixed << std::setprecision(6) << vpmass << "_" << chimass << "_" << kappa << "_" << alpha << ".root";
	std::string sdmfile = dmfile.str();
    if(stat(sdmfile.c_str(), &buffer) == 0)
    {
        dm_cache = new TFile(sdmfile.c_str());
        dm_cache->GetObject("nth", dm_energy);
        dm_cache->GetObject("theta", theta_avg);
        dm_cache->GetObject("xsection", xsection);
        loadedDM = true;
    }
    else
    {
        dm_cache = new TFile(sdmfile.c_str(), "NEW");
        dm_energy = new TH1D("nth","nth",20,0.3,2);
	    theta_avg = new TProfile("theta","theta",20,0.3,2);
	    xsection = new TVectorD(1);
	    xsection->operator()(0) = 0;
        loadedDM = false;
    }
}

int SensitivityAnalysis::Analyze(TFile* file, TBranch* branch) {
    if(loadedDM) return 0;
    TVectorD *v;
	file->GetObject("xsection", v);
	if(v)
	    xsection->operator()(0) += v->operator()(0);
	else
		std::cout << "Root couldn't find crosssection object in file.\n";
	
	iterateDarkmatter((int)dm_detector_scale, branch, vpmass, chimass, kappa, alpha,
    [=](Particle&, double, double){},
    [=](Particle& darkmatter, Particle&, Particle& electron){
        double E = electron.E;
        E += Random::Gauss(smear_mean, smear_sigma / std::sqrt(E));
        if(E>0){
            dm_energy->Fill(E);
            theta_avg->Fill(E, electron.theta);
        }
    });
    
    return 0;
}

void SensitivityAnalysis::UnInit() {
    if(!loadedDM)
    {
        dm_cache->cd();
        xsection->Write();
        dm_cache->Write();
    }
    double desired_pot = 6e20;
    double dm_pot = std::pow(files.size(),2)/(12*xsection->operator()(0)*1e-40*100*6.022140857e23);
    double nu_pot = 25000000;
    
    dm_energy->Scale(desired_pot/dm_pot/dm_detector_scale);
    nu_energy->Scale(desired_pot/nu_pot);

    std::vector<double> N_th, N_ex, cs_mean;
    int nbins = dm_energy->GetSize()-1;
    for(int i = 1; i<nbins; ++i){
        N_th.push_back(dm_energy->GetBinContent(i)+nu_energy->GetBinContent(i));
        N_ex.push_back(nu_energy->GetBinContent(i));
        cs_mean.push_back(theta_avg->GetBinContent(i));
    }
    /*for(auto i = N_th.begin(); i!=N_th.end(); ++i)
        std::cout << *i << " ";
    std::cout << std::endl << std::endl;
    for(auto i = N_ex.begin(); i!=N_ex.end(); ++i)
        std::cout << *i << " ";
    std::cout << std::endl << std::endl;
    for(auto i = cs_mean.begin(); i!=cs_mean.end(); ++i)
        std::cout << *i << " ";
    std::cout << std::endl << std::endl;*/
    chisqr = chisq_pullfunc(N_th, N_ex, cs_mean);
    std::cout << "Chisqr is " << chisqr << std::endl;
    
    delete dm_energy;
    delete nu_energy;
    delete theta_avg;
    delete nu_cache;
    delete dm_cache;
}

SensitivityScan::SensitivityScan() {
}

int SensitivityScan::Process(const std::vector<std::string>& folders){
    std::vector<double> mixings, chisqr, masses, chisqr1d, masses1d;
    double chiscan;
    double ivpmass, ichimass, ialpha, ikappa;
    double vpmass, chimass, alpha, kappa;
    if(DMParameters(folders[0], ivpmass, ichimass, ikappa, ialpha)){
        std::cout << "Couldn't read DM parameters from " << folders[0] << std::endl;
        return 1;
    }
    
    for(auto i = folders.begin(); i!=folders.end(); ++i){
        struct stat st;
        stat(i->c_str(), &st);
        if(S_ISDIR(st.st_mode))
        {
            std::cout << *i << std::endl;
            if(DMParameters(*i, vpmass, chimass, kappa, alpha)){
                std::cout << "Couldn't read DM parameters!\n";
                continue;
            }
            SensitivityAnalysis analysis;
            if(analysis.Process(*i)){
                std::cout << "Processing failed!\n";
                continue;
            }
            if(analysis.chisqr!=analysis.chisqr)
            {
                std::cout << "Calculated a nan chisqr!\n";
                continue;
            }
            mixings.push_back(std::pow(kappa,2)*alpha*std::pow(chimass/vpmass,4));
            chisqr.push_back(analysis.chisqr);
            masses.push_back(chimass);
            if(ivpmass==vpmass && ialpha==alpha && ikappa==kappa){
                chisqr1d.push_back(analysis.chisqr);
                masses1d.push_back(chimass);
            }
        }
    }
    
    if(mixings.size()==0){
        std::cout << "No sensitivities were scanned...\n";
        return 1;
    }
    
    std::stringstream sstitle1, sstitle2;
    sstitle1 << std::setprecision(6);
    sstitle1 << "Sensitivity M_{vp}=" << ivpmass << " #kappa=" << ikappa << " #alpha=" << ialpha << ";M_{#chi} (GeV);#Chi^{2}";
    
    sstitle2 << std::setprecision(6);
    sstitle2 << "Sensitivity #alpha=" << ialpha << ";m_{#Chi} (GeV);y=#epsilon^{2}#alpha(m_{#Chi}/m_{A'})^{4}";
    
    std::string title1 = sstitle1.str();
    std::string title2 = sstitle2.str();
    
    TCanvas* canvas = new TCanvas("c1","canvas",1024,576);
    /*TGraph* graph = new TGraph(masses1d.size(),masses1d.data(),chisqr1d.data());
    graph->SetTitle(title1.c_str());
    graph->Draw("");
    canvas->Update();
    canvas->SaveAs("Sensitivity.png");*/
    canvas->SetLogx();
    canvas->SetLogy();
    canvas->SetLogz();
    
    gStyle->SetPalette(55);
    //TH2D hist("blah",title2.c_str(),1,0.1,3,1,1e-7,2e-5);
    //hist.SetStats(false);
    TGraph2D* graph2 = new TGraph2D(masses.size(),masses.data(),mixings.data(),chisqr.data());
    
    //hist.Draw("");
    graph2->Draw("CONT1Z");
    for(auto i = masses.begin(); i!=masses.end(); ++i)
    {
        for(auto o = mixings.begin(); o!=mixings.end(); ++o)
        {
            TMarker* m = new TMarker(*i, *o, 0);
            m->Draw("");
        }
    }
    canvas->Update();
    //graph2->GetZaxis()->SetRangeUser(1e-4, 1e3);
    canvas->SaveAs("SensitivityContours.png");
    
    double level = 0.9;
    graph2->GetHistogram()->SetContour(1, &level);
    //hist.Draw("");
    graph2->Draw("CONT1Z");
    for(auto i = masses.begin(); i!=masses.end(); ++i)
    {
        for(auto o = mixings.begin(); o!=mixings.end(); ++o)
        {
            TMarker* m = new TMarker(*i, *o, 0);
            m->Draw("");
        }
    }
    canvas->Update();
    //graph2->GetZaxis()->SetRangeUser(1e-4, 1e3);
    canvas->SaveAs("Sensitivity90.png");
    
    delete graph2;
    //delete graph;
    delete canvas;
}

