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
#include <TGaxis.h>
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
#include "dk2nu.h"

// Detector sensitivity headers
#include "Particle.h"
#include "DUNEdet.h"
#include "DMElscattering.h"
#include "Random.h"

const double emass = 0.000511;
//double desired_pot = 1.47e21 * 3.5 * 0.5; //80GeV, 3.5 years, 50% Duty factor
double desired_pot = 1.1e21 * 3.5 * 0.5; //120GeV, 3.5 years, 50% Duty factor
double nu_pot = 50000000;

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

inline double calcDMPOT(double totalevents, double xsection){
    return totalevents/(xsection*1e-36*12*100*6.022140857e23);
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

void iterateNeutrino(
    std::function<void(Particle&, double)> fproduction,
    std::function<void(Particle&, double)> fdetector,
    std::function<void(Particle&, Particle&, Particle&, double, double)> fscatter
){
    std::string neutrino_path = options["neutrinos"].as<std::string>();
    struct stat buffer;
    if(stat(neutrino_path.c_str(), &buffer))
    {
        std::cout << "Failed to open neutrino directory. \"" << neutrino_path << "\".\n";
        return;
    }

    DUNEDetector det;
    DMscattering scatter;

    Particle neutrino(0);
    Particle electron(emass);

    std::vector<std::string> neutrino_files;
    getFolderFiles(neutrino_path, neutrino_files);
    for (auto filen = neutrino_files.begin(); filen != neutrino_files.end(); ++filen) {
        TFile *neutrinos = new TFile(filen->c_str());
        std::cout << neutrinos->GetPath() << std::endl;
        TTree *neutrino_tree = (TTree *) neutrinos->Get("dk2nuTree");
        if(!neutrino_tree){
            delete neutrinos;
            continue;
        }
        bsim::Dk2Nu* dk2nu = new bsim::Dk2Nu;
        neutrino_tree->SetBranchAddress("dk2nu",&dk2nu);

        long long neutrino_entries = neutrino_tree->GetEntries();
        // int total_attempts = 0, success_attempts = 0;
        for (long long i = 0; i < neutrino_entries; ++i) {
            neutrino_tree->GetEvent(i);

            if(dk2nu->decay.ntype!=12) continue;

            double nimpwt = dk2nu->decay.nimpwt;

            bsim::NuRay prod_ray = dk2nu->nuray[0];
            neutrino.FourMomentum(prod_ray.px, prod_ray.py, prod_ray.pz, prod_ray.E);
            neutrino.calcOptionalKinematics();
            fproduction(neutrino, nimpwt * prod_ray.wgt);

            bsim::NuRay det_ray = dk2nu->nuray[1];
            neutrino.FourMomentum(det_ray.px, det_ray.py, det_ray.pz, det_ray.E);
            neutrino.calcOptionalKinematics();

            bsim::Ancestor& parent_info = dk2nu->ancestor[1];
            bsim::Ancestor& neutrino_info = dk2nu->ancestor[dk2nu->ancestor.size()-1];
            double tmin, tmax;
            double nu_x = neutrino_info.startx/100.0, nu_y = neutrino_info.starty/100.0, nu_z = neutrino_info.startz/100.0;
            if(det.intersect(nu_x, nu_y, nu_z, neutrino.normpx, neutrino.normpy, neutrino.normpz, tmin, tmax))
            {
                fdetector(neutrino, nimpwt * det_ray.wgt);
                if (scatter.probscatterNeutrino(neutrino, tmax - tmin)) {
                    double distance = Random::Flat(tmin,tmax);
                    double time = (neutrino_info.startt - parent_info.startt) + distance / ( neutrino.norm / neutrino.E * 0.299792458 ); //neutrino creation time + travel time
                    Particle neutrino2 = neutrino;
                    scatter.scattereventNeutrino(neutrino2, electron);
                    neutrino2.calcOptionalKinematics();
                    electron.calcOptionalKinematics();
                    fscatter(neutrino, neutrino2, electron, time, nimpwt * det_ray.wgt * scatter.pMax0);
                }
            }
        }
        std::cout << std::endl;

        delete dk2nu;
        delete neutrino_tree;
        delete neutrinos;
    }
}

void iterateDarkmatter(TBranch* branch, double vpmass, double chimass, double kappa, double alpha,
    std::function<void(Particle&, double, double)> fdetector,
    std::function<void(Particle&, Particle&, Particle&, double, double)> fscatter
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
    auto doIterations = [&] (int index) {
        for(Int_t i = 0; i < nentries; ++i)
        {
            branch->GetEntry(i);
            TRootLHEFParticle* lhefdarkmatter = (TRootLHEFParticle*)array->At(index);
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
                if(scatter.probscatter(vpmass,chimass,kappa,alpha,darkmatter,pathlength)) {
                    double time = Random::Flat(tmin,tmax) / ( darkmatter.norm / darkmatter.E * 0.299792458 );
                    Particle darkmatter2 = darkmatter;
                    scatter.scatterevent(vpmass, chimass, kappa, alpha, darkmatter2, electron);
                    darkmatter2.calcOptionalKinematics();
                    electron.calcOptionalKinematics();
                    fscatter(darkmatter, darkmatter2, electron, time, scatter.pMax0);
                }
            }
        }
    };
    doIterations(darkmatter_index);
    doIterations(antidarkmatter_index);
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

    std::cout << std::endl;
    bool first = true;
    for(int i = 0; i<files.size(); ++i) {
        TFile file(files[i].c_str());
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
            if(Init(&file, branch))
                return 1;
        }
        if(Analyze(&file, branch))
            continue;
        std::cout << "\r(" << (i+1) << " / " << files.size() << ") completed..." << std::flush;
    }
    std::cout << std::endl;

    UnInit();
    return 0;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
DetectorAnalysis::DetectorAnalysis(){
    xsection = 0;
    totalevents = 0;
}

DetectorAnalysis::~DetectorAnalysis(){
}

int DetectorAnalysis::Init(TFile* file, TBranch* branch) {
    //if(gOptions[OPT_DETECTOR])
    //  detectorType = gOptions[OPT_DETECTOR].last()->arg;
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

    std::stringstream ssrunname;
    ssrunname << std::fixed << std::setprecision(6) << "DM_" << vpmass << "_" << chimass << "_" << kappa << "_" << alpha;
    std::string runname = ssrunname.str();

    std::string neutrino_path = options["neutrinos"].as<std::string>();
    struct stat buffer;
    if(stat(neutrino_path.c_str(), &buffer))
    {
        std::cout << "Failed to open neutrino directory. \"" << neutrino_path << "\".\n";
        return 1;
    }
    std::string neutrino_root = neutrino_path+"/neutrinos.root";
    bool neutrino_exists = (stat(neutrino_root.c_str(), &buffer) == 0);

    mkdir("histograms", S_IRWXU | S_IRWXG | S_IRWXO);
    mkdir("root", S_IRWXU | S_IRWXG | S_IRWXO);
    mkdir(("histograms/"+runname).c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    TFile* dm_output = new TFile(("root/"+runname+".root").c_str(), "RECREATE");
    TFile* nu_output;
    if(neutrino_exists)
        nu_output = new TFile(neutrino_root.c_str(), "READ");
    else
        nu_output = new TFile(neutrino_root.c_str(), "CREATE");

    hists = new DMHistograms(runname, dm_output, nu_output, neutrino_exists);

    if(!neutrino_exists){
        iterateNeutrino(
        [=](Particle& neutrino, double weight){
            hists->AddProductionNu(neutrino.E, neutrino.pz, neutrino.pt, neutrino.theta, /*ntype,*/ weight);
        },
        [=](Particle& neutrino, double weight){
            hists->AddDetectorNu(neutrino.E, neutrino.pz, neutrino.pt, neutrino.theta, weight);
        },
        [=](Particle& neutrino, Particle&, Particle& electron, double time, double weight){
            hists->AddScatterNu(neutrino.E, neutrino.pz, neutrino.pt, neutrino.theta, weight);
            hists->AddScatterBgElectron(electron.E, electron.pz, electron.pt, electron.theta, time, weight);
        });
        hists->SaveNeutrinos();
    }

    return 0;
}

int DetectorAnalysis::Analyze(TFile* file, TBranch* branch) {
    {
        TVectorD *v;
        TH1D* E, *Px, *Py, *Pz, *Th, *Phi;
        file->GetObject("energies", E);
        file->GetObject("xmomentum", Px);
        file->GetObject("ymomentum", Py);
        file->GetObject("zmomentum", Pz);
        file->GetObject("theta", Th);
        file->GetObject("phi", Phi);
        if(E) hists->AddProductionDM(E, Px, Py, Pz, Th, Phi);
        else std::cout << "Root couldn't find histogram objects in file.\n";

        file->GetObject("xsection", v);
        if(v){
            xsection += v->operator()(0) / (double)files.size(); //Divide by # of files to average
            totalevents += E->GetEntries() / 2.0; //Divide by two because anti-darkmatter is included and we want # of events
        }
        else
            std::cout << "Root couldn't find crosssection object in file.\n";
    }

    iterateDarkmatter(branch, vpmass, chimass, kappa, alpha,
    [=](Particle& darkmatter, double tmin, double tmax){
        hists->AddDetectorDM(darkmatter.E, darkmatter.normpx*tmin, darkmatter.normpy*tmin, darkmatter.px, darkmatter.py, darkmatter.pz, darkmatter.pt, darkmatter.theta, darkmatter.phi);
    },
    [=](Particle& darkmatter, Particle&, Particle& electron, double time, double weight){
        hists->AddScatterDM(darkmatter.E, darkmatter.pz, darkmatter.pt, darkmatter.theta, weight);
        hists->AddScatterSigElectron(electron.E, electron.pz, electron.pt, electron.theta, time, weight);
    });
    return 0;
}

void DetectorAnalysis::UnInit() {
    if(totalevents!=0)
    {
        if(normalize_histos){
            hists->NormalizeHistograms();
        }
        else{
            if(xsection!=0){
                double dm_pot = calcDMPOT(totalevents, xsection);
                double dm_scale = desired_pot/dm_pot;
                std::cout << "DM Scale: " << dm_scale << std::endl;
                hists->ScaleDarkmatter(dm_scale);
            }
            double nu_scale = desired_pot/nu_pot;
            hists->ScaleNeutrinos(nu_scale);
        }
        hists->SaveHistograms();
    }

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
    detector_efficiency = 0.9;
    nbins = 20;
    minbin = 0.3;
    maxbin = 15;
    binstart_max = 0;
    binend_max = 0;
    loadedDM = false;
}

SensitivityAnalysis::~SensitivityAnalysis(){
}

int SensitivityAnalysis::Init(TFile* file, TBranch* branch) {
    mkdir("detector", S_IRWXU | S_IRWXG | S_IRWXO);

    std::string neutrino_path = options["neutrinos"].as<std::string>();
    struct stat buffer;
    if(stat(neutrino_path.c_str(), &buffer))
    {
        std::cout << "Failed to open neutrino directory. \"" << neutrino_path << "\".\n";
        return 1;
    }
    std::string neutrino_root = neutrino_path+"/neutrinos_detector.root";
    bool neutrino_exists = (stat(neutrino_root.c_str(), &buffer) == 0);

    if(neutrino_exists)
    {
        nu_cache = new TFile(neutrino_root.c_str());
        nu_cache->GetObject("nuenergy", nu_energy);
    }
    else
    {
        nu_cache = new TFile(neutrino_root.c_str(), "NEW");
        nu_energy = new TH1D("nuenergy","nuenergy",nbins,minbin,maxbin);
        iterateNeutrino(
        [=](Particle&, double){},
        [=](Particle&, double){},
        [=](Particle&, Particle& neutrino, Particle& electron, double time, double weight){
            double E = electron.E;
            E += Random::Gauss(smear_mean, smear_sigma / std::sqrt(E));
            if(E>0){
                nu_energy->Fill(E, weight);
            }
        });
        nu_cache->Write();
    }

    std::stringstream dmfile;
    dmfile << "detector/DM_" << std::fixed << std::setprecision(6) << vpmass << "_" << chimass << "_" << kappa << "_" << alpha << ".root";
    std::string sdmfile = dmfile.str();
    if(stat(sdmfile.c_str(), &buffer) == 0)
    {
        std::cout << "Loading Precached: " << sdmfile << std::endl;
        dm_cache = new TFile(sdmfile.c_str());
        dm_cache->GetObject("nth", dm_energy);
        dm_cache->GetObject("theta", theta_avg);
        dm_cache->GetObject("TVectorT<double>", xsection);
        loadedDM = true;
    }
    else
    {
        std::cout << "Creating: " << sdmfile << std::endl;
        dm_cache = new TFile(sdmfile.c_str(), "NEW");
        dm_energy = new TH1D("nth","nth",nbins,minbin,maxbin);
        theta_avg = new TProfile("theta","theta",nbins,minbin,maxbin);
        xsection = new TVectorD(2);
        xsection->operator()(0) = 0;
        xsection->operator()(1) = 0;
        loadedDM = false;
    }

    return 0;
}

int SensitivityAnalysis::Analyze(TFile* file, TBranch* branch) {
    if(loadedDM) return 0;
    TVectorD *v;
    TH1D* E;
    file->GetObject("energies", E);
    file->GetObject("xsection", v);
    if(v){
        xsection->operator()(0) += v->operator()(0) / (double)files.size(); //Divide by # of files to average
        xsection->operator()(1) += E->GetEntries() / 2.0; //Divide by two because anti-darkmatter is included and we want # of events
    }
    else
        std::cout << "Root couldn't find crosssection object in file.\n";

    iterateDarkmatter(branch, vpmass, chimass, kappa, alpha,
    [=](Particle&, double, double){},
    [=](Particle& darkmatter, Particle&, Particle& electron, double time, double weight){
        double E = electron.E;
        E += Random::Gauss(smear_mean, smear_sigma / std::sqrt(E));
        if(E>0){
            dm_energy->Fill(E, weight);
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
    double sigma = xsection->operator()(0);
    double Nx = xsection->operator()(1);
    double dm_pot = calcDMPOT(Nx, sigma);

    dm_energy->Scale(desired_pot/dm_pot);
    nu_energy->Scale(desired_pot/nu_pot);

    int nbins = dm_energy->GetSize()-1;

    int binstart = 1, binend = nbins;
    // for(int binstart = 1; binstart < nbins; ++binstart){
        // for(int binend = binstart+1; binend<=nbins; ++binend){
            std::vector<double> N_th, N_ex, cs_mean;
            for(int i = binstart; i<binend; ++i){
                double N_x = dm_energy->GetBinContent(i);
                double N_v = nu_energy->GetBinContent(i);
                if(N_v > 0){
                    N_th.push_back((N_x+N_v) * detector_efficiency);
                    N_ex.push_back(N_v * detector_efficiency);
                    cs_mean.push_back(theta_avg->GetBinContent(i));
                }
            }
            double newchisqr = chisq_pullfunc(N_th, N_ex, cs_mean);
            // if(newchisqr > chisqr){
                chisqr = newchisqr;
                binstart_max = dm_energy->GetBinLowEdge(binstart);
                binend_max = dm_energy->GetBinLowEdge(binend)+dm_energy->GetBinWidth(binend);
            // }
        // }
    // }

    std::cout << std::scientific << "Chisqr is " << chisqr << ". Cut from " << std::fixed << binstart_max << "(GeV) to " << binend_max << "(GeV)" << std::endl;

    delete xsection;
    delete dm_energy;
    delete nu_energy;
    delete theta_avg;
    delete nu_cache;
    delete dm_cache;
}

SensitivityScan::SensitivityScan() {
}

int SensitivityScan::Process(const std::vector<std::string>& folders){
    std::vector<double> masses, kappas, mixings, sigmas, chisqr;
    double vpmass, chimass, alpha, kappa;
    std::ofstream sensitivity_output("sensitivity.txt");

    for(auto i = folders.begin(); i!=folders.end(); ++i){
        struct stat st;
        stat(i->c_str(), &st);
        if(S_ISDIR(st.st_mode))
        {
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
            masses.push_back(chimass);
            kappas.push_back(kappa);
            chisqr.push_back(analysis.chisqr);

            double mixing = std::pow(kappa,2)*alpha*std::pow(chimass/vpmass,4);
            mixings.push_back(mixing);
            double sigma = 3.0*mixing/(299792458.0 * 137.035999074 * chimass * chimass);
            sigmas.push_back(sigma);

            sensitivity_output << vpmass << "," << chimass << "," << kappa << "," << alpha << "," << mixing << "," << sigma << "," << analysis.chisqr << "," << analysis.binstart_max << "," << analysis.binend_max << std::endl;
        }
    }

    if(masses.size()==0){
        std::cout << "No sensitivities were scanned...\n";
        return 1;
    }

    {
        std::vector<TMarker*> markers;
        TCanvas* canvas = new TCanvas("c1","canvas",1024,576);
        gStyle->SetPalette(55);
        TGraph2D* graph = new TGraph2D(masses.size());

        for(int i = 0; i<masses.size(); ++i)
        {
            double mass = std::log10(masses[i]), kappa = std::log10(kappas[i]);
            graph->SetPoint(i, mass, kappa, chisqr[i]);
        }

        graph->Draw("CONT1Z");
        canvas->Update();
        {
            std::ofstream curve("curve.txt");
            TList* list = graph->GetContourList(0.9);
            if(!list) std::cout << "Contour list is NULL!" << std::endl;
            TIter next(list);
            TGraph* obj;
            while(obj = (TGraph*)next()){
                bool waslast = false;
                double used_length = 0;
                double lastx, lasty;
                for(int i = 0; i<obj->GetN(); ++i){
                    double x, y;
                    obj->GetPoint(i, x, y);
                    if(waslast){
                        const double spacing = 0.09;

                        double dirx = x - lastx, diry = y - lasty;
                        double length = std::sqrt(dirx*dirx + diry*diry);
                        double normx = dirx / length;
                        double normy = diry / length;

                        while(used_length < length)
                        {
                            double newx = lastx + normx*used_length;
                            double newy = lasty + normy*used_length;
                            double perpx1 = newx - normy*spacing;
                            double perpy1 = newy + normx*spacing;
                            double perpx2 = newx + normy*spacing;
                            double perpy2 = newy - normx*spacing;

                            TMarker* m;
                            m = new TMarker(newx, newy, 0);
                            m->SetMarkerStyle(7);
                            markers.push_back(m);

                            m = new TMarker(perpx1, perpy1, 0);
                            m->SetMarkerStyle(7);
                            markers.push_back(m);

                            m = new TMarker(perpx2, perpy2, 0);
                            m->SetMarkerStyle(7);
                            markers.push_back(m);

                            curve << std::pow(10, newx) << " " << std::pow(10, newy) << std::endl;
                            curve << std::pow(10, perpx1) << " " << std::pow(10, perpy1) << std::endl;
                            curve << std::pow(10, perpx2) << " " << std::pow(10, perpy2) << std::endl;

                            used_length += spacing;
                        }
                        used_length -= length;
                    }
                    else
                        waslast = true;
                    lastx = x; lasty = y;
                }
            }
        }
        double level = 0.9;
        graph->GetHistogram()->SetContour(1, &level);
        graph->Draw("CONT1Z");
        for(auto i = markers.begin(); i!=markers.end(); ++i)
            (*i)->Draw();
        canvas->SaveAs("Curve.png");

        for(auto i = markers.begin(); i!=markers.end(); ++i)
            delete *i;
        delete graph;
        delete canvas;
    }
    {
        std::vector<TMarker*> markers;
        std::stringstream sstitle;
        sstitle << std::setprecision(6);
        sstitle << "90% Confidence Level - 1.925#times10^{21} POT at 120GeV, m_{V} = 3m_{#chi}, #alpha=" << alpha << ";m_{#chi} (GeV);Y=#kappa^{2}#alpha(m_{#chi}/m_{V})^{4}";
        std::string title = sstitle.str();

        TCanvas* canvas = new TCanvas("c1","canvas",1024,576);
        // canvas->SetFrameLineColor(0);
        canvas->SetLogz();
        gStyle->SetPalette(55);
        //gStyle->SetNumberContours(100);

        TGraph2D* graph = new TGraph2D(masses.size());
        graph->SetTitle(title.c_str());

        for(int i = 0; i<masses.size(); ++i)
        {
            double mass = std::log10(masses[i]), mixing = std::log10(mixings[i]);
            graph->SetPoint(i, mass, mixing, chisqr[i]);

            TMarker* m = new TMarker(mass, mixing, 0);
            m->SetMarkerStyle(7);
            markers.push_back(m);
        }

        TGaxis* xLogaxis;
        TGaxis* yLogaxis;
        {
            double minx = graph->GetXmin();
            double maxx = graph->GetXmax();
            double miny = graph->GetYmin();
            double maxy = graph->GetYmax();
            xLogaxis = new TGaxis(minx, miny, maxx, miny, *std::min_element(masses.begin(), masses.end()), *std::max_element(masses.begin(), masses.end()), 510, "G");
            yLogaxis = new TGaxis(minx, miny, minx, maxy, *std::min_element(mixings.begin(), mixings.end()), *std::max_element(mixings.begin(), mixings.end()), 510, "G");
        }
        //graph->SetNpx(500);
        //graph->SetNpy(500);

        //graph->Draw("col box scat");
        //graph->Draw("SAME CONT1Z");

        graph->Draw("CONT1Z");
        xLogaxis->Draw();
        yLogaxis->Draw();

        for(auto i = markers.begin(); i!=markers.end(); ++i)
            (*i)->Draw();

        graph->GetHistogram()->GetXaxis()->SetLabelSize(0);
        graph->GetHistogram()->GetXaxis()->SetTickLength(0);
        graph->GetHistogram()->GetYaxis()->SetLabelSize(0);
        graph->GetHistogram()->GetYaxis()->SetTickLength(0);
        graph->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);

        canvas->Update();
        canvas->SaveAs("SensitivityMixing.png");

        double level = 0.9;
        graph->GetHistogram()->SetContour(1, &level);

        //graph->Draw("col box scat");
        //graph->Draw("SAME CONT1Z");

        graph->Draw("CONT1Z");
        xLogaxis->Draw();
        yLogaxis->Draw();

        // for(auto i = markers.begin(); i!=markers.end(); ++i)
            // (*i)->Draw();

        graph->GetHistogram()->GetZaxis()->SetLabelOffset(999);
        graph->GetHistogram()->GetZaxis()->SetLabelSize(0);
        graph->GetHistogram()->GetZaxis()->SetTickLength(0);

        canvas->Update();
        canvas->SaveAs("SensitivityMixing90.png");

        for(auto i = markers.begin(); i!=markers.end(); ++i)
            delete *i;
        delete xLogaxis;
        delete yLogaxis;
        delete graph;
        delete canvas;
    }
    {
        std::vector<TMarker*> markers;
        std::stringstream sstitle;
        sstitle << std::setprecision(6);
        sstitle << "90% Confidence Level - 1.925#times10^{21} POT at 120GeV, m_{V} = 3m_{#chi}, #alpha=" << alpha << ";m_{#chi} (GeV);#sigma=3Y/(c#alpha_{e}m_{#chi}^{2})";
        std::string title = sstitle.str();

        TCanvas* canvas = new TCanvas("c1","canvas",1024,576);
        // canvas->SetFrameLineColor(0);
        canvas->SetLogz();
        gStyle->SetPalette(55);
        //gStyle->SetNumberContours(100);

        TGraph2D* graph = new TGraph2D(masses.size());
        graph->SetTitle(title.c_str());

        for(int i = 0; i<masses.size(); ++i)
        {
            double mass = std::log10(masses[i]), sigma = std::log10(sigmas[i]);
            graph->SetPoint(i, mass, sigma, chisqr[i]);

            TMarker* m = new TMarker(mass, sigma, 0);
            m->SetMarkerStyle(7);
            markers.push_back(m);
        }

        TGaxis* xLogaxis;
        TGaxis* yLogaxis;
        {
            double minx = graph->GetXmin();
            double maxx = graph->GetXmax();
            double miny = graph->GetYmin();
            double maxy = graph->GetYmax();
            xLogaxis = new TGaxis(minx, miny, maxx, miny, *std::min_element(masses.begin(), masses.end()), *std::max_element(masses.begin(), masses.end()), 510, "G");
            yLogaxis = new TGaxis(minx, miny, minx, maxy, *std::min_element(sigmas.begin(), sigmas.end()), *std::max_element(sigmas.begin(), sigmas.end()), 510, "G");
        }
        //graph->SetNpx(500);
        //graph->SetNpy(500);

        //graph->Draw("col box scat");
        //graph->Draw("SAME CONT1Z");

        graph->Draw("CONT1Z");
        xLogaxis->Draw();
        yLogaxis->Draw();

        for(auto i = markers.begin(); i!=markers.end(); ++i)
            (*i)->Draw();

        graph->GetHistogram()->GetXaxis()->SetLabelSize(0);
        graph->GetHistogram()->GetXaxis()->SetTickLength(0);
        graph->GetHistogram()->GetYaxis()->SetLabelSize(0);
        graph->GetHistogram()->GetYaxis()->SetTickLength(0);
        graph->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);

        canvas->Update();
        canvas->SaveAs("SensitivitySigma.png");

        double level = 0.9;
        graph->GetHistogram()->SetContour(1, &level);

        //graph->Draw("col box scat");
        //graph->Draw("SAME CONT1Z");

        graph->Draw("CONT1Z");
        xLogaxis->Draw();
        yLogaxis->Draw();

        // for(auto i = markers.begin(); i!=markers.end(); ++i)
            // (*i)->Draw();

        graph->GetHistogram()->GetZaxis()->SetLabelOffset(999);
        graph->GetHistogram()->GetZaxis()->SetLabelSize(0);
        graph->GetHistogram()->GetZaxis()->SetTickLength(0);

        canvas->Update();
        canvas->SaveAs("SensitivitySigma90.png");

        for(auto i = markers.begin(); i!=markers.end(); ++i)
            delete *i;
        delete xLogaxis;
        delete yLogaxis;
        delete graph;
        delete canvas;
    }
}

