// main program
#include <iostream>
/*#include <iomanip> 
#include <fstream>
#include <string>
#include <cctype>
#include <sstream>
#include <cmath>
#include <time.h>      
#include "Particle.h"      
#include "Random.h"
#include "Kinematics.h"
#include "sanfordwang.h"
#include "parameter.h"
#include "decay.h"
#include "detector.h"
#include "DMscattering.h"
*/
#include <sys/stat.h>
#include <vector>
#include <sstream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLegend.h>

// Header file for the classes stored in the TTree if any.
#include "ExRootClasses.h"
#include "TClonesArray.h"
#include "TObject.h"

int main (int argc, char** argv)
{
	//Read LHEF Data
	struct datavectors {
		std::vector<double> px, py, pz, theta, e;
	};
	datavectors particledata[5];
		
	{
		if(argc<2) {std::cout << "Expected input root file.\n"; return 1;}
		
		TFile *f = new TFile(argv[1]);
		if (!f || !f->IsOpen()) {
			std::cout << "Root couldn't find input file name.\n";
			return 1;
		}
		
		TTree *fChain = nullptr;
		f->GetObject("LHEF",fChain);
		if(!fChain) {
			std::cout << "Root couldn't find LHEF tree in the input root file.\n";
			return 1;
		}
		Long64_t nentries = fChain->GetEntries();	
		
		TBranch *b_Particles = fChain->GetBranch("Particle");
		if(!b_Particles) {
			std::cout << "Root couldn't find Particle branch in the LHEF tree.\n";
			return 1;
		}
	
		TClonesArray* array = new TClonesArray("TRootLHEFParticle", 5);
		b_Particles->SetAddress(&array);
		
		std::cout<<"Number of Entries " << nentries << std::endl;
		
		if(argc>2) {
			//Output histograms to folder
			mkdir("images", S_IRWXU | S_IRWXG | S_IRWXO);
			std::string params(argv[2]);
			std::string outputfilen(argv[3]);

			TCanvas* canvas = new TCanvas("output", "", 1600,900);
			canvas->SetFillColor(33);
			canvas->SetFrameFillColor(17);
			canvas->SetGrid();
			TLegend *legend=new TLegend(0.85,0.65,0.95,0.85);
			legend->SetTextFont(72);
			legend->SetTextSize(0.02);
			struct hists {
				TH1D *px, *py, *pz, *theta, *e;
			};
			hists particle_hists[5];
			std::string particle_names[5] = {"Pqd","Pqd1","Chi","Chi bar","V"}; 
			for(int i = 2; i<4; ++i)
			{
				std::string pxname = "X Momentum" + params;
				std::string pyname = "Y Momentum " + params;
				std::string pzname = "Z momentum " + params;
				std::string thetaname = "Theta " + params;
				std::string ename = "Energy " + params;
				
				std::stringstream numstr;
				numstr << i;
				
				particle_hists[i].px = new TH1D((pxname+numstr.str()).c_str(), pxname.c_str(), 100, -2.5, 2.5);
				particle_hists[i].py = new TH1D((pyname+numstr.str()).c_str(), pyname.c_str(), 100, -2.5, 2.5);
				particle_hists[i].pz = new TH1D((pzname+numstr.str()).c_str(), pzname.c_str(), 100, -70.0, 0.0);
				particle_hists[i].theta = new TH1D((thetaname+numstr.str()).c_str(), thetaname.c_str(), 100, -3.2, 3.2);
				particle_hists[i].e = new TH1D((ename+numstr.str()).c_str(), ename.c_str(), 100, 0.0, 60.0);
				
				particle_hists[i].px->SetLineColor(i);
				particle_hists[i].px->SetMarkerSize(0.8);
				particle_hists[i].px->SetStats(0);
				particle_hists[i].py->SetLineColor(i);
				particle_hists[i].py->SetMarkerSize(0.8);
				particle_hists[i].py->SetStats(0);
				particle_hists[i].pz->SetLineColor(i);
				particle_hists[i].pz->SetMarkerSize(0.8);
				particle_hists[i].pz->SetStats(0);
				particle_hists[i].theta->SetLineColor(i);
				particle_hists[i].theta->SetMarkerSize(0.8);
				particle_hists[i].theta->SetStats(0);
				particle_hists[i].e->SetLineColor(i);
				particle_hists[i].e->SetMarkerSize(0.8);
				particle_hists[i].e->SetStats(0);
			}
			
			for(Int_t i = 0; i < nentries; ++i)
			{
				fChain->GetEntry(i);
				for(Int_t j = 0; j < array->GetEntries(); ++j)
				{
					TRootLHEFParticle* particle = (TRootLHEFParticle*)array->At(j);
					hists* histograms;
					switch(particle->PID)
					{
						case 33:
						histograms = particle_hists+2;
						break;
						case -33:
						histograms = particle_hists+3;
						break;
						//case 32:
						//histograms = particle_hists+4;
						//break;
						default:
						histograms = 0;
					}
					if(histograms)
					{
						histograms->px->Fill(particle->Px);
						histograms->py->Fill(particle->Py);
						histograms->pz->Fill(particle->Pz);
						histograms->theta->Fill(particle->Phi);
						histograms->e->Fill(particle->E);
					}
				}
			}
			
			auto saveimage = [&](TH1D* histo1, TH1D* histo2, const std::string& name){
				std::stringstream filename;
				filename << "images/" << name << "_" << outputfilen << ".png";
				histo1->Draw();
				histo2->Draw("same");
				legend->Clear();
				legend->AddEntry(histo1,"Chi");
				legend->AddEntry(histo2,"Chi Bar");
				legend->Draw();
				canvas->Update();
				canvas->SaveAs(filename.str().c_str());
			};
			/*saveimage(particle_hists[0].px, "pxd");
			saveimage(particle_hists[0].py, "pyd");
			saveimage(particle_hists[0].pz, "pzd");
			saveimage(particle_hists[0].theta, "thetad");
			saveimage(particle_hists[0].e, "ed");
			saveimage(particle_hists[1].px, "pxdbar");
			saveimage(particle_hists[1].py, "pydbar");
			saveimage(particle_hists[1].pz, "pzdbar");
			saveimage(particle_hists[1].theta, "thetadbar");
			saveimage(particle_hists[1].e, "edbar");*/
			saveimage(particle_hists[2].px, particle_hists[3].px, "pxchi");
			saveimage(particle_hists[2].py, particle_hists[3].py, "pychi");
			saveimage(particle_hists[2].pz, particle_hists[3].pz, "pzchi");
			saveimage(particle_hists[2].theta, particle_hists[3].theta, "thetachi");
			saveimage(particle_hists[2].e, particle_hists[3].e, "echi");
			/*saveimage(particle_hists[4].px, "pxv");
			saveimage(particle_hists[4].py, "pyv");
			saveimage(particle_hists[4].pz, "pzv");
			saveimage(particle_hists[4].theta, "thetav");
			saveimage(particle_hists[4].e, "ev");*/
			
			for(int i = 2; i<4; ++i)
			{
				delete particle_hists[i].px;
				delete particle_hists[i].py;
				delete particle_hists[i].pz;
				delete particle_hists[i].theta;
				delete particle_hists[i].e;
			}
			
			delete canvas;
		}
		else
		{
			for(Int_t i = 0; i < nentries; ++i)
			{
				fChain->GetEntry(i);
				for(Int_t j = 0; j < array->GetEntries(); ++j)
				{
					TRootLHEFParticle* particle = (TRootLHEFParticle*)array->At(j);
					std::cout << i << " " << particle->PID << std::endl;
				}
			}
		}
		
		delete array;
		delete f;
	}
	
	
	//-------------------
	// Read in parameters
	//-------------------
	/*parameter par;
	int NMC; 
	double MV, MX, alphaD, kappa;
	double Pi, Mpi, Me, alphaEM;
	Pi = par.pi();
	NMC = par.NMCtrials();
	Mpi = par.MassPion();
	Me = par.MassElectron();
	MV = par.MassDP();
	MX = par.MassDM();
	alphaEM = par.alEM();
	alphaD = par.alD();
	kappa = par.kap();

	std::cout << "--------------------" << std::endl;	
	std::cout << "Run Parameters:" << std::endl;	
	std::cout << "--------------------" << std::endl;	
	//std::cout << "Pi = " << Pi << std::endl;	
	std::cout << "Minimum # of DM intersecting detector = " << NMC  << std::endl;	
	//std::cout << "Pion Mass = " << Mpi << " GeV" << std::endl;	
	//std::cout << "Electron Mass = " << Me << " GeV" << std::endl;	
	std::cout << "Dark Photon Mass = " << MV << " GeV" << std::endl;	
	std::cout << "Dark Matter Mass = " << MX << " GeV" << std::endl;	
	//std::cout << "alphaEM = " << alphaEM << std::endl;	
	std::cout << "alphaD = " << alphaD << std::endl;	
	std::cout << "kappa = " << kappa << std::endl;	

	//-------------
	// Declarations
	//-------------

	sanfordwang SW;	
	Kinematics kin;	
	decay dec;
	detector det;
	DMscattering scatter;	


	double ppix, ppiy, ppiz, Epi;  
	double thetaX1, thetaX2, thetaX;
	double thetacut = 0.0111113;	

	double probMax = 10e-15;	
	double ne = 5.1e+23;
	int Nscatter;	

	double EeMin, EeMax;
	double xe, ye, Thetae, Phie, Ee;
	double pe, pex, pey, pez;
	double dsig, sig, psig;
	double dsigMax, psigMax;
	double probe, Re;
	int Nelectron;
	int eswitch;

	int DMswitch1;
	int DMswitch2;

	int Npion = 0;		
	int NDM = 0;		
	double acc,fracDMscat,probAv;	
	acc = 0.0;	
	fracDMscat = 0.0;
	probAv = 0.0;

	// 13231	
	Random::Random(35459);	

	ofstream myfileout ("events.dat");

	// check if files are open	
	if (myfileout.is_open())
	{
		for (int i = 1; i < NMC; i++) 
		{
			Particle pion(Mpi);
			Particle darkphoton(MV);
			Particle darkmatter1(MX);
			Particle darkmatter2(MX);
			Particle electron1(Me);
			Particle electron2(Me);

			DMswitch1 = 0;
			DMswitch2 = 0;
			while (DMswitch1 == 0 && DMswitch2 == 0) 
			{
				//-------------------
				// Generate pion
				//-------------------

				Npion = Npion+1;	
				S/Users/sshahsav/Downloads/detector 2/main.cppW.pionGen(ppix,ppiy,ppiz,Epi);
				pion.FourMomentum(ppix,ppiy,ppiz,Epi);

				//--------------------------
				// decay pion to dark matter
				//--------------------------

				decay dec;
				dec.DecayDM(darkmatter1,darkmatter2,darkphoton,pion);		

				//--------------------	
				//	intersect detector
				//--------------------

				det.intersect(DMswitch1,NDM,darkmatter1);
				det.intersect(DMswitch2,NDM,darkmatter2);

			}	
			//-------------------------------
			// check if DM scatters
			// ------------------------------

			scatter.probscatter(DMswitch1,Nscatter,probMax,MV,MX,kappa,alphaD,darkmatter1);	
			scatter.probscatter(DMswitch2,Nscatter,probMax,MV,MX,kappa,alphaD,darkmatter2);	

			//--------
			// Scatter
			//--------	

			scatter.scatterevent(DMswitch1,Nelectron,MV,MX,kappa,alphaD,darkmatter1,electron1);
			scatter.scatterevent(DMswitch2,Nelectron,MV,MX,kappa,alphaD,darkmatter2,electron2);

			//-------
			// output
			//-------	
			if (DMswitch1 == 2 || DMswitch2 == 2) 
			{
				myfileout << "event "  << i << std::endl;	
				myfileout << "pion    "  << setw(15) << pion.px << setw(15) << pion.py << setw(15) << pion.pz << setw(15) << pion.E << std::endl;
				if (DMswitch1 == 2) 
				{
					myfileout << "DM      "  << setw(15) << darkmatter1.px << setw(15) << darkmatter1.py << setw(15) << darkmatter1.pz << setw(15) << darkmatter1.E << std::endl;
					myfileout << "electron"  << setw(15) << electron1.px << setw(15) << electron1.py << setw(15) << electron1.pz << setw(15) << electron1.E << std::endl;
				}
				if (DMswitch2 == 2) 
				{
					myfileout << "DM      "  << setw(15) << darkmatter2.px << setw(15) << darkmatter2.py << setw(15) << darkmatter2.pz << setw(15) << darkmatter2.E << std::endl;
					myfileout << "electron"  << setw(15) << electron2.px << setw(15) << electron2.py << setw(15) << electron2.pz << setw(15) << electron2.E << std::endl;
				}
				myfileout << " " << std::endl;
			}
		}
		myfileout.close(); //close output files
	}	
	else std::cout << "Unable to open file";	

	std::cout << "--------------------" << std::endl;	
	std::cout << "Run Output:" << std::endl;	
	std::cout << "--------------------" << std::endl;	

	acc = (double)NDM/Npion;
	fracDMscat = (double)Nelectron/NDM;
	probAv = fracDMscat*probMax;

	std::cout << "Number of pions = " << Npion << std::endl;	
	std::cout << "Number of DM intersecting detector = " << NDM << std::endl;	
	std::cout << "Number of electrons = " << Nelectron << std::endl;
	std::cout << "acceptance = " << acc << std::endl;	
	std::cout << "Maximum scattering probability = "  << probMax << std::endl;
	std::cout << "Average scattering probability = "  << probAv << std::endl;

	std::cout << "--------------------" << std::endl;	
	std::cout << "Events stored in file ``events.dat'' " << std::endl;	
	std::cout << "--------------------" << std::endl;
	*/

	return 0;
}

