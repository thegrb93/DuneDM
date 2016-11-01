#include "DuneDMImporter.h"

void DuneDMImporter::ReadTree()
{
	if (!fChain) return;

	Long64_t nentries = fChain->GetEntriesFast();
	for (Long64_t jentry=0; jentry<nentries; ++jentry) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		
		//---------------- Looping over all particle ----------------------------------

		std::cout<<"Number of Particles"<< Particle_ <<std::endl;

		for(Int_t g4 =0; g4<Particle_; g4++)
		{
			if(Particle_PID[g4] ==33)
			{


			}
		}
	}
}


DuneDMImporter::DuneDMImporter(const char* filen) : fChain(0) 
{
}

DuneDMImporter::~DuneDMImporter()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DuneDMImporter::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DuneDMImporter::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

int DuneDMImporter::Init(const char* filen)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).
	TFile *f = new TFile(filen);
	if (!f || !f->IsOpen()) {
		std::cout << "Root couldn't find input file name.\n";
		return 1;
	}
	f->GetObject("LHEF",fChain);

	
	if(!fChain) {
		std::cout << "Root couldn't find LHEF branch in the input root file.\n";
		return 1;
	}
	
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("Event", &Event_, &b_Event_);
	fChain->SetBranchAddress("Event.fUniqueID", Event_fUniqueID, &b_Event_fUniqueID);
	fChain->SetBranchAddress("Event.fBits", Event_fBits, &b_Event_fBits);
	fChain->SetBranchAddress("Event.Number", Event_Number, &b_Event_Number);
	fChain->SetBranchAddress("Event.Nparticles", Event_Nparticles, &b_Event_Nparticles);
	fChain->SetBranchAddress("Event.ProcessID", Event_ProcessID, &b_Event_ProcessID);
	fChain->SetBranchAddress("Event.Weight", Event_Weight, &b_Event_Weight);
	fChain->SetBranchAddress("Event.ScalePDF", Event_ScalePDF, &b_Event_ScalePDF);
	fChain->SetBranchAddress("Event.CouplingQED", Event_CouplingQED, &b_Event_CouplingQED);
	fChain->SetBranchAddress("Event.CouplingQCD", Event_CouplingQCD, &b_Event_CouplingQCD);
	fChain->SetBranchAddress("Event_size", &Event_size, &b_Event_size);
	fChain->SetBranchAddress("Rwgt", &Rwgt_, &b_Rwgt_);
	fChain->SetBranchAddress("Rwgt.fUniqueID", &Rwgt_fUniqueID, &b_Rwgt_fUniqueID);
	fChain->SetBranchAddress("Rwgt.fBits", &Rwgt_fBits, &b_Rwgt_fBits);
	fChain->SetBranchAddress("Rwgt.Weight", &Rwgt_Weight, &b_Rwgt_Weight);
	fChain->SetBranchAddress("Rwgt_size", &Rwgt_size, &b_Rwgt_size);
	fChain->SetBranchAddress("Particle", &Particle_, &b_Particle_);
	fChain->SetBranchAddress("Particle.fUniqueID", Particle_fUniqueID, &b_Particle_fUniqueID);
	fChain->SetBranchAddress("Particle.fBits", Particle_fBits, &b_Particle_fBits);
	fChain->SetBranchAddress("Particle.PID", Particle_PID, &b_Particle_PID);
	fChain->SetBranchAddress("Particle.Status", Particle_Status, &b_Particle_Status);
	fChain->SetBranchAddress("Particle.Mother1", Particle_Mother1, &b_Particle_Mother1);
	fChain->SetBranchAddress("Particle.Mother2", Particle_Mother2, &b_Particle_Mother2);
	fChain->SetBranchAddress("Particle.ColorLine1", Particle_ColorLine1, &b_Particle_ColorLine1);
	fChain->SetBranchAddress("Particle.ColorLine2", Particle_ColorLine2, &b_Particle_ColorLine2);
	fChain->SetBranchAddress("Particle.Px", Particle_Px, &b_Particle_Px);
	fChain->SetBranchAddress("Particle.Py", Particle_Py, &b_Particle_Py);
	fChain->SetBranchAddress("Particle.Pz", Particle_Pz, &b_Particle_Pz);
	fChain->SetBranchAddress("Particle.E", Particle_E, &b_Particle_E);
	fChain->SetBranchAddress("Particle.M", Particle_M, &b_Particle_M);
	fChain->SetBranchAddress("Particle.PT", Particle_PT, &b_Particle_PT);
	fChain->SetBranchAddress("Particle.Eta", Particle_Eta, &b_Particle_Eta);
	fChain->SetBranchAddress("Particle.Phi", Particle_Phi, &b_Particle_Phi);
	fChain->SetBranchAddress("Particle.Rapidity", Particle_Rapidity, &b_Particle_Rapidity);
	fChain->SetBranchAddress("Particle.LifeTime", Particle_LifeTime, &b_Particle_LifeTime);
	fChain->SetBranchAddress("Particle.Spin", Particle_Spin, &b_Particle_Spin);
	fChain->SetBranchAddress("Particle_size", &Particle_size, &b_Particle_size);
}

