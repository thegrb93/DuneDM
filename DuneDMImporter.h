#pragma once

#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"

class DuneDMImporter {
public :
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain

	// Fixed size dimensions of array or collections stored in the TTree if any.
	//  const Int_t kMaxEvent = 1;
	//  const Int_t kMaxRwgt = 1;


	// Declaration of leaf types
	Int_t           Event_;
	UInt_t          Event_fUniqueID[1];   //[Event_]
	UInt_t          Event_fBits[1];   //[Event_]
	Long64_t        Event_Number[1];   //[Event_]
	Int_t           Event_Nparticles[1];   //[Event_]
	Int_t           Event_ProcessID[1];   //[Event_]
	Double_t        Event_Weight[1];   //[Event_]
	Double_t        Event_ScalePDF[1];   //[Event_]
	Double_t        Event_CouplingQED[1];   //[Event_]
	Double_t        Event_CouplingQCD[1];   //[Event_]
	Int_t           Event_size;
	Int_t           Rwgt_;
	UInt_t          Rwgt_fUniqueID[1];   //[Rwgt_]
	UInt_t          Rwgt_fBits[1];   //[Rwgt_]
	Double_t        Rwgt_Weight[1];   //[Rwgt_]
	Int_t           Rwgt_size;
	Int_t           Particle_;
	UInt_t          Particle_fUniqueID[5];   //[Particle_]
	UInt_t          Particle_fBits[5];   //[Particle_]
	Int_t           Particle_PID[5];   //[Particle_]
	Int_t           Particle_Status[5];   //[Particle_]
	Int_t           Particle_Mother1[5];   //[Particle_]
	Int_t           Particle_Mother2[5];   //[Particle_]
	Int_t           Particle_ColorLine1[5];   //[Particle_]
	Int_t           Particle_ColorLine2[5];   //[Particle_]
	Double_t        Particle_Px[5];   //[Particle_]
	Double_t        Particle_Py[5];   //[Particle_]
	Double_t        Particle_Pz[5];   //[Particle_]
	Double_t        Particle_E[5];   //[Particle_]
	Double_t        Particle_M[5];   //[Particle_]
	Double_t        Particle_PT[5];   //[Particle_]
	Double_t        Particle_Eta[5];   //[Particle_]
	Double_t        Particle_Phi[5];   //[Particle_]
	Double_t        Particle_Rapidity[5];   //[Particle_]
	Double_t        Particle_LifeTime[5];   //[Particle_]
	Double_t        Particle_Spin[5];   //[Particle_]
	Int_t           Particle_size;

	// List of branches
	TBranch        *b_Event_;   //!
	TBranch        *b_Event_fUniqueID;   //!
	TBranch        *b_Event_fBits;   //!
	TBranch        *b_Event_Number;   //!
	TBranch        *b_Event_Nparticles;   //!
	TBranch        *b_Event_ProcessID;   //!
	TBranch        *b_Event_Weight;   //!
	TBranch        *b_Event_ScalePDF;   //!
	TBranch        *b_Event_CouplingQED;   //!
	TBranch        *b_Event_CouplingQCD;   //!
	TBranch        *b_Event_size;   //!
	TBranch        *b_Rwgt_;   //!
	TBranch        *b_Rwgt_fUniqueID;   //!
	TBranch        *b_Rwgt_fBits;   //!
	TBranch        *b_Rwgt_Weight;   //!
	TBranch        *b_Rwgt_size;   //!
	TBranch        *b_Particle_;   //!
	TBranch        *b_Particle_fUniqueID;   //!
	TBranch        *b_Particle_fBits;   //!
	TBranch        *b_Particle_PID;   //!
	TBranch        *b_Particle_Status;   //!
	TBranch        *b_Particle_Mother1;   //!
	TBranch        *b_Particle_Mother2;   //!
	TBranch        *b_Particle_ColorLine1;   //!
	TBranch        *b_Particle_ColorLine2;   //!
	TBranch        *b_Particle_Px;   //!
	TBranch        *b_Particle_Py;   //!
	TBranch        *b_Particle_Pz;   //!
	TBranch        *b_Particle_E;   //!
	TBranch        *b_Particle_M;   //!
	TBranch        *b_Particle_PT;   //!
	TBranch        *b_Particle_Eta;   //!
	TBranch        *b_Particle_Phi;   //!
	TBranch        *b_Particle_Rapidity;   //!
	TBranch        *b_Particle_LifeTime;   //!
	TBranch        *b_Particle_Spin;   //!
	TBranch        *b_Particle_size;   //!

	std::vector<double> v_Particle_Px, v_Particle_Py, v_Particle_Pz;
	std::vector<double> v_Particle_E;
	std::vector<double> v_Particle_Phi;

	DuneDMImporter();
	~DuneDMImporter();
	void Init(const char* filen);
	void ReadTree();
};

