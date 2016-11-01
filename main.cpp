// main program
#include <iostream>
#include <iomanip> 
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
using namespace std;

int main ()
{
	
//-------------------
// Read in parameters
//-------------------
parameter par;
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
	
cout << "--------------------" << endl;	
cout << "Run Parameters:" << endl;	
cout << "--------------------" << endl;	
//cout << "Pi = " << Pi << endl;	
cout << "Minimum # of DM intersecting detector = " << NMC  << endl;	
//cout << "Pion Mass = " << Mpi << " GeV" << endl;	
//cout << "Electron Mass = " << Me << " GeV" << endl;	
cout << "Dark Photon Mass = " << MV << " GeV" << endl;	
cout << "Dark Matter Mass = " << MX << " GeV" << endl;	
//cout << "alphaEM = " << alphaEM << endl;	
cout << "alphaD = " << alphaD << endl;	
cout << "kappa = " << kappa << endl;	

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
double ne = 5.1e+23;<
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

		myfileout << "event "  << i << endl;	
		myfileout << "pion    "  << setw(15) << pion.px << setw(15) << pion.py << setw(15) << pion.pz << setw(15) << pion.E << endl;
			if (DMswitch1 == 2) 
			{
				myfileout << "DM      "  << setw(15) << darkmatter1.px << setw(15) << darkmatter1.py << setw(15) << darkmatter1.pz << setw(15) << darkmatter1.E << endl;
				myfileout << "electron"  << setw(15) << electron1.px << setw(15) << electron1.py << setw(15) << electron1.pz << setw(15) << electron1.E << endl;
			}
		if (DMswitch2 == 2) 
		{
			myfileout << "DM      "  << setw(15) << darkmatter2.px << setw(15) << darkmatter2.py << setw(15) << darkmatter2.pz << setw(15) << darkmatter2.E << endl;
			myfileout << "electron"  << setw(15) << electron2.px << setw(15) << electron2.py << setw(15) << electron2.pz << setw(15) << electron2.E << endl;
		}
		myfileout << " " << endl;
	}

}
	
myfileout.close(); //close output files
}	
else cout << "Unable to open file";	

	cout << "--------------------" << endl;	
	cout << "Run Output:" << endl;	
	cout << "--------------------" << endl;	
	
	acc = (double)NDM/Npion;
	fracDMscat = (double)Nelectron/NDM;
	probAv = fracDMscat*probMax;
		
	cout << "Number of pions = " << Npion << endl;	
	cout << "Number of DM intersecting detector = " << NDM << endl;	
	cout << "Number of electrons = " << Nelectron << endl;
	cout << "acceptance = " << acc << endl;	
	cout << "Maximum scattering probability = "  << probMax << endl;
	cout << "Average scattering probability = "  << probAv << endl;
	
	cout << "--------------------" << endl;	
	cout << "Events stored in file ``events.dat'' " << endl;	
	cout << "--------------------" << endl;	
	
return 0;
}
