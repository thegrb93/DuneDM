
//############ Main programme               ############//
//############ Written By Animesh Chatterjee#########//
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <string>
#include <cctype>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <time.h>      
#include "Particle.h"      
#include "Random.h"
#include "Kinematics.h"
#include "DUNEdet.h"
#include "DMElscattering.h"
using namespace std;

int main()
{
	//-------------------
	// Read in parameters
	//-------------------

	double MV, MX, alphaD, kappa;
	double Pi, Me, alphaEM;
	double a,b,c,d,e;

	int DMswitch1;
	int DMswitch2;

	MV= 4.00;
	MX=2.25;
	Me=0.000511;
	Pi = 3.14159265;
	alphaD = 0.1;
	kappa = 0.001;

	double probMax = 10e-15;	
	double ne = 5.1e+23;
	int Nscatter;	
	int Nelectron;

	//-------------
	// Declarations
	//-------------


	Kinematics kin;	
	detector det;
	DMscattering scatter;	


	//---------------
	// Input file
	//---------------

	ifstream in;
	in.open("outputfile.dat");

	ofstream out,out1;
	out.open("ELXSoutput1.dat");
	out1.open("ELXSoutput2.dat");
	int NDM = 0;
	int n=0;
	int num=0;
	while(!in.eof())
	{
		in>>a>>b>>c>>d>>e;

		//---------- Particle mass is defined here---------------------//

		Particle darkphoton(MV);
		Particle darkmatter1(MX);
		Particle darkmatter2(MX);
		Particle electron1(Me);
		Particle electron2(Me);

		DMswitch1 = 0;
		DMswitch2 = 0;

		if (d<0)
		{
			darkmatter1.FourMomentum(b,c,-1*d,e);
			darkmatter2.FourMomentum(b,c,-1*d,e);




			//--------------------	
			//   intersect detector1
			//--------------------

			det.intersect(DMswitch1,NDM,darkmatter1);
			//det.intersect(DMswitch2,NDM,darkmatter2);


			//-------------------------------
			// check if DM scatters
			// ------------------------------

			scatter.probscatter(DMswitch1,Nscatter,probMax,MV,MX,kappa,alphaD,darkmatter1);	
			//scatter.probscatter(DMswitch2,Nscatter,probMax,MV,MX,kappa,alphaD,darkmatter2);	


			//cout<<"DMswitch"<<DMswitch1<<endl;
			//---------------------------------
			// Scattering of DM-Electron
			//---------------------------------	

			scatter.scatterevent(DMswitch1,Nelectron,MV,MX,kappa,alphaD,darkmatter1,electron1);
			//scatter.scatterevent(DMswitch2,Nelectron,MV,MX,kappa,alphaD,darkmatter2,electron2);

			//-----------------------------------------------
			//---------------- Output File -------------------
			//------------------------------------------------

			if (DMswitch1 == 2 ) 
			{


				out<<darkmatter1.pz<<"\t"<<darkmatter1.E<<"\t"<<electron1.pz<<"\t"<<electron1.E<<endl;         
				num++;

			}
		} 
		n++;
	}

	cout<<"number of DM particle"<<n<<"# DM scatters"<<NDM<<"Number of electron interacted"<<num<<endl;
	in.close();
	out.close();
	out1.close();
	return 0;
}

