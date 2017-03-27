#include <cmath>

#include "DUNEdet.h"
#include "Particle.h"
#include "Kinematics.h"

// polar angle measured from center of detector at which DM enters
double MiniBoonDetector::thetaenter (double theta) {
	double rthetaenter;
	double Rdet = 6.4;
	double ddet = 500.0;
	double Adet, Bdet, Cdet, rad, root;
	Adet = Rdet*Rdet;
	Bdet = -2*ddet*Rdet*sin(theta)*sin(theta);
	Cdet = ddet*ddet*sin(theta)*sin(theta)-Rdet*Rdet*cos(theta)*cos(theta);
	rad  = Bdet*Bdet-4*Adet*Cdet;
        if(rad >0)
        {
	root = (-Bdet+sqrt(rad))/(2*Adet);
	rthetaenter = acos(root);}
        else
        rthetaenter=0.0; 
        
	return(rthetaenter);
}
// polar angle measured from center of detector at which DM exits
double MiniBoonDetector::thetaexit (double theta) {
	double rthetaexit;
	double Rdet = 6.4;
	double ddet = 500.0;
	double Adet, Bdet, Cdet, rad, root;
	Adet = Rdet*Rdet;
	Bdet = -2*ddet*Rdet*sin(theta)*sin(theta);
	Cdet = ddet*ddet*sin(theta)*sin(theta)-Rdet*Rdet*cos(theta)*cos(theta);
	rad  = Bdet*Bdet-4*Adet*Cdet;
        if(rad>0)
        {
	root = (-Bdet-sqrt(rad))/(2*Adet);
	rthetaexit = acos(root);
        }
        else
        rthetaexit=0.0;
	return(rthetaexit);
}
// distance from target to DM entrance point of detector
double MiniBoonDetector::Lenter (double theta) {
	double rLenter;
	double Rdet = 6.4;
	double ddet = 500.0;
	double thenter = thetaenter(theta);
	rLenter = Rdet*sin(thenter)/sin(theta); 
	return(rLenter);
}
// distance from target DM exit point of detector
double MiniBoonDetector::Lexit (double theta) {
	double rLexit;
	double Rdet = 6.4;
	double ddet = 500.0;
	double thexit = thetaexit(theta);
	rLexit = Rdet*sin(thexit)/sin(theta); 
	return(rLexit);
}
// distance DM travels through detector
double MiniBoonDetector::Ldet (double theta) {
	double rLdet;
	double Ldetenter, Ldetexit;
	Ldetenter = Lenter(theta); 
	Ldetexit = Lexit(theta); 
        //std::cout<<"Lenter"<<Ldetenter<<"Lexit"<<Ldetexit<<std::endl;
	rLdet = Ldetexit - Ldetenter;
	return(rLdet);
}
//
void MiniBoonDetector::intersect(int &dswitch, int &Ndm, Particle &DM){
	// compute polar angle
	double thetacut = 0.0128;	
	double polarX;
	Kinematics kin;
        
	polarX = kin.theta(DM.px,DM.py,DM.pz,DM.E);
       // std::cout<<"Polar X is"<<polarX<<std::endl;
	// check if dark matter intersects detector
	if (polarX < thetacut) 
	{
		dswitch = 1;
		Ndm = Ndm+1;
	}		
}
