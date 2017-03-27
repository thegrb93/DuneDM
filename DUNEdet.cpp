#include <cmath>

#include "DUNEdet.h"
#include "Particle.h"
#include "Kinematics.h"

bool intersectAABB(float x, float y, float z, float dirx, float diry, float dirz, float minsx, float minsy, float minsz, float maxsx, float maxsy, float maxsz, float &t) const
{
	float invdirx = 1/dirx, invdiry = 1/diry, invdirz = 1/dirz;
	int sign[3];
	sign[0] = (invdirx < 0);
	sign[1] = (invdiry < 0);
	sign[2] = (invdirz < 0);
	
	float tmin, tmax, tymin, tymax, tzmin, tzmax;

	tmin = (bounds[r.sign[0]].x - r.orig.x) * r.invdir.x;
	tmax = (bounds[1-r.sign[0]].x - r.orig.x) * r.invdir.x;
	tymin = (bounds[r.sign[1]].y - r.orig.y) * r.invdir.y;
	tymax = (bounds[1-r.sign[1]].y - r.orig.y) * r.invdir.y;

	if ((tmin > tymax) || (tymin > tmax))
	return false;

	if (tymin > tmin)
	tmin = tymin;
	if (tymax < tmax)
	tmax = tymax;

	tzmin = (bounds[r.sign[2]].z - r.orig.z) * r.invdir.z;
	tzmax = (bounds[1-r.sign[2]].z - r.orig.z) * r.invdir.z;

	if ((tmin > tzmax) || (tzmin > tmax))
	return false;

	if (tzmin > tmin)
	tmin = tzmin;
	if (tzmax < tmax)
	tmax = tzmax;

	t = tmin;

	if (t < 0) {
	t = tmax;
	if (t < 0) return false;
	}

	return true;
} 


// polar angle measured from center of detector at which DM enters
double DUNEDetector::thetaenter (double theta) {
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
double DUNEDetector::thetaexit (double theta) {
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
double DUNEDetector::Lenter (double theta) {
	double rLenter;
	double Rdet = 6.4;
	double ddet = 500.0;
	double thenter = thetaenter(theta);
	rLenter = Rdet*sin(thenter)/sin(theta); 
	return(rLenter);
}
// distance from target DM exit point of detector
double DUNEDetector::Lexit (double theta) {
	double rLexit;
	double Rdet = 6.4;
	double ddet = 500.0;
	double thexit = thetaexit(theta);
	rLexit = Rdet*sin(thexit)/sin(theta); 
	return(rLexit);
}
// distance DM travels through detector
double DUNEDetector::Ldet (double theta) {
	double rLdet;
	double Ldetenter, Ldetexit;
	Ldetenter = Lenter(theta); 
	Ldetexit = Lexit(theta); 
        //std::cout<<"Lenter"<<Ldetenter<<"Lexit"<<Ldetexit<<std::endl;
	rLdet = Ldetexit - Ldetenter;
	return(rLdet);
}
//
void DUNEDetector::intersect(int &dswitch, int &Ndm, Particle &DM){
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
