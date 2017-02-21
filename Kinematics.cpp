#include <cmath>

#include "Kinematics.h"
//
double Kinematics::pAbs (double px, double py, double pz, double p0) {
	double rpAbs;
	rpAbs = sqrt(px*px+py*py+pz*pz);
	return(rpAbs);
}
//
double Kinematics::pt (double px, double py, double pz, double p0) {
	double rpt;
	rpt = sqrt(px*px+py*py);
	return(rpt);
}
//
double Kinematics::theta (double px, double py, double pz, double p0) {
	double rtheta;
	rtheta = acos(pz/pAbs(px,py,pz,p0));
	return(rtheta);
}
//
double Kinematics::phi (double px, double py, double pz, double p0) {
	double rphi;
	if (px > 0 && py > 0) 
	{
		rphi=atan(fabs(py/px));
	}
	else if (px < 0 && py > 0) 
	{
		rphi=3.141592653589793-atan(fabs(py/px));
	}	
	else if (px < 0 && py < 0) 
	{
		rphi=3.141592653589793+atan(fabs(py/px));
	}	
	else
	{
		rphi = 2*3.141592653589793 - atan(fabs(py/px));
	}	
	return(rphi);
}
//
double Kinematics::invariantmass (double px, double py, double pz, double p0) {
	double rinvariantmass;
	rinvariantmass = sqrt(p0*p0-py*py-pz*pz-py*py);
	return(rinvariantmass);
}
//
double Kinematics::lambda (double a, double b, double c) {
	double rlam;
	rlam = a*a+b*b+c*c-2*a*b-2*a*c-2*b*c;
	return(rlam);
}
//
