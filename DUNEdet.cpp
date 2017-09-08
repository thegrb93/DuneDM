#include <cmath>
#include <algorithm>

#include "DUNEdet.h"
#include "Particle.h"

static bool RaySlabIntersect(double slabmin, double slabmax, double raystart, double raydir, double& tbenter, double& tbexit)
{
    // ray parallel to the slab
    if (fabs(raydir) < 1.0E-9f)
    {
        // ray parallel to the slab, but ray not inside the slab planes
        if(raystart < slabmin || raystart > slabmax)
        {
            return false;
        }
            // ray parallel to the slab, but ray inside the slab planes
        else
        {
            return true;
        }
    }

    // slab's enter and exit parameters
    double tsenter = (slabmin - raystart) / raydir;
    double tsexit = (slabmax - raystart) / raydir;

    // order the enter / exit values.
    if(tsenter > tsexit)
    {
        double temp = tsexit;
        tsexit = tsenter;
        tsenter = temp;
    }

    // make sure the slab interval and the current box intersection interval overlap
    if (tbenter > tsexit || tsenter > tbexit)
    {
        // nope. Ray missed the box.
        return false;
    }
        // yep, the slab and current intersection interval overlap
    else
    {
        // update the intersection interval
        tbenter = std::max(tbenter, tsenter);
        tbexit = std::min(tbexit, tsexit);
        return true;
    }
}

static bool intersectAABB(double x, double y, double z,
                   double dirx, double diry, double dirz,
                   double minsx, double minsy, double minsz,
                   double maxsx, double maxsy, double maxsz, double &tenter, double& texit)
{
    tenter = 0;
    texit = std::numeric_limits<double>::max();
    if (!RaySlabIntersect(minsx, maxsx, x, dirx, tenter, texit))
        return false;
    if (!RaySlabIntersect(minsy, maxsy, y, diry, tenter, texit))
        return false;
    if (!RaySlabIntersect(minsz, maxsz, z, dirz, tenter, texit))
        return false;
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
double DUNEDetector::Ldet (Particle &DM) {
    double pLen = sqrt(DM.px*DM.px+DM.py*DM.py+DM.pz*DM.pz);
    double tmin = 0, tmax = 0;
    if(intersectAABB(0, 0, 0, DM.px/pLen, DM.py/pLen, DM.pz/pLen, minsx, minsy, minsz, maxsx, maxsy, maxsz, tmin, tmax))
    {
        return tmax - tmin;
    }
}
//
bool DUNEDetector::intersect(double dx, double dy, double dz, double& tmin, double& tmax){
    return intersectAABB(0, 0, 0, dx, dy, dz, minsx, minsy, minsz, maxsx, maxsy, maxsz, tmin, tmax);
}

DUNEDetector::DUNEDetector():
posx(0), posy(0), posz(570)
{
    double width = 3.5, height = 3.5, depth = 6.4;
    minsx = posx - width/2;
    maxsx = posx + width/2;
    minsy = posy - height/2;
    maxsy = posy + height/2;
    minsz = posz - depth/2;
    maxsz = posz + depth/2;
}
