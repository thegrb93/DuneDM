#include <cmath>

#include "Particle.h"
#include <algorithm>

Particle::Particle(double mass){
	m = mass;
	px = 0.0;
	py = 0.0;
	pz = 0.0;
	E = mass;
}

void Particle::FourMomentum(double PX, double PY, double PZ, double P0){

	px = PX;
	py = PY;
	pz = PZ;
	E = P0;
}

void Particle::Lorentz(Particle parent){
	
	double M1 = parent.m;	
	double E1 = parent.E;
	double p1x = parent.px;
	double p1y = parent.py;
	double p1z = parent.pz;
	double E2 = E;	
	double p2x = px;
	double p2y = py;
	double p2z = pz;
	double E3, p3x, p3y, p3z;		
	double p1perp, p1, ctheta1, stheta1, cphi1, sphi1;
	double gamma, v;	
	double Lam11, Lam12, Lam13, Lam14;	
	double Lam21, Lam22, Lam23, Lam24;	
	double Lam31, Lam32, Lam33, Lam34;
	double Lam41, Lam42, Lam43, Lam44;
	p1 = sqrt(p1x*p1x+p1y*p1y+p1z*p1z);
	p1perp = sqrt(p1x*p1x+p1y*p1y);
	ctheta1 = p1z/p1;
	stheta1 = p1perp/p1;
	cphi1 = p1x/p1/stheta1;
	sphi1 = p1y/p1/stheta1;
	gamma = E1/M1;
	v = p1/E1;
	Lam11 = gamma;
	Lam12 = 0;
	Lam13 = 0;
	Lam14 = gamma*v;
	Lam21 = gamma*v*stheta1*cphi1;
	Lam22 = ctheta1*cphi1;
	Lam23 = -sphi1;
	Lam24 = gamma*stheta1*cphi1;
	Lam31 = gamma*v*stheta1*sphi1;
	Lam32 = ctheta1*sphi1;
	Lam33 = cphi1;
	Lam34 = gamma*stheta1*sphi1;
	Lam41 = gamma*v*ctheta1;
	Lam42 = -stheta1;
	Lam43 = 0;
	Lam44 = gamma*ctheta1;
	E3   = Lam11*E2+Lam12*p2x+Lam13*p2y+Lam14*p2z;	
	p3x  = Lam21*E2+Lam22*p2x+Lam23*p2y+Lam24*p2z;	
	p3y  = Lam31*E2+Lam32*p2x+Lam33*p2y+Lam34*p2z;	
	p3z  = Lam41*E2+Lam42*p2x+Lam43*p2y+Lam44*p2z;
		
	FourMomentum(p3x, p3y, p3z, E3);
}

void Particle::calcOptionalKinematics() {
    norm = sqrt(px*px+py*py+pz*pz);
    normpx = px/norm; normpy = py/norm; normpz = pz/norm;
    pt = sqrt(px*px+py*py);
    theta = acos(normpz);
    phi = atan2(py, px);
}

