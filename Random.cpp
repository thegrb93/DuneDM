#include "Random.h"
#include <time.h>
#include <iostream>
#include <cmath>

//------------------------------------------------------------------------------
// The constructors
// Set seed from the system timer or initiate seed yourself
//
Random::Random()        { srand(time(NULL)); }
Random::Random(int seed){ srand(seed); }

//------------------------------------------------------------------------------
// Flat (uniform) distribution
// Returns a uniformly distributed random real number between [min,max]
//
double Random::Flat(double min, double max){
	return ( min + (max-min)*rand()/RAND_MAX ) ;
}
//------------------------------------------------------------------------------
// Integer (uniform) distribution
// Returns a uniformly distributed random integer number between [min,max]
//
int Random::Integer(int min, int max){
	return ( min + int(max*Flat()) ) ;
}
//------------------------------------------------------------------------------
// Gaussian Distribution
// Returns a normally distributed deviate with mean and sigma
// The routine uses the Box-Muller transformation of uniform deviates.
//
double Random::Gauss(double mean, double sigma){
	double x, y, z;
	
	while(1){
        x = 2.0 * Flat() - 1.0;
        y = 2.0 * Flat() - 1.0;
        z = x*x + y*y;
        if( z <= 1.0 ) break;
	}
	return ( mean + sigma*x*sqrt(-2.0*log(z)/z) );
}
//------------------------------------------------------------------------------
// Exponential (decay) distribution
// Returns a random number between times t1 and t2 
// according to f(t) = exp (-t/tau)
//
double Random::Exponential(double tau, double tmin, double tmax){
	double r1 =  exp(-tmin/tau);
	double r2 =  exp(-tmax/tau);
	double ed = -tau*log(r2 + Flat() * (r1-r2) );
	
	return (ed);
}
//------------------------------------------------------------------------------
// Breit-Wigner Distribution
// Returns a random number from a Breit-Wigner distribution 
// for center mean Full Width Half Maximum fwhm
//
double Random::BreitWigner(double mean, double fwhm){
	double x, y, z;
	
	while(1){
        x = 2.0 * Flat() - 1.0;
        y = 2.0 * Flat() - 1.0;
        z = x*x + y*y;
        if( z <= 1.0 ) break;
	}
	return ( mean + 0.5*fwhm*x/y );
}

