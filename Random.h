#pragma once
/*  Random.h
 
 Random Class  (Aug 2008)
 Returns random numbers from a special distributions of popular algorithms
 */

class Random{
public:
	Random();
	Random(int);
	static double Flat       (double =0.0, double =1.0);
	static double Gauss      (double =0.0, double =1.0);
	static double Exponential(double, double =0.0, double =1.0e100);
	static double BreitWigner(double =0.0, double =1.0);
	static int    Integer    (int =0, int =100);
};

