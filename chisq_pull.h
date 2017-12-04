#ifndef _CHISQ_PULL_TILT_H__
#define _CHISQ_PULL_TILT_H__

#include <vector>

std::vector<double> matrix_inv(std::vector<std::vector<double> > &, 
                               std::vector<double> &);

double chisq_pullfunc(std::vector<double> &, std::vector<double> &,
		      std::vector<double> &,std::vector<double> &);
double chisq_pullfunc(std::vector<double> &,
		      std::vector<double> &,std::vector<double> &);

#endif
