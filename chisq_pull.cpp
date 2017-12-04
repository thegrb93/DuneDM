/*
 * chisq_pull.cc
 * This file is part of LDM code
 *
 * Copyright (C) 2013 - Animesh Chatterjee
 *
 * atmosnu is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * abc_odesolver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with abc_odesolver. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>
#define ARMA_NO_DEBUG

#include "chisq_pull.h"

using namespace std;
using namespace arma;
static const unsigned int row = 5;
static const unsigned int col = 5;

static double flux_norm = 0.20;
static double tilt_fact = 0.05;
static double zen_param = 0.05;
static double xsec_err  = 0.1;
static double sys_err   = 0.05;

inline double sqr(double x){return x*x;}

vector<double> matrix_inv(vector<vector<double> > &A, vector<double> &B)
{
	const unsigned int row1 = A.size();
	const unsigned int col1 = A.front().size();
	vec D(B.size());
	mat E(row1, col1);
	
	for (unsigned int i = 0; i < D.n_rows; i++)
	{
		D(i) = B.at(i);

	}

	for (unsigned int i = 0; i < E.n_rows; ++i)
	{
		for (unsigned int j = 0; j < A.at(i).size(); ++j)
		{
			E(i, j) = A.at(i).at(j);

		}
	}

	vec x = solve(E, D);
	vector<double> r;
	for (unsigned int i = 0; i < x.n_rows; i++)
	{
		r.push_back(x(i));

	}

	return r;

}


double chisq_pullfunc(vector<double> &N_th, vector<double> &N_ex,
		      vector<double> &cs_mean, vector<double> &tilt)
{
	double  x_n = 0, s_n = 0, x_k = 0, row1, st_n=0;

	// Errors on different parameter
        row1=N_th.size();

	vector<vector<double> > f(row1, vector<double> (col));

	for (unsigned int i = 0; i < f.size(); i++)
	{
		f.at(i).at(0) = N_th.at(i) * flux_norm;		// flux normalization
		f.at(i).at(1) = N_th.at(i) * tilt.at(i);	// tilt
		f.at(i).at(2) = cs_mean.at(i) * zen_param;	//zenith unc
		f.at(i).at(3) = N_th.at(i) * xsec_err;		// cross section error
		f.at(i).at(4) = N_th.at(i) * sys_err;		// over all systematic
	}

	vector<double> b(row);

	for (unsigned int i = 0; i < b.size(); i++)
	{
		b.at(i) = 0;
		for (unsigned int k = 0; k < N_th.size(); k++)
		{
			
			b.at(i)  += ( f.at(k).at(i) * (N_ex.at(k) - N_th.at(k)))
			             / N_ex.at(k);

		}
	}

	vector<vector<double> > A(row, vector<double> (col));
	for (unsigned int i = 0; i < A.size(); ++i)
	{
		for (unsigned int j = 0; j < A.at(i).size(); ++j)
		{
			if (i == j)
				A.at(i).at(j) = 1;
			else
				A.at(i).at(j) = 0;
		}
	}

	for (unsigned int i = 0; i < A.size(); ++i)
	{
		for (unsigned int j = 0; j < A.at(i).size(); ++j)
		{
			for (unsigned int k = 0; k < N_th.size(); k++)
			{
				A.at(i).at(j) += ((f.at(k).at(i) * f.at(k).at(j))/N_ex.at(k));
			}
		}
	}

	// Minimization of chisq with respect to xi
	vector<double> xi = matrix_inv(A, b);

	//chisq(N_ex-N_th)
	for (unsigned int k = 0; k < N_th.size(); k++)
	{
		s_n = 0;
		for (unsigned int i = 0; i < xi.size(); i++)
		{
			s_n += xi.at(i) * f.at(k).at(i);
		}

		x_n +=
		    (sqr(N_th.at(k) + s_n - N_ex.at(k)) / N_ex.at(k));

	}

	for(unsigned int k = 0; k < N_th.size(); k++)
	{
		st_n += sqr(N_ex.at(k) - N_th.at(k)) / N_ex.at(k);
	}

	// chisq(xi)
	for (unsigned int i = 0; i < xi.size(); i++)
	{
		x_k += sqr(xi.at(i));

	}

	double chisq_pull = x_n + x_k;
	return chisq_pull;
}


double chisq_pullfunc(vector<double> &N_th, vector<double> &N_ex,
		      vector<double> &cs_mean)
{
	double  x_n = 0, s_n = 0, x_k = 0, row1,st_n=0;

	// Errors on different parameter
        row1=N_th.size();

	vector<vector<double> > f(row1, vector<double> (col));

	for (unsigned int i = 0; i < f.size(); i++)
	{
		f.at(i).at(0) = N_th.at(i) * flux_norm;		// flux normalization
		f.at(i).at(1) = N_th.at(i) * tilt_fact;		// tilt
		f.at(i).at(2) = cs_mean.at(i) * zen_param;	//zenith unc
		f.at(i).at(3) = N_th.at(i) * xsec_err;		// cross section error
		f.at(i).at(4) = N_th.at(i) * sys_err;		// over all systematic
	}

	vector<double> b(row);

	for (unsigned int i = 0; i < b.size(); i++)
	{
		b.at(i) = 0;
		for (unsigned int k = 0; k < N_th.size(); k++)
		{
			
			b.at(i)  += ( f.at(k).at(i) * (N_ex.at(k) - N_th.at(k)))
			             / N_ex.at(k);

		}
	}

	vector<vector<double> > A(row, vector<double> (col));
	for (unsigned int i = 0; i < A.size(); ++i)
	{
		for (unsigned int j = 0; j < A.at(i).size(); ++j)
		{
			if (i == j)
				A.at(i).at(j) = 1;
			else
				A.at(i).at(j) = 0;
		}
	}

	for (unsigned int i = 0; i < A.size(); ++i)
	{
		for (unsigned int j = 0; j < A.at(i).size(); ++j)
		{
			for (unsigned int k = 0; k < N_th.size(); k++)
			{
				A.at(i).at(j) += ((f.at(k).at(i) * f.at(k).at(j))/N_ex.at(k));
			}
		}
	}

	// Minimization of chisq with respect to xi
	vector<double> xi = matrix_inv(A, b);

	//chisq(N_ex-N_th)
	for (unsigned int k = 0; k < N_th.size(); k++)
	{
		s_n = 0;
		for (unsigned int i = 0; i < xi.size(); i++)
		{
			s_n += xi.at(i) * f.at(k).at(i);
		}

		x_n +=
		    (sqr(N_th.at(k) + s_n - N_ex.at(k)) / N_ex.at(k));

	}

	for(unsigned int k = 0; k < N_th.size(); k++)
	{
		st_n += sqr(N_ex.at(k) - N_th.at(k)) / N_ex.at(k);
	}

	// chisq(xi)
	for (unsigned int i = 0; i < xi.size(); i++)
	{
		x_k += sqr(xi.at(i));

	}

	double chisq_pull = x_n + x_k;
	return chisq_pull;
}
