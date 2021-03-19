#ifndef GMRES_H
#define GMRES_H

#define MKL_Complex16 std::complex<double>

#include "Utils.h"
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
#include <mkl.h>
using namespace std;

class GMRes
{
  public:
	const int krylovDemension;
	const int mRestart;
	const int unknows;
	const double tolerance;


	double norm_b,norm_residual;

	const double ONE = 1.0;
	const double NEGONE = -1.0;
	const double ZERO = 0.0;
	double *residual, *H, *v, *q, *y, *x, *beta, *res_n, **givens;
	std::vector<double*>Q;


	GMRes(int krylovDemension, int mRestart, int unknows, double tolerance);
	~GMRes();
	double Solve(double* A, int* IA, int* JA, double* x0, double* b);
	double RestartGMRes(double *A, int* IA, int* JA, double *x0, double *b);
	double InnerLoop(double *A, int* IA, int* JA, double *x0, double *b);
	void Update(double* H, double* x, double* beta, std::vector<double*> Q, int iteration);
};

#endif // GMRES_H