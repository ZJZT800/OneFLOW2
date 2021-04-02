#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include "Utils.h"
#include "BasicComputation.h"


class CGS
{
public:
	BasicCompute basic;
	int unknowns;
	int entrynum;
	double Tol;
	int maxIter;

	double res;
	const double ONE = 1.0;
	const double NEGONE = -1.0;
	const double ZERO = 0.0;
	
public:
	CGS(double* A, int* IA, int* JA, double* b, double* x0, double Tol, int maxIter, int unknowns, int entrynum);
	CGS(double* sp, double** aii, double* rhs, double* x, double Tol, int maxIter, int unknowns);
	~CGS();

public:
	void Solve(double* A, int* IA, int* JA, double* b, double* x);
	void SolveDirect(double* sp, double** aii, double* rhs, double* x);
	void PrecondILU(double* sp, double** aii, double* apTmp);
	void PrecondRU(double** aii, double* apTmp, double* b);
};