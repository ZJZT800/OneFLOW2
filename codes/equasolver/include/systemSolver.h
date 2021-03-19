#pragma once

#include "HXDefine.h"
#include "UCom.h"
#include "FaceTopo.h"

BeginNameSpace(ONEFLOW)

class SolveEqua 
{
public:
	enum SolverType
	{
		GMRES = 1,
		BiCGStab = 2,
		CGS = 3,
	};

public:
	SolveEqua();
	SolveEqua(RealField& sp, RealField2D& ai, RealField& rhs, RealField& x, Real res, std::string TypeNum, int MaxIter, double Tol, int Restarts);
	~SolveEqua();

public:
	std::string TypeNum;
	int restarts = 10;
	int MaxIter = 500;
	double Tol = 1e-8;
	double ONE = 1.0;
	double ZERO = 0.0;
	double ress;
	int unknowns = 0;
	int entryNum = 0;

	double* x0;
	double* b;
	double* A;
	int* IA;
	int* JA;

	void solveGMRES();
	void CSR(RealField& sp, RealField2D& ai, RealField& rhs);
	void Init();
	void Deallocate();
	void vecCopy(double* x0, RealField& x);
};

EndNameSpace
