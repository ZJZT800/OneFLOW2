#pragma once

#include "HXDefine.h"
#include "UCom.h"
#include "FaceTopo.h"
//#include "pGAMG.h"

BeginNameSpace(ONEFLOW)

class SolveEqua 
{

public:
	SolveEqua();
	SolveEqua(RealField& sp, RealField2D& ai, RealField& rhs, RealField& x, Real res, std::string TypeNum, int MaxIter, double Tol, int Restarts, bool ifPrecond);
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
	int groupsize = 2;
	int max_sweep = 30;

	double* x0;
	double* b;
	double* A;
	int* IA;
	int* JA;

	double* sp;
	double** ai;
	double* bb;
	double* xx;

public:
	void solveGMRES();
	void CSR(RealField& sp, RealField2D& ai, RealField& rhs);
	void CSRsp(RealField& sp, RealField2D& ai, RealField& rhs);
	void Init();
	void Deallocate();
	void Assignment(RealField& sp, RealField2D& ai, RealField& rhs, RealField& x);
	void vecCopy(double* x0, RealField& x);
	void InitCGS();
	void DeallocateCGS();
	void InitAMG();

};

EndNameSpace
