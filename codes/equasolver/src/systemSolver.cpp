#include "GMRes.h"
#include "CGS.h"
#include "BiCGStab.h"
#include "Utils.h"
#include "systemSolver.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>


BeginNameSpace(ONEFLOW)


SolveEqua::SolveEqua()
{
}

SolveEqua::SolveEqua(RealField& sp, RealField2D& ai, RealField& rhs, RealField& x, Real res, std::string TypeNum, int MaxIter, double Tol, int Restarts, bool ifPrecond)
	: MaxIter(MaxIter), Tol(Tol), restarts(Restarts), TypeNum(TypeNum)
{	
	unknowns = ONEFLOW::ug.nCell;
	entryNum = (ONEFLOW::ug.nFace - ONEFLOW::ug.nBFace) * 2 + ONEFLOW::ug.nCell;

	if (TypeNum == "GMRES")
	{
		this->CSR(sp, ai, rhs);
		this->solveGMRES();
		this->vecCopy(x0, x);
		this->Deallocate();

	}
	else if (TypeNum == "CGS")
	{
		Assignment(sp, ai, rhs, x);
		CGS(this->sp, this->ai, this->bb, this->xx, Tol, MaxIter, unknowns);
		this->vecCopy(this->xx, x);
		this->DeallocateCGS();
	}
	else if (TypeNum == "PGAMG")
	{
		this->CSRsp(sp, ai, rhs);
		system("pause");
		/*typedef void (_stdcall*PGAMG)(int nrows, int ncols, int nnonzeros, int* IA, int* JA, double* d, double* A, double* u, double* b, int groupsize, int max_sweep, double res);
		HINSTANCE hDLL;
		hDLL = LoadLibrary(TEXT("Dll3_2.dll"));
		if (hDLL == NULL)
		{
			std::cout << "failed load Dll3.dll" << std::endl;
		}
		PGAMG pGAMG = (PGAMG)GetProcAddress(hDLL, "PGAMG");
		if (pGAMG == NULL)
		{
			std::cout << "failed get procaddress" << std::endl;
		}
		double a = GetLastError();
		system("pause");*/
		//PGAMG(unknowns, unknowns, entryNum, IA, JA, this->sp, A, this->x0, this->b, groupsize, max_sweep, ress);
		
		this->vecCopy(this->x0, x);
		this->Deallocate();
		//FreeLibrary(hDLL);
	}
	else if (TypeNum == "BiCGStab")
	{
		RealField2D aii;
		aii.resize(ug.nFace, 2);
		for (int fId = 0; fId < ug.nFace; ++fId)
		{
			aii[fId][0] = ai[0][fId];
			aii[fId][1] = ai[1][fId];
		}
		BiCGStab(sp, aii, rhs, x, unknowns, MaxIter, Tol, ifPrecond);
	}
	res = ress;
}

SolveEqua::~SolveEqua()
{

}

void SolveEqua::solveGMRES()
{

	GMRes solver(MaxIter, restarts, unknowns, Tol);

	ress = solver.Solve(A, IA, JA, x0, b);
}

void SolveEqua::CSR(RealField& sp, RealField2D& ai, RealField& rhs)
{
	RealField dj;
	unknowns = ug.nCell;
	dj.resize(ug.nCell);
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		dj[cId] = (*ug.c2f)[cId].size();
		for (int iFace = 0; iFace < (*ug.c2f)[cId].size(); ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];
			if (fId < ug.nBFace)
			{
				dj[cId] -= 1;
			}
		}
		entryNum += dj[cId];
	}
	entryNum += ug.nCell;
	this->Init();

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		IA[0] = 0;
		int n = IA[cId];
		int fn = (*ug.c2f)[cId].size();
		IA[cId + 1] = IA[cId] + dj[cId] + 1;
		int tempCout = 0;
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];
			int lc = (*ug.lcf)[fId];

			if (fId > ug.nBFace - 1)
			{
				int rc = (*ug.rcf)[fId];
				if (cId == lc)
				{
					A[n + tempCout] = -ai[0][fId];
					JA[n + tempCout] = rc;
					tempCout += 1;
				}
				else if (cId == rc)
				{
					A[n + tempCout] = -ai[1][fId];
					JA[n + tempCout] = lc;
					tempCout += 1;
				}
			}
			else
			{
				continue;
			}
		}

		int fj = dj[cId];
		A[n + fj] = sp[cId];
		JA[n + fj] = cId;
	}
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		b[cId] = rhs[cId];
	}
	dj.resize(0);
}

void SolveEqua::CSRsp(RealField& sp, RealField2D& ai, RealField& rhs)
{
	RealField dj;
	int* Counttmp = ArrayUtils<int>::onetensor(unknowns);
	int lc, rc, k;
	unknowns = ug.nCell;
	dj.resize(ug.nCell);
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		dj[cId] = (*ug.c2f)[cId].size();
		for (int iFace = 0; iFace < (*ug.c2f)[cId].size(); ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];
			if (fId < ug.nBFace)
			{
				dj[cId] -= 1;
			}
		}
		entryNum += dj[cId];
	}

	this->InitAMG();

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		this->b[cId] = rhs[cId];
		this->sp[cId] = sp[cId];
		IA[0] = 1;
		IA[cId + 1] = IA[cId] + dj[cId];
	}
	std::ofstream param("param.txt");
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		lc = (*ug.lcf)[fId];
		rc = (*ug.rcf)[fId];

		k = IA[lc] + Counttmp[lc];
		A[k] = -ai[fId][0];
		JA[k] = rc + 1;
		Counttmp[lc] += 1;
		param << IA[lc] << '\t' << JA[k] << '\t' << A[k] << std::endl;
		k = IA[rc] + Counttmp[rc];
		A[k] = -ai[fId][1];
		JA[k] = lc + 1;
		Counttmp[rc] += 1;
		param << IA[rc] << '\t' << JA[k] << '\t' << A[k] << std::endl;
	}
	param.close();
	ArrayUtils<int>::delonetensor(Counttmp);
}

void SolveEqua::Init()
{
	b = ArrayUtils<double>::onetensor(unknowns);
	x0 = ArrayUtils<double>::onetensor(unknowns);
	A = ArrayUtils<double>::onetensor(entryNum);
	IA = ArrayUtils<int>::onetensor(unknowns + 1);
	JA = ArrayUtils<int>::onetensor(entryNum);
	sp = ArrayUtils<double>::onetensor(unknowns);
}

void SolveEqua::Deallocate()
{
	ArrayUtils<double>::delonetensor(b);
	ArrayUtils<double>::delonetensor(x0);
	ArrayUtils<double>::delonetensor(A);
	ArrayUtils<int>::delonetensor(IA);
	ArrayUtils<int>::delonetensor(JA);
	ArrayUtils<double>::delonetensor(sp);
	b = NULL;
	x0 = NULL;
	A = NULL;
	IA = NULL;
	JA = NULL;
	sp = NULL;
}

void SolveEqua::Assignment(RealField& spp, RealField2D& aii, RealField& rhs, RealField& xx)
{
	this->InitCGS();
	for (int cId = 0; cId < ONEFLOW::ug.nCell; cId++)
	{
		this->sp[cId] = spp[cId];
		this->xx[cId] = xx[cId];
		this->bb[cId] = rhs[cId];
	}
	for (int fId = ONEFLOW::ug.nBFace; fId < ONEFLOW::ug.nFace; ++fId)
	{
		this->ai[fId][0] = aii[0][fId];
		this->ai[fId][1] = aii[1][fId];
	}
}

void SolveEqua::vecCopy(double* x0, RealField& x)
{
	int k = x.size();
	for (int i = 0; i < k; i++)
	{
		x[i] = x0[i];
	}
}

void SolveEqua::InitCGS()
{
	sp = ArrayUtils<double>::onetensor(unknowns);
	ai = ArrayUtils<double>::twotensor(ONEFLOW::ug.nFace, 2);
	bb = ArrayUtils<double>::onetensor(unknowns);
	xx = ArrayUtils<double>::onetensor(unknowns);
}

void SolveEqua::DeallocateCGS()
{
	ArrayUtils<double>::delonetensor(sp);
	ArrayUtils<double>::deltwotensor(ai);
	ArrayUtils<double>::delonetensor(bb);
	ArrayUtils<double>::delonetensor(xx);
}

void SolveEqua::InitAMG()
{
	b = ArrayUtils<double>::onetensor(unknowns);
	x0 = ArrayUtils<double>::onetensor(unknowns + 1);
	A = ArrayUtils<double>::onetensor(entryNum);
	IA = ArrayUtils<int>::onetensor(unknowns + 1);
	JA = ArrayUtils<int>::onetensor(entryNum);
	sp = ArrayUtils<double>::onetensor(unknowns);
}



EndNameSpace