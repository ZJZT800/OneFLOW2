#include "GMRes.h"
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

SolveEqua::SolveEqua(RealField& sp, RealField2D& ai, RealField& rhs, RealField& x, Real res, std::string TypeNum, int MaxIter, double Tol, int Restarts)
	: MaxIter(MaxIter), Tol(Tol), restarts(Restarts), TypeNum(TypeNum)
{
	this->CSR(sp, ai, rhs);
	if (TypeNum == "GMRES")
	{
		this->solveGMRES();
	}
	this->vecCopy(x0, x);

	res = ress;
	this->Deallocate();
}

SolveEqua::~SolveEqua()
{

}

void SolveEqua::solveGMRES()
{
	//double dur;
	//clock_t start, end;

	GMRes solver(MaxIter, restarts, unknowns, Tol);

	//start = clock();

	ress = solver.Solve(A, IA, JA, x0, b);

	//end = clock();
	//dur = (double)(end - start);

	//std::cout << dur << "ms" << std::endl;
	//printf("Use Time:%f\n", (dur / CLOCKS_PER_SEC));

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

	//ofstream file("CoeMatrix.txt", ios::app);
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
					//file << cId + 1 << "\t" << rc + 1 << "\t" << setprecision(18) << A[n + tempCout] << std::endl;
					tempCout += 1;
				}
				else if (cId == rc)
				{
					A[n + tempCout] = -ai[1][fId];
					JA[n + tempCout] = lc;
					//file << cId + 1 << "\t" << lc + 1 << "\t" << setprecision(18) << A[n + tempCout] << std::endl;
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
		//file << cId + 1 << "\t" << cId + 1 << "\t" << setprecision(18) << A[n + fj] << std::endl;
	}
	//file.close();
	//ofstream RHS("rhs.txt", ios::app);
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		b[cId] = rhs[cId];
		//RHS << setprecision(18) << b[cId] << endl;
	}
	//RHS.close();

}

void SolveEqua::Init()
{
	b = ArrayUtils<double>::onetensor(unknowns);
	x0 = ArrayUtils<double>::onetensor(unknowns);
	A = ArrayUtils<double>::onetensor(entryNum);
	IA = ArrayUtils<int>::onetensor(unknowns + 1);
	JA = ArrayUtils<int>::onetensor(entryNum);
}

void SolveEqua::Deallocate()
{
	ArrayUtils<double>::delonetensor(b);
	ArrayUtils<double>::delonetensor(x0);
	ArrayUtils<double>::delonetensor(A);
	ArrayUtils<int>::delonetensor(IA);
	ArrayUtils<int>::delonetensor(JA);


}

void SolveEqua::vecCopy(double* x0, RealField& x)
{
	int k = x.size();
	for (int i = 0; i < k; i++)
	{
		x[i] = x0[i];
	}
}

EndNameSpace