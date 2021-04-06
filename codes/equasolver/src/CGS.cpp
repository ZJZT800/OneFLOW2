#include "CGS.h"
#include <iomanip>

BeginNameSpace(ONEFLOW)

CGS::CGS(double* A, int* IA, int* JA, double* b, double* x0, double Tol, int maxIter, int unknowns, int entrynum)
	:unknowns(unknowns), entrynum(entrynum), Tol(Tol), maxIter(maxIter)
{
	this->Solve(A, IA, JA, b, x0);
}

CGS::CGS(double* sp, double** aii, double* rhs, double* x, double Tol, int maxIter, int unknowns)
	: unknowns(unknowns), Tol(Tol), maxIter(maxIter)
{
	this->SolveDirect(sp, aii, rhs, x);
}

CGS::~CGS()
{
}

void CGS::Solve(double* A, int* IA, int* JA, double* b, double* x)
{
	double* residual = ArrayUtils<double>::onetensor(unknowns);
	double* u = ArrayUtils<double>::onetensor(unknowns);
	double* r_hat = ArrayUtils<double>::onetensor(unknowns);
	double* p = ArrayUtils<double>::onetensor(unknowns);
	double* q = ArrayUtils<double>::onetensor(unknowns);
	double* r_old = ArrayUtils<double>::onetensor(unknowns);
	double alpha, beta, initial_residual;
	int iter = 0;

	basic.CSRMVs(A, IA, JA, x, residual, unknowns);

	basic.vecPlus(b, residual, unknowns, -1);        //get the residual vector

	initial_residual = basic.norm(residual, unknowns);

	res = initial_residual;

	basic.vecCopy(p,  residual, unknowns);
	basic.vecCopy(u,  residual, unknowns);
	basic.vecCopy(r_hat,  residual, unknowns);

	while (res > Tol && iter < maxIter)
	{
		double* temp = ArrayUtils<double>::onetensor(unknowns);
		double* temp2 = ArrayUtils<double>::onetensor(unknowns);
		double* temp3 = ArrayUtils<double>::onetensor(unknowns);

		alpha = basic.dot(residual, r_hat, unknowns);
		basic.CSRMVs(A, IA, JA, p, temp, unknowns);    //temp = A*p
		alpha = alpha / basic.dot(temp, r_hat, unknowns);

		basic.vecPlusb(temp3, temp, unknowns, NEGONE * alpha);    //temp3 = alpha*A*p
		basic.vecPlus(u, temp3, q, unknowns, ONE);     //q=u-alpha*A*p

		basic.vecPlus(q, u, temp, unknowns, ONE);
		basic.vecPlusb(x, temp, unknowns, alpha);    //update x

		basic.CSRMVs(A, IA, JA, temp, temp2, unknowns);    //temp2 = A*temp
		basic.vecCopy(r_old, residual, unknowns);

		basic.vecPlusb(residual, temp2, unknowns, NEGONE*alpha);  //update residual

		beta = basic.dot(residual, r_hat, unknowns);
		beta = beta / basic.dot(r_old, r_hat, unknowns);
		basic.vecPlus(residual, q, u, unknowns, beta);    //update u

		basic.vecPlus(u, p, unknowns, beta * beta);
		
		basic.vecPlusb(p, q, unknowns, beta);
		
		iter++;

		ArrayUtils<double>::delonetensor(temp);
		ArrayUtils<double>::delonetensor(temp2);
		ArrayUtils<double>::delonetensor(temp3);
		temp = NULL;
		temp2 = NULL;
		temp3 = NULL;
	}

	ArrayUtils<double>::delonetensor(u);
	ArrayUtils<double>::delonetensor(r_hat);
	ArrayUtils<double>::delonetensor(p);
	ArrayUtils<double>::delonetensor(q);
	ArrayUtils<double>::delonetensor(r_old);
	ArrayUtils<double>::delonetensor(residual);
	u = NULL;
	r_hat = NULL;
	p = NULL;
	q = NULL;
	residual = NULL;
	r_old = NULL;
}

void CGS::SolveDirect(double* sp, double** aii, double* rhs, double* x)
{
	double* residual = ArrayUtils<double>::onetensor(unknowns);
	double* u = ArrayUtils<double>::onetensor(unknowns);
	double* r_hat = ArrayUtils<double>::onetensor(unknowns);
	double* p = ArrayUtils<double>::onetensor(unknowns);
	double* q = ArrayUtils<double>::onetensor(unknowns);
	double* r_old = ArrayUtils<double>::onetensor(unknowns);
	double* apTmp = ArrayUtils<double>::onetensor(unknowns);
	double alpha, beta, initial_residual, res0;
	int iter = 0;

	PrecondILU(sp, aii, apTmp);
	PrecondRU(aii, apTmp, rhs);

	for (int i = 0; i < unknowns; ++i)
	{
		residual[i] = rhs[i];
		p[i] = residual[i];
		u[i] = residual[i];
		r_hat[i] = residual[i];
	}

	res0 = basic.dot(rhs, rhs, unknowns);
	res = res0;
	//basic.vecCopy(p, residual, unknowns);
	//basic.vecCopy(u, residual, unknowns);
	//basic.vecCopy(r_hat, residual, unknowns);

	while (res > Tol && iter < maxIter)
	{
		double* temp = ArrayUtils<double>::onetensor(unknowns);
		double* temp2 = ArrayUtils<double>::onetensor(unknowns);
		double* temp3 = ArrayUtils<double>::onetensor(unknowns);

		alpha = basic.dot(residual, r_hat, unknowns);
		basic.MVs(sp, aii, p, temp);

		PrecondRU(aii, apTmp, temp);

		alpha = alpha / basic.dot(temp, r_hat, unknowns);

		for (int i = 0; i < unknowns; ++i)
		{
			q[i] = u[i] - alpha * temp[i];
		}
		//basic.vecPlusb(temp3, temp, unknowns, NEGONE * alpha);    //temp3 = alpha*A*p
		//basic.vecPlus(u, temp3, q, unknowns, ONE);     //q=u-alpha*A*p

		basic.vecPlus(q, u, temp, unknowns, ONE);
		basic.vecPlusb(x, temp, unknowns, alpha);    //update x

		basic.MVs(sp, aii, temp, temp2);
		PrecondRU(aii, apTmp, temp2);

		for (int i = 0; i < unknowns; ++i)
		{
			r_old[i] = residual[i];
			residual[i] = residual[i] - alpha * temp2[i];
		}
		//basic.vecCopy(r_old, residual, unknowns);
		//basic.vecPlusb(residual, temp2, unknowns, NEGONE * alpha);  //update residual

		res = basic.dot(residual, residual, unknowns);
		res = sqrt(abs(res / res0));
		
		if (res > Tol)
		{
			beta = basic.dot(residual, r_hat, unknowns);
			beta = beta / basic.dot(r_old, r_hat, unknowns);
			double betas = beta * beta;
			for (int i = 0; i < unknowns; ++i)
			{
				u[i] = residual[i] + beta * q[i];
				p[i] = u[i] + betas * p[i] + beta * q[i];
			}

			//basic.vecPlus(residual, q, u, unknowns, beta);    //update u

			//basic.vecPlus(u, p, unknowns, beta * beta);

			//basic.vecPlusb(p, q, unknowns, beta);
		}
		else
		{
			ArrayUtils<double>::delonetensor(temp);
			ArrayUtils<double>::delonetensor(temp2);
			ArrayUtils<double>::delonetensor(temp3);
			temp = NULL;
			temp2 = NULL;
			temp3 = NULL;
			continue;
		}

		iter++;
		
		ArrayUtils<double>::delonetensor(temp);
		ArrayUtils<double>::delonetensor(temp2);
		ArrayUtils<double>::delonetensor(temp3);
		temp = NULL;
		temp2 = NULL;
		temp3 = NULL;
	}

	ArrayUtils<double>::delonetensor(u);
	ArrayUtils<double>::delonetensor(r_hat);
	ArrayUtils<double>::delonetensor(p);
	ArrayUtils<double>::delonetensor(q);
	ArrayUtils<double>::delonetensor(r_old);
	ArrayUtils<double>::delonetensor(residual);
	ArrayUtils<double>::delonetensor(apTmp);
	u = NULL;
	r_hat = NULL;
	p = NULL;
	q = NULL;
	residual = NULL;
	r_old = NULL;
	apTmp = NULL;
}

void CGS::PrecondILU(double* sp, double** aii, double* apTmp)
{
	for (int cId = 0; cId < ONEFLOW::ug.nCell; ++cId)
	{
		apTmp[cId] = sp[cId];
	}
	for (int fId = ONEFLOW::ug.nBFace; fId < ONEFLOW::ug.nFace; ++fId)
	{
		int lc = (*ONEFLOW::ug.lcf)[fId];
		int rc = (*ONEFLOW::ug.rcf)[fId];
		apTmp[rc] = apTmp[rc] - aii[fId][0] * aii[fId][1] / apTmp[lc];
	}
	for (int cId = 0; cId < ONEFLOW::ug.nCell; ++cId)
	{
		apTmp[cId] = 1 / apTmp[cId];
	}
}

void CGS::PrecondRU(double** aii, double* apTmp, double* b)
{
	double* bTmp = ArrayUtils<double>::onetensor(unknowns);

	//forward
	for (int cId = 0; cId < ONEFLOW::ug.nCell; ++cId)
	{
		bTmp[cId] =  b[cId];
	}
	for (int fId = ONEFLOW::ug.nBFace; fId < ONEFLOW::ug.nFace; ++fId)
	{
		int lc = (*ONEFLOW::ug.lcf)[fId];
		int rc = (*ONEFLOW::ug.rcf)[fId];
		bTmp[rc] = bTmp[rc] + aii[fId][1] * apTmp[lc] * bTmp[lc];
	}
	//backforward
	for (int cId = 0; cId < ONEFLOW::ug.nCell; ++cId)
	{
		b[cId] = bTmp[cId] * apTmp[cId];
	}
	for (int fId = ONEFLOW::ug.nFace - 1; fId >= ONEFLOW::ug.nBFace; --fId)
	{
		int lc = (*ONEFLOW::ug.lcf)[fId];
		int rc = (*ONEFLOW::ug.rcf)[fId];
		b[lc] = b[lc] + aii[fId][0] * apTmp[lc] * b[rc];
	}
	ArrayUtils<double>::delonetensor(bTmp);
	bTmp = NULL;
}

EndNameSpace