#include "CGS.h"

CGS::CGS(double* A, int* IA, int* JA, double* b, double* x0, double Tol, int maxIter, int unknowns, int entrynum)
	:unknowns(unknowns), entrynum(entrynum), Tol(Tol), maxIter(maxIter)
{
	this->Solve(A, IA, JA, b, x0);
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

double* CGS::precond(const double* pre)
{
	double* Temp = ArrayUtils<double>::onetensor(unknowns);



	return(Temp);
}
