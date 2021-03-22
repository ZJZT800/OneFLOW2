#include "GMRes.h"

using namespace std;

GMRes::GMRes(int krylovDemension, int mRestart, int unknows, double tolerance)
	: krylovDemension(krylovDemension), mRestart(mRestart), unknows(unknows), tolerance(tolerance)
{
	res_n = ArrayUtils<double>::onetensor(unknows);
	residual = ArrayUtils<double>::onetensor(unknows);
	//Q = ArrayUtils<double>::twotensor((krylovDemension + 1), unknows);
	H = ArrayUtils<double>::twotensor(krylovDemension, krylovDemension + 1);    //ColMajor
	v = ArrayUtils<double>::onetensor(unknows);
	q = ArrayUtils<double>::onetensor(unknows);
	y = ArrayUtils<double>::onetensor(krylovDemension + 1);
	x = ArrayUtils<double>::onetensor(unknows);
	beta = ArrayUtils<double>::onetensor(krylovDemension + 1);
	givens = ArrayUtils<double>::twotensor(krylovDemension, 2);
	Q.resize(krylovDemension + 1);
}

GMRes::~GMRes()
{
}

int GMRes::Solve(double* A, int* IA, int* JA, double* x0, double* b)
{
	double error = RestartGMRes(A, IA, JA, x0, b);
	this->Deallocate();
	return error;
}

double GMRes::RestartGMRes(double* A, int* IA, int* JA, double* x0, double* b)
{
	double error = 1;
	int outloops = 0;
	norm_b = basic.norm(b, unknows);
	error = norm_b;
	if (norm_b < 1e-5)
	{
		norm_b = 1.0;
	}
	while (outloops < mRestart && error > tolerance * norm_b)
	{
		outloops++;
		error = InnerLoop(A, IA, JA, x0, b);
	}
	return error;
}

double GMRes::InnerLoop(double* A, int* IA, int* JA, double* x0, double* b)
{

	basic.CSRMVs(A, IA, JA, x0, residual, unknows);

	basic.vecPlus(b, residual, unknows, -1);        //get the residual vector

	norm_residual = basic.norm(residual, unknows);

	beta[0] = norm_residual;

	double one_norm = ONE / norm_residual;

	basic.VCs(residual, q, one_norm, unknows);

	Q[0] = ArrayUtils<double>::onetensor(unknows);

	basic.vecCopy(Q[0], q, unknows);

	for (int k = 0; k < krylovDemension; k++)     //k represents iter
	{
		Q[k + 1] = ArrayUtils<double>::onetensor(unknows);

		basic.CSRMVs(A, IA, JA, q, v, unknows);

		for (int j = 0; j < k + 1; j++)
		{
			H[k][j] = basic.dot(Q[j], v, unknows);      //colMajor

			double coe = NEGONE * H[k][j];
			
			basic.vecPlusb(v, Q[j], unknows, coe);
		}
		H[k][k + 1] = basic.norm(v, unknows);

		double coe = ONE / H[k][k + 1];

		basic.VCs(v, q, coe, unknows);
	
		basic.vecCopy(Q[k + 1], q, unknows);

		//givens rotation for single rhs
		for (int j = 0; j < k; j++)
		{
			basic.rot(&H[k][j], &H[k][j + 1], &givens[j][0], &givens[j][1]);
		}

		basic.rotg(&H[k][k], &H[k][k + 1], &givens[k][0], &givens[k][1]);
		
		//givens rotation for residual vector
		basic.rot(&beta[k], &beta[k + 1], &givens[k][0], &givens[k][1]);
		
		//if convergence?
		if (abs(beta[k + 1]) < tolerance * norm_b)
		{
			/*std::cout << "tolerance: " << tolerance << std::endl;
			std::cout << "beta[k + 1]: " << beta[k + 1] << std::endl;
			std::cout << "norm_b: " << norm_b << std::endl;*/

			Update(H, x, beta, Q, k + 1);
			basic.vecCopy(x0, x, unknows);
			for (int i = 0; i < k + 2; ++i)
			{
				ArrayUtils<double>::delonetensor(Q[i]);
			}
			return(beta[k + 1]);
		}
	}

	Update(H, x, beta, Q, krylovDemension);
	basic.vecCopy(x0, x, unknows);
	for (int i = 0; i < krylovDemension + 1; ++i)
	{
		ArrayUtils<double>::delonetensor(Q[i]);
	}
	return(beta[krylovDemension - 1]);
}

void GMRes::Update(double** H, double* x, double* beta, std::vector<double*> Q, int iteration)
{
	int lupe;
	for (lupe = iteration - 1; lupe >= 0; --lupe)
	{
		beta[lupe] = beta[lupe] / H[lupe][lupe];

		for (int innerlupe = lupe - 1; innerlupe >= 0; --innerlupe)
		{
			beta[innerlupe] -= beta[lupe] * H[lupe][innerlupe];
		}
	}
	for (lupe = 0; lupe < iteration; lupe++)
	{
		basic.VCs(Q[lupe], Q[lupe], beta[lupe], unknows);
		basic.vecPlus(Q[lupe], x, unknows);
	}
}

void GMRes::Deallocate()
{
	ArrayUtils<double>::delonetensor(res_n);
	ArrayUtils<double>::delonetensor(residual);
	std::vector<double*>().swap(Q);
	ArrayUtils<double>::deltwotensor(H);
	ArrayUtils<double>::delonetensor(v);
	ArrayUtils<double>::delonetensor(q);
	ArrayUtils<double>::delonetensor(y);
	ArrayUtils<double>::delonetensor(x);
	ArrayUtils<double>::delonetensor(beta);
	ArrayUtils<double>::deltwotensor(givens);
	res_n = NULL;
	residual = NULL;
	H = NULL;
	v = NULL;
	q = NULL;
	y = NULL;
	x = NULL;
	beta = NULL;
	givens = NULL;
}

