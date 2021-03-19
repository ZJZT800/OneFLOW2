#include "GMRes.h"

using namespace std;

GMRes::GMRes(int krylovDemension, int mRestart, int unknows, double tolerance)
	: krylovDemension(krylovDemension), mRestart(mRestart), unknows(unknows), tolerance(tolerance)
{
	res_n = ArrayUtils<double>::onetensor(unknows);
	residual = ArrayUtils<double>::onetensor(unknows);
	H = ArrayUtils<double>::onetensor((krylovDemension + 1) * krylovDemension);
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
	ArrayUtils<double>::delonetensor(res_n);
	ArrayUtils<double>::delonetensor(residual);
	ArrayUtils<double>::delonetensor(H);
	ArrayUtils<double>::delonetensor(v);
	ArrayUtils<double>::delonetensor(q);
	ArrayUtils<double>::delonetensor(y);
	ArrayUtils<double>::delonetensor(x);
	ArrayUtils<double>::delonetensor(beta);
	ArrayUtils<double>::deltwotensor(givens);
}

double GMRes::Solve(double* A, int* IA, int* JA, double* x0, double* b)
{
	double error = RestartGMRes(A, IA, JA, x0, b);
	return error;
}

double GMRes::RestartGMRes(double* A, int* IA, int* JA, double* x0, double* b)
{
	double error = 1;
	int outloops = 0;
	norm_b = cblas_dnrm2(unknows, b, 1);
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
	char transa = 'N';
	mkl_cspblas_dcsrgemv(&transa, &unknows, A, IA, JA, x0, residual);

	cblas_daxpby(unknows, ONE, b, 1, NEGONE, residual, 1);    //get the residual vector

	norm_residual = cblas_dnrm2(unknows, residual, 1);

	beta[0] = norm_residual;
	double one_norm = ONE / norm_residual;

	cblas_daxpby(unknows, one_norm, residual, 1, ZERO, q, 1);    //constant*vector

	Q[0] = ArrayUtils<double>::onetensor(unknows);

	cblas_dcopy(unknows, q, 1, Q[0], 1);

	for (int k = 0; k < krylovDemension; k++)     //k represents iter
	{
		Q[k + 1] = ArrayUtils<double>::onetensor(unknows);

		mkl_cspblas_dcsrgemv(&transa, &unknows, A, IA, JA, q, v);

		for (int j = 0; j < k + 1; j++)
		{
			H[j + k * (krylovDemension + 1)] = cblas_ddot(unknows, Q[j], 1, v, 1);

			double coe = NEGONE * H[j + k * (krylovDemension + 1)];
			cblas_daxpy(unknows, coe, Q[j], 1, v, 1);
		}
		H[(k + 1) + k * (krylovDemension + 1)] = cblas_dnrm2(unknows, v, 1);

		double coe = ONE / H[(k + 1) + k * (krylovDemension + 1)];

		cblas_daxpby(unknows, coe, v, 1, ZERO, q, 1);                 //Normalize

		cblas_dcopy(unknows, q, 1, Q[k + 1], 1);

		//givens rotation for single rhs
		for (int j = 0; j < k; j++)
		{
			cblas_drot(1, &H[j + k * (krylovDemension + 1)], 1, &H[(j + 1) + k * (krylovDemension + 1)], 1, givens[j][0], givens[j][1]);
		}
		cblas_drotg(&H[k + k * (krylovDemension + 1)], &H[(k + 1) + k * (krylovDemension + 1)], &givens[k][0], &givens[k][1]);
		H[(k + 1) + k * (krylovDemension + 1)] = 0.0;

		//givens rotation for residual vector
		cblas_drot(1, &beta[k], 1, &beta[k + 1], 1, givens[k][0], givens[k][1]);

		//if convergence?
		if (abs(beta[k + 1]) < tolerance * norm_b)
		{
			Update(H, x, beta, Q, k + 1);
			cblas_dcopy(unknows, x, 1, x0, 1);

			for (int i = 0; i < k + 1; i++)
			{
				ArrayUtils<double>::delonetensor(Q[i]);
			}
			return(beta[k + 1]);
		}
	}

	Update(H, x, beta, Q, krylovDemension);
	cblas_dcopy(unknows, x, 1, x0, 1);
	for (int i = 0; i < krylovDemension - 1; i++)
	{
		ArrayUtils<double>::delonetensor(Q[i]);
	}

	return(beta[krylovDemension - 1]);
}

void GMRes::Update(double* H, double* x, double* beta, std::vector<double*> Q, int iteration)
{
	int* ipiv = ArrayUtils<int>::onetensor(iteration * iteration);
	double* tmp = ArrayUtils<double>::onetensor(unknows * iteration);
	double* tmpx = ArrayUtils<double>::onetensor(unknows * iteration);
	double* tempH = ArrayUtils<double>::onetensor(iteration * iteration);
	for (int i = 0; i < iteration; i++)
	{
		for (int j = 0; j < unknows; j++)
		{
			tmpx[i * unknows + j] = Q[i][j];
		}
	}

	for (int i = 0; i < iteration; i++)
	{
		for (int j = 0; j < iteration; j++)
		{
			tempH[i * iteration + j] = H[i * (krylovDemension + 1) + j];
		}
	}

	LAPACKE_dtrtri(LAPACK_COL_MAJOR, 'U', 'N', iteration, tempH, iteration);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, unknows, iteration, iteration, 1, tmpx, unknows, tempH, iteration, 0, tmp, unknows);
	cblas_dgemv(CblasColMajor, CblasNoTrans, unknows, iteration, 1, tmp, unknows, beta, 1, 1, x, 1);

	ArrayUtils<int>::delonetensor(ipiv);
	ArrayUtils<double>::delonetensor(tmp);
	ArrayUtils<double>::delonetensor(tmpx);
	ArrayUtils<double>::delonetensor(tempH);
}

