#include "BasicComputation.h"
#include <iomanip>

BeginNameSpace(ONEFLOW)

BasicCompute::BasicCompute()
{
}

BasicCompute::~BasicCompute()
{
}

void BasicCompute::MVs(double* sp, double** aii, double* b, double* ans)
{
	for (int cId = 0; cId < ONEFLOW::ug.nCell; ++cId)
	{
		ans[cId] = 0.0;
		ans[cId] = sp[cId] * b[cId];
	}
	for (int fId = ONEFLOW::ug.nBFace; fId < ONEFLOW::ug.nFace; ++fId)
	{
		int lc = (*ONEFLOW::ug.lcf)[fId];
		int rc = (*ONEFLOW::ug.rcf)[fId];
		ans[lc] -= b[rc] * aii[fId][0];
		ans[rc] -= b[lc] * aii[fId][1];
	}
}

void BasicCompute::MVs(RealField& sp, RealField2D& aii, RealField& b, RealField& ans)
{
	for (int cId = 0; cId < ONEFLOW::ug.nCell; ++cId)
	{
		ans[cId] = 0.0;
		ans[cId] = sp[cId] * b[cId];
	}
	for (int fId = ONEFLOW::ug.nBFace; fId < ONEFLOW::ug.nFace; ++fId)
	{
		int lc = (*ONEFLOW::ug.lcf)[fId];
		int rc = (*ONEFLOW::ug.rcf)[fId];
		ans[lc] -= b[rc] * aii[fId][0];
		ans[rc] -= b[lc] * aii[fId][1];
	}
}

void BasicCompute::MVs(ONEFLOW::RealField* sp, ONEFLOW::RealField2D* aii, double* b, double* ans)
{
	for (int cId = 0; cId < ONEFLOW::ug.nCell; ++cId)
	{
		ans[cId] = 0.0;
		ans[cId] = (*sp)[cId] * b[cId];
	}
	for (int fId = ONEFLOW::ug.nBFace; fId < ONEFLOW::ug.nFace; ++fId)
	{
		int lc = (*ONEFLOW::ug.lcf)[fId];
		int rc = (*ONEFLOW::ug.rcf)[fId];
		ans[lc] -= b[rc] * (*aii)[fId][0];
		ans[rc] -= b[lc] * (*aii)[fId][1];
	}
}

void BasicCompute::CSRMVs(double* A, int* IA, int* JA, double* b, double* ans, int rows)
{
	int cols, rowIndex, colIndex, frontEntry;
	double* temp = ArrayUtils<double>::onetensor(rows);
	for (int i = 1; i < rows + 1; ++i)
	{
		rowIndex = i - 1;
		cols = IA[i] - IA[i - 1];
		frontEntry = IA[i - 1];
		for (int j = 0; j < cols; ++j)
		{
			colIndex = JA[frontEntry + j];
			temp[rowIndex] += A[frontEntry + j] * b[colIndex];
		}
		ans[rowIndex] = temp[rowIndex];
	}
	ArrayUtils<double>::delonetensor(temp);
	temp = NULL;
}

void BasicCompute::apxy(double* b, double* c, int rows, double realnum)
{
	for (int i = 0; i < rows; ++i)
	{
		b[i] = realnum * b[i] + b[i] * c[i];
	}
}

void BasicCompute::apxy(double* b, double* c, int rows)
{
	for (int i = 0; i < rows; ++i)
	{
		b[i] += b[i] * c[i];
	}
}

double BasicCompute::dot(double* b, double* c, int rows)
{
	double ans = 0.0;
	for (int i = 0; i < rows; ++i)
	{
		ans += b[i] * c[i];
	}
	return(ans);
}

Real BasicCompute::dot(RealField& b, RealField& c, int rows)
{
	Real ans = 0.0;
	for (int i = 0; i < rows; ++i)
	{
		ans += b[i] * c[i];
	}
	return(ans);
}

double BasicCompute::norm(double* b, int rows)
{
	double ans = 0.0;
	for (int i = 0; i < rows; ++i)
	{
		ans += b[i] * b[i];
	}
	return(sqrt(ans));
}

double BasicCompute::norm(ONEFLOW::RealField* b, int rows)
{
	double ans = 0.0;
	for (int i = 0; i < rows; ++i)
	{
		ans += (*b)[i] * (*b)[i];
	}
	return(sqrt(ans));
}

void BasicCompute::VCs(double* b, double* c, double realnum, int rows)
{
	for (int i = 0; i < rows; ++i)
	{
		c[i] = b[i] * realnum;
	}
}

void BasicCompute::vecPlus(double* b, double* c, int rows)
{
	for (int i = 0; i < rows; ++i)
	{
		c[i] = b[i] + c[i];
	}
}

void BasicCompute::vecPlus(double* b, ONEFLOW::RealField* c, int rows)
{
	for (int i = 0; i < rows; ++i)
	{
		(*c)[i] = b[i] + (*c)[i];
	}
}

void BasicCompute::vecPlus(double* b, double* c, int rows, double realnum)
{
	for (int i = 0; i < rows; ++i)
	{
		c[i] = b[i] + realnum * c[i];
	}
}

void BasicCompute::vecPlus(ONEFLOW::RealField* b, double* c, int rows, double realnum)
{
	for (int i = 0; i < rows; ++i)
	{
		c[i] = (*b)[i] + realnum * c[i];
	}
}

void BasicCompute::vecPlus(double* b, double* c, double* d, int rows, double realnum)
{
	for (int i = 0; i < rows; ++i)
	{
		d[i] = b[i] + realnum * c[i];
	}
}

void BasicCompute::vecPlusb(double* b, double* c, int rows, double realnum)
{
	for (int i = 0; i < rows; ++i)
	{
		b[i] = b[i] + realnum * c[i];
	}
}

void BasicCompute::vecPlusb(ONEFLOW::RealField* b, double* c, int rows, double realnum)
{
	for (int i = 0; i < rows; ++i)
	{
		(*b)[i] = (*b)[i] + realnum * c[i];
	}
}

void BasicCompute::vecCopy(double* b, double* c, int rows)
{
	for (int i = 0; i < rows; ++i)
	{
		b[i] = c[i];
	}
}

void BasicCompute::rotg(double* X, double* Y, double* c, double* s)
{
	double temp;
	if (*Y == 0.0)
	{
		*c = 1.0;
		*s = 0.0;
	}
	else if (abs(*X) > abs(*Y))
	{
		temp = *Y / *X;
		*c = 1.0 / sqrt(1.0 + temp * temp);
		*s = temp * (*c);
	}
	else if (abs(*X) < abs(*Y))
	{
		temp = *X / *Y;
		*s = 1.0 / sqrt(1.0 + temp * temp);
		*c = temp * (*s);
	}
	double tmp = *X * (*c) + *Y * (*s);
	*Y = -(*s) * (*X) + (*c) * (*Y);
	*X = tmp;
}

void BasicCompute::rot(double* X, double* Y, double* c, double* s)
{
	double tmp = (*X) * (*c) + (*Y) * (*s);
	(*Y) = -(*s) * (*X) + (*c) * (*Y);
	(*X) = tmp;
}

EndNameSpace
