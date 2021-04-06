#include "BiCGStab.h"

BeginNameSpace(ONEFLOW)
BiCGStab::BiCGStab(RealField& sp, RealField2D& aii, RealField& rhs, RealField& x, int unknowns, int Iter, Real Tol, bool ifPrecond)
{
	if (ifPrecond)
	{
		SolvePrecond(sp, aii, rhs, x, unknowns, Iter, Tol);
	}
	else
	{
		Solve(sp, aii, rhs, x, unknowns, Iter, Tol);
	}
}

BiCGStab::~BiCGStab()
{
}

void BiCGStab::InitPre(int unknowns)
{
	r.resize(unknowns);
	p.resize(unknowns);
	v.resize(unknowns);
	t.resize(unknowns);
	s.resize(unknowns);
	rHat.resize(unknowns);
	r_temp.resize(unknowns);
	apTmp.resize(unknowns);
}

void BiCGStab::Init(int unknowns)
{
	r.resize(unknowns);
	p.resize(unknowns);
	v.resize(unknowns);
	t.resize(unknowns);
	s.resize(unknowns);
	rHat.resize(unknowns);
	r_temp.resize(unknowns);
}

void BiCGStab::Solve(RealField& sp, RealField2D& aii, RealField& rhs, RealField& x, int unknowns, int MaxIter, Real Tol)
{
	int i, j;
	Real alpha, rho1, beta, omega;
	Init(unknowns);

	r = rhs;
	p = r;
	rHat = r;

	rho1 = basic.dot(rhs, rhs, unknowns);
	Real res = sqrt(rho1);
	int iter = 0;
	while (res > Tol && iter < MaxIter)
	{
		basic.MVs(sp, aii, p, v);
		alpha = basic.dot(rHat, r, unknowns) / basic.dot(rHat, v, unknowns);

		for (int i = 0; i < unknowns; ++i)
		{
			s[i] = r[i] - alpha * v[i];
		}
		basic.MVs(sp, aii, s, t);
		omega = basic.dot(t, s, unknowns) / basic.dot(t, t, unknowns);

		for (int i = 0; i < unknowns; ++i)
		{
			x[i] = x[i] + alpha * p[i] + omega * s[i];
			r_temp[i] = r[i];
			r[i] = s[i] - omega * t[i];
		}

		res = basic.dot(r, r, unknowns);
		res = sqrt(abs(res / rho1));

		if (res > Tol)
		{
			beta = basic.dot(rHat, r, unknowns);
			beta = beta / basic.dot(rHat, r_temp, unknowns);
			beta = beta * (alpha / omega);
			for (int i = 0; i < unknowns; ++i)
			{
				p[i] = r[i] + beta * (p[i] - omega * s[i]);
			}
		}
		else
		{
			continue;
		}

		iter++;
	}
}

void BiCGStab::SolvePrecond(RealField& sp, RealField2D& aii, RealField& rhs, RealField& x, int unknowns, int MaxIter, Real Tol)
{
	int i, j;
	Real alpha, rho1, beta, omega;
	InitPre(unknowns);

	PrecondILU(sp, aii, apTmp, unknowns);
	PrecondRU(aii, apTmp, rhs, unknowns);

	r = rhs;
	p = r;
	rHat = r;

	rho1 = basic.dot(rhs, rhs, unknowns);
	Real res = sqrt(rho1);
	int iter = 0;
	while (res > Tol && iter < MaxIter)
	{
		basic.MVs(sp, aii, p, v);
		PrecondRU(aii, apTmp, v, unknowns);
		alpha = basic.dot(rHat, r, unknowns) / basic.dot(rHat, v, unknowns);

		for (int i = 0; i < unknowns; ++i)
		{
			s[i] = r[i] - alpha * v[i];
		}
		basic.MVs(sp, aii, s, t);
		PrecondRU(aii, apTmp, t, unknowns);
		omega = basic.dot(t, s, unknowns) / basic.dot(t, t, unknowns);

		for (int i = 0; i < unknowns; ++i)
		{
			x[i] = x[i] + alpha * p[i] + omega * s[i];
			r_temp[i] = r[i];
			r[i] = s[i] - omega * t[i];
		}

		res = basic.dot(r, r, unknowns);
		res = sqrt(abs(res / rho1));

		if (res > Tol)
		{
			beta = basic.dot(rHat, r, unknowns);
			beta = beta / basic.dot(rHat, r_temp, unknowns);
			beta = beta * (alpha / omega);
			for (int i = 0; i < unknowns; ++i)
			{
				p[i] = r[i] + beta * (p[i] - omega * s[i]);
			}
		}
		else
		{
			continue;
		}

		iter++;
	}
}

void BiCGStab::PrecondILU(RealField& sp, RealField2D& aii, RealField& apTmp, int unknowns)
{
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		apTmp[cId] = sp[cId];
	}
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];
		apTmp[rc] = apTmp[rc] - aii[fId][0] * aii[fId][1] / apTmp[lc];
	}
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		apTmp[cId] = 1 / apTmp[cId];
	}
}

void BiCGStab::PrecondRU(RealField2D& aii, RealField& apTmp, RealField& b, int unknowns)
{
	RealField bTmp;
	bTmp.resize(unknowns);

	//forward
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		bTmp[cId] = b[cId];
	}
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];
		bTmp[rc] = bTmp[rc] + aii[fId][1] * apTmp[lc] * bTmp[lc];
	}
	//backforward
	for (int cId = 0; cId < ONEFLOW::ug.nCell; ++cId)
	{
		b[cId] = bTmp[cId] * apTmp[cId];
	}
	for (int fId = ug.nFace - 1; fId >= ug.nBFace; --fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];
		b[lc] = b[lc] + aii[fId][0] * apTmp[lc] * b[rc];
	}
}

EndNameSpace