#pragma once
#include "HXDefine.h"
BeginNameSpace(ONEFLOW)

class TimeSpan;

class INsInv
{
public:
	INsInv();
	~INsInv();
public:
	void FluxInit();
	void MomPreInit();
	void DeleteMomPreVar();
	void PressCorInit();
	void DeletePressCorVar();

public:
	RealField flux;
	RealField spu, bu, bv, bw;
	RealField rb,ub, vb, wb, pb;
	RealField spp, bp;
	RealField dun, dup;
	RealField  ppf, pp;
	RealField2D ai;//ajp;

public:
	Real remax_up, remax_vp, remax_wp, remax_pp, res_u, res_v, res_w, res_p;
};

extern INsInv iinv;

class SIMPLEC
{
public:
	SIMPLEC();
	//~SIMPLEC();
public:
	void Run();
public:
	static bool Converge();
public:
	void SolveInnerIter();

protected:
	void InnerProcess();
	void OuterProcess(TimeSpan * timeSpan);
};

void SIMPLECSolve();

EndNameSpace