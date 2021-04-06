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
	void OldValueInit();

public:
	RealField flux;
	RealField spu, bu, bv, bw;
	RealField ub, vb, wb, pb, fvis_cof, fvisb_cof;
	RealField u_old, v_old, w_old;
	RealField spp, bp;
	RealField dun, dup;
	RealField  ppf, pp;
	RealField2D ai;

public:
	Real remax_u, remax_v, remax_w, remax_pp, res_u, res_v, res_w, res_p;
};

extern INsInv iinv;

class SIMPLEC
{
public:
	SIMPLEC();
	~SIMPLEC();
public:
	void Run();
public:
	void ConveResInit();
	static bool Converge();
public:
	void SolveInnerIter();

protected:
	void InnerProcess();
	void OuterProcess(TimeSpan * timeSpan);
	void SaveOldTimeValue();
};

void SIMPLECSolve();

EndNameSpace