#pragma once
#include "HXDefine.h"
BeginNameSpace(ONEFLOW)

class TimeSpan;

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