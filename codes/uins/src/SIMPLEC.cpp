
#include "SIMPLEC.h"
#include "INsInvterm.h"
//#include "UINsInvterm.h"
#include "Mesh.h"
#include "Ctrl.h"
#include "Solver.h"
#include "InnerIter.h"
#include "Iteration.h"
#include "DataBase.h"
#include "UINsSolver.h"
#include "Zone.h"
#include "TimeIntegral.h"
#include "Parallel.h"
#include "SolverState.h"
#include "GridState.h"
#include "CmxTask.h"
#include "BgField.h"
#include "TimeSpan.h"
#include "UINsSolver.h"
#include <iostream>
#include <Rhs.h>

using namespace std;


BeginNameSpace(ONEFLOW)

SIMPLEC::SIMPLEC()
{
	;
}

/*SIMPLEC::~SIMPLEC()
{
	;
}*/


void SIMPLECSolve()
{
	SIMPLEC * Simplec = new SIMPLEC();
	Simplec->Run();
	delete Simplec;
}

void SIMPLEC::Run()
{
	/*double rhs_u = 1e-8;
	double rhs_v = 1e-8;
	double rhs_w = 1e-8;

	iinv.remax_up = 1;
	iinv.remax_vp = 1;
	iinv.remax_wp = 1;*/

	/*int transt = ONEFLOW::GetDataValue< int >("transt");
	if (transt != 0)
	{
		TimeSpan * timeSpan = new TimeSpan();
		while (SimuIterState::Running())
		{
			Iteration::outerSteps++;
			ctrl.currTime += ctrl.pdt;

			int maxIterSteps = GetDataValue< int >("maxIterSteps");
			while (iinv.remax_up > rhs_u || iinv.remax_vp > rhs_v || iinv.remax_wp > rhs_w)
			{
				if (Iteration::innerSteps >= maxIterSteps) break;

				Iteration::innerSteps++;

				this->SolveInnerIter();

			}
			this->OuterProcess(timeSpan);
		}
		delete timeSpan;
	}*/
	//else
	//{
	//iinv.remax_up > rhs_u || iinv.remax_vp > rhs_v || iinv.remax_wp > rhs_w
		//int maxIterSteps = GetDataValue< int >("maxIterSteps");

	    Iteration::innerSteps = 0;
		while (!SIMPLEC::Converge())
		{
			Iteration::innerSteps++;

			this->SolveInnerIter();

		}
	//}
}

bool SIMPLEC::Converge()
{
	int maxIterSteps = GetDataValue< int >("maxIterSteps");
	if (Iteration::innerSteps == 0) return false;
	if (Iteration::innerSteps < maxIterSteps) return false;
	//if (ctrl.idualtime == 0) return true;
	bool flag = true;

	return flag;
}

void SIMPLEC::SolveInnerIter()
{
	Inner * uINsSolver = new Inner;
	uINsSolver->SolveFlow();
	delete uINsSolver;

	this->InnerProcess();
}

void SIMPLEC::InnerProcess()
{
	ONEFLOW::MsMgTask("POST_PROCESS");
}

void SIMPLEC::OuterProcess(TimeSpan * timeSpan)
{
	if (Iteration::outerSteps % Iteration::nFieldSave == 0)
	{
		if (Parallel::IsServer())
		{
			cout << "dumping field...";
			cout << "  finished " << endl;
			timeSpan->ShowTimeSpan();
		}
	}
}

EndNameSpace