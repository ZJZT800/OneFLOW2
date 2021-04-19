
#include "SIMPLEC.h"
#include "INsCom.h"
#include "UCom.h"
#include "UINsCom.h"
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

INsInv iinv;

INsInv::INsInv()
{
	;
}

INsInv::~INsInv()
{
	;
}

void INsInv::FieldInit()
{
	flux.resize(ug.nFace);
	ub.resize(ug.nBFace);
	vb.resize(ug.nBFace);
	wb.resize(ug.nBFace);
	pb.resize(ug.nBFace);
	dun.resize(ug.nFace);
	fvisb_cof.resize(ug.nBFace);
}

void INsInv::MomPreInit()
{
	bu.resize(ug.nCell);
	bv.resize(ug.nCell);
	bw.resize(ug.nCell);
	spu.resize(ug.nCell);
	ai.resize(2, ug.nFace);
}

void INsInv::DeleteMomPreVar()
{
	bu.clear();
	bv.clear();
	bw.clear();
}

void INsInv::PressCorInit()
{
	bp.resize(ug.nCell);
	spp.resize(ug.nCell);
	pp.resize(ug.nCell);
	ppf.resize(ug.nFace);
	dup.resize(ug.nCell);

	iinv.spp = 0;
	iinv.bp = 0;
}

void INsInv::DeletePressCorVar()
{
	bp.clear();
	spp.clear();
	pp.clear();
	ppf.clear();
	dup.clear();
	spu.clear();
	ai.clear();
}

void INsInv::OldValueInit()
{
	u_old.resize(ug.nCell);
	v_old.resize(ug.nCell);
	w_old.resize(ug.nCell);
	
	u_old = inscom.inflow[1];
	v_old = inscom.inflow[2];
	w_old = inscom.inflow[3];

	int solve_energy = GetDataValue< int >("solve_energy");
	if (solve_energy == 1)
	{
		;
	}
}

void INsInv::DeleteOldValue()
{
	u_old.clear();
	v_old.clear();
	w_old.clear();
}


SIMPLEC::SIMPLEC()
{
	;
}

SIMPLEC::~SIMPLEC()
{
	;
}


void SIMPLECSolve()
{
	SIMPLEC * Simplec = new SIMPLEC();
	Simplec->Run();
	delete Simplec;
}

void SIMPLEC::Run()
{
	int transt = ONEFLOW::GetDataValue< int >("transt");
	if (transt != 0)
	{
		TimeSpan * timeSpan = new TimeSpan();
		// outer loop(Unsteady loop)
		while (!SimuIterState::Running())
		{
			Iteration::outerSteps++;
			ctrl.currTime += ctrl.pdt;

			ConveResInit();
			//Inner loop(Steady loop)
			Iteration::innerSteps = 0;
			while (!SIMPLEC::Converge())
			{
				Iteration::innerSteps++;

				this->SolveInnerIter();
			}

			this->SaveOldTimeValue();
			this->OuterProcess(timeSpan);
		}

		    iinv.DeleteOldValue();
		    delete timeSpan;
	}
   else
   {
		TimeSpan * timeSpan = new TimeSpan();
	    //Inner loop(Steady loop)
		ConveResInit();
	    Iteration::innerSteps = 0;
	    while (!SIMPLEC::Converge())
	    {
		   Iteration::innerSteps++;

		   this->SolveInnerIter();
	    }
		this->OuterProcess(timeSpan);
		delete timeSpan;
   }

}

void SIMPLEC::ConveResInit()
{
	iinv.remax_u = 1;
	iinv.remax_v = 1;
	iinv.remax_w = 1;
	iinv.remax_pp = 1;
}

bool SIMPLEC::Converge()
{
	int maxIterSteps = GetDataValue< int >("maxIterSteps");
	Real ConRes_u = GetDataValue< Real >("ConRes_u");
	Real ConRes_v = GetDataValue< Real >("ConRes_v");
	Real ConRes_w = GetDataValue< Real >("ConRes_w");
	Real ConRes_pp = GetDataValue< Real >("ConRes_pp");

	if (Iteration::innerSteps == 0) return false;
	if (Iteration::innerSteps < maxIterSteps && (iinv.remax_u > ConRes_u ||iinv.remax_v >ConRes_v|| iinv.remax_w > ConRes_w|| iinv.remax_pp > ConRes_pp))
	    return false;

	//if (ctrl.idualtime == 0) return true;
	bool flag = true;

	return flag;
}

void SIMPLEC::SolveInnerIter()
{
	int solve_flow = GetDataValue< int >("solve_flow");
	int solve_energy = GetDataValue< int >("solve_energy");
	int solve_turb = GetDataValue< int >("solve_turb");
	int solve_multicomponent = GetDataValue< int >("solve_multicomponent");

	Inner * uINsSolver = new Inner;

	if (solve_flow ==1)
	{
		uINsSolver->SolveFlow();
	}
	if (solve_energy == 1)
	{
		uINsSolver->SolveEnergy();
	}
	if (solve_turb ==1)
	{
		uINsSolver->SolveTurb();
	}
	if (solve_multicomponent == 1)
	{
		uINsSolver->SolveMultiComp();
	}
		
	delete uINsSolver;

	this->InnerProcess();
}

void SIMPLEC::InnerProcess()
{
	ONEFLOW::MsMgTask("POST_PROCESS");
}

void SIMPLEC::SaveOldTimeValue()
{
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		iinv.u_old[cId] = (*uinsf.u)[0][cId];
		iinv.v_old[cId] = (*uinsf.v)[0][cId];
		iinv.w_old[cId] = (*uinsf.w)[0][cId];
	}
}

void SIMPLEC::OuterProcess(TimeSpan * timeSpan)
{
	int transt = ONEFLOW::GetDataValue< int >("transt");

	if (transt == 0)
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
	else
	{
		int  maxSteps = ONEFLOW::GetDataValue< int >("maxSteps");

		if (Iteration::outerSteps == maxSteps)
		{
			if (Parallel::IsServer())
			{
				cout << "dumping field...";
				cout << "  finished " << endl;
				timeSpan->ShowTimeSpan();
			}
		}
	}
}

EndNameSpace