/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "InnerIter.h"
#include "UINsInitField.h"
#include "UNsBcSolver.h"
#include "UINsConvTerm.h"
#include "UINsDiffusTerm.h"
#include "UINsBcTerm.h"
#include "UINsTranst.h"
#include "UINsSrcTerm.h"
#include "UINsPrimeEqu.h"
#include "UINsRelax.h"
#include "UINsSolveVar.h"
#include "UINsCmpFlux.h"
#include "UINsPressEqu.h"
#include "UINsMomCorrect.h"
#include "UINsRes.h"
#include "Zone.h"
#include "DataBase.h"
#include "Iteration.h"
#include "NsCom.h"
#include "UCom.h"
#include "UNsCom.h"
#include "UnsGrid.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "UNsInvFlux.h"
#include "UNsVisFlux.h"
#include "UNsUnsteady.h"
#include "Ctrl.h"
#include "INsCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "INsIdx.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

Inner::Inner()
{
    ;
}

Inner::~Inner()
{
    ;
}


/*void INSCmpGamaT(int flag)
{
	UnsGrid * grid = Zone::GetUnsGrid();

	ug.Init();
	uinsf.Init();
	ug.SetStEd(flag);

	if (inscom.chemModel == 1)
	{
	}
	else
	{
		//if ( ZoneState::zid == 0 )
		//{
		//    cout << " ug.ist = " << ug.ist  << " ug.ied = " << ug.ied << "\n";
		//    int kkk = 1;
		//}
		Real oamw = one;
		for (int cId = ug.ist; cId < ug.ied; ++cId)
		{
			Real & density = ( * uinsf.q )[ IIDX::IIR ][ cId ];
			Real & pressure = ( * uinsf.q )[ IIDX::IIP ][ cId ];
			//( * uinsf.gama )[ 0 ][ cId ] = inscom.gama_ref;
			//( * uinsf.tempr )[ IIDX::IITT ][ cId ] = pressure / ( inscom.statecoef * density * oamw );
			//(*uinsf.tempr)[IIDX::IITT][cId] = 0;
		}
	}
}*/

void Inner::FieldInit()
{
	InitField();      //流场变量初始化
}

void Inner::SolveFlow()
{
	Mom_pre();

	Pres_cor();

	UpdateRes();
}

void Inner::SolveEnergy()
{
	;
}

void Inner::SolveTurb()
{
	;
}

void Inner::SolveMultiComp()
{
	;
}

void Mom_pre()
{
	iinv.MomPreInit();

	CmpMomConv(); //计算对流项

	CmpMomDiffus(); //计算扩散项

	CmpMomBc();

	int transt = ONEFLOW::GetDataValue< int >("transt");
	if (transt != 0) MomTranst();  //瞬态项

	CmpMomSrc(); //计算压力梯度和动量方程系数

	MomPrimeEqu();     //和对流扩散项平级

	MomRelaxation();

	SolveMom(); //求解动量方程

	CmpFaceflux(); //计算界面流量

	iinv.DeleteMomPreVar();
}

void Pres_cor()
{
	iinv.PressCorInit();

	PresEqu(); //计算压力修正方程系数

	SolvePress();  //需要解压力修正方程组，增设单元修正压力未知量

	UpdateCorSpeed();  //需要先增设界面修正速度未知量并进行求解,更新单元速度和压力

	UpdateFaceflux();   //更新界面流量

	iinv.DeletePressCorVar();
}

void InitField()
{
	UINsInitField * uINsInitField = new UINsInitField();
	uINsInitField->InitInputField();
	delete uINsInitField;
}

void CmpMomConv()
{
	string vary = "mom";
	string conv_ischeme = ONEFLOW::GetDataValue< string >("conv_ischeme");
	UINsConvTerm * uINsConvTerm = new UINsConvTerm();
	uINsConvTerm->CmpConvTerm(vary, conv_ischeme);
	delete uINsConvTerm;
}

void CmpMomDiffus()
{
	string vary = "mom";
	string diffus_ischeme = ONEFLOW::GetDataValue< string >("diffus_ischeme");
	UINsDiffusTerm * uINsDiffusTerm = new UINsDiffusTerm;
	uINsDiffusTerm->CmpDiffusTerm(vary, diffus_ischeme);
	delete uINsDiffusTerm;
}

void CmpMomBc()
{
	string vary = "mom";
	UINsBcTerm * uINsBcTerm = new UINsBcTerm;
	uINsBcTerm->CmpMomBcTerm(vary);
	delete uINsBcTerm;
}

void MomTranst()
{
	string vary = "mom";
	UINsTranst * uINsTranst = new UINsTranst();
	uINsTranst->CmpTranstTerm(vary);
	delete uINsTranst;
}

void CmpMomSrc()
{
	string vary = "mom";
	UINsSrcTerm * uINsSrcTerm = new UINsSrcTerm();
	uINsSrcTerm->CmpMomSrcTerm(vary);
	delete uINsSrcTerm;
}

void MomPrimeEqu()
{
	string vary = "mom";
	UINsPrimeEqu * uINsPrimeEqu = new UINsPrimeEqu();
	uINsPrimeEqu->PrimeEqu(vary);
	delete uINsPrimeEqu;
}

void MomRelaxation()
{
	string vary = "mom";
	Real mom_relax = GetDataValue< Real >("mom_relax");
	UINsRelax * uINsRelax = new UINsRelax();
	uINsRelax->Relax(vary,mom_relax);
	delete uINsRelax;
}

void SolveMom()
{
	string vary = "mom";
	UINsSolveVar * uINsSolveVar = new UINsSolveVar();
	uINsSolveVar->SolveVar(vary);
	delete uINsSolveVar;
}

void CmpFaceflux()
{
	UINsCmpFlux * uINsCmpFlux = new UINsCmpFlux();
	uINsCmpFlux->CmpFlux();
	delete uINsCmpFlux;
}

void PresEqu()
{
	UINsPressEqu * uINsPressEqu = new UINsPressEqu();
	uINsPressEqu->PressEqu();
	delete uINsPressEqu;
}

void SolvePress()
{
	string vary = "press";
	UINsSolveVar * uINsSolveVar = new UINsSolveVar();
	uINsSolveVar->SolveVar(vary);
	delete uINsSolveVar;
}

void UpdateCorSpeed()
{
	UINsMomCorrect * uINsMomCorrect = new UINsMomCorrect();
	uINsMomCorrect->UpdateCorrectSpeed();
	delete uINsMomCorrect;
}

void UpdateFaceflux()
{
	UINsMomCorrect * uINsMomCorrect = new UINsMomCorrect();
	uINsMomCorrect->UpdateCorrectFaceflux();
	delete uINsMomCorrect;
}

void UpdateRes()
{
	UINsRes * uINsRes = new UINsRes();
	uINsRes->UpdateIterRes();
	delete uINsRes;
}

//void INsCorrectSpeed()
//{
//	UINsInvterm * uINsInvterm = new UINsInvterm();
//	uINsInvterm->CmpCorrectSpeed();
//	delete uINsInvterm;
//}



/*void INsCmpChemSrc()
{
	;
}*/

/*void INsCmpTurbEnergy()
{
	;
}*/

//void INsCmpDualTimeStepSrc()
//{
//	UINsUnsteady * uinsUnsteady = new UINsUnsteady();
//	uinsUnsteady->CmpDualTimeSrc();
//	delete uinsUnsteady;
//}


EndNameSpace