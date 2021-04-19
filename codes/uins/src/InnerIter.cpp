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
#include "EquaDiscrete.h"
#include "SIMPLEC.h"
#include "systemSolver.h"
#include "UpdateVary.h"
#include "FaceValue.h"
#include "UINsInitField.h"
#include "UNsBcSolver.h"
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
	//Discrete format,mom relax coefficient,transt, Parameter Setting of Solution Equations
	string conv_ischeme = ONEFLOW::GetDataValue< string >("conv_ischeme");
	string diffus_ischeme = ONEFLOW::GetDataValue< string >("diffus_ischeme");
	Real mom_relax = GetDataValue< Real >("mom_relax");
	int transt = ONEFLOW::GetDataValue< int >("transt");

	int MaxIter = GetDataValue< int >("EquaMomIter");
	Real Tol = GetDataValue< Real >("EquaMomTol");
	int mRestarts = GetDataValue< int >("EquaMomRestarts");
	string gt = GetDataValue< string >("EquaMomMethod");
	bool ifPrecond = true;

	RealField diff_value;
	diff_value.resize(ug.nCell);
	RealField fdiffus_cof;
	fdiffus_cof.resize(ug.nFace);

	iinv.MomPreInit(); //allocation  mom storage

	string EquaVary = "vary_u";
	FvisCoeff(EquaVary, (*uinsf.vis_coef)[0], iinv.fvisb_cof, fdiffus_cof);
	EquaDisc(EquaVary, (*uinsf.u)[0], iinv.ub, iinv.u_old, iinv.flux, (*uinsf.p)[0], iinv.pb, fdiffus_cof, conv_ischeme, diffus_ischeme, mom_relax, transt, iinv.spu, iinv.ai, iinv.bu, iinv.remax_u);
	SolveEqua(iinv.spu, iinv.ai, iinv.bu, diff_value, iinv.res_u, gt, MaxIter, Tol, mRestarts, ifPrecond);
	UpateVar((*uinsf.u)[0], diff_value);
	diff_value = 0;
	iinv.spu = 0;
	iinv.ai[0] = 0;
	iinv.ai[1] = 0;

	EquaVary = "vary_v";
	EquaDisc(EquaVary, (*uinsf.v)[0], iinv.vb, iinv.v_old, iinv.flux, (*uinsf.p)[0], iinv.pb, fdiffus_cof, conv_ischeme, diffus_ischeme, mom_relax, transt, iinv.spu, iinv.ai, iinv.bv, iinv.remax_v);
	SolveEqua(iinv.spu, iinv.ai, iinv.bv, diff_value, iinv.res_v, gt, MaxIter, Tol, mRestarts, ifPrecond);
	UpateVar((*uinsf.v)[0], diff_value);
	diff_value = 0;
	iinv.spu = 0;
	iinv.ai[0] = 0;
	iinv.ai[1] = 0;

	EquaVary = "vary_w";
	EquaDisc(EquaVary, (*uinsf.w)[0], iinv.wb, iinv.w_old, iinv.flux, (*uinsf.p)[0], iinv.pb, fdiffus_cof, conv_ischeme, diffus_ischeme, mom_relax, transt, iinv.spu, iinv.ai, iinv.bw, iinv.remax_w);
	SolveEqua(iinv.spu, iinv.ai, iinv.bw, diff_value, iinv.res_w, gt, MaxIter, Tol, mRestarts, ifPrecond);
	UpateVar((*uinsf.w)[0], diff_value);

	iinv.DeleteMomPreVar(); //free mom storage

	UINsCmpFlux((*uinsf.rho)[0], (*uinsf.u)[0], iinv.ub, (*uinsf.v)[0], iinv.vb, (*uinsf.w)[0], iinv.wb, (*uinsf.p)[0], iinv.pb, iinv.spu, iinv.flux, iinv.dun);
}

void Pres_cor()
{
	iinv.PressCorInit();

	UINsPressEqu((*uinsf.rho)[0], iinv.flux, iinv.spu, iinv.ai, iinv.dup, iinv.spp, iinv.bp, iinv.remax_pp);

	int MaxIter = GetDataValue< int >("EquaPressIter");
	Real Tol = GetDataValue< Real >("EquaPressTol");
	int mRestarts = GetDataValue< int >("EquaPressRestarts");
	string gt = GetDataValue< string >("EquaPressMethod");
	bool ifPrecond = true;
	SolveEqua(iinv.spp, iinv.ai, iinv.bp, iinv.pp, iinv.res_p, gt, MaxIter, Tol, mRestarts, ifPrecond);

	UpdateCorPressAndSpeed();  //Update correction press and speed 

	UpdateFaceflux();   //Update face flux

	iinv.DeletePressCorVar();
}


void InitField()
{
	UINsInitField * uINsInitField = new UINsInitField();
	uINsInitField->InitInputField();
	delete uINsInitField;
}

void UpdateCorPressAndSpeed()
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

EndNameSpace