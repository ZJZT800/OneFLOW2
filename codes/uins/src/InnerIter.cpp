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
#include "UNsBcSolver.h"
#include "UINsPressCorrect.h"
#include "UINsMomPre.h"
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

//#include "UINsCorrectPress.h"
//#include "UINsCorrectSpeed.h"
#include "INsCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "INsIdx.h"
//#include "UINsInvterm.h"
//#include "UINsVisterm.h"
//#include "UINsUnsteady.h"
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
	INsPreflux();      //����������ʼ��
}

void Inner::SolveFlow()
{
	Mom_pre();

	Pres_cor();

	INsUpdateRes();
}

void Mom_pre()
{
	INsCmpConv(); //���������

	INsCmpDiffus(); //������ɢ��

	int transt = ONEFLOW::GetDataValue< int >("transt");
	if (transt != 0) INsTranst();  //˲̬��

	INsCmpSrc(); //����ѹ���ݶȺͶ�������ϵ��

	MomEqu();     //�Ͷ�����ɢ��ƽ��

	Relaxation();

	SolveMom(); //��⶯������

	INsCmpFaceflux(); //�����������
}

void Pres_cor()
{
	PresEqu(); //����ѹ����������ϵ��

	INsCmpPressCorrectEquandUpdatePress();  //��Ҫ��ѹ�����������飬���赥Ԫ����ѹ��δ֪��

	INsCmpSpeedCorrectandUpdateSpeed();  //��Ҫ��������������ٶ�δ֪�����������,���µ�Ԫ�ٶȺ�ѹ��

	INsUpdateFaceflux();   //���½�������
}

void INsPreflux()
{
	UINsMomPre * uINsMomPre = new UINsMomPre();
	uINsMomPre->CmpINsPreflux();
	delete uINsMomPre;
}

void INsCmpConv()
{
	UINsMomPre * uINsMomPre = new UINsMomPre();
	uINsMomPre->CmpConv();
	delete uINsMomPre;
}

void INsCmpDiffus()
{
	UINsMomPre * uINsMomPre = new UINsMomPre();
	uINsMomPre->CmpDiffus();
	delete uINsMomPre;
}

void INsTranst()
{
	UINsMomPre * uINsMomPre = new UINsMomPre();
	uINsMomPre->CmpTranst();
	delete uINsMomPre;
}

void INsCmpSrc()
{
	UINsMomPre * uINsMomPre = new UINsMomPre();
	uINsMomPre->CmpSrc();
	delete uINsMomPre;
}

void MomEqu()
{
	UINsMomPre * uINsMomPre = new UINsMomPre();
	uINsMomPre->MomEquCoeff();
	delete uINsMomPre;
}

void Relaxation()
{
	UINsMomPre * uINsMomPre = new UINsMomPre();
	uINsMomPre->RelaxMom(0.8);
	delete uINsMomPre;
}

void SolveMom()
{
	UINsMomPre * uINsMomPre = new UINsMomPre();
	uINsMomPre->Solveuvw();
	delete uINsMomPre;
}

void INsCmpFaceflux()
{
	UINsMomPre * uINsMomPre = new UINsMomPre();
	uINsMomPre->CmpFaceflux();
	delete uINsMomPre;
}

void PresEqu()
{
	UINsPressCorrect * uINsPressCorrect = new UINsPressCorrect();
	uINsPressCorrect->PresEquCoeff();
	delete uINsPressCorrect;
}

void INsCmpPressCorrectEquandUpdatePress()
{
	UINsPressCorrect * uINsPressCorrect = new UINsPressCorrect();
	uINsPressCorrect->CmpPressCorrectEqu();
	delete uINsPressCorrect;
}

void INsUpdateFaceflux()
{
	UINsPressCorrect * uINsPressCorrect = new UINsPressCorrect();
	uINsPressCorrect->UpdateFaceflux();
	delete uINsPressCorrect;
}

void INsCmpSpeedCorrectandUpdateSpeed()
{
	UINsPressCorrect * uINsPressCorrect = new UINsPressCorrect();
	uINsPressCorrect->UpdateSpeed();
	delete uINsPressCorrect;
}

void INsUpdateRes()
{
	UINsRes * uINsRes = new UINsRes();
	uINsRes->UpdateINsRes();
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