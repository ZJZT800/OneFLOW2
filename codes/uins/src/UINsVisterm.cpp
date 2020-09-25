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

#include "UINsVisterm.h"
#include "INsInvterm.h"
#include "UINsInvterm.h"
#include "INsVisterm.h"
#include "Iteration.h"
#include "HeatFlux.h"
#include "UGrad.h"
#include "Zone.h"
#include "BcData.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "INsCtrl.h"
#include "UCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "VisGrad.h"
#include "UINsGrad.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "ULimiter.h"
#include "UINsLimiter.h"
#include "FieldImp.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )


UINsVisterm::UINsVisterm()
{
	;
}

UINsVisterm::~UINsVisterm()
{
    ;
}


void UINsVisterm::CmpViscoff()
{
    if ( vis_model.vismodel == 0 ) return;
    ug.Init();
    uinsf.Init();
    visQ.Init( inscom.nEqu );

    Alloc();

    this->PrepareField();
    this->CmpVisterm();

    DeAlloc();
}

void UINsVisterm::Alloc()
{
	uinsf.qf = new MRField( inscom.nEqu, ug.nFace );
}

void UINsVisterm::DeAlloc()
{
    delete uinsf.qf;
}

void UINsVisterm::PrepareField()
{

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		(*uinsf.qf)[IIDX::IIU][ug.fId] = iinv.uf[ug.fId];
		(*uinsf.qf)[IIDX::IIV][ug.fId] = iinv.vf[ug.fId];
		(*uinsf.qf)[IIDX::IIW][ug.fId] = iinv.wf[ug.fId];
	}
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real dxl = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dyl = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dzl = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real dxr = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
		Real dyr = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
		Real dzr = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

		Real delt1 = DIST(dxl, dyl, dzl);
		Real delt2 = DIST(dxr, dyr, dzr);
		Real delta = 1.0 / (delt1 + delt2);

		Real cl = delt2 * delta;
		Real cr = delt1 * delta;

		(*uinsf.qf)[IIDX::IIU][ug.fId] = cl * (*uinsf.q)[IIDX::IIU][ug.lc] + cr * (*uinsf.q)[IIDX::IIU][ug.rc];
		(*uinsf.qf)[IIDX::IIV][ug.fId] = cl * (*uinsf.q)[IIDX::IIV][ug.lc] + cr * (*uinsf.q)[IIDX::IIV][ug.rc];
		(*uinsf.qf)[IIDX::IIW][ug.fId] = cl * (*uinsf.q)[IIDX::IIW][ug.lc] + cr * (*uinsf.q)[IIDX::IIW][ug.rc];

	}

	ONEFLOW::CmpINsGrad((*uinsf.qf)[IIDX::IIU], (*uinsf.dqdx)[IIDX::IIU], (*uinsf.dqdy)[IIDX::IIU], (*uinsf.dqdz)[IIDX::IIU]);
	ONEFLOW::CmpINsGrad((*uinsf.qf)[IIDX::IIV], (*uinsf.dqdx)[IIDX::IIV], (*uinsf.dqdy)[IIDX::IIV], (*uinsf.dqdz)[IIDX::IIV]);
	ONEFLOW::CmpINsGrad((*uinsf.qf)[IIDX::IIW], (*uinsf.dqdx)[IIDX::IIW], (*uinsf.dqdy)[IIDX::IIW], (*uinsf.dqdz)[IIDX::IIW]);

}

void UINsVisterm::CmpPreandVisGrad()
{
	(*uinsf.dqdx)[IIDX::IIR] = 0;
	(*uinsf.dqdy)[IIDX::IIR] = 0;
	(*uinsf.dqdz)[IIDX::IIR] = 0;

	(*uinsf.dqdx)[IIDX::IIU] = 0;
	(*uinsf.dqdy)[IIDX::IIU] = 0;
	(*uinsf.dqdz)[IIDX::IIU] = 0;

	(*uinsf.dqdx)[IIDX::IIV] = 0;
	(*uinsf.dqdy)[IIDX::IIV] = 0;
	(*uinsf.dqdz)[IIDX::IIV] = 0;

	(*uinsf.dqdx)[IIDX::IIW] = 0;
	(*uinsf.dqdy)[IIDX::IIW] = 0;
	(*uinsf.dqdz)[IIDX::IIW] = 0;

	(*uinsf.dqdx)[IIDX::IIP] = 0;
	(*uinsf.dqdy)[IIDX::IIP] = 0;
	(*uinsf.dqdz)[IIDX::IIP] = 0;


	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		if (fId == 432)
		{
			int kkk = 1;
		}
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real dxl = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dyl = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dzl = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real dxr = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
		Real dyr = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
		Real dzr = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

		Real delt1 = DIST(dxl, dyl, dzl);
		Real delt2 = DIST(dxr, dyr, dzr);
		Real delta = 1.0 / (delt1 + delt2);

		Real cl = delt2 * delta;
		Real cr = delt1 * delta;

		Real value1 = cl * (*uinsf.q)[IIDX::IIU][ug.lc] + cr * (*uinsf.q)[IIDX::IIU][ug.rc];
		Real value2 = cl * (*uinsf.q)[IIDX::IIV][ug.lc] + cr * (*uinsf.q)[IIDX::IIV][ug.rc];
		Real value3 = cl * (*uinsf.q)[IIDX::IIW][ug.lc] + cr * (*uinsf.q)[IIDX::IIW][ug.rc];
		Real value4 = cl * (*uinsf.q)[IIDX::IIP][ug.lc] + cr * (*uinsf.q)[IIDX::IIP][ug.rc];

		/*Real value1 = iinv.uf[ug.fId];
		Real value2 = iinv.vf[ug.fId];
		Real value3 = iinv.wf[ug.fId];
		Real value4 = iinv.pf[ug.fId];*/

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		(*uinsf.dqdx)[IIDX::IIU][ug.lc] += fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.lc] += fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.lc] += fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.lc] += fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.lc] += fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.lc] += fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.lc] += fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.lc] += fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.lc] += fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.lc] += fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.lc] += fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.lc] += fnza * value4;

		//if (ug.fId < ug.nBFace) continue;

		(*uinsf.dqdx)[IIDX::IIU][ug.rc] += -fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.rc] += -fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.rc] += -fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.rc] += -fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.rc] += -fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.rc] += -fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.rc] += -fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.rc] += -fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.rc] += -fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] += -fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.rc] += -fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.rc] += -fnza * value4;
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		if (fId == 432)
		{
			int kkk = 1;
		}
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real value1 = iinv.uf[ug.fId];
		Real value2 = iinv.vf[ug.fId];
		Real value3 = iinv.wf[ug.fId];
		Real value4 = iinv.pf[ug.fId];

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		(*uinsf.dqdx)[IIDX::IIU][ug.lc] += fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.lc] += fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.lc] += fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.lc] += fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.lc] += fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.lc] += fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.lc] += fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.lc] += fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.lc] += fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.lc] += fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.lc] += fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.lc] += fnza * value4;
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		Real ovol = one / (*ug.cvol)[ug.cId];
		(*uinsf.dqdx)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIP][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIP][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIP][ug.cId] *= ovol;
	}

}


void UINsVisterm::CmpVisterm()
{
    for ( int fId = ug.nBFace; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;

        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        if ( fId == 147489 )
        {
            int kkk = 1;
        }

        this->CmpFaceVisterm();  //要改动

    }
	ug.nRegion = ug.bcRecord->bcInfo->bcType.size();
	BcInfo * bcInfo = ug.bcRecord->bcInfo;

	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
		{
			ug.bcfId = ibc;

			BcInfo * bcInfo = ug.bcRecord->bcInfo;

			ug.fId = bcInfo->bcFace[ug.ir][ibc];
			ug.bcr = bcInfo->bcRegion[ug.ir][ibc];

			ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			inscom.bcdtkey = 0;
			if (ug.bcr == -1) return; //interface
			int dd = bcdata.r2d[ug.bcr];
			if (dd != -1)
			{
				inscom.bcdtkey = 1;
				inscom.bcflow = &bcdata.dataList[dd];
			}

			this->CmpBcFaceVisterm();
		}
	}

}

void UINsVisterm::CmpFaceVisterm()
{

	iinv.l2rdx = (*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc];  //界面左右单元中心距
	iinv.l2rdy = (*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc];
	iinv.l2rdz = (*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc];

	iinv.c2d = sqrt(iinv.l2rdx * iinv.l2rdx + iinv.l2rdy * iinv.l2rdy + iinv.l2rdz * iinv.l2rdz);

	iinv.chalength = GetDataValue< Real >("chara_length");

	iinv.vis = 1* iinv.chalength / inscom.reynolds;  //动力粘度

	iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	iinv.Fn[ug.fId] = iinv.vis * (*ug.farea)[ug.fId]  / (iinv.dist);

	iinv.Puf = (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*(*ug.zfn)[ug.fId];  //▽q*n

	iinv.Pvf = (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*(*ug.zfn)[ug.fId];
	
	iinv.Pwf = (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId];
	
	iinv.Pdu = -((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*iinv.l2rdx+
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*iinv.l2rdy+
	 	                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*iinv.l2rdz) / (iinv.dist);

	iinv.Pdv = -((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*iinv.l2rdx  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*iinv.l2rdy  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*iinv.l2rdz)/ (iinv.dist);

	iinv.Pdw = -((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*iinv.l2rdx  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*iinv.l2rdy  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*iinv.l2rdz )  / (iinv.dist);

	
	iinv.Ftu1 = iinv.Puf *(*ug.farea)[ug.fId]* iinv.vis;   //扩散项中归入源项的部分1
    iinv.Ftv1 = iinv.Pvf *(*ug.farea)[ug.fId]* iinv.vis;
	iinv.Ftw1 = iinv.Pwf *(*ug.farea)[ug.fId]* iinv.vis;

	iinv.Ftu2 = iinv.Pdu*(*ug.farea)[ug.fId] * iinv.vis;   //扩散项中归入源项的部分2
	iinv.Ftv2 = iinv.Pdv*(*ug.farea)[ug.fId] * iinv.vis;
	iinv.Ftw2 = iinv.Pdw*(*ug.farea)[ug.fId] * iinv.vis;



	iinv.PufT = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.zfn)[ug.fId]);

	iinv.PvfT = ((iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId]);

	iinv.PwfT = ((iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId]);



	iinv.Pud = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.xfn)[ug.fId];

    iinv.Pvd = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.yfn)[ug.fId];

	iinv.Pwd =((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.zfn)[ug.fId];

	iinv.FtuT1 = iinv.PufT * (*ug.farea)[ug.fId] * iinv.vis;  //Г(▽V)T，表面源项
	iinv.FtvT1 = iinv.PvfT * (*ug.farea)[ug.fId] * iinv.vis;
	iinv.FtwT1 = iinv.PwfT * (*ug.farea)[ug.fId] * iinv.vis;

	iinv.FtuT2 = iinv.Pud*(*ug.farea)[ug.fId]*(-2/3)*iinv.vis; //λ(∇∙V ⃗)I∙(n_n ) ⃗dS 表面源项
	iinv.FtvT2 = iinv.Pvd*(*ug.farea)[ug.fId]*(-2/3)*iinv.vis;
	iinv.FtwT2 = iinv.Pwd*(*ug.farea)[ug.fId]*(-2/3)*iinv.vis;

	iinv.ai[ug.fId][0] += iinv.Fn[ug.fId];
	iinv.ai[ug.fId][1] += iinv.Fn[ug.fId];

	iinv.spc[ug.lc] += iinv.Fn[ug.fId];
	iinv.spc[ug.rc] += iinv.Fn[ug.fId];
	
	iinv.buc[ug.lc] += iinv.Ftu1 + iinv.Ftu2;
	iinv.buc[ug.rc] += -iinv.Ftu1 - iinv.Ftu2;

	iinv.bvc[ug.lc] += iinv.Ftv1 + iinv.Ftv2;
	iinv.bvc[ug.rc] += -iinv.Ftv1 - iinv.Ftv2;

	iinv.bwc[ug.lc] += iinv.Ftw1 + iinv.Ftw2;
	iinv.bwc[ug.rc] += -iinv.Ftw1 - iinv.Ftw2;
}

void UINsVisterm::CmpBcFaceVisterm()
{
	iinv.l2rdx = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];  //界面左右单元中心距
	iinv.l2rdy = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	iinv.l2rdz = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	iinv.c2d = sqrt(iinv.l2rdx * iinv.l2rdx + iinv.l2rdy * iinv.l2rdy + iinv.l2rdz * iinv.l2rdz);

	iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);

	iinv.chalength = GetDataValue< Real >("chara_length");

	iinv.vis = 1 * iinv.chalength / inscom.reynolds;  //动力粘度

	iinv.Fn[ug.fId] = iinv.vis * (*ug.farea)[ug.fId]/ iinv.dist;

	iinv.Fbu = iinv.vis * (*ug.farea)[ug.fId] * iinv.uf[ug.fId]  / (iinv.dist);
	iinv.Fbv = iinv.vis * (*ug.farea)[ug.fId] * iinv.vf[ug.fId]  / (iinv.dist);
	iinv.Fbw = iinv.vis * (*ug.farea)[ug.fId] * iinv.wf[ug.fId] / (iinv.dist);


	/*if (ug.bctype == BC::OUTFLOW)
	{
		iinv.Puf = ((iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.dist;

		iinv.Pvf = ((iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.dist;

		iinv.Pwf = ((iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.dist;

		iinv.Pdu = -((iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] * iinv.l2rdx + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.yfn)[ug.fId] * iinv.l2rdy + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.zfn)[ug.fId] * iinv.l2rdz) / (iinv.dist*iinv.dist);

		iinv.Pdv = -((iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.xfn)[ug.fId] * iinv.l2rdx + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] * iinv.l2rdy + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.zfn)[ug.fId] * iinv.l2rdz) / (iinv.dist*iinv.dist);

		iinv.Pdw = -((iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.xfn)[ug.fId] * iinv.l2rdx + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.yfn)[ug.fId] * iinv.l2rdy + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId] * iinv.l2rdz) / (iinv.dist*iinv.dist);

	}
	else if (ug.bctype == BC::INFLOW)
	{
		iinv.Puf = ((iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.dist;

		iinv.Pvf = ((iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.dist;

		iinv.Pwf = ((iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.dist;

		iinv.Pdu = -((iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] * iinv.l2rdx + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.yfn)[ug.fId] * iinv.l2rdy + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.zfn)[ug.fId] * iinv.l2rdz) / (iinv.dist*iinv.dist);

		iinv.Pdv = -((iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.xfn)[ug.fId] * iinv.l2rdx + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] * iinv.l2rdy + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.zfn)[ug.fId] * iinv.l2rdz) / (iinv.dist*iinv.dist);

		iinv.Pdw = -((iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.xfn)[ug.fId] * iinv.l2rdx + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.yfn)[ug.fId] * iinv.l2rdy + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId] * iinv.l2rdz) / (iinv.dist*iinv.dist);

	}
	else
	{
		iinv.Puf = ((iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.dist;

		iinv.Pvf = ((iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.dist;

		iinv.Pwf = ((iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) / iinv.dist;

		iinv.Pdu = -((iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] * iinv.l2rdx + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.yfn)[ug.fId] * iinv.l2rdy + (iinv.uf[ug.fId] - (*uinsf.q)[IIDX::IIU][ug.lc])*(*ug.zfn)[ug.fId] * iinv.l2rdz) / (iinv.dist*iinv.dist);

		iinv.Pdv = -((iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.xfn)[ug.fId] * iinv.l2rdx + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] * iinv.l2rdy + (iinv.vf[ug.fId] - (*uinsf.q)[IIDX::IIV][ug.lc])*(*ug.zfn)[ug.fId] * iinv.l2rdz) / (iinv.dist*iinv.dist);
		
		iinv.Pdw = -((iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.xfn)[ug.fId] * iinv.l2rdx + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.yfn)[ug.fId] * iinv.l2rdy + (iinv.wf[ug.fId] - (*uinsf.q)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId] * iinv.l2rdz) / (iinv.dist*iinv.dist);

	
	}*/
	
	iinv.Puf = ((*uinsf.dqdx)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] + ((*uinsf.dqdy)[IIDX::IIU][ug.lc])*(*ug.yfn)[ug.fId] + ((*uinsf.dqdz)[IIDX::IIU][ug.lc])*(*ug.zfn)[ug.fId];  //▽q*n

	iinv.Pvf = ((*uinsf.dqdx)[IIDX::IIV][ug.lc])*(*ug.xfn)[ug.fId] + ((*uinsf.dqdy)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] + ((*uinsf.dqdz)[IIDX::IIV][ug.lc])*(*ug.zfn)[ug.fId];

	iinv.Pwf = ((*uinsf.dqdx)[IIDX::IIW][ug.lc])*(*ug.xfn)[ug.fId] + ((*uinsf.dqdy)[IIDX::IIW][ug.lc])*(*ug.yfn)[ug.fId] + ((*uinsf.dqdz)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId];

	iinv.Pdu = -((*uinsf.dqdx)[IIDX::IIU][ug.lc] * iinv.l2rdx + (*uinsf.dqdy)[IIDX::IIU][ug.lc] * iinv.l2rdy + (*uinsf.dqdz)[IIDX::IIU][ug.lc] * iinv.l2rdz) / iinv.dist;

	iinv.Pdv = -((*uinsf.dqdx)[IIDX::IIV][ug.lc] * iinv.l2rdx + (*uinsf.dqdy)[IIDX::IIV][ug.lc] * iinv.l2rdy + (*uinsf.dqdz)[IIDX::IIV][ug.lc] * iinv.l2rdz) / iinv.dist;

	iinv.Pdw = -((*uinsf.dqdx)[IIDX::IIW][ug.lc] * iinv.l2rdx + (*uinsf.dqdy)[IIDX::IIW][ug.lc] * iinv.l2rdy + (*uinsf.dqdz)[IIDX::IIW][ug.lc] * iinv.l2rdz) / iinv.dist;


	iinv.Ftu1 = iinv.Puf *(*ug.farea)[ug.fId] * iinv.vis;   //扩散项中归入源项的部分1
	iinv.Ftv1 = iinv.Pvf *(*ug.farea)[ug.fId] * iinv.vis;
	iinv.Ftw1 = iinv.Pwf *(*ug.farea)[ug.fId] * iinv.vis;

	iinv.Ftu2 = iinv.Pdu*(*ug.farea)[ug.fId] * iinv.vis;   //扩散项中归入源项的部分2
	iinv.Ftv2 = iinv.Pdv*(*ug.farea)[ug.fId] * iinv.vis;
	iinv.Ftw2 = iinv.Pdw*(*ug.farea)[ug.fId] * iinv.vis;


	iinv.PufT = ((*uinsf.dqdx)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] +((*uinsf.dqdx)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] +
		        ((*uinsf.dqdx)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId];

	iinv.PvfT = ((*uinsf.dqdy)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] +((*uinsf.dqdy)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] +
		        ((*uinsf.dqdy)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId];

	iinv.PwfT = ((*uinsf.dqdz)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] +((*uinsf.dqdz)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] +
		        ((*uinsf.dqdz)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId];


	iinv.Pud = ((*uinsf.dqdx)[IIDX::IIU][ug.lc]+(*uinsf.dqdy)[IIDX::IIV][ug.lc]+(*uinsf.dqdz)[IIDX::IIW][ug.lc])*(*ug.xfn)[ug.fId];

	iinv.Pvd = ((*uinsf.dqdx)[IIDX::IIU][ug.lc] +(*uinsf.dqdy)[IIDX::IIV][ug.lc]+(*uinsf.dqdz)[IIDX::IIW][ug.lc])*(*ug.yfn)[ug.fId];

	iinv.Pwd = ((*uinsf.dqdx)[IIDX::IIU][ug.lc]+(*uinsf.dqdy)[IIDX::IIV][ug.lc] +(*uinsf.dqdz)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId];

	iinv.FtuT1 = iinv.PufT * (*ug.farea)[ug.fId] * iinv.vis;  //Г(▽V)T，表面源项
	iinv.FtvT1 = iinv.PvfT * (*ug.farea)[ug.fId] * iinv.vis;
	iinv.FtwT1 = iinv.PwfT * (*ug.farea)[ug.fId] * iinv.vis;

	iinv.FtuT2 = iinv.Pud*(*ug.farea)[ug.fId] * (-2 / 3)*iinv.vis; //λ(∇∙V ⃗)I∙(n_n ) ⃗dS 表面源项
	iinv.FtvT2 = iinv.Pvd*(*ug.farea)[ug.fId] * (-2 / 3)*iinv.vis;
	iinv.FtwT2 = iinv.Pwd*(*ug.farea)[ug.fId] * (-2 / 3)*iinv.vis;

	iinv.ai[ug.fId][0] += iinv.Fn[ug.fId];
	iinv.ai[ug.fId][1] += iinv.Fn[ug.fId];

	iinv.spc[ug.lc] += iinv.Fn[ug.fId];
	
	iinv.buc[ug.lc] += iinv.Fbu +iinv.Ftu1 + iinv.Ftu2;

	iinv.bvc[ug.lc] += iinv.Fbv+ iinv.Ftv1 + iinv.Ftv2;

	iinv.bwc[ug.lc] += iinv.Fbw+ iinv.Ftw1 + iinv.Ftw2;
}


void UINsVisterm::CmpUnsteadcoff()
{
	iinv.timestep = GetDataValue< Real >("global_dt");

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.spt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]/ iinv.timestep;  //矩阵对角线元素的非稳态项
		
		if (ctrl.currTime == 0.0005 && Iteration::innerSteps == 1)
		{
			iinv.up[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];
			iinv.vp[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
			iinv.wp[ug.cId] = (*uinsf.q)[IIDX::IIW][ug.cId];

			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]* iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]* iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]* iinv.wp[ug.cId]/ iinv.timestep;
		}
		else
		{
			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.wp[ug.cId] / iinv.timestep;
		}
	}

}



void UINsVisterm::CmpINsSrc()
{
	Alloc();

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		(*uinsf.qf)[IIDX::IIP][ug.fId] = iinv.pf[ug.fId];
	}
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real dxl = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dyl = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dzl = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real dxr = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
		Real dyr = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
		Real dzr = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

		Real delt1 = DIST(dxl, dyl, dzl);
		Real delt2 = DIST(dxr, dyr, dzr);
		Real delta = 1.0 / (delt1 + delt2);

		Real cl = delt2 * delta;
		Real cr = delt1 * delta;

		(*uinsf.qf)[IIDX::IIP][ug.fId] = cl * (*uinsf.q)[IIDX::IIP][ug.lc] + cr * (*uinsf.q)[IIDX::IIP][ug.rc];

	}

	ONEFLOW::CmpINsGrad((*uinsf.qf)[IIDX::IIP], (*uinsf.dqdx)[IIDX::IIP], (*uinsf.dqdy)[IIDX::IIP], (*uinsf.dqdz)[IIDX::IIP]);

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.spc[ug.cId] += iinv.spt[ug.cId];

		iinv.buc[ug.cId] += iinv.but[ug.cId]- (*ug.cvol)[ug.cId] * (*uinsf.dqdx)[IIDX::IIP][ug.cId];
		iinv.bvc[ug.cId] += iinv.bvt[ug.cId] -(*ug.cvol)[ug.cId] * (*uinsf.dqdy)[IIDX::IIP][ug.cId];
		iinv.bwc[ug.cId] += iinv.bwt[ug.cId] - (*ug.cvol)[ug.cId] * (*uinsf.dqdz)[IIDX::IIP][ug.cId];

		int fn = (*ug.c2f)[ug.cId].size();
		iinv.dj[ug.cId] = fn;

		if (ctrl.currTime == iinv.timestep && Iteration::innerSteps == 1)
		{
			iinv.sj.resize(ug.nCell, fn);
			iinv.sd.resize(ug.nCell, fn);
		}
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			if (fId > ug.nBFace - 1)
			{
				if (ug.cId == ug.lc)
				{
					iinv.sj[ug.cId][iFace] = -iinv.ai[ug.fId][0];  //矩阵非零系数，动量方程中与主单元相邻的单元面通量
					iinv.sd[ug.cId][iFace] = ug.rc;
				}
				else if (ug.cId == ug.rc)
				{
					iinv.sj[ug.cId][iFace] = -iinv.ai[ug.fId][1];  //矩阵非零系数，动量方程中与主单元相邻的单元面通量
					iinv.sd[ug.cId][iFace] = ug.lc;
				}
			}
			else
			{
				iinv.dj[ug.cId] -= 1;
			}

		}


	}

	DeAlloc();
}



EndNameSpace

