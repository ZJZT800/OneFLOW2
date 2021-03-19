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
#include "ZoneState.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "INsCtrl.h"
#include "UCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "VisGrad.h"
//#include "UINsGrad.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "BcData.h"
#include "ULimiter.h"
#include "UINsLimiter.h"
#include "FieldImp.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include <iostream>
using namespace std;

BeginNameSpace(ONEFLOW)


UINsVisterm::UINsVisterm()
{
	;
}

UINsVisterm::~UINsVisterm()
{
	;
}


void UINsVisterm::CmpDiffus()
{
	this->CmpDiffusTerm();
}

void UINsVisterm::Alloc()
{
	uinsf.qf = new MRField(inscom.nEqu, ug.nFace);
}

void UINsVisterm::DeAlloc()
{
	delete uinsf.qf;
}

/*void UINsVisterm::PrepareField()
{

	ONEFLOW::CmpINsGrad((*uinsf.q)[IIDX::IIU], (*uinsf.dqdx)[IIDX::IIU], (*uinsf.dqdy)[IIDX::IIU], (*uinsf.dqdz)[IIDX::IIU]);
	ONEFLOW::CmpINsGrad((*uinsf.q)[IIDX::IIV], (*uinsf.dqdx)[IIDX::IIV], (*uinsf.dqdy)[IIDX::IIV], (*uinsf.dqdz)[IIDX::IIV]);
	ONEFLOW::CmpINsGrad((*uinsf.q)[IIDX::IIW], (*uinsf.dqdx)[IIDX::IIW], (*uinsf.dqdy)[IIDX::IIW], (*uinsf.dqdz)[IIDX::IIW]);

}*/

/*void UINsVisterm::CmpPreandVisGrad()
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

}*/


void UINsVisterm::CmpDiffusTerm()
{
	RealField dudx, dudy, dudz;
	RealField dvdx, dvdy, dvdz;
	RealField dwdx, dwdy, dwdz;

	RealField ub, vb, wb;
	RealField uf, vf, wf;

	dudx.resize(ug.nCell);
	dudy.resize(ug.nCell);
	dudz.resize(ug.nCell);
	dvdx.resize(ug.nCell);
	dvdy.resize(ug.nCell);
	dvdz.resize(ug.nCell);
	dwdx.resize(ug.nCell);
	dwdy.resize(ug.nCell);
	dwdz.resize(ug.nCell);

	ub.resize(ug.nBFace);
	vb.resize(ug.nBFace);
	wb.resize(ug.nBFace);

	uf.resize(ug.nFace);
	vf.resize(ug.nFace);
	wf.resize(ug.nFace);

	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->BcVelocity(ub, vb, wb);
	uINsInvterm->FaceVelocity(ub, vb, wb, uf, vf, wf);


	ONEFLOW::CmpINsGrad(uf, dudx, dudy, dudz);
	ONEFLOW::CmpINsGrad(vf, dvdx, dvdy, dvdz);
	ONEFLOW::CmpINsGrad(wf, dwdx, dwdy, dwdz);

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		this->InDiffusCoff(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz,fId);
	}

	//Direchlet Boundary Condition
	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		Real ub1 = ub[fId];
		Real vb1 = vb[fId];
		Real wb1 = wb[fId];

		this->BcDiffusCoff(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, ub1,vb1,wb1,fId);
	}

}

void UINsVisterm::InDiffusCoff(RealField& dudx, RealField& dudy, RealField& dudz, RealField& dvdx, RealField& dvdy, RealField& dvdz, RealField& dwdx, RealField& dwdy, RealField& dwdz, int& fId)
{
	int lc = (*ug.lcf)[fId];
	int rc = (*ug.rcf)[fId];

	Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];  
	Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
	Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];
	 
	Real vis = GetDataValue< Real >("k_vis");

	Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

	Real Fn = (*ug.a1)[fId] * (*ug.a1)[fId] + (*ug.a2)[fId] * (*ug.a2)[fId] + (*ug.a3)[fId] * (*ug.a3)[fId];

	Fn = Fn / dist;

	Real T1 = (*ug.a1)[fId] - l2rdx * Fn;
	Real T2 = (*ug.a2)[fId] - l2rdy * Fn;
	Real T3 = (*ug.a3)[fId] - l2rdz * Fn;

	Real fdudx = (*ug.fl)[fId] * dudx[lc] + (*ug.fr)[fId] * dudx[rc];
	Real fdudy = (*ug.fl)[fId] * dudy[lc] + (*ug.fr)[fId] * dudy[rc];
	Real fdudz = (*ug.fl)[fId] * dudz[lc] + (*ug.fr)[fId] * dudz[rc];
	Real fdvdx = (*ug.fl)[fId] * dvdx[lc] + (*ug.fr)[fId] * dvdx[rc];
	Real fdvdy = (*ug.fl)[fId] * dvdy[lc] + (*ug.fr)[fId] * dvdy[rc];
	Real fdvdz = (*ug.fl)[fId] * dvdz[lc] + (*ug.fr)[fId] * dvdz[rc];
	Real fdwdx = (*ug.fl)[fId] * dwdx[lc] + (*ug.fr)[fId] * dwdx[rc];
	Real fdwdy = (*ug.fl)[fId] * dwdy[lc] + (*ug.fr)[fId] * dwdy[rc];
	Real fdwdz = (*ug.fl)[fId] * dwdz[lc] + (*ug.fr)[fId] * dwdz[rc];

	iinv.ai[0][fId] += vis * Fn;
	iinv.ai[1][fId] += vis * Fn;

	iinv.bu[lc] += vis * (fdudx * T1 + fdudy * T2 + fdudz * T3);
	iinv.bu[rc] -= vis * (fdudx * T1 + fdudy * T2 + fdudz * T3);

	iinv.bv[lc] += vis * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);
	iinv.bv[rc] -= vis * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);

	iinv.bw[lc] += vis * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);
	iinv.bw[rc] -= vis * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);
}

void UINsVisterm::BcDiffusCoff(RealField& dudx, RealField& dudy, RealField& dudz, RealField& dvdx, RealField& dvdy, RealField& dvdz, RealField& dwdx, RealField& dwdy, RealField& dwdz,Real& ub1, Real& vb1, Real& wb1, int& fId)
{
	int lc = (*ug.lcf)[fId];

	Real l2rdx = (*ug.xfc)[fId] - (*ug.xcc)[lc];  
	Real l2rdy = (*ug.yfc)[fId] - (*ug.ycc)[lc];
	Real l2rdz = (*ug.zfc)[fId] - (*ug.zcc)[lc];

	Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

	Real Fn = (*ug.a1)[fId] * (*ug.a1)[fId] + (*ug.a2)[fId] * (*ug.a2)[fId] + (*ug.a3)[fId] * (*ug.a3)[fId];

	Fn = Fn / dist;

	Real vis = GetDataValue< Real >("k_vis");

	Real T1 = (*ug.a1)[fId] - l2rdx * Fn;
	Real T2 = (*ug.a2)[fId] - l2rdy * Fn;
	Real T3 = (*ug.a3)[fId] - l2rdz * Fn;

	Real fdudx = dudx[lc];
	Real fdudy = dudy[lc];
	Real fdudz = dudz[lc];
	Real fdvdx = dvdx[lc];
	Real fdvdy = dvdy[lc];
	Real fdvdz = dvdz[lc];
	Real fdwdx = dwdx[lc];
	Real fdwdy = dwdy[lc];
	Real fdwdz = dwdz[lc];

	iinv.spu[lc] += vis * Fn;

	iinv.bu[lc] += vis * Fn * ub1 + vis * (fdudx * T1 + fdudy * T2 + fdudz * T3);
	iinv.bv[lc] += vis * Fn * vb1 + vis * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);
	iinv.bw[lc] += vis * Fn * wb1 + vis * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);

}


void UINsVisterm::CmpTranst()
{
	 Real timestep = GetDataValue< Real >("global_dt");

	for (int cId = 0; cId < ug.nCell; ++cId)
	{

		iinv.spu[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][cId] / timestep;  //矩阵对角线元素的非稳态项

		iinv.bu[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * (*uinsf.q)[IIDX::IIU][cId] / timestep; //源项的非稳态项
		iinv.bv[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * (*uinsf.q)[IIDX::IIV][cId] / timestep;
		iinv.bw[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * (*uinsf.q)[IIDX::IIW][cId] / timestep;
	}

}


void UINsVisterm::CmpSrc()
{
	RealField dpdx, dpdy, dpdz;
	RealField pb, pf;

	dpdx.resize(ug.nCell);
	dpdy.resize(ug.nCell);
	dpdz.resize(ug.nCell);

	pb.resize(ug.nBFace);
	pf.resize(ug.nFace);

	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->BcPressure(pb);
	uINsInvterm->FacePressure(pb, pf);

	ONEFLOW::CmpINsGrad(pf, dpdx, dpdy, dpdz);

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		Real vol = (*ug.cvol)[cId];

		iinv.bu[cId] += -vol * dpdx[cId];
		iinv.bv[cId] += -vol * dpdy[cId];
		iinv.bw[cId] += -vol * dpdz[cId];
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		iinv.spu[lc] += iinv.ai[0][fId];
		iinv.spu[rc] += iinv.ai[1][fId];
	}
}

void UINsVisterm::MomEquCoeff()
{
	iinv.remax_up = 0;
	iinv.remax_vp = 0;
	iinv.remax_wp = 0;

	for (int fId = ug.nBFace; fId < ug.nFace; fId++)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		iinv.bu[lc] += iinv.ai[0][fId] * (*uinsf.q)[IIDX::IIU][rc];
		iinv.bv[lc] += iinv.ai[0][fId] * (*uinsf.q)[IIDX::IIV][rc];
		iinv.bw[lc] += iinv.ai[0][fId] * (*uinsf.q)[IIDX::IIW][rc];

		iinv.bu[rc] += iinv.ai[1][fId] * (*uinsf.q)[IIDX::IIU][lc];
		iinv.bv[rc] += iinv.ai[1][fId] * (*uinsf.q)[IIDX::IIV][lc];
		iinv.bw[rc] += iinv.ai[1][fId] * (*uinsf.q)[IIDX::IIW][lc];
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		iinv.bu[cId] -= iinv.spu[cId] * (*uinsf.q)[IIDX::IIU][cId];
		iinv.bv[cId] -= iinv.spu[cId] * (*uinsf.q)[IIDX::IIV][cId];
		iinv.bw[cId] -= iinv.spu[cId] * (*uinsf.q)[IIDX::IIW][cId];

		/*iinv.remax_up = MAX(abs(iinv.remax_up), abs(iinv.buc[cId]));
		iinv.remax_vp = MAX(abs(iinv.remax_vp), abs(iinv.bvc[cId]));
		iinv.remax_wp = MAX(abs(iinv.remax_wp), abs(iinv.bwc[cId]));*/

		iinv.remax_up += abs(iinv.bu[cId]);
		iinv.remax_vp += abs(iinv.bv[cId]);
		iinv.remax_wp += abs(iinv.bw[cId]);
	}
}

void UINsVisterm::RelaxMom(Real a)
{
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.spu[cId] = iinv.spu[cId] * (1 + a);
	}

}

EndNameSpace

