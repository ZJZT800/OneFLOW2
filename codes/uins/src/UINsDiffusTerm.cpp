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

#include "UINsDiffusTerm.h"
#include "FaceValue.h"
#include "CmpGrad.h"
#include "BcData.h"
#include "Zone.h"
#include "Atmosphere.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "UCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Multigrid.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "FieldImp.h"
#include "Iteration.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "Ctrl.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace(ONEFLOW)

UINsDiffusTerm::UINsDiffusTerm()
{
	;
}

UINsDiffusTerm::~UINsDiffusTerm()
{
	;
}


void UINsDiffusTerm::CmpDiffusTerm(string &Equa_vary, string &ischeme)
{
	RealField dudx, dudy, dudz;
	RealField dvdx, dvdy, dvdz;
	RealField dwdx, dwdy, dwdz;

	RealField uf, vf, wf, fvis_cof;

	dudx.resize(ug.nCell);
	dudy.resize(ug.nCell);
	dudz.resize(ug.nCell);
	dvdx.resize(ug.nCell);
	dvdy.resize(ug.nCell);
	dvdz.resize(ug.nCell);
	dwdx.resize(ug.nCell);
	dwdy.resize(ug.nCell);
	dwdz.resize(ug.nCell);

	uf.resize(ug.nFace);
	vf.resize(ug.nFace);
	wf.resize(ug.nFace);

	fvis_cof.resize(ug.nFace);

	FaceValue(iinv.ub, uf, (*uinsf.u)[0]);
	FaceValue(iinv.vb, vf, (*uinsf.v)[0]);
	FaceValue(iinv.wb, wf, (*uinsf.w)[0]);

	CmpUnsGrad(uf, dudx, dudy, dudz);
	CmpUnsGrad(vf, dvdx, dvdy, dvdz);
	CmpUnsGrad(wf, dwdx, dwdy, dwdz);

	Initfvis_cof(fvis_cof,Equa_vary);

	if (Equa_vary == "mom")
	{
		if (ischeme == "CENTRAL")
		{
			for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
			{
				int lc = (*ug.lcf)[fId];
				int rc = (*ug.rcf)[fId];

				Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
				Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
				Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];

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

				iinv.ai[0][fId] += fvis_cof[fId] * Fn;
				iinv.ai[1][fId] += fvis_cof[fId] * Fn;

				iinv.bu[lc] += fvis_cof[fId] * (fdudx * T1 + fdudy * T2 + fdudz * T3);
				iinv.bu[rc] -= fvis_cof[fId] * (fdudx * T1 + fdudy * T2 + fdudz * T3);

				iinv.bv[lc] += fvis_cof[fId] * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);
				iinv.bv[rc] -= fvis_cof[fId] * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);

				iinv.bw[lc] += fvis_cof[fId] * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);
				iinv.bw[rc] -= fvis_cof[fId] * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);
			}
		}
	}

	else if (Equa_vary == "energy")
	{
		if (ischeme == "CENTRAL")
		{
			;
		}
	}

}

void UINsDiffusTerm::Initfvis_cof(RealField&fvis_cof,string &Equa_vary)
{
	for (int fId = ug.nBFace; fId < ug.nFace; fId++)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		if (Equa_vary == "mom")
		{
			fvis_cof[fId] = (*ug.fl)[fId] * (*uinsf.vis_coef)[0][lc] + (*ug.fr)[fId] * (*uinsf.vis_coef)[0][rc];
		}
		else if(Equa_vary == "energy")
		{
			;
		}
	}

	ug.nRegion = ug.bcRecord->bcInfo->bcType.size();
	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
		{
			BcInfo* bcInfo = ug.bcRecord->bcInfo;
			int fId = bcInfo->bcFace[ug.ir][ibc];
			int lc = (*ug.lcf)[fId];

			if (ug.bctype == BC::SOLID_SURFACE)
			{
				if (Equa_vary == "mom")
				{
					iinv.fvisb_cof[fId] = (*uinsf.vis_coef)[0][lc];
					fvis_cof[fId] = iinv.fvisb_cof[fId];
				}
				else if (Equa_vary == "energy")
				{
					;
				}
			}
			else if (ug.bctype == BC::INFLOW)
			{
				if (Equa_vary == "mom")
				{
					iinv.fvisb_cof[fId] = (*uinsf.vis_coef)[0][lc];
					fvis_cof[fId] = iinv.fvisb_cof[fId];
				}
				else if (Equa_vary == "energy")
				{
					;
				}
			}
			else if (ug.bctype == BC::OUTFLOW)
			{
				if (Equa_vary == "mom")
				{
					iinv.fvisb_cof[fId] = (*uinsf.vis_coef)[0][lc];
					fvis_cof[fId] = iinv.fvisb_cof[fId];
				}
				else if (Equa_vary == "energy")
				{
					;
				}
			}
		}
	}
}

EndNameSpace