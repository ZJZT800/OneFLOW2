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

#include "UINsMomCorrect.h"
#include "CmpGrad.h"
#include "UGrad.h"
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

UINsMomCorrect::UINsMomCorrect()
{
	;
}

UINsMomCorrect::~UINsMomCorrect()
{
	;
}

void UINsMomCorrect::UpdateCorrectSpeed()
{
	RealField dppdx, dppdy, dppdz;
	dppdx.resize(ug.nCell);
	dppdy.resize(ug.nCell);
	dppdz.resize(ug.nCell);
	 
	CmpUnsGrad(iinv.ppf, dppdx, dppdy, dppdz);

	//Update pressure
	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		int bcType = ug.bcRecord->bcType[fId];

		if (bcType == BC::OUTFLOW)
		{
			iinv.ppf[fId] = 0;
		}

		else if (bcType == BC::SOLID_SURFACE)
		{
			iinv.ppf[fId] = iinv.pp[lc];
		}

		else if (bcType == BC::INFLOW)
		{
			iinv.ppf[fId] = iinv.pp[lc];
		}

		else if (ug.bctype == BC::SYMMETRY)
		{
			iinv.ppf[fId] = iinv.pp[lc];
		}
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		iinv.ppf[fId] = (*ug.fl)[fId] * iinv.pp[lc] + (*ug.fr)[fId] * iinv.pp[rc];
	}

	Real press_relax = GetDataValue< Real >("press_relax");

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		(*uinsf.p)[0][cId] = (*uinsf.p)[0][cId] + press_relax * (iinv.pp[cId]);
	}


	for (int fId = 0; fId < ug.nBFace; fId++)
	{
		int bcType = ug.bcRecord->bcType[fId];
		int lc = (*ug.lcf)[fId];

		if (bcType == BC::SOLID_SURFACE)
		{
			iinv.pb[fId] = (*uinsf.p)[0][lc];
		}
		else if (bcType == BC::INFLOW)
		{
			iinv.pb[fId] = (*uinsf.p)[0][lc];
		}
		else if (bcType == BC::OUTFLOW)
		{
			iinv.pb[fId] += 0;
		}
		else if (ug.bctype == BC::SYMMETRY)
		{
			iinv.pb[fId] = (*uinsf.p)[0][lc];
		}
	}

	for (int fId = 0; fId < ug.nBFace; fId++)
	{
		int rc = (*ug.rcf)[fId];
		(*uinsf.p)[0][rc] = iinv.pb[fId];
	}

	// Update speed
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		(*uinsf.u)[0][cId] -= (*ug.cvol)[cId] / iinv.dup[cId] * dppdx[cId];
		(*uinsf.v)[0][cId] -= (*ug.cvol)[cId] / iinv.dup[cId] * dppdy[cId];
		(*uinsf.w)[0][cId] -= (*ug.cvol)[cId] / iinv.dup[cId] * dppdz[cId];
	}

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
			int rc = (*ug.rcf)[fId];

			if (ug.bctype == BC::INFLOW)
			{
				;
			}

			else if (ug.bctype == BC::SOLID_SURFACE)
			{
				;
			}

			else if (ug.bctype == BC::OUTFLOW)
			{
				iinv.ub[fId] = iinv.ub[fId] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdx[lc];

				iinv.vb[fId] = iinv.vb[fId] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdy[lc];

				iinv.wb[fId] = iinv.wb[fId] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdz[lc];
			}

			else if (ug.bctype == BC::SYMMETRY)
			{
				Real vnRelative = (*uinsf.u)[0][lc] * (*ug.xfn)[fId] + (*uinsf.v)[0][lc] * (*ug.yfn)[fId] + (*uinsf.w)[0][lc] * (*ug.zfn)[fId];

				iinv.ub[fId] = (*uinsf.u)[0][lc] - (*ug.xfn)[fId] * vnRelative;

				iinv.vb[fId] = (*uinsf.v)[0][lc] - (*ug.yfn)[fId] * vnRelative;

				iinv.wb[fId] = (*uinsf.u)[0][lc] - (*ug.zfn)[fId] * vnRelative;
			}

			(*uinsf.u)[0][rc] = iinv.ub[fId];
			(*uinsf.v)[0][rc] = iinv.vb[fId];
			(*uinsf.w)[0][rc] = iinv.wb[fId];

		}
	}
}

void UINsMomCorrect::UpdateCorrectFaceflux()
{
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		Real dupf;
		dupf = (*ug.fl)[fId] * ((*ug.cvol)[lc] / iinv.dup[lc]) + (*ug.fr)[fId] * ((*ug.cvol)[rc] / iinv.dup[rc]);
		Real Df1 = dupf * (*ug.a1)[fId];
		Real Df2 = dupf * (*ug.a2)[fId];
		Real Df3 = dupf * (*ug.a3)[fId];

		Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
		Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
		Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];

		Real Df = Df1 * (*ug.a1)[fId] + Df2 * (*ug.a2)[fId] + Df3 * (*ug.a3)[fId];

		Real dist = l2rdx * (*ug.a1)[fId] + l2rdy * (*ug.a2)[fId] + l2rdz * (*ug.a3)[fId];

		Real fux = (*uinsf.rho)[0][lc] * Df / dist * (iinv.pp[lc] - iinv.pp[rc]);
		iinv.flux[fId] = iinv.flux[fId] + fux;

		Real uf1 = (*ug.fl)[fId] * (*uinsf.u)[0][lc] + (*ug.fr)[fId] * (*uinsf.u)[0][rc];
		Real vf1 = (*ug.fl)[fId] * (*uinsf.v)[0][lc] + (*ug.fr)[fId] * (*uinsf.v)[0][rc];
		Real wf1 = (*ug.fl)[fId] * (*uinsf.w)[0][lc] + (*ug.fr)[fId] * (*uinsf.w)[0][rc];

		Real un = uf1 * (*ug.a1)[fId] + vf1 * (*ug.a2)[fId] + wf1 * (*ug.a3)[fId];

		iinv.dun[fId] = iinv.flux[fId] / ((*uinsf.rho)[0][lc] + SMALL) - un;
	}

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
				;
			}

			else if (ug.bctype == BC::INFLOW)
			{
				;
			}

			else if (ug.bctype == BC::OUTFLOW)
			{
				Real dupf, dvpf, dwpf;
				dupf = (*ug.cvol)[lc] / iinv.dup[lc];
				dvpf = (*ug.cvol)[lc] / iinv.dup[lc];
				dwpf = (*ug.cvol)[lc] / iinv.dup[lc];
				Real Df1 = dupf * (*ug.a1)[fId];
				Real Df2 = dvpf * (*ug.a2)[fId];
				Real Df3 = dwpf * (*ug.a3)[fId];

				Real l2rdx = (*ug.xfc)[fId] - (*ug.xcc)[lc];
				Real l2rdy = (*ug.yfc)[fId] - (*ug.ycc)[lc];
				Real l2rdz = (*ug.zfc)[fId] - (*ug.zcc)[lc];

				Real Df = Df1 * (*ug.a1)[fId] + Df2 * (*ug.a2)[fId] + Df3 * (*ug.a3)[fId];

				Real dist = l2rdx * (*ug.a1)[fId] + l2rdy * (*ug.a2)[fId] + l2rdz * (*ug.a3)[fId];

				Real fux = (*uinsf.rho)[0][lc] * Df / dist * (iinv.pp[lc] - iinv.ppf[fId]);
				iinv.flux[fId] = iinv.flux[fId] + fux;
			}

			else if (ug.bctype == BC::SYMMETRY)
			{
				;
			}
		}
	}
}

EndNameSpace