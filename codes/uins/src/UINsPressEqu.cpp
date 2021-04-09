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

#include "UINsPressEqu.h"
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

UINsPressEqu::UINsPressEqu()
{
	;
}

UINsPressEqu::~UINsPressEqu()
{
	;
}


void UINsPressEqu::PressEqu()
{
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.dup[cId] = iinv.spu[cId];
	}
	for (int fId = ug.nBFace; fId < ug.nFace; fId++)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];
		iinv.dup[lc] = iinv.dup[lc] - iinv.ai[0][fId];
		iinv.dup[rc] = iinv.dup[rc] - iinv.ai[1][fId];
	}

	iinv.ai[0] = 0;
	iinv.ai[1] = 0;

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		Real duf = (*ug.fl)[fId] * ((*ug.cvol)[lc] / iinv.dup[lc]) + (*ug.fr)[fId] * ((*ug.cvol)[rc] / iinv.dup[rc]);
		Real Sf1 = duf * (*ug.a1)[fId];
		Real Sf2 = duf * (*ug.a2)[fId];
		Real Sf3 = duf * (*ug.a3)[fId];

		Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
		Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
		Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];

		Real dist = l2rdx * (*ug.a1)[fId] + l2rdy * (*ug.a2)[fId] + l2rdz * (*ug.a3)[fId];

		Real Sfarea = Sf1 * (*ug.a1)[fId] + Sf2 * (*ug.a2)[fId] + Sf3 * (*ug.a3)[fId];

		//iinv.rf = (*ug.fl)[ug.fId] * (*uinsf.q)[IIDX::IIR][ug.lc] + ((*ug.fr)[ug.fId]) * (*uinsf.q)[IIDX::IIR][ug.rc];

		iinv.spp[lc] += (*uinsf.rho)[0][lc] * Sfarea / dist;
		iinv.spp[rc] += (*uinsf.rho)[0][lc] * Sfarea / dist;
		iinv.ai[0][fId] = (*uinsf.rho)[0][lc] * Sfarea / dist;
		iinv.ai[1][fId] = (*uinsf.rho)[0][lc] * Sfarea / dist;

		iinv.bp[lc] -= iinv.flux[fId];
		iinv.bp[rc] += iinv.flux[fId];
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];

		Real duf = (*ug.cvol)[lc] / iinv.dup[lc];
		Real Sf1 = duf * (*ug.a1)[fId];
		Real Sf2 = duf * (*ug.a2)[fId];
		Real Sf3 = duf * (*ug.a3)[fId];

		Real l2rdx = (*ug.xfc)[fId] - (*ug.xcc)[lc];
		Real l2rdy = (*ug.yfc)[fId] - (*ug.ycc)[lc];
		Real l2rdz = (*ug.zfc)[fId] - (*ug.zcc)[lc];

		Real dist = l2rdx * (*ug.a1)[fId] + l2rdy * (*ug.a2)[fId] + l2rdz * (*ug.a3)[fId];

		Real Sfarea = Sf1 * (*ug.a1)[fId] + Sf2 * (*ug.a2)[fId] + Sf3 * (*ug.a3)[fId];

		int bcType = ug.bcRecord->bcType[fId];

		if (bcType == BC::OUTFLOW)
		{
			iinv.spp[lc] += (*uinsf.rho)[0][lc] * Sfarea / dist;
		}

		else if (bcType == BC::SOLID_SURFACE)
		{
			;
		}

		else if (bcType == BC::INFLOW)
		{
			;
		}

		else if (ug.bctype == BC::SYMMETRY)
		{
			;
		}

		iinv.bp[lc] -= iinv.flux[fId];
	}

	iinv.remax_pp = 0;
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.remax_pp += pow(iinv.bp[cId], 2); 
	}

	    iinv.remax_pp = sqrt(iinv.remax_pp);
}


EndNameSpace