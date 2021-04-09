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

#include "UINsPrimeEqu.h"
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

UINsPrimeEqu::UINsPrimeEqu()
{
	;
}

UINsPrimeEqu::~UINsPrimeEqu()
{
	;
}


void UINsPrimeEqu::PrimeEqu(string &Equa_vary)
{
	if (Equa_vary == "mom")
	{
		iinv.remax_u = 0;
		iinv.remax_v = 0;
		iinv.remax_w = 0;

		for (int fId = ug.nBFace; fId < ug.nFace; fId++)
		{
			int lc = (*ug.lcf)[fId];
			int rc = (*ug.rcf)[fId];

			iinv.spu[lc] += iinv.ai[0][fId];
			iinv.spu[rc] += iinv.ai[1][fId];

			iinv.bu[lc] += iinv.ai[0][fId] * (*uinsf.u)[0][rc];
			iinv.bv[lc] += iinv.ai[0][fId] * (*uinsf.v)[0][rc];
			iinv.bw[lc] += iinv.ai[0][fId] * (*uinsf.w)[0][rc];

			iinv.bu[rc] += iinv.ai[1][fId] * (*uinsf.u)[0][lc];
			iinv.bv[rc] += iinv.ai[1][fId] * (*uinsf.v)[0][lc];
			iinv.bw[rc] += iinv.ai[1][fId] * (*uinsf.w)[0][lc];
		}

		for (int cId = 0; cId < ug.nCell; ++cId)
		{
			iinv.bu[cId] -= iinv.spu[cId] * (*uinsf.u)[0][cId];
			iinv.bv[cId] -= iinv.spu[cId] * (*uinsf.v)[0][cId];
			iinv.bw[cId] -= iinv.spu[cId] * (*uinsf.w)[0][cId];

			iinv.remax_u += pow(iinv.bu[cId], 2);
			iinv.remax_v += pow(iinv.bv[cId], 2);
			iinv.remax_w += pow(iinv.bw[cId], 2);
		}

		iinv.remax_u = sqrt(iinv.remax_u);
		iinv.remax_v = sqrt(iinv.remax_v);
		iinv.remax_w = sqrt(iinv.remax_w);
	}
	else if (Equa_vary == "energy")
	{
		;
	}
}

EndNameSpace