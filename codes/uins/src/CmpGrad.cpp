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

#include "CmpGrad.h"
//#include "INsInvterm.h"
//#include "UINsVisterm.h"
//#include "UINsGrad.h"
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

void CmpUnsGrad(RealField & q, RealField & dqdx, RealField & dqdy, RealField & dqdz)
{
	dqdx = 0;
	dqdy = 0;
	dqdz = 0;

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		dqdx[lc] += (*ug.a1)[fId] * q[fId];
		dqdy[lc] += (*ug.a2)[fId] * q[fId];
		dqdz[lc] += (*ug.a3)[fId] * q[fId];
	}
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		dqdx[lc] += (*ug.a1)[fId] * q[fId];
		dqdy[lc] += (*ug.a2)[fId] * q[fId];
		dqdz[lc] += (*ug.a3)[fId] * q[fId];

		dqdx[rc] -= (*ug.a1)[fId] * q[fId];
		dqdy[rc] -= (*ug.a2)[fId] * q[fId];
		dqdz[rc] -= (*ug.a3)[fId] * q[fId];
	}
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		Real vol = (*ug.cvol)[cId];
		dqdx[cId] /= vol;
		dqdy[cId] /= vol;
		dqdz[cId] /= vol;
	}
}

EndNameSpace