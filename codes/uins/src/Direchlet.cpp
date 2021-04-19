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

#include "Direchlet.h"
//#include "INsInvterm.h"
//#include "UINsVisterm.h"
//#include "UINsGrad.h"
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

void DirechletBc(string&EquaVary,RealField& dqdx, RealField& dqdy, RealField& dqdz, Real& qb, RealField &fdiffus_cof, string &diffus_ischeme,int& fId, RealField &spu,RealField &bu)
{
	if (EquaVary == "vary_h")
	{
	;
	}
	else
	{
		if (diffus_ischeme == "CENTRAL")
		{
			int lc = (*ug.lcf)[fId];

			Real l2rdx = (*ug.xfc)[fId] - (*ug.xcc)[lc];
			Real l2rdy = (*ug.yfc)[fId] - (*ug.ycc)[lc];
			Real l2rdz = (*ug.zfc)[fId] - (*ug.zcc)[lc];

			Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

			Real Fn = (*ug.a1)[fId] * (*ug.a1)[fId] + (*ug.a2)[fId] * (*ug.a2)[fId] + (*ug.a3)[fId] * (*ug.a3)[fId];

			Fn = Fn / dist;

			Real T1 = (*ug.a1)[fId] - l2rdx * Fn;
			Real T2 = (*ug.a2)[fId] - l2rdy * Fn;
			Real T3 = (*ug.a3)[fId] - l2rdz * Fn;

			Real fdqdx = dqdx[lc];
			Real fdqdy = dqdy[lc];
			Real fdqdz = dqdz[lc];

			spu[lc] += fdiffus_cof[fId] * Fn;
			bu[lc] += fdiffus_cof[fId] * Fn * qb + fdiffus_cof[fId] * (fdqdx * T1 + fdqdy * T2 + fdqdz * T3);
		}
	}
}

EndNameSpace