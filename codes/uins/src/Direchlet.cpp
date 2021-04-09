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

void DirechletMom(RealField& dudx, RealField& dudy, RealField& dudz, RealField& dvdx, RealField& dvdy, RealField& dvdz, RealField& dwdx, RealField& dwdy, RealField& dwdz, Real& ub1, Real& vb1, Real& wb1, int& fId)
{
	int lc = (*ug.lcf)[fId];

	Real l2rdx = (*ug.xfc)[fId] - (*ug.xcc)[lc];
	Real l2rdy = (*ug.yfc)[fId] - (*ug.ycc)[lc];
	Real l2rdz = (*ug.zfc)[fId] - (*ug.zcc)[lc];

	Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

	Real Fn = (*ug.a1)[fId] * (*ug.a1)[fId] + (*ug.a2)[fId] * (*ug.a2)[fId] + (*ug.a3)[fId] * (*ug.a3)[fId];

	Fn = Fn / dist;

	//Real vis_cof = GetDataValue< Real >("vis_coef");

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

	iinv.spu[lc] += iinv.fvisb_cof[fId] * Fn;

	iinv.bu[lc] += iinv.fvisb_cof[fId] * Fn * ub1 + iinv.fvisb_cof[fId] * (fdudx * T1 + fdudy * T2 + fdudz * T3);
	iinv.bv[lc] += iinv.fvisb_cof[fId] * Fn * vb1 + iinv.fvisb_cof[fId] * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);
	iinv.bw[lc] += iinv.fvisb_cof[fId] * Fn * wb1 + iinv.fvisb_cof[fId] * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);
}

EndNameSpace