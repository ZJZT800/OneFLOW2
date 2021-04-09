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

#include "UINsCmpFlux.h"
#include "FaceValue.h"
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

UINsCmpFlux::UINsCmpFlux()
{
	;
}

UINsCmpFlux::~UINsCmpFlux()
{
	;
}


void UINsCmpFlux::CmpFlux()
{
	RealField dpdx, dpdy, dpdz;
	RealField pb, pf;

	dpdx.resize(ug.nCell);
	dpdy.resize(ug.nCell);
	dpdz.resize(ug.nCell);

	pf.resize(ug.nFace);

	FaceValue(iinv.pb, pf, (*uinsf.p)[0]);
	CmpUnsGrad(pf, dpdx, dpdy, dpdz);

	Real mom_relax = GetDataValue< Real >("mom_relax");
	Real rurf = mom_relax / (1 + mom_relax);

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		Real dpdx1 = dpdx[lc];
		Real dpdx2 = dpdx[rc];
		Real dpdy1 = dpdy[lc];
		Real dpdy2 = dpdy[rc];
		Real dpdz1 = dpdz[lc];
		Real dpdz2 = dpdz[rc];

		Real rl, ul, vl, wl, pl;
		Real rr, ur, vr, wr, pr;

		Real VdU1, VdU2, Vdvu;

		Real uf, vf, wf;

		Real vnflow;

		rl = (*uinsf.rho)[0][lc];
		ul = (*uinsf.u)[0][lc];
		vl = (*uinsf.v)[0][lc];
		wl = (*uinsf.w)[0][lc];
		pl = (*uinsf.p)[0][lc];

		ur = (*uinsf.u)[0][rc];
		vr = (*uinsf.v)[0][rc];
		wr = (*uinsf.w)[0][rc];
		pr = (*uinsf.p)[0][rc];

		Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
		Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
		Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];

		VdU1 = (*ug.cvol)[lc] / iinv.spu[lc];
		VdU2 = (*ug.cvol)[rc] / iinv.spu[rc];

		Vdvu = (*ug.fl)[fId] * VdU1 + (*ug.fr)[fId] * VdU2;

		Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

		Real Df1 = Vdvu * (*ug.a1)[fId] / dist;
		Real Df2 = Vdvu * (*ug.a2)[fId] / dist;
		Real Df3 = Vdvu * (*ug.a3)[fId] / dist;

		Real dx1 = (*ug.xfc)[fId] - (*ug.xcc)[lc];
		Real dy1 = (*ug.yfc)[fId] - (*ug.ycc)[lc];
		Real dz1 = (*ug.zfc)[fId] - (*ug.zcc)[lc];

		Real dx2 = (*ug.xcc)[rc] - (*ug.xfc)[fId];
		Real dy2 = (*ug.ycc)[rc] - (*ug.yfc)[fId];
		Real dz2 = (*ug.zcc)[rc] - (*ug.zfc)[fId];

		Real fdpdx = dpdx1 * dx1 + dpdx2 * dx2 - (pr - pl);
		Real fdpdy = dpdy1 * dy1 + dpdy2 * dy2 - (pr - pl);
		Real fdpdz = dpdz1 * dz1 + dpdz2 * dz2 - (pr - pl);

		uf = ul * (*ug.fl)[fId] + ur * (*ug.fr)[fId];
		vf = vl * (*ug.fl)[fId] + vr * (*ug.fr)[fId];
		wf = wl * (*ug.fl)[fId] + wr * (*ug.fr)[fId];

		vnflow = (*ug.a1)[fId] * (uf + fdpdx * Df1) + (*ug.a2)[fId] * (vf + fdpdy * Df2) + (*ug.a3)[fId] * (wf + fdpdz * Df3) + rurf * iinv.dun[fId];
		iinv.flux[fId] = rl * vnflow;
	}

	ug.nRegion = ug.bcRecord->bcInfo->bcType.size();
	BcInfo* bcInfo = ug.bcRecord->bcInfo;

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

			Real dpdx1 = dpdx[lc];
			Real dpdy1 = dpdy[lc];
			Real dpdz1 = dpdz[lc];

			Real pb1 = iinv.pb[fId];

			Real rl, ul, vl, wl, pl;

			Real VdU1, Vdvu;

			Real uf, vf, wf;

			Real vnflow;

			if (ug.bctype == BC::SOLID_SURFACE)
			{
				rl = (*uinsf.rho)[0][lc];

				vnflow = (*ug.a1)[fId] * iinv.ub[fId] + (*ug.a2)[fId] * iinv.vb[fId] + (*ug.a3)[fId] * iinv.wb[fId];

				iinv.flux[fId] = rl * vnflow;
			}

			else if (ug.bctype == BC::INFLOW)
			{
				vnflow = (*ug.a1)[fId] * iinv.ub[fId] + (*ug.a2)[fId] * iinv.vb[fId] + (*ug.a3)[fId] * iinv.wb[fId];

				iinv.flux[fId] = inscom.inflow[IIDX::IIR] * vnflow;
			}

			else if (ug.bctype == BC::OUTFLOW)
			{
				rl = (*uinsf.rho)[0][lc];
				ul = (*uinsf.u)[0][lc];
				vl = (*uinsf.v)[0][lc];
				wl = (*uinsf.w)[0][lc];
				pl = (*uinsf.p)[0][lc];

				Real l2rdx = (*ug.xfc)[fId] - (*ug.xcc)[ug.lc];
				Real l2rdy = (*ug.yfc)[fId] - (*ug.ycc)[ug.lc];
				Real l2rdz = (*ug.zfc)[fId] - (*ug.zcc)[ug.lc];

				VdU1 = (*ug.cvol)[lc] / iinv.spu[lc];

				Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

				Real Df1 = VdU1 * (*ug.a1)[fId] / dist;
				Real Df2 = VdU1 * (*ug.a2)[fId] / dist;
				Real Df3 = VdU1 * (*ug.a3)[fId] / dist;

				Real dx1 = (*ug.xfc)[fId] - (*ug.xcc)[lc];
				Real dy1 = (*ug.yfc)[fId] - (*ug.ycc)[lc];
				Real dz1 = (*ug.zfc)[fId] - (*ug.zcc)[lc];

				Real fdpdx = dpdx1 * dx1 - (pb1 - pl);
				Real fdpdy = dpdy1 * dy1 - (pb1 - pl);
				Real fdpdz = dpdz1 * dz1 - (pb1 - pl);

				iinv.ub[fId] = ul + fdpdx * Df1;
				iinv.vb[fId] = vl + fdpdy * Df2;
				iinv.wb[fId] = wl + fdpdz * Df3;

				vnflow = (*ug.a1)[fId] * iinv.ub[fId] + (*ug.a2)[fId] * iinv.vb[fId] + (*ug.a3)[fId] * iinv.wb[fId] + rurf * iinv.dun[fId];

				iinv.flux[fId] = rl * vnflow;
			}

			else if (ug.bctype == BC::SYMMETRY)
			{
				iinv.flux[fId] = 0;
			}
		}
	}

	/*RealField massflux = 0;
massflux.resize(ug.nCell);
for (int cId = 0; cId < ug.nCell; cId++)
{
	int fn = (*ug.c2f)[cId].size();
	for (int iFace = 0; iFace < fn; iFace++)
	{
		int fId = (*ug.c2f)[cId][iFace];
		massflux[cId] += iinv.fq[fId];
	}
	std::cout << "cId: " << cId << ", massflux[cId]: " << massflux[cId] << std::endl;
}*/
}

EndNameSpace