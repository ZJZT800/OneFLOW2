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

#include "UINsInitField.h"
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

UINsInitField::UINsInitField()
{
	;
}

UINsInitField::~UINsInitField()
{
	;
}


void UINsInitField::InitInputField()
{
		iinv.FieldInit();
		ug.Init();
		uinsf.Init();
		int transt = ONEFLOW::GetDataValue< int >("transt");
		if(transt != 0) iinv.OldValueInit();

		Real  rl, ul, vl, wl, pl;
		Real  rr, ur, vr, wr, pw;

		int solve_energy = GetDataValue< int >("solve_energy");

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
				ug.bcr = bcInfo->bcRegion[ug.ir][ibc];
				ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];
				int lc = (*ug.lcf)[fId];

				if (ug.bcr == -1) return; //interface
				int dd = bcdata.r2d[ug.bcr];
				if (dd != -1)
				{
					ug.bcdtkey = 1;
					inscom.bcflow = &bcdata.dataList[dd];
				}

				if (ug.bctype == BC::SOLID_SURFACE)
				{

					if (ug.bcdtkey == 0)
					{

						iinv.ub[fId] = (*ug.vfx)[fId];

						iinv.vb[fId] = (*ug.vfy)[fId];

						iinv.wb[fId] = (*ug.vfz)[fId];

						iinv.pb[fId] = (*uinsf.p)[0][lc];

					}
					else
					{
						iinv.ub[fId] = (*inscom.bcflow)[1];

						iinv.vb[fId] = (*inscom.bcflow)[2];

						iinv.wb[fId] = (*inscom.bcflow)[3];

						iinv.pb[fId] = (*uinsf.p)[0][lc];
					}

					if (solve_energy == 1)
					{
						;
					}
					iinv.flux[fId] = 0;
					iinv.fvisb_cof[fId] = (*uinsf.vis_coef)[0][lc];
				}

				else if (ug.bctype == BC::OUTFLOW)
				{
					iinv.ub[fId] = (*uinsf.u)[0][lc];
					iinv.vb[fId] = (*uinsf.v)[0][lc];
					iinv.wb[fId] = (*uinsf.w)[0][lc];
					iinv.pb[fId] = (*uinsf.p)[0][lc];

					if (solve_energy == 1)
					{
						;
					}
					iinv.flux[fId] = (*uinsf.rho)[0][lc] * ((*ug.a1)[fId] * iinv.ub[fId] + (*ug.a2)[fId] * iinv.vb[fId] + (*ug.a3)[fId] * iinv.wb[fId]);
					iinv.fvisb_cof[fId] = (*uinsf.vis_coef)[0][lc];
				}

				else if (ug.bctype == BC::INFLOW)
				{
					iinv.ub[fId] = inscom.inflow[IIDX::IIU];

					iinv.vb[fId] = inscom.inflow[IIDX::IIV];

					iinv.wb[fId] = inscom.inflow[IIDX::IIW];

					iinv.pb[fId] = inscom.inflow[IIDX::IIP];

					if (solve_energy == 1)
					{
						;
					}
					iinv.flux[fId] = inscom.inflow[0] * ((*ug.a1)[fId] * iinv.ub[fId] + (*ug.a2)[fId] * iinv.vb[fId] + (*ug.a3)[fId] * iinv.wb[fId]);
					iinv.fvisb_cof[fId] = (*uinsf.vis_coef)[0][lc];
				}

				else if (ug.bctype == BC::SYMMETRY)
				{
					;
				}

			}
		}

		Real uf, vf, wf;

		for (int fId = ug.nBFace; fId < ug.nFace; fId++)
		{
			int lc = (*ug.lcf)[fId];
			int rc = (*ug.rcf)[fId];

			rl = (*uinsf.rho)[0][lc];
			ul = (*uinsf.u)[0][lc];
			vl = (*uinsf.v)[0][lc];
			wl = (*uinsf.w)[0][lc];
			ur = (*uinsf.u)[0][rc];
			vr = (*uinsf.v)[0][rc];
			wr = (*uinsf.w)[0][rc];

			uf = ul * (*ug.fl)[fId] + ur * (*ug.fr)[fId];

			vf = vl * (*ug.fl)[fId] + vr * (*ug.fr)[fId];

			wf = wl * (*ug.fl)[fId] + wr * (*ug.fr)[fId];

			iinv.flux[fId] = rl * ((*ug.a1)[fId] * uf + (*ug.a2)[fId] * vf + (*ug.a3)[fId] * wf);
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