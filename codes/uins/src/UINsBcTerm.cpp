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

#include "UINsBcTerm.h"
#include "FaceValue.h"
#include "CmpGrad.h"
#include "Direchlet.h"
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

UINsBcTerm::UINsBcTerm()
{
	;
}

UINsBcTerm::~UINsBcTerm()
{
	;
}


void UINsBcTerm::CmpMomBcTerm(string &Equa_vary)
{
	if (Equa_vary == "mom")
	{
		RealField dudx, dudy, dudz;
		RealField dvdx, dvdy, dvdz;
		RealField dwdx, dwdy, dwdz;

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

		uf.resize(ug.nFace);
		vf.resize(ug.nFace);
		wf.resize(ug.nFace);

		FaceValue(iinv.ub, uf, (*uinsf.u)[0]);
		FaceValue(iinv.vb, vf, (*uinsf.v)[0]);
		FaceValue(iinv.wb, wf, (*uinsf.w)[0]);

		CmpUnsGrad(uf, dudx, dudy, dudz);
		CmpUnsGrad(vf, dvdx, dvdy, dvdz);
		CmpUnsGrad(wf, dwdx, dwdy, dwdz);

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

				Real clr = MAX(0, iinv.flux[fId]);
				Real crl = clr - iinv.flux[fId];

				if (ug.bctype == BC::SOLID_SURFACE)
				{
					Real ub1 = iinv.ub[fId];
					Real vb1 = iinv.vb[fId];
					Real wb1 = iinv.wb[fId];

					DirechletMom(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, ub1, vb1, wb1, fId);
				}
				else if (ug.bctype == BC::INFLOW)
				{
					if (iinv.flux[fId] < 0)
					{
						Real crl = clr - iinv.flux[fId];
						iinv.spu[lc] += crl;
						iinv.bu[lc] += crl * iinv.ub[fId];
						iinv.bv[lc] += crl * iinv.vb[fId];
						iinv.bw[lc] += crl * iinv.wb[fId];
					}

					Real ub1 = iinv.ub[fId];
					Real vb1 = iinv.vb[fId];
					Real wb1 = iinv.wb[fId];

					DirechletMom(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, ub1, vb1, wb1, fId);
				}
				else if (ug.bctype == BC::OUTFLOW)
				{
					if (iinv.flux[fId] < 0)
					{
						Real crl = clr - iinv.flux[fId];
						iinv.spu[lc] += crl;
						iinv.bu[lc] += crl * iinv.ub[fId];
						iinv.bv[lc] += crl * iinv.vb[fId];
						iinv.bw[lc] += crl * iinv.wb[fId];
					}

					Real ub1 = iinv.ub[fId];
					Real vb1 = iinv.vb[fId];
					Real wb1 = iinv.wb[fId];

					DirechletMom(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, ub1, vb1, wb1, fId);
				}
			}
		}
	}

	else if (Equa_vary == "energy")
	{
		;
	}
}

EndNameSpace