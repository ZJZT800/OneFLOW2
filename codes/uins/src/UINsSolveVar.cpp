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

#include "UINsSolveVar.h"
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

UINsSolveVar::UINsSolveVar()
{
	;
}

UINsSolveVar::~UINsSolveVar()
{
	;
}


void UINsSolveVar::SolveVar(string &Equa_vary)
{
	if (Equa_vary == "mom")
	{
		RealField uCorrect, vCorrect, wCorrect;
		uCorrect.resize(ug.nCell);
		vCorrect.resize(ug.nCell);
		wCorrect.resize(ug.nCell);

		int MaxIter = GetDataValue< int >("EquaMomIter");
		Real Tol = GetDataValue< Real >("EquaMomTol");
		int mRestarts = GetDataValue< int >("EquaMomRestarts");
		string gt = GetDataValue< string >("EquaMomMethod");
		bool ifPrecond = true;
		SolveEqua(iinv.spu, iinv.ai, iinv.bu, uCorrect, iinv.res_u, gt, MaxIter, Tol, mRestarts, ifPrecond);
		SolveEqua(iinv.spu, iinv.ai, iinv.bv, vCorrect, iinv.res_v, gt, MaxIter, Tol, mRestarts, ifPrecond);
		SolveEqua(iinv.spu, iinv.ai, iinv.bw, wCorrect, iinv.res_w, gt, MaxIter, Tol, mRestarts, ifPrecond);

		for (int cId = 0; cId < ug.nCell; cId++)
		{
			(*uinsf.u)[0][cId] += uCorrect[cId];
			(*uinsf.v)[0][cId] += vCorrect[cId];
			(*uinsf.w)[0][cId] += wCorrect[cId];
		}
	}

	else if (Equa_vary == "press")
	{
		int MaxIter = GetDataValue< int >("EquaPressIter");
		Real Tol = GetDataValue< Real >("EquaPressTol");
		int mRestarts = GetDataValue< int >("EquaPressRestarts");
		string gt = GetDataValue< string >("EquaPressMethod");
		bool ifPrecond = true;
		SolveEqua(iinv.spp, iinv.ai, iinv.bp, iinv.pp, iinv.res_p, gt, MaxIter, Tol, mRestarts, ifPrecond);
		//std::cout << iinv.res_p << std::endl;

		//�߽絥Ԫ
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

	}

	else if (Equa_vary == "energy")
	{
		;
	}
}

EndNameSpace