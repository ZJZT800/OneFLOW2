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

#include "UINsTranst.h"
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

UINsTranst::UINsTranst()
{
	;
}

UINsTranst::~UINsTranst()
{
	;
}


void UINsTranst::CmpTranstTerm(string &Equa_vary)
{
	if (Equa_vary == "mom")
	{
		Real timestep = GetDataValue< Real >("global_dt");

		for (int cId = 0; cId < ug.nCell; ++cId)
		{
			iinv.spu[cId] += (*ug.cvol)[cId] * (*uinsf.rho)[0][cId] / timestep;  //矩阵对角线元素的非稳态项

			iinv.bu[cId] += (*ug.cvol)[cId] * (*uinsf.rho)[0][cId] * iinv.u_old[cId] / timestep; //源项的非稳态项
			iinv.bv[cId] += (*ug.cvol)[cId] * (*uinsf.rho)[0][cId] * iinv.v_old[cId] / timestep;
			iinv.bw[cId] += (*ug.cvol)[cId] * (*uinsf.rho)[0][cId] * iinv.w_old[cId] / timestep;
		}
	}
	else if (Equa_vary == "energy")
	{
		;
	}
}

void SaveOldValue()
{
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		iinv.u_old[cId] = (*uinsf.u)[0][cId];
		iinv.v_old[cId] = (*uinsf.v)[0][cId];
		iinv.w_old[cId] = (*uinsf.w)[0][cId];

		int solve_energy = GetDataValue< int >("solve_energy");
		if (solve_energy == 1)
		{
			;
		}
	}
}

EndNameSpace