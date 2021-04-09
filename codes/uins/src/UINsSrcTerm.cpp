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

#include "UINsSrcTerm.h"
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

UINsSrcTerm::UINsSrcTerm()
{
	;
}

UINsSrcTerm::~UINsSrcTerm()
{
	;
}


void UINsSrcTerm::CmpMomSrcTerm(string &Equa_vary)
{
	if (Equa_vary == "mom")
	{
		RealField dpdx, dpdy, dpdz;
		RealField pb, pf;

		dpdx.resize(ug.nCell);
		dpdy.resize(ug.nCell);
		dpdz.resize(ug.nCell);

		pf.resize(ug.nFace);

		FaceValue(iinv.pb, pf, (*uinsf.p)[0]);
        CmpUnsGrad(pf, dpdx, dpdy, dpdz);

		for (int cId = 0; cId < ug.nCell; ++cId)
		{
			Real vol = (*ug.cvol)[cId];

			iinv.bu[cId] += -vol * dpdx[cId];
			iinv.bv[cId] += -vol * dpdy[cId];
			iinv.bw[cId] += -vol * dpdz[cId];
		}
	}
	else if (Equa_vary == "energy")
	{
		;
	}
}

EndNameSpace