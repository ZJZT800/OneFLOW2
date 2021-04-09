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

#include "UINsConvTerm.h"
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

UINsConvTerm::UINsConvTerm()
{
	;
}

UINsConvTerm::~UINsConvTerm()
{
	;
}


void UINsConvTerm::CmpConvTerm(string &Equa_vary, string &ischeme)
{
	if (Equa_vary == "mom")
	{
		if (ischeme == "FOU")
		{
			for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
			{
				Real clr = MAX(0, iinv.flux[fId]);
				Real crl = clr - iinv.flux[fId];

				iinv.ai[0][fId] = crl;
				iinv.ai[1][fId] = clr;
			}
		}

		else if (ischeme == "CENTRAL")
		{
			;
		}
		else if (ischeme == "SOU")
		{
			;
		}
	}
	else if (Equa_vary == "energy")
	{
		if (ischeme == "FOU")
		{
			for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
			{
				Real clr = MAX(0, iinv.flux[fId]);
				Real crl = clr - iinv.flux[fId];

				iinv.ai[0][fId] = crl;
				iinv.ai[1][fId] = clr;
			}
		}
		else if (ischeme == "CENTRAL")
		{
			;
		}
		else if (ischeme == "SOU")
		{
			;
		}

	}
}

EndNameSpace