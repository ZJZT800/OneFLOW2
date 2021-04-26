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

#include "FaceValue.h"
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

void FaceValue(RealField& phi,RealField& phib,RealField& phif)
{
	phif = 0;

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		phif[fId] = phib[fId];
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		phif[fId] = (*ug.fl)[fId] * phi[lc] + (*ug.fr)[fId] * phi[rc];
	}
}

void FvisCoeff(string &Equa_vary, RealField &diffus_coef, RealField &fdiffusb_coef,RealField&fvis_cof)
{
	for (int fId = ug.nBFace; fId < ug.nFace; fId++)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		if (Equa_vary == "vary_u"|| Equa_vary == "vary_v"|| Equa_vary == "vary_w")
		{
			fvis_cof[fId] = (*ug.fl)[fId] * diffus_coef[lc] + (*ug.fr)[fId] * diffus_coef[rc];
		}
		else if (Equa_vary == "vary_h")
		{
			;
		}
	}

	for (int fId = 0; fId < ug.nBFace; fId++)
	{
		int lc = (*ug.lcf)[fId];

		if (Equa_vary == "vary_u" || Equa_vary == "vary_v" || Equa_vary == "vary_w")
		{
			fdiffusb_coef[fId] =  diffus_coef[lc];

			fvis_cof[fId] = fdiffusb_coef[fId];
		}
		else if (Equa_vary == "vary_h")
		{
			;
		}
	}
}

EndNameSpace