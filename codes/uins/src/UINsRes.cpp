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

#include "UINsRes.h"
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

UINsRes::UINsRes()
{
	;
}

UINsRes::~UINsRes()
{
	;
}


void UINsRes::UpdateIterRes()
{

	std::cout << "iinv.remax_u:" << iinv.remax_u << std::endl;
	std::cout << "iinv.remax_v:" << iinv.remax_v << std::endl;
	std::cout << "iinv.remax_w:" << iinv.remax_w << std::endl;
	std::cout << "iinv.remax_pp:" << iinv.remax_pp << std::endl;

	ofstream fileres_up("residual_u.txt", ios::app);
	//fileres_p << "residual_p:" <<residual_p << endl;
	fileres_up << iinv.remax_u << endl;
	fileres_up.close();

	ofstream fileres_vp("residual_v.txt", ios::app);
	//fileres_p << "residual_p:" <<residual_p << endl;
	fileres_vp << iinv.remax_v << endl;
	fileres_vp.close();

	ofstream fileres_wp("residual_w.txt", ios::app);
	//fileres_p << "residual_p:" <<residual_p << endl;
	fileres_wp << iinv.remax_w << endl;
	fileres_wp.close();

	ofstream fileres_pp("residual_pp.txt", ios::app);
	//fileres_p << "residual_p:" <<residual_p << endl;
	fileres_pp << iinv.remax_pp << endl;
	fileres_pp.close();
}


EndNameSpace