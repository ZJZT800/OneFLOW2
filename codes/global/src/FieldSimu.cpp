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
#include "FieldSimu.h"
#include "Iteration.h"
#include "Ctrl.h"
#include "NsCom.h"
#include "UsdData.h"
#include "DataBase.h"
#include "MultiBlock.h"
#include "SolverMap.h"
#include "SIMPLEC.h"
#include "CmxTask.h"
#include "Multigrid.h"
#include "BcData.h"

BeginNameSpace( ONEFLOW )

void FieldSimu()
{
    InitFlowSimuGlobal();
    MultiBlock::LoadGridAndBuildLink();
    MultiBlock::ProcessFlowWallDist();
    SolverMap::CreateSolvers();
    InitializeSolver();
	
	int compressible = ONEFLOW::GetDataValue< int >("compressible");

	if (compressible == 1)
	{
		SIMPLECSolve();
	}
	else
	{
		MultigridSolve();
	}
}

void InitFlowSimuGlobal()
{
    vis_model.Init();
    ctrl.Init();
    Iteration::Init();
    usd.InitBasic();
    bcdata.Init();
}

void InitializeSolver()
{
    ONEFLOW::MsMgTask( "INIT_FLOWFIELD" );
}

EndNameSpace