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


#pragma once
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

class Inner
{
public:
    Inner ();
    ~Inner();
public:
    void UpdateNsResiduals();
	void UpdateINsResiduals();
    void FieldInit();
    void SolveFlow();
	void SolveEnergy();
	void SolveTurb();
	void SolveMultiComp();
};

void INsCmpBc();
void INSCmpGamaT(int flag);
void INsCmpChemSrc();
void INsCmpTurbEnergy();
void PresEqu();
void InitField();
void Mom_pre();
void CmpMomConv();
void CmpMomDiffus();
void CmpMomBc();
void MomTranst();
void CmpMomSrc();
void MomPrimeEqu();
void MomRelaxation();
void SolveMom();
void CmpFaceflux();
void Pres_cor();
void PresEqu();
void SolvePress();
void UpdateCorSpeed();
void UpdateFaceflux();
void UpdateRes();


EndNameSpace