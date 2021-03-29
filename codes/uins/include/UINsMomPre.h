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
#include "SIMPLEC.h"
#include "systemSolver.h"

BeginNameSpace(ONEFLOW)

class UINsFField;
class Limiter;
class LimField;
class SolveEqua;

class UINsMomPre : public INsInv
{
public:
	UINsMomPre();
    ~UINsMomPre();
public:
	void CmpINsPreflux();
    void CmpConv();
	void InitMomCoeff();
    void CmpConvTerm();
	void InConvCoff(int&fId);
	void BcConvCoff(Real &ub1, Real &vb1, Real &wb1, int&fId);
	void BcVelocity(RealField& ub, RealField& vb, RealField& wb);
	void CmpDiffus();
	void CmpDiffusTerm();
	void InDiffusCoff(RealField& dudx, RealField& dudy, RealField& dudz, RealField& dvdx, RealField& dvdy, RealField& dvdz, RealField& dwdx, RealField& dwdy, RealField& dwdz, int& fId);
	void BcDiffusCoff(RealField& dudx, RealField& dudy, RealField& dudz, RealField& dvdx, RealField& dvdy, RealField& dvdz, RealField& dwdx, RealField& dwdy, RealField& dwdz, Real& ub1, Real& vb1, Real& wb1, int& fId);
	void BcPressure(RealField& pb);
	void FaceVelocity(RealField& ub, RealField& vb, RealField& wb, RealField& uf, RealField& vf, RealField& wf);
	void FacePressure(RealField& pb, RealField& pf);
	void CmpTranst();
	void CmpSrc();
	void MomEquCoeff();
	void RelaxMom(Real a);
	void CmpINsMomRes();
	void Solveuvw();
	void CmpFaceflux();
	void CmpINsFaceflux(Real & dpdx1, Real & dpdx2, Real & dpdy1, Real & dpdy2, Real & dpdz1, Real & dpdz2, Real & rurf, int& fId);
	void CmpINsBcFaceflux(Real& dpdx1, Real& dpdy1, Real& dpdz1, Real& pb1, Real & rurf, int& fId);
};
//void PrimToQ(RealField & prim, Real gama, RealField & q);
//extern UINsMomPre NonZero;
EndNameSpace