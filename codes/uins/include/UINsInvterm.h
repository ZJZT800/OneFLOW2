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
#include "INsInvterm.h"
#include "systemSolver.h"
#include "poisson.h"

BeginNameSpace(ONEFLOW)

class UINsFField;
class Limiter;
class LimField;
class SolveMRhs;


class UINsInvterm : public INsInvterm
{
public:
    UINsInvterm();
    ~UINsInvterm();
public:
    void Alloc();
    void DeAlloc();
	void CmpINsTimestep();
	void CmpINsPreflux();
	void Init();
    void CmpConv();
	void InitMomCoeff();
    void CmpConvTerm();
	void BcVelocity(RealField& ub, RealField& vb, RealField& wb);
	void BcPressure(RealField& pb);
	void FaceVelocity(RealField& ub, RealField& vb, RealField& wb, RealField& uf, RealField& vf, RealField& wf);
	void FacePressure(RealField& pb, RealField& pf);
    void CmpInvFace();
    void CmpLimiter();
	void SolveEquation(RealField& sp, RealField2D& ai, RealField& b, RealField& x, Real res);
	void CmpFaceflux();
	void CmpFaceVelocityValue(RealField& uf, RealField& vf, RealField& wf);
	void CmpINsMomRes();
	void CmpINsPreRes();
	void PresEquCoeff();
	void InitPressCoeff();
	void maxmin(RealField& a, Real& max_a, Real& min_a);
	void CmpPressCorrectEqu();
	void UpdateFaceflux();
	void CmpUpdateINsFaceflux(int& fId);
    void CmpDun(int& fId);
	void CmpUpdateINsBcFaceflux(int& fId);
	void UpdateSpeed();
	void UpdateINsRes();
    void AddFlux();
    void PrepareFaceValue();
	void PrepareProFaceValue();
	void CmpPreGrad();
	//void CmpINsinvTerm();
    //void UpdateFaceInvFlux();
    void ReadTmp();
public:
    void GetQlQrField();
    void ReconstructFaceValueField();
    void BoundaryQlQrFixField();
    //void Init();
    void Solveuvw();
public:
    Limiter* limiter;
    LimField* limf;
    MRField* iinvflux;
public:
    Real Number;
};
//void PrimToQ(RealField & prim, Real gama, RealField & q);
extern UINsInvterm NonZero;
EndNameSpace