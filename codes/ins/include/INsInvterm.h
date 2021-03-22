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
#include "HXArray.h"
BeginNameSpace( ONEFLOW )


class INsInv
{
public:
    INsInv();
    ~INsInv();
public:
    void Init();
/*public:
    Real gama;
    Real gama1;
    Real gama2;*/ 
public:
    //RealField prim, prim1, prim2;
    //RealField q, q1, q2;
    //RealField dq;
	RealField flux; //flux1, flux2;
	//RealField uf, vf, wf, pf;
	//RealField fq, VdU, Vdvu;
	RealField spu, bu, bv, bw;
	RealField spp, bp;
	RealField dun, dup;
	RealField  ppf, pp;
	//RealField  dppbdx, dppbdy, dppbdz;
	RealField2D ai;//ajp;

public:
    /*Real rl, ul, vl, wl, pl, hl, el;
    Real rr, ur, vr, wr, pr, hr, er;
    Real rm, um, vm, wm, pm, hm, em;
	Real rf, vnflow, fux;*/
	Real remax_up, remax_vp, remax_wp, remax_pp, res_u, res_v, res_w, res_p;

};

extern INsInv iinv;

//class INsInvterm
//{
//public:
//	INsInvterm();
//    ~INsInvterm();
//public:
    //typedef void (INsInvterm:: * InvtermPointer )();
//public:
    //void Solve();
//public:
    //void SetPointer( int schemeIndex );
	//InvtermPointer InvtermPointer;
//public:
	//void CmpINsinvFlux();
	//void CmpINsBcinvFlux();
	/*void InConvCoff(int&fId);
	void BcConvCoff(Real &ub1, Real &vb1, Real &wb1,int&fId);
	void CmpINsFaceflux(Real & dpdx1, Real & dpdx2, Real & dpdy1, Real & dpdy2, Real & dpdz1, Real & dpdz2, int& fId );
	void CmpINsBcFaceflux(Real& dpdx1, Real& dpdy1, Real& dpdz1, Real& pb1, int& fId);
	void CmpInPressCoeff(int& fId);
	void CmpBcPressCoeff(int& fId);*/
	//void Roe      ();
    //void RoeOld   (){};
    //void HybridRoe(){};
    //void Vanleer  ();
    //void Steger   ();
    //void Hlle     ();
    //void LaxFriedrichs();
    //void Ausmp    ();
    //void AusmpUp  ();
    //void Ausmdv   ();
    //void Ausmw    ();
    //void Ausmpw   ();
    //void Slau2();
//public:
   // void ModifyAbsoluteEigenvalue();
//};

/*void INsCmpEnthalpy( RealField & prim, Real gama, Real & enthalpy );
void INsCmpTotalEnthalpyChange( RealField & prim, Real & gama, RealField & dq, Real & dh );
void INsPrimToQ( RealField & prim, Real gama, RealField & q );
void INsQToPrim( RealField & q, Real gama, RealField & prim, RealField & temp );
void INsCmpInternalEnergy( RealField & prim, Real gama, Real & em );*/

EndNameSpace