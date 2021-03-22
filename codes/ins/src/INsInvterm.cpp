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

//#include "UINsCorrectSpeed.h"
#include "INsInvterm.h"
#include "Iteration.h"
#include "UINsCom.h"
#include "Zone.h"
#include "DataBase.h"
#include "UCom.h"
#include "UGrad.h"
#include "Com.h"
#include "INsCom.h"
#include "INsIDX.h"
#include "HXMath.h"
#include "Ctrl.h"
#include "Boundary.h"
#include "BcRecord.h"

BeginNameSpace(ONEFLOW)

INsInv iinv;



INsInv::INsInv()
{
	;
}

INsInv::~INsInv()
{
	;
}

void INsInv::Init()
{
	bu.resize(ug.nCell);
	bv.resize(ug.nCell);
	bw.resize(ug.nCell);
	bp.resize(ug.nCell);
	flux.resize(ug.nFace);
	spu.resize(ug.nCell);
	ai.resize(2, ug.nFace);
	spp.resize(ug.nCell);
	pp.resize(ug.nCell);
	ppf.resize(ug.nFace);
	dup.resize(ug.nCell);
	dun.resize(ug.nFace);
}


/*void INsInvterm::Solve()
{
}*/


/*void INsInvterm::InConvCoff(int&fId)
{
	Real clr = MAX(0, iinv.flux[fId]);   
	Real crl = clr - iinv.flux[fId];

	iinv.ai[0][fId] = crl;
	iinv.ai[1][fId] = clr;
}*/

/*void INsInvterm::BcConvCoff(Real &ub1, Real &vb1, Real &wb1, int&fId)
{
	int lc = (*ug.lcf)[fId];
	int bctype = ug.bcRecord->bcType[fId];

	Real clr = MAX(0, iinv.flux[fId]);
	Real crl = clr-iinv.flux[fId];

	if (bctype == BC::SOLID_SURFACE|| bctype == BC::SYMMETRY)
	{
		;
	}
	else if (bctype == BC::INFLOW)
	{
		if (iinv.flux[fId] < 0)
		{
			Real crl = clr - iinv.flux[fId];
			iinv.spu[lc] += crl;
			iinv.bu[lc] += crl * ub1;
			iinv.bv[lc] += crl * vb1;
			iinv.bw[lc] += crl * wb1;
		}
	}
	else if (bctype == BC::OUTFLOW)
	{
		if (iinv.flux[fId] < 0)
		{
			Real crl = clr - iinv.flux[fId];
			iinv.spu[lc] += crl;
			iinv.bu[lc] += crl * ub1;
			iinv.bv[lc] += crl * vb1;
			iinv.bw[lc] += crl * wb1;
		}
	}

}*/

/*void INsInvterm::CmpINsFaceflux(Real & dpdx1, Real & dpdx2, Real & dpdy1, Real & dpdy2, Real & dpdz1, Real & dpdz2, int& fId)
{
	int lc = (*ug.lcf)[fId];
	int rc = (*ug.rcf)[fId];
	
	Real rl, ul, vl, wl, pl;
	Real rr, ur, vr, wr, pr;

	Real VdU1, VdU2, Vdvu;

	Real uf, vf, wf;

	Real vnflow;

	INsExtractl(*uinsf.q, rl, ul, vl, wl, pl,lc);
	INsExtractr(*uinsf.q, rr, ur, vr, wr, pr,rc);

	Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
	Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
	Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];
	Real rurf = 0.8/(1+0.8);

	VdU1 = (*ug.cvol)[lc] / iinv.spu[lc];
	VdU2 = (*ug.cvol)[rc] / iinv.spu[rc];

	Vdvu = (*ug.fl)[fId] * VdU1 + (*ug.fr)[fId] * VdU2;

	Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

	Real Df1 = Vdvu *(*ug.a1)[fId] / dist;
	Real Df2 = Vdvu *(*ug.a2)[fId] / dist;
	Real Df3 = Vdvu *(*ug.a3)[fId] / dist;

	Real dx1 = (*ug.xfc)[fId] - (*ug.xcc)[lc];
	Real dy1 = (*ug.yfc)[fId] - (*ug.ycc)[lc];
	Real dz1 = (*ug.zfc)[fId] - (*ug.zcc)[lc];

	Real dx2 = (*ug.xcc)[rc] - (*ug.xfc)[fId];
	Real dy2 = (*ug.ycc)[rc] - (*ug.yfc)[fId];
	Real dz2 = (*ug.zcc)[rc] - (*ug.zfc)[fId];

	Real fdpdx = dpdx1 * dx1 + dpdx2 * dx2 - (pr - pl);
	Real fdpdy = dpdy1 * dy1 + dpdy2 * dy2 - (pr - pl);
	Real fdpdz = dpdz1 * dz1 + dpdz2 * dz2 - (pr - pl);

	uf = ul * (*ug.fl)[fId] + ur * (*ug.fr)[fId];
	vf = vl * (*ug.fl)[fId] + vr * (*ug.fr)[fId];
	wf = wl * (*ug.fl)[fId] + wr * (*ug.fr)[fId];

	vnflow = (*ug.a1)[fId] * (uf + fdpdx * Df1) + (*ug.a2)[fId] * (vf + fdpdy * Df2) + (*ug.a3)[fId] * (wf + fdpdz * Df3) +rurf*iinv.dun[fId];
	iinv.flux[fId] = rl * vnflow;  

}*/


/*void INsInvterm::CmpINsBcFaceflux(Real& dpdx1, Real& dpdy1, Real& dpdz1, Real& pb1, int& fId)
{
	int lc = (*ug.lcf)[fId];
	
	Real rl, ul, vl, wl, pl;

	Real VdU1, Vdvu;

	Real uf, vf, wf;

	Real vnflow;

	INsExtractl(*uinsf.q, rl, ul, vl, wl, pl, lc);

	if (ug.bctype == BC::SOLID_SURFACE)
	{
		if (ug.bcdtkey == 0)
		{

			uf = (*ug.vfx)[fId];

			vf = (*ug.vfy)[fId];

			wf = (*ug.vfz)[fId];

			vnflow = (*ug.a1)[fId] * uf + (*ug.a2)[fId] * vf + (*ug.a3)[fId] * wf;

			iinv.flux[fId] = rl * vnflow;
		}
		else
		{

			uf = (*inscom.bcflow)[IIDX::IIU];

			vf = (*inscom.bcflow)[IIDX::IIV];

			wf = (*inscom.bcflow)[IIDX::IIW];

			vnflow = (*ug.a1)[fId] * uf + (*ug.a2)[fId] * vf + (*ug.a3)[fId] * wf;

			iinv.flux[fId] = rl * vnflow;
		}

	}

	else if (ug.bctype == BC::INFLOW)
	{    
		uf = inscom.inflow[IIDX::IIU];

		vf = inscom.inflow[IIDX::IIV];

		wf = inscom.inflow[IIDX::IIW];

		vnflow = (*ug.a1)[fId] * uf + (*ug.a2)[fId] * vf + (*ug.a3)[fId] * wf;

		iinv.flux[fId] = inscom.inflow[IIDX::IIR] * vnflow;
	}

	else if (ug.bctype == BC::OUTFLOW)
	{

		Real l2rdx = (*ug.xfc)[fId] - (*ug.xcc)[ug.lc];
		Real l2rdy = (*ug.yfc)[fId] - (*ug.ycc)[ug.lc];
		Real l2rdz = (*ug.zfc)[fId] - (*ug.zcc)[ug.lc];
		Real rurf = 0.8 / (1 + 0.8);

		VdU1 = (*ug.cvol)[lc] / iinv.spu[lc];

		Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

		Real Df1 = VdU1 * (*ug.a1)[fId] / dist;
		Real Df2 = VdU1 * (*ug.a2)[fId] / dist;
		Real Df3 = VdU1 * (*ug.a3)[fId] / dist;

		Real dx1 = (*ug.xfc)[fId] - (*ug.xcc)[lc];
		Real dy1 = (*ug.yfc)[fId] - (*ug.ycc)[lc];
		Real dz1 = (*ug.zfc)[fId] - (*ug.zcc)[lc];

		Real fdpdx = dpdx1 * dx1 - (pb1 - pl);
		Real fdpdy = dpdy1 * dy1 - (pb1 - pl);
		Real fdpdz = dpdz1 * dz1 - (pb1 - pl);

		uf = ul+fdpdx * Df1;
		vf = vl+fdpdy * Df2;
		wf = wl+fdpdz * Df3;

		vnflow = (*ug.a1)[fId] * uf + (*ug.a2)[fId] * vf + (*ug.a3)[fId] * wf + rurf * iinv.dun[fId];
		iinv.flux[fId] = rl * vnflow;
	}

	else if (ug.bctype == BC::SYMMETRY)
	{
		uf = 0;

		vf = 0;

		wf = 0;

		iinv.flux[fId] = 0;
	}

}*/

/*void INsInvterm::CmpInPressCoeff(int& fId)
{
	   int lc = (*ug.lcf)[fId];
	   int rc = (*ug.rcf)[fId];
		
		Real duf = (*ug.fl)[fId] * ((*ug.cvol)[lc] / iinv.dup[lc]) + (*ug.fr)[fId] * ((*ug.cvol)[rc] / iinv.dup[rc]);
		Real Sf1 = duf * (*ug.a1)[fId];
		Real Sf2 = duf * (*ug.a2)[fId];
		Real Sf3 = duf * (*ug.a3)[fId];

		Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
		Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
		Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];

		Real dist = l2rdx * (*ug.a1)[fId] + l2rdy * (*ug.a2)[fId] + l2rdz * (*ug.a3)[fId];

		Real Sfarea = Sf1 * (*ug.a1)[fId] + Sf2 * (*ug.a2)[fId] + Sf3 * (*ug.a3)[fId];

		//iinv.rf = (*ug.fl)[ug.fId] * (*uinsf.q)[IIDX::IIR][ug.lc] + ((*ug.fr)[ug.fId]) * (*uinsf.q)[IIDX::IIR][ug.rc];

		iinv.spp[lc] += (*uinsf.q)[IIDX::IIR][lc] * Sfarea / dist;
		iinv.spp[rc] += (*uinsf.q)[IIDX::IIR][lc] * Sfarea / dist;
		iinv.ai[0][fId] = (*uinsf.q)[IIDX::IIR][lc] * Sfarea / dist;
		iinv.ai[1][fId] = (*uinsf.q)[IIDX::IIR][lc] * Sfarea / dist;

		iinv.bp[lc] -= iinv.flux[fId];
		iinv.bp[rc] += iinv.flux[fId];
}*/

/*void INsInvterm::CmpBcPressCoeff(int& fId)
{
	int lc = (*ug.lcf)[fId];

	Real duf = (*ug.cvol)[lc] / iinv.dup[lc];
	Real Sf1 = duf * (*ug.a1)[fId];
	Real Sf2 = duf * (*ug.a2)[fId];
	Real Sf3 = duf * (*ug.a3)[fId];

	Real l2rdx = (*ug.xfc)[fId] - (*ug.xcc)[lc];
	Real l2rdy = (*ug.yfc)[fId] - (*ug.ycc)[lc];
	Real l2rdz = (*ug.zfc)[fId] - (*ug.zcc)[lc];

	Real dist = l2rdx * (*ug.a1)[fId] + l2rdy * (*ug.a2)[fId] + l2rdz * (*ug.a3)[fId];

	Real Sfarea = Sf1 * (*ug.a1)[fId] + Sf2 * (*ug.a2)[fId] + Sf3 * (*ug.a3)[fId];

	int bcType = ug.bcRecord->bcType[fId];

	if (bcType == BC::OUTFLOW)
	{
		iinv.spp[lc] += (*uinsf.q)[IIDX::IIR][lc] * Sfarea / dist;
	}

	else if (bcType == BC::SOLID_SURFACE)
	{
		;
	}

	else if (bcType == BC::INFLOW)
	{
		;
	}
	
	else if (ug.bctype == BC::SYMMETRY)
	{
		;
	}

	iinv.bp[lc] -= iinv.flux[fId];
}*/

EndNameSpace