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
#include "INsVisterm.h"
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
	int nEqu = inscom.nEqu;

	q.resize(nEqu);
}

INsInvterm::INsInvterm()
{
	;
}

INsInvterm::~INsInvterm()
{
	;
}

void INsInvterm::Solve()
{
}


void INsInvterm::CmpINsinvTerm()
{
	Real clr = MAX(0, iinv.fq[ug.fId]);   
	Real crl = clr - iinv.fq[ug.fId];

	iinv.ai[ug.fId][0] = crl;
	iinv.ai[ug.fId][1] = clr;

}

void INsInvterm::CmpINsBcinvTerm()
{

	Real clr = MAX(0, iinv.fq[ug.fId]);
	Real crl = clr-iinv.fq[ug.fId];

	iinv.spc[ug.lc] += crl;

	iinv.buc[ug.lc] += crl * iinv.uf[ug.fId];

	iinv.bvc[ug.lc] += crl * iinv.vf[ug.fId];

	iinv.bwc[ug.lc] += crl * iinv.wf[ug.fId];
}

void INsInvterm::CmpINsFaceflux(RealField & dpdx, RealField & dpdy, RealField & dpdz)
{
	INsExtractl(*uinsf.q, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtractr(*uinsf.q, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	Real l2rdx = (*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc];
	Real l2rdy = (*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc];
	Real l2rdz = (*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc];
	Real rurf = 0.8/(1+0.8);

	iinv.VdU[ug.lc] = (*ug.cvol)[ug.lc] / iinv.spc[ug.lc];
	iinv.VdU[ug.rc] = (*ug.cvol)[ug.rc] / iinv.spc[ug.rc];

	iinv.Vdvu[ug.fId] = (*ug.fl)[ug.fId] * iinv.VdU[ug.lc] + (*ug.fr)[ug.fId] * iinv.VdU[ug.rc];

	Real dist = (*ug.a1)[ug.fId] * l2rdx + (*ug.a2)[ug.fId] * l2rdy + (*ug.a3)[ug.fId] * l2rdz;

	Real Df1 = iinv.Vdvu[ug.fId] *(*ug.a1)[ug.fId] / dist;
	Real Df2 = iinv.Vdvu[ug.fId] *(*ug.a2)[ug.fId] / dist;
	Real Df3 = iinv.Vdvu[ug.fId] *(*ug.a3)[ug.fId] / dist;

	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dx2 = (*ug.xcc)[ug.rc] - (*ug.xfc)[ug.fId];
	Real dy2 = (*ug.ycc)[ug.rc] - (*ug.yfc)[ug.fId];
	Real dz2 = (*ug.zcc)[ug.rc] - (*ug.zfc)[ug.fId];

	Real fdpdx = dpdx[ug.lc] * dx1 + dpdx[ug.rc] * dx2 - (iinv.pr - iinv.pl);
	Real fdpdy = dpdy[ug.lc] * dy1 + dpdy[ug.rc] * dy2 - (iinv.pr - iinv.pl);
	Real fdpdz = dpdz[ug.lc] * dz1 + dpdz[ug.rc] * dz2 - (iinv.pr - iinv.pl);

	iinv.uf[ug.fId] = iinv.ul * (*ug.fl)[ug.fId] + iinv.ur * (*ug.fr)[ug.fId];
	iinv.vf[ug.fId] = iinv.vl * (*ug.fl)[ug.fId] + iinv.vr * (*ug.fr)[ug.fId];
	iinv.wf[ug.fId] = iinv.wl * (*ug.fl)[ug.fId] + iinv.wr * (*ug.fr)[ug.fId];

	iinv.rf = (*ug.fl)[ug.fId] *iinv.rl + (*ug.fr)[ug.fId]*iinv.rr;
	iinv.vnflow = (*ug.a1)[ug.fId] * (iinv.uf[ug.fId] + fdpdx * Df1) + (*ug.a2)[ug.fId] * (iinv.vf[ug.fId] + fdpdy * Df2) + (*ug.a3)[ug.fId] * (iinv.wf[ug.fId] + fdpdz * Df3) +rurf*iinv.dun[ug.fId];
	iinv.fq[ug.fId] = iinv.rf * iinv.vnflow;  

}


void INsInvterm::CmpINsBcFaceflux(RealField& dpdx, RealField& dpdy, RealField& dpdz)
{
	INsExtractl(*uinsf.q, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);

	if (ug.bctype == BC::SOLID_SURFACE)
	{
		if (ug.bcdtkey == 0)
		{
			iinv.rf = iinv.rl;    

			iinv.uf[ug.fId] = (*ug.vfx)[ug.fId];

			iinv.vf[ug.fId] = (*ug.vfy)[ug.fId];

			iinv.wf[ug.fId] = (*ug.vfz)[ug.fId];

			iinv.vnflow = (*ug.a1)[ug.fId] * iinv.uf[ug.fId] + (*ug.a2)[ug.fId] * iinv.vf[ug.fId] + (*ug.a3)[ug.fId] * iinv.wf[ug.fId];

			iinv.fq[ug.fId] = iinv.rf * iinv.vnflow;
		}
		else
		{
			iinv.rf = iinv.rl;    

			iinv.uf[ug.fId] = (*inscom.bcflow)[IIDX::IIU];

			iinv.vf[ug.fId] = (*inscom.bcflow)[IIDX::IIV];

			iinv.wf[ug.fId] = (*inscom.bcflow)[IIDX::IIW];

			iinv.vnflow = (*ug.a1)[ug.fId] * iinv.uf[ug.fId] + (*ug.a2)[ug.fId] * iinv.vf[ug.fId] + (*ug.a3)[ug.fId] * iinv.wf[ug.fId];

			iinv.fq[ug.fId] = iinv.rf * iinv.vnflow;
		}

	}

	else if (ug.bctype == BC::INFLOW)
	{
		iinv.rf = inscom.inflow[IIDX::IIR];    

		iinv.uf[ug.fId] = inscom.inflow[IIDX::IIU];

		iinv.vf[ug.fId] = inscom.inflow[IIDX::IIV];

		iinv.wf[ug.fId] = inscom.inflow[IIDX::IIW];

		iinv.vnflow = (*ug.a1)[ug.fId] * iinv.uf[ug.fId] + (*ug.a2)[ug.fId] * iinv.vf[ug.fId] + (*ug.a3)[ug.fId] * iinv.wf[ug.fId];

		iinv.fq[ug.fId] = iinv.rf * iinv.vnflow;
	}

	else if (ug.bctype == BC::OUTFLOW)
	{

		Real l2rdx = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real l2rdy = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real l2rdz = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];
		Real rurf = 0.8 / (1 + 0.8);

		iinv.VdU[ug.lc] = (*ug.cvol)[ug.lc] / iinv.spc[ug.lc];

		iinv.Vdvu[ug.fId] = iinv.VdU[ug.lc];

		Real dist = (*ug.a1)[ug.fId] * l2rdx + (*ug.a2)[ug.fId] * l2rdy + (*ug.a3)[ug.fId] * l2rdz;

		Real Df1 = iinv.Vdvu[ug.fId] * (*ug.a1)[ug.fId] / dist;
		Real Df2 = iinv.Vdvu[ug.fId] * (*ug.a2)[ug.fId] / dist;
		Real Df3 = iinv.Vdvu[ug.fId] * (*ug.a3)[ug.fId] / dist;

		Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real fdpdx = dpdx[ug.lc] * dx1 - (iinv.pf[ug.fId] - iinv.pl);
		Real fdpdy = dpdy[ug.lc] * dy1 - (iinv.pf[ug.fId] - iinv.pl);
		Real fdpdz = dpdz[ug.lc] * dz1 - (iinv.pf[ug.fId] - iinv.pl);

		iinv.uf[ug.fId] += fdpdx * Df1;
		iinv.vf[ug.fId] += fdpdy * Df2;
		iinv.wf[ug.fId] += fdpdz * Df3;

		iinv.rf = iinv.rl;
		iinv.vnflow = (*ug.a1)[ug.fId] * (iinv.uf[ug.fId]) + (*ug.a2)[ug.fId] * (iinv.vf[ug.fId]) + (*ug.a3)[ug.fId] * (iinv.wf[ug.fId]) + rurf * iinv.dun[ug.fId];
		iinv.fq[ug.fId] = iinv.rf * iinv.vnflow;
	}

	else if (ug.bctype == BC::SYMMETRY)
	{

		iinv.rf = iinv.rl;

		iinv.uf[ug.fId] = 0;

		iinv.vf[ug.fId] = 0;

		iinv.wf[ug.fId] = 0;

		iinv.vnflow = (*ug.a1)[ug.fId] * iinv.uf[ug.fId] + (*ug.a2)[ug.fId] * iinv.vf[ug.fId] + (*ug.a3)[ug.fId] * iinv.wf[ug.fId];

		iinv.fq[ug.fId] = iinv.rf * iinv.vnflow;
	}

}

void INsInvterm::CmpINsFaceCorrectPresscoef()
{
		
		Real duf = (*ug.fl)[ug.fId] * ((*ug.cvol)[ug.lc] / iinv.dup[ug.lc]) + (*ug.fr)[ug.fId] * ((*ug.cvol)[ug.rc] / iinv.dup[ug.rc]);
		Real Sf1 = duf * (*ug.a1)[ug.fId];
		Real Sf2 = duf * (*ug.a2)[ug.fId];
		Real Sf3 = duf * (*ug.a3)[ug.fId];

		Real l2rdx = (*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc];
		Real l2rdy = (*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc];
		Real l2rdz = (*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc];

		Real dist = l2rdx * (*ug.a1)[ug.fId] + l2rdy * (*ug.a2)[ug.fId] + l2rdz * (*ug.a3)[ug.fId];

		Real Sfarea = Sf1 * (*ug.a1)[ug.fId] + Sf2 * (*ug.a2)[ug.fId] + Sf3 * (*ug.a3)[ug.fId];

		iinv.rf = (*ug.fl)[ug.fId] * (*uinsf.q)[IIDX::IIR][ug.lc] + ((*ug.fr)[ug.fId]) * (*uinsf.q)[IIDX::IIR][ug.rc];

		iinv.spp[ug.lc] += iinv.rf * Sfarea / dist;
		iinv.spp[ug.rc] += iinv.rf * Sfarea / dist;
		iinv.ajp[ug.fId][0] = iinv.rf * Sfarea / dist;
		iinv.ajp[ug.fId][1] = iinv.rf * Sfarea / dist;

		iinv.bp[ug.lc] -= iinv.fq[ug.fId];
		iinv.bp[ug.rc] += iinv.fq[ug.fId];
}

void INsInvterm::CmpINsBcFaceCorrectPresscoef()
{

	Real duf = (*ug.cvol)[ug.lc] / iinv.dup[ug.lc];
	Real Sf1 = duf * (*ug.a1)[ug.fId];
	Real Sf2 = duf * (*ug.a2)[ug.fId];
	Real Sf3 = duf * (*ug.a3)[ug.fId];

	Real l2rdx = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real l2rdy = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real l2rdz = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dist = l2rdx * (*ug.a1)[ug.fId] + l2rdy * (*ug.a2)[ug.fId] + l2rdz * (*ug.a3)[ug.fId];

	Real Sfarea = Sf1 * (*ug.a1)[ug.fId] + Sf2 * (*ug.a2)[ug.fId] + Sf3 * (*ug.a3)[ug.fId];

	iinv.rf = (*uinsf.q)[IIDX::IIR][ug.lc];

	int bcType = ug.bcRecord->bcType[ug.fId];

	if (bcType == BC::OUTFLOW)
	{
		iinv.spp[ug.lc] += iinv.rf * Sfarea / dist;
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

	iinv.bp[ug.lc] -= iinv.fq[ug.fId];
}

EndNameSpace