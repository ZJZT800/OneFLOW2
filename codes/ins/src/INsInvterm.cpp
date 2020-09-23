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
#include "Com.h"
#include "INsCom.h"
#include "INsIDX.h"
#include "HXMath.h"
#include "Ctrl.h"
#include "Boundary.h"
#include "BcRecord.h"

BeginNameSpace( ONEFLOW )

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
    prim.resize(nEqu);
    prim1.resize( nEqu );
    prim2.resize( nEqu );

    q.resize( nEqu );
    q1.resize( nEqu );
    q2.resize( nEqu );

    dq.resize( nEqu );

    flux.resize( nEqu );
    flux1.resize( nEqu );
    flux2.resize( nEqu );

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

void INsInvterm::CmpINsinvFlux()
{

		INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);

		INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

		Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real dx2 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
		Real dy2 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
		Real dz2 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

		Real de1 = DIST(dx1, dy1, dz1);
		Real de2 = DIST(dx2, dy2, dz2);
		Real de = 1.0 / (de1 + de2);

		iinv.f1[ug.fId] = de2 * de;  //左单元权重
		iinv.f2[ug.fId] = de1 * de;  //右单元权重

		iinv.rf[ug.fId] = iinv.f1[ug.fId]*iinv.rl + iinv.f2[ug.fId]*iinv.rr;    //初始界面上的值（u、v、w ）

		iinv.uf[ug.fId] = iinv.f1[ug.fId] * iinv.ul + iinv.f2[ug.fId] * iinv.ur;

		iinv.vf[ug.fId] = iinv.f1[ug.fId] * iinv.vl + iinv.f2[ug.fId] * iinv.vr;

		iinv.wf[ug.fId] = iinv.f1[ug.fId] * iinv.wl + iinv.f2[ug.fId] * iinv.wr;

		iinv.pf[ug.fId] = iinv.f1[ug.fId] * iinv.pl + iinv.f2[ug.fId] * iinv.pr;

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;  //初始界面上 V*n

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea; //初始界面上的质量通量

}

void INsInvterm::CmpINsBcinvFlux()
{

	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);

	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dx2 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
	Real dy2 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
	Real dz2 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

	Real de1 = DIST(dx1, dy1, dz1);
	Real de2 = DIST(dx2, dy2, dz2);
	Real de = 1.0 / (de1 + de2);

	iinv.f1[ug.fId] = de2 * de;  //左单元权重
	iinv.f2[ug.fId] = de1 * de;  //右单元权重

	if(ug.bctype == BC::SOLID_SURFACE)
	{
		if (inscom.bcdtkey == 0)
		{
			iinv.rf[ug.fId] = iinv.rl;    //初始界面上的值（u、v、w ）

			iinv.uf[ug.fId] = (*ug.vfx)[ug.fId];

			iinv.vf[ug.fId] = (*ug.vfy)[ug.fId];

			iinv.wf[ug.fId] = (*ug.vfz)[ug.fId];

			iinv.pf[ug.fId] = iinv.pl;

			iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

			iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;

		}
		else
		{
			iinv.rf[ug.fId] = iinv.rl;    //初始界面上的值（u、v、w ）

			iinv.uf[ug.fId] = (*inscom.bcflow)[IIDX::IIU];

			iinv.vf[ug.fId] = (*inscom.bcflow)[IIDX::IIV];

			iinv.wf[ug.fId] = (*inscom.bcflow)[IIDX::IIW];

			iinv.pf[ug.fId] = iinv.pl;

			iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

			iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;
		}

	}
	else if (ug.bctype == BC::OUTFLOW)
	{

		iinv.rf[ug.fId] = iinv.rl;    //初始界面上的值（u、v、w ）

		iinv.uf[ug.fId] = iinv.ul;

		iinv.vf[ug.fId] = iinv.vl;

		iinv.wf[ug.fId] = iinv.wl;

		iinv.pf[ug.fId] = iinv.pl;

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId]* iinv.vnflow * gcom.farea;
	}

	else if(ug.bctype == BC::INFLOW)
	{
		iinv.rf[ug.fId] = inscom.inflow[IIDX::IIR];    //初始界面上的值（u、v、w ）

		iinv.uf[ug.fId] = inscom.inflow[IIDX::IIU];

		iinv.vf[ug.fId] = inscom.inflow[IIDX::IIV];

		iinv.wf[ug.fId] = inscom.inflow[IIDX::IIW];

		iinv.pf[ug.fId] = inscom.inflow[IIDX::IIP];

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea; //初始界面上的质量通量
	}
}

void INsInvterm::CmpINsinvTerm()
{ 
		Real clr = MAX(0, iinv.fq[ug.fId]);  //从界面左侧单元流入右侧单元的质量流量

		Real crl = clr - iinv.fq[ug.fId];   //从界面右侧单元流入左侧单元的质量流量
		
		iinv.ai[ug.fId][0] = crl;//非对角线元素
		iinv.ai[ug.fId][1] = clr;//非对角线元素

		iinv.spc[ug.lc] += clr;
		iinv.spc[ug.rc] += crl;
}

void INsInvterm::CmpINsBcinvTerm()
{
	
		Real clr = MAX(0, iinv.fq[ug.fId]);  //从界面左侧单元流入右侧单元的初始质量流量
		Real crl = clr - iinv.fq[ug.fId];   //从界面右侧单元流入左侧单元的初始质量流量


		iinv.ai[ug.fId][0] = crl;//clr;
		iinv.ai[ug.fId][1] = clr;//crl;

		int bcType = ug.bcRecord->bcType[ug.fId];

		if (bcType == BC::SOLID_SURFACE)
		{
			;
		}

		if (bcType == BC::INFLOW)
		{
			if (iinv.fq[ug.fId] < 0)
			{
				iinv.spc[ug.lc] += crl;

				iinv.buc[ug.lc] += crl * iinv.uf[ug.fId];
				iinv.bvc[ug.lc] += crl * iinv.vf[ug.fId];
				iinv.bwc[ug.lc] += crl * iinv.wf[ug.fId];
			}
		}
		else if (bcType == BC::OUTFLOW)
		{
				iinv.spc[ug.lc] += clr;

				iinv.buc[ug.lc] += crl * iinv.uf[ug.fId];
				iinv.bvc[ug.lc] += crl * iinv.vf[ug.fId];
				iinv.bwc[ug.lc] += crl * iinv.wf[ug.fId];
		}
}

void INsInvterm::CmpINsFaceflux()
{

	INsExtractl(*uinsf.q, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtractr(*uinsf.q, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	// Rie-Chow插值
	Real l2rdx = (*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc];
	Real l2rdy = (*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc];
	Real l2rdz = (*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc];

	iinv.VdU[ug.lc] = (*ug.cvol1)[ug.lc] / iinv.spc[ug.lc];
	iinv.VdU[ug.rc] = (*ug.cvol1)[ug.rc] / iinv.spc[ug.rc];
	iinv.VdV[ug.lc] = (*ug.cvol1)[ug.lc] / iinv.spc[ug.lc];
	iinv.VdV[ug.rc] = (*ug.cvol1)[ug.rc] / iinv.spc[ug.rc];
	iinv.VdW[ug.lc] = (*ug.cvol1)[ug.lc] / iinv.spc[ug.lc];
	iinv.VdW[ug.rc] = (*ug.cvol1)[ug.rc] / iinv.spc[ug.rc];

	iinv.Vdvu[ug.fId] = (*ug.fl)[ug.fId] * iinv.VdU[ug.lc] + (*ug.fr)[ug.fId] * iinv.VdU[ug.rc];
	iinv.Vdvv[ug.fId] = (*ug.fl)[ug.fId] * iinv.VdV[ug.lc] + (*ug.fr)[ug.fId] * iinv.VdV[ug.rc];
	iinv.Vdvw[ug.fId] = (*ug.fl)[ug.fId] * iinv.VdW[ug.lc] + (*ug.fr)[ug.fId] * iinv.VdW[ug.rc];

	Real dist = (*ug.a1)[ug.fId] * l2rdx + (*ug.a2)[ug.fId] * l2rdy + (*ug.a3)[ug.fId] * l2rdz;

	Real Df1 = iinv.Vdvu[ug.fId] * (*ug.a1)[ug.fId] / dist;
	Real Df2 = iinv.Vdvv[ug.fId] * (*ug.a2)[ug.fId] / dist;
	Real Df3 = iinv.Vdvw[ug.fId] * (*ug.a3)[ug.fId] / dist;

	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dx2 = (*ug.xcc)[ug.rc] - (*ug.xfc)[ug.fId];
	Real dy2 = (*ug.ycc)[ug.rc] - (*ug.yfc)[ug.fId];
	Real dz2 = (*ug.zcc)[ug.rc] - (*ug.zfc)[ug.fId];

	Real fdpdx = iinv.dpdx[ug.lc] * dx1 + iinv.dpdx[ug.rc] * dx2 - (iinv.pr - iinv.pl);
	Real fdpdy = iinv.dpdy[ug.lc] * dy1 + iinv.dpdy[ug.rc] * dy2 - (iinv.pr - iinv.pl);
	Real fdpdz = iinv.dpdz[ug.lc] * dz1 + iinv.dpdz[ug.rc] * dz2 - (iinv.pr - iinv.pl);

	iinv.rf[ug.fId] = half * (iinv.rl + iinv.rr);
	iinv.uf[ug.fId] = iinv.ul * (*ug.fl)[ug.fId] + iinv.ur * (*ug.fr)[ug.fId];
	iinv.vf[ug.fId] = iinv.vl * (*ug.fl)[ug.fId] + iinv.vr * (*ug.fr)[ug.fId];
	iinv.wf[ug.fId] = iinv.wl * (*ug.fl)[ug.fId] + iinv.wr * (*ug.fr)[ug.fId];


	iinv.vnflow = (*ug.a1)[ug.fId] * (iinv.uf[ug.fId]+ fdpdx * Df1) + (*ug.a2)[ug.fId] * (iinv.vf[ug.fId]+ fdpdy * Df2) + (*ug.a3)[ug.fId] *(iinv.wf[ug.fId]+ fdpdz * Df3) - (*ug.vfn)[ug.fId] - iinv.dun[ug.fId];
	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow;


	/*INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	
	iinv.Vau = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * gcom.xfn / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * gcom.xfn / (iinv.spc[ug.rc])); //Df*n，分子
	iinv.Vav = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * gcom.yfn / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * gcom.yfn / (iinv.spc[ug.rc]));
	iinv.Vaw = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * gcom.zfn / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * gcom.zfn / (iinv.spc[ug.rc]));

	iinv.dist = gcom.xfn * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + gcom.yfn * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + gcom.zfn * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	iinv.Deun = iinv.Vau / iinv.dist;   //Df*n/e*n
	iinv.Devn = iinv.Vav / iinv.dist;
	iinv.Dewn = iinv.Vaw / iinv.dist;


	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dx2 = (*ug.xcc)[ug.rc] - (*ug.xfc)[ug.fId];
	Real dy2 = (*ug.ycc)[ug.rc] - (*ug.yfc)[ug.fId];
	Real dz2 = (*ug.zcc)[ug.rc] - (*ug.zfc)[ug.fId];

	iinv.Bpe = (*uinsf.dqdx)[IIDX::IIP][ug.lc] * (dx1) + (*uinsf.dqdy)[IIDX::IIP][ug.lc] * (dy1) + (*uinsf.dqdz)[IIDX::IIP][ug.lc] * (dz1) +
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] * (dx2) + (*uinsf.dqdy)[IIDX::IIP][ug.rc] * (dy2) + (*uinsf.dqdz)[IIDX::IIP][ug.rc] * (dz2) -
		(iinv.pr - iinv.pl);

	iinv.rf[ug.fId] = half * (iinv.rl+ iinv.rr);

	iinv.uf[ug.fId] = iinv.f1[ug.fId] * iinv.ul + iinv.f2[ug.fId] * iinv.ur;  //下一时刻的界面预测速度
	iinv.vf[ug.fId] = iinv.f1[ug.fId] * iinv.vl + iinv.f2[ug.fId] * iinv.vr;
	iinv.wf[ug.fId] = iinv.f1[ug.fId] * iinv.wl + iinv.f2[ug.fId] * iinv.wr;

	iinv.vnflow = (*ug.xfn)[ug.fId] * (iinv.uf[ug.fId]+iinv.Deun * iinv.Bpe) + (*ug.yfn)[ug.fId] * (iinv.vf[ug.fId]+iinv.Devn * iinv.Bpe) + (*ug.zfn)[ug.fId] * (iinv.wf[ug.fId]+iinv.Dewn * iinv.Bpe);

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * (*ug.farea)[ug.fId];  //下一时刻界面预测通量*/

}


void INsInvterm::CmpINsBcFaceflux()
{

	INsExtractl(*uinsf.q, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);

	

	if (ug.bctype == BC::SOLID_SURFACE)
	{
		if (inscom.bcdtkey == 0)
		{
			iinv.rf[ug.fId] = iinv.rl;

			iinv.uf[ug.fId] = (*ug.vfx)[ug.fId];

			iinv.vf[ug.fId] = (*ug.vfy)[ug.fId];

			iinv.wf[ug.fId] = (*ug.vfz)[ug.fId];

			iinv.vnflow = (*ug.a1)[ug.fId] * iinv.uf[ug.fId] + (*ug.a2)[ug.fId] * iinv.vf[ug.fId] + (*ug.a3)[ug.fId] * iinv.wf[ug.fId] - (*ug.vfn)[ug.fId] - iinv.dun[ug.fId];

			iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow;
		}
		else
		{
			iinv.rf[ug.fId] = iinv.rl;

			iinv.uf[ug.fId] = (*inscom.bcflow)[IIDX::IIU];

			iinv.vf[ug.fId] = (*inscom.bcflow)[IIDX::IIV];

			iinv.wf[ug.fId] = (*inscom.bcflow)[IIDX::IIW];

			iinv.vnflow = (*ug.a1)[ug.fId] * iinv.uf[ug.fId] + (*ug.a2)[ug.fId] * iinv.vf[ug.fId] + (*ug.a3)[ug.fId] * iinv.wf[ug.fId] - (*ug.vfn)[ug.fId] - iinv.dun[ug.fId];

			iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow;
		}

	}

	else if (ug.bctype == BC::INFLOW)
	{
		iinv.rf[ug.fId] = inscom.inflow[IIDX::IIR];

		iinv.uf[ug.fId] = inscom.inflow[IIDX::IIU];

		iinv.vf[ug.fId] = inscom.inflow[IIDX::IIV];

		iinv.wf[ug.fId] = inscom.inflow[IIDX::IIW];

		iinv.vnflow = (*ug.a1)[ug.fId] * iinv.uf[ug.fId] + (*ug.a2)[ug.fId] * iinv.vf[ug.fId] + (*ug.a3)[ug.fId] * iinv.wf[ug.fId] - (*ug.vfn)[ug.fId] - iinv.dun[ug.fId];

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow;
	}

	else if (ug.bctype == BC::OUTFLOW)
	{
		Real l2rdx = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real l2rdy = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real l2rdz = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		iinv.VdU[ug.lc] = (*ug.cvol1)[ug.lc] / iinv.spc[ug.lc];
		iinv.VdV[ug.lc] = (*ug.cvol1)[ug.lc] / iinv.spc[ug.lc];
		iinv.VdW[ug.lc] = (*ug.cvol1)[ug.lc] / iinv.spc[ug.lc];

		iinv.Vdvu[ug.fId] = iinv.VdU[ug.lc];
		iinv.Vdvv[ug.fId] = iinv.VdV[ug.lc];
		iinv.Vdvw[ug.fId] = iinv.VdW[ug.lc];

		Real dist = (*ug.a1)[ug.fId] * l2rdx + (*ug.a2)[ug.fId] * l2rdy + (*ug.a3)[ug.fId] * l2rdz;

		Real Df1 = iinv.Vdvu[ug.fId] * (*ug.a1)[ug.fId] / dist;
		Real Df2 = iinv.Vdvv[ug.fId] * (*ug.a2)[ug.fId] / dist;
		Real Df3 = iinv.Vdvw[ug.fId] * (*ug.a3)[ug.fId] / dist;

		Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real fdpdx = iinv.dpdx[ug.lc] * dx1 - (iinv.pr - iinv.pl);
		Real fdpdy = iinv.dpdy[ug.lc] * dy1 - (iinv.pr - iinv.pl);
		Real fdpdz = iinv.dpdz[ug.lc] * dz1 - (iinv.pr - iinv.pl);

		iinv.rf[ug.fId] = iinv.rl;

		iinv.uf[ug.cId] = iinv.ul + Df1;

		iinv.vf[ug.cId] = iinv.vl + Df2;

		iinv.wf[ug.cId] = iinv.wl + Df3;

		iinv.vnflow = (*ug.a1)[ug.fId] * (iinv.ul+ fdpdx * Df1) + (*ug.a2)[ug.fId] * (iinv.vf[ug.fId]+ fdpdy * Df2) + (*ug.a3)[ug.fId] * (iinv.wf[ug.fId]+ fdpdz * Df3) - (*ug.vfn)[ug.fId] - iinv.dun[ug.fId];

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow;
	}


	/*INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	iinv.Vau = ((*ug.cvol1)[ug.lc] * gcom.xfn / (iinv.spc[ug.lc])); //Df*n，分子
	iinv.Vav = ((*ug.cvol1)[ug.lc] * gcom.yfn / (iinv.spc[ug.lc]));
	iinv.Vaw = ((*ug.cvol1)[ug.lc] * gcom.zfn / (iinv.spc[ug.lc]));

	iinv.dist = gcom.xfn * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + gcom.yfn * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) + gcom.zfn * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);

	iinv.Deun = iinv.Vau / iinv.dist;   //Df*n/e*n
	iinv.Devn = iinv.Vav / iinv.dist;
	iinv.Dewn = iinv.Vaw / iinv.dist;
	
	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	iinv.Bpe = (*uinsf.dqdx)[IIDX::IIP][ug.lc] * (dx1)+(*uinsf.dqdy)[IIDX::IIP][ug.lc] * (dy1)+(*uinsf.dqdz)[IIDX::IIP][ug.lc] * (dz1)-
		(iinv.pf[ug.fId] - iinv.pl);

	if (ug.bctype == BC::SOLID_SURFACE)
	{
		if (inscom.bcdtkey == 0)
		{
			iinv.rf[ug.fId] = iinv.rl;    //初始界面上的值（u、v、w ）

			iinv.uf[ug.fId] = (*ug.vfx)[ug.fId];

			iinv.vf[ug.fId] = (*ug.vfy)[ug.fId];

			iinv.wf[ug.fId] = (*ug.vfz)[ug.fId];

			iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

			iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;

		}
		else
		{
			iinv.rf[ug.fId] = iinv.rl;    //初始界面上的值（u、v、w ）

			iinv.uf[ug.fId] = (*inscom.bcflow)[IIDX::IIU];

			iinv.vf[ug.fId] = (*inscom.bcflow)[IIDX::IIV];

			iinv.wf[ug.fId] = (*inscom.bcflow)[IIDX::IIW];


			iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

			iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;
			
		}

	}

	else if (ug.bctype == BC::INFLOW)
	{
		iinv.rf[ug.fId] =  inscom.inflow[IIDX::IIR];    //初始界面上的值（u、v、w ）

		iinv.uf[ug.fId] =  inscom.inflow[IIDX::IIU];

		iinv.vf[ug.fId] =  inscom.inflow[IIDX::IIV];

		iinv.wf[ug.fId] =  inscom.inflow[IIDX::IIW];

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea; //初始界面上的质量通量

	}

	else if (ug.bctype == BC::OUTFLOW)
	{

		iinv.rf[ug.fId] = iinv.rl;    //初始界面上的值（u、v、w ）

		iinv.uf[ug.fId] = iinv.ul + iinv.Deun;

		iinv.vf[ug.fId] = iinv.vl + iinv.Devn;

		iinv.wf[ug.fId] = iinv.wl + iinv.Dewn;

		iinv.vnflow = gcom.xfn * (iinv.ul + iinv.Deun * iinv.Bpe) + gcom.yfn * (iinv.vl + iinv.Devn * iinv.Bpe) + gcom.zfn * (iinv.wl + iinv.Dewn * iinv.Bpe) - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;

	}*/
									  
}

void INsInvterm::CmpINsFaceCorrectPresscoef()
{

	Real Sf1 = iinv.Vdvu[ug.fId] * (*ug.a1)[ug.fId];
	Real Sf2 = iinv.Vdvv[ug.fId] * (*ug.a2)[ug.fId];
	Real Sf3 = iinv.Vdvw[ug.fId] * (*ug.a3)[ug.fId];

	Real r2ldx = (*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc];
	Real r2ldy = (*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc];
	Real r2ldz = (*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc];

	Real dist = r2ldx * (*ug.a1)[ug.fId] + r2ldy * (*ug.a2)[ug.fId] + r2ldz * (*ug.a3)[ug.fId];

	Real Sfarea = Sf1 * (*ug.a1)[ug.fId] + Sf2 * (*ug.a2)[ug.fId] + Sf3 * (*ug.a3)[ug.fId];

	iinv.spp[ug.lc] += iinv.rf[ug.fId] * Sfarea / dist;
	iinv.spp[ug.rc] += iinv.rf[ug.fId] * Sfarea / dist;

	iinv.ajp[ug.fId] += iinv.rf[ug.fId] * Sfarea / dist;

	iinv.bp[ug.lc] -= iinv.fq[ug.fId];
	iinv.bp[ug.rc] += iinv.fq[ug.fId];


	/*iinv.Vdvu[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] /(iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc]  / (iinv.spc[ug.rc]));  // -Mf*n，用于求面速度修正量
	iinv.Vdvv[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));
	iinv.Vdvw[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc]  / (iinv.spc[ug.rc]));
	
	
	iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	iinv.ajp[ug.fId] = iinv.rf[ug.fId] * (iinv.Vdvu[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + iinv.Vdvv[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Vdvw[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) * (*ug.farea)[ug.fId] / iinv.dist;*/
}

void INsInvterm::CmpINsBcFaceCorrectPresscoef()
{
	
	int bcType = ug.bcRecord->bcType[ug.fId];

	Real Sf1 = iinv.VdU[ug.lc] * (*ug.a1)[ug.fId];
	Real Sf2 = iinv.VdV[ug.lc] * (*ug.a2)[ug.fId];
	Real Sf3 = iinv.VdW[ug.lc] * (*ug.a3)[ug.fId];

	Real r2ldx = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real r2ldy = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real r2ldz = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dist = r2ldx * (*ug.a1)[ug.fId] + r2ldy * (*ug.a2)[ug.fId] + r2ldz * (*ug.a3)[ug.fId];

	Real Sfarea = Sf1 * (*ug.a1)[ug.fId] + Sf2 * (*ug.a2)[ug.fId] + Sf3 * (*ug.a3)[ug.fId];

	iinv.spp[ug.lc] += iinv.rf[ug.fId] * Sfarea / dist;
	iinv.ajp[ug.fId] = 0;

	iinv.bp[ug.lc] += iinv.fq[ug.fId];

	
	/*int bcType = ug.bcRecord->bcType[ug.fId];

	iinv.Vdvu[ug.fId] = (*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc]);
	iinv.Vdvv[ug.fId] = (*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc]);
	iinv.Vdvw[ug.fId] = (*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc]);

	iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);

	iinv.ajp[ug.fId] = 0;*/

}


  














EndNameSpace