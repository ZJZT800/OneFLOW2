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



		iinv.rf[ug.fId] = (iinv.rl + iinv.rr) * half;    //初始界面上的值（u、v、w ）

		iinv.uf[ug.fId] = (iinv.ul + iinv.ur) * half;

		iinv.vf[ug.fId] = (iinv.vl + iinv.vr) * half;

		iinv.wf[ug.fId] = (iinv.wl + iinv.wr) * half;

		iinv.pf[ug.fId] = (iinv.pl + iinv.pr) * half;

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;  //初始界面上 V*n

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea; //初始界面上的质量通量

}

void INsInvterm::CmpINsBcinvFlux()
{

	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);

	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	int bcType = ug.bcRecord->bcType[ug.fId];

	if(bcType == BC::SOLID_SURFACE)
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
	else if (bcType == BC::OUTFLOW)
	{

		iinv.rf[ug.fId] = iinv.rl;    //初始界面上的值（u、v、w ）

		iinv.uf[ug.fId] = iinv.ul;

		iinv.vf[ug.fId] = iinv.vl;

		iinv.wf[ug.fId] = iinv.wl;

		iinv.pf[ug.fId] = iinv.pl;

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId]* iinv.vnflow * gcom.farea;
	}

	else if(bcType == BC::INFLOW)
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
		
		iinv.ai[ug.fId][0] = clr;
		iinv.ai[ug.fId][1] = crl;

		iinv.biu[ug.fId][0] = 0;
		iinv.biu[ug.fId][1] = 0;

		iinv.biv[ug.fId][0] = 0;
		iinv.biv[ug.fId][1] = 0;

		iinv.biw[ug.fId][0] = 0;
		iinv.biw[ug.fId][1] = 0;

		/*if (ug.fId < ug.nBFace)
		{
			int bcType = ug.bcRecord->bcType[ug.fId];

			if (bcType == BC::INFLOW)
			{
				iinv.Tqu += iinv.rf[ug.fId] * (*ug.xfn)[ug.fId] * iinv.uf[ug.fId] * (*ug.farea)[ug.fId];
				iinv.Tqv += iinv.rf[ug.fId] * (*ug.yfn)[ug.fId] * iinv.vf[ug.fId] * (*ug.farea)[ug.fId];
				iinv.Tqw += iinv.rf[ug.fId] * (*ug.zfn)[ug.fId] * iinv.wf[ug.fId] * (*ug.farea)[ug.fId];

				iinv.Tq += iinv.rf[ug.fId] * ((*ug.xfn)[ug.fId] * iinv.uf[ug.fId]+ (*ug.yfn)[ug.fId] * iinv.vf[ug.fId]+ (*ug.zfn)[ug.fId] * iinv.wf[ug.fId])* (*ug.farea)[ug.fId];
			}
		}*/
		/*iinv.ai[0][ug.fId] = clr;
		iinv.ai[1][ug.fId] = crl;*/
}

void INsInvterm::CmpINsBcinvTerm()
{
	
		Real clr = MAX(0, iinv.fq[ug.fId]);  //从界面左侧单元流入右侧单元的初始质量流量
		Real crl = clr - iinv.fq[ug.fId];   //从界面右侧单元流入左侧单元的初始质量流量


		iinv.ai[ug.fId][0] = clr;
		iinv.ai[ug.fId][1] = crl;

		iinv.biu[ug.fId][0] = crl*iinv.uf[ug.fId];
		iinv.biu[ug.fId][1] = 0;

		iinv.biv[ug.fId][0] = crl * iinv.vf[ug.fId];
		iinv.biv[ug.fId][1] = 0;

		iinv.biw[ug.fId][0] = crl * iinv.wf[ug.fId];
		iinv.biw[ug.fId][1] = 0;
}

void INsInvterm::CmpINsFaceflux()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
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

	iinv.dlf = sqrt(dx1*dx1 + dy1 * dy1 + dz1 * dz1);
	iinv.dfr = sqrt(dx2*dx2 + dy2 * dy2 + dz2 * dz2);

	iinv.Bpe = (*uinsf.dqdx)[IIDX::IIP][ug.lc] * (dx1) + (*uinsf.dqdy)[IIDX::IIP][ug.lc] * (dy1) + (*uinsf.dqdz)[IIDX::IIP][ug.lc] * (dz1) +
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] * (dx2) + (*uinsf.dqdy)[IIDX::IIP][ug.rc] * (dy2) + (*uinsf.dqdz)[IIDX::IIP][ug.rc] * (dz2) -
		(iinv.pr - iinv.pl);

	iinv.rf[ug.fId] = half * (iinv.rl+ iinv.rr);
	iinv.uf[ug.fId] = (iinv.f1[ug.fId] * iinv.ul + iinv.f2[ug.fId] * iinv.ur)+iinv.Deun * iinv.Bpe;  //下一时刻的界面预测速度
	iinv.vf[ug.fId] = (iinv.f1[ug.fId] * iinv.vl + iinv.f2[ug.fId] * iinv.vr)+iinv.Devn * iinv.Bpe;
	iinv.wf[ug.fId] = (iinv.f1[ug.fId] * iinv.wl + iinv.f2[ug.fId] * iinv.wr)+iinv.Dewn * iinv.Bpe;
	

	iinv.vnflow = (*ug.xfn)[ug.fId] * iinv.uf[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vf[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wf[ug.fId];

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * (*ug.farea)[ug.fId];  //下一时刻界面预测通量


	Real clr = MAX(0, iinv.fq[ug.fId]);  //从界面左侧单元流入右侧单元的初始质量流量

	Real crl = clr - iinv.fq[ug.fId];   //从界面右侧单元流入左侧单元的初始质量流量

	iinv.ai[ug.fId][0] = clr;
	iinv.ai[ug.fId][1] = crl;

}


void INsInvterm::CmpINsBcFaceflux()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	iinv.Vau =  (*ug.cvol1)[ug.lc] * gcom.xfn / (iinv.spc[ug.lc]); //Df*n，分子
	iinv.Vav =  (*ug.cvol1)[ug.lc] * gcom.yfn / (iinv.spc[ug.lc]);
	iinv.Vaw =  (*ug.cvol1)[ug.lc] * gcom.zfn / (iinv.spc[ug.lc]);

	iinv.dist = gcom.xfn * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + gcom.yfn * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) + gcom.zfn * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);

	iinv.Deun = iinv.Vau / iinv.dist;   //Df*n/e*n
	iinv.Devn = iinv.Vav / iinv.dist;
	iinv.Dewn = iinv.Vaw / iinv.dist;
	
	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	iinv.dlf = sqrt(dx1*dx1 + dy1 * dy1 + dz1 * dz1);

	iinv.Bpe = (*uinsf.dqdx)[IIDX::IIP][ug.lc] * (dx1)+(*uinsf.dqdy)[IIDX::IIP][ug.lc] * (dy1)+(*uinsf.dqdz)[IIDX::IIP][ug.lc] * (dz1)-
		(iinv.pf[ug.fId] - iinv.pl);

	int bcType = ug.bcRecord->bcType[ug.fId];

	if (bcType == BC::SOLID_SURFACE)
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

	else if (bcType == BC::INFLOW)
	{
		iinv.rf[ug.fId] = inscom.inflow[IIDX::IIR];    //初始界面上的值（u、v、w ）

		iinv.uf[ug.fId] = inscom.inflow[IIDX::IIU];

		iinv.vf[ug.fId] = inscom.inflow[IIDX::IIV];

		iinv.wf[ug.fId] = inscom.inflow[IIDX::IIW];

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea; //初始界面上的质量通量
	}

	else if (bcType == BC::OUTFLOW)
	{

		iinv.rf[ug.fId] = iinv.rl;    //初始界面上的值（u、v、w ）

		iinv.uf[ug.fId] = iinv.ul+iinv.Deun * iinv.Bpe;

		iinv.vf[ug.fId] = iinv.vl+ iinv.Devn * iinv.Bpe;

		iinv.wf[ug.fId] = iinv.wl+ iinv.Dewn * iinv.Bpe;

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;
	}

	



	Real clr = MAX(0, iinv.fq[ug.fId]);  //从界面左侧单元流入右侧单元的初始质量流量

	Real crl = clr - iinv.fq[ug.fId];   //从界面右侧单元流入左侧单元的初始质量流量

	iinv.ai[ug.fId][0] = clr;
	iinv.ai[ug.fId][1] = crl;
									  
}

void INsInvterm::CmpINsFaceCorrectPresscoef()
{

	iinv.Vdvu[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] /(iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc]  / (iinv.spc[ug.rc]));  // -Mf*n，用于求面速度修正量
	iinv.Vdvv[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));
	iinv.Vdvw[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc]  / (iinv.spc[ug.rc]));
	
	
	iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	iinv.ajp[ug.fId] = iinv.rf[ug.fId] * (iinv.Vdvu[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + iinv.Vdvv[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Vdvw[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) * (*ug.farea)[ug.fId] / iinv.dist;
}

void INsInvterm::CmpINsBcFaceCorrectPresscoef()
{
	int bcType = ug.bcRecord->bcType[ug.fId];

	iinv.Vdvu[ug.fId] = (*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc]);
	iinv.Vdvv[ug.fId] = (*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc]);
	iinv.Vdvw[ug.fId] = (*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc]);
	iinv.ajp[ug.fId] = 0;

	/*if (bcType == BC::SOLID_SURFACE)
	{
		iinv.Vdvu[ug.fId] = 0;
		iinv.Vdvv[ug.fId] = 0;
		iinv.Vdvw[ug.fId] = 0;
		iinv.ajp[ug.fId] = 0;
	}
	else if (bcType == BC::INFLOW)
	{
		iinv.Vdvu[ug.fId] = 0;
		iinv.Vdvv[ug.fId] = 0;
		iinv.Vdvw[ug.fId] = 0;
		iinv.ajp[ug.fId] = 0;
	}
	else if (bcType == BC::OUTFLOW)
	{
		iinv.Vdvu[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));  // -Mf*n，用于求面速度修正量
		iinv.Vdvv[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));
		iinv.Vdvw[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));

		iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

		iinv.ajp[ug.fId] = iinv.rf[ug.fId] * (iinv.Vdvu[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + iinv.Vdvv[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Vdvw[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) * (*ug.farea)[ug.fId] / iinv.dist;
	}*/


		/*iinv.Vdvu[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));  // -Mf*n，用于求面速度修正量
		iinv.Vdvv[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));
		iinv.Vdvw[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));

		iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

		iinv.ajp[ug.fId] = iinv.rf[ug.fId] * (iinv.Vdvu[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + iinv.Vdvv[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Vdvw[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) * (*ug.farea)[ug.fId] / iinv.dist;*/

}


  














EndNameSpace