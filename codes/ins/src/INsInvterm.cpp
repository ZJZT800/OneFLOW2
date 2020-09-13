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
	prim.resize(nEqu);
	prim1.resize(nEqu);
	prim2.resize(nEqu);

	q.resize(nEqu);
	q1.resize(nEqu);
	q2.resize(nEqu);

	dq.resize(nEqu);

	flux.resize(nEqu);
	flux1.resize(nEqu);
	flux2.resize(nEqu);

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



	iinv.rf[ug.fId] = (iinv.rl + iinv.rr) * half;    //��ʼ�����ϵ�ֵ��u��v��w ��

	iinv.uf[ug.fId] = (iinv.ul + iinv.ur) * half;

	iinv.vf[ug.fId] = (iinv.vl + iinv.vr) * half;

	iinv.wf[ug.fId] = (iinv.wl + iinv.wr) * half;

	iinv.pf[ug.fId] = (iinv.pl + iinv.pr) * half;

	iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;  //��ʼ������ V*n

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea; //��ʼ�����ϵ�����ͨ��

}

void INsInvterm::CmpINsBcinvFlux()
{

	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);

	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	if (ug.bctype == BC::SOLID_SURFACE)
	{
		if (inscom.bcdtkey == 0)     //��ֹ����״̬ʱ�̱ڱ߽�����ٶ�Ӧ��Ϊ��
		{
			iinv.rf[ug.fId] = iinv.rl;    //��ʼ�����ϵ�ֵ��u��v��w ��

			iinv.uf[ug.fId] = (*ug.vfx)[ug.fId];

			iinv.vf[ug.fId] = (*ug.vfy)[ug.fId];

			iinv.wf[ug.fId] = (*ug.vfz)[ug.fId];

			iinv.pf[ug.fId] = iinv.pl;

			iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

			iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;

		}
		else
		{
			iinv.rf[ug.fId] = iinv.rl;    //��ʼ�����ϵ�ֵ��u��v��w ��

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

		iinv.rf[ug.fId] = iinv.rl;    //��ʼ�����ϵ�ֵ��u��v��w ��

		iinv.uf[ug.fId] = iinv.ul;

		iinv.vf[ug.fId] = iinv.vl;

		iinv.wf[ug.fId] = iinv.wl;

		iinv.pf[ug.fId] = iinv.pl;

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;
	}

	else if (ug.bctype == BC::INFLOW)
	{
		iinv.rf[ug.fId] = inscom.inflow[IIDX::IIR];    //��ʼ�����ϵ�ֵ��u��v��w ��

		iinv.uf[ug.fId] = inscom.inflow[IIDX::IIU];

		iinv.vf[ug.fId] = inscom.inflow[IIDX::IIV];

		iinv.wf[ug.fId] = inscom.inflow[IIDX::IIW];

		iinv.pf[ug.fId] = inscom.inflow[IIDX::IIP];

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea; //��ʼ�����ϵ�����ͨ��
	}
}

void INsInvterm::CmpINsinvTerm()
{
	Real clr = MAX(0, iinv.fq[ug.fId]);  //�ӽ�����൥Ԫ�����Ҳ൥Ԫ����������

	Real crl = clr - iinv.fq[ug.fId];   //�ӽ����Ҳ൥Ԫ������൥Ԫ����������

	ug.lc = (*ug.lcf)[ug.fId];
	ug.rc = (*ug.rcf)[ug.fId];
	iinv.ai[ug.fId][0] += crl;
	iinv.ai[ug.fId][1] += clr;

	//diagonal coefficient
	iinv.spc[ug.lc] += clr;
	iinv.spc[ug.rc] += crl;

	iinv.buc[ug.lc] = 0;
	iinv.buc[ug.rc] = 0;

	iinv.bvc[ug.lc] = 0;
	iinv.bvc[ug.rc] = 0;

	iinv.bwc[ug.lc] = 0;
	iinv.bwc[ug.rc] = 0;
}

void INsInvterm::CmpINsBcinvTerm()
{

	Real clr = MAX(0, iinv.fq[ug.fId]);  //�ӽ�����൥Ԫ�����Ҳ൥Ԫ�ĳ�ʼ��������
	Real crl = clr - iinv.fq[ug.fId];   //�ӽ����Ҳ൥Ԫ������൥Ԫ�ĳ�ʼ��������


	//iinv.ai[ug.fId][0] = clr;           //�߽�����Ӧ��ֻ�����Խ���ϵ������Ӱ��
	//iinv.ai[ug.fId][1] = 0;

	iinv.spc[ug.lc] += clr;

	iinv.buc[ug.lc] = crl * iinv.uf[ug.fId];
	iinv.buc[ug.rc] = 0;

	iinv.bvc[ug.lc] = crl * iinv.vf[ug.fId];
	iinv.bvc[ug.rc] = 0;

	iinv.bwc[ug.lc] = crl * iinv.wf[ug.fId];
	iinv.bwc[ug.rc] = 0;
}

void INsInvterm::CmpINsFaceflux()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	// Rie-Chow��ֵ
	iinv.Vau = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * gcom.xfn / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * gcom.xfn / (iinv.spc[ug.rc])); //Df*n������
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

	iinv.Bpe = (*uinsf.dqdx)[IIDX::IIP][ug.lc] * (dx1)+(*uinsf.dqdy)[IIDX::IIP][ug.lc] * (dy1)+(*uinsf.dqdz)[IIDX::IIP][ug.lc] * (dz1)+
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] * (dx2)+(*uinsf.dqdy)[IIDX::IIP][ug.rc] * (dy2)+(*uinsf.dqdz)[IIDX::IIP][ug.rc] * (dz2)-
		(iinv.pr - iinv.pl);      //ѹ����ƽ���ݶ�*d

	iinv.rf[ug.fId] = half * (iinv.rl + iinv.rr);
	iinv.uf[ug.fId] = (iinv.f1[ug.fId] * iinv.ul + iinv.f2[ug.fId] * iinv.ur) + iinv.Deun * iinv.Bpe;  //��һʱ�̵Ľ���Ԥ���ٶ�
	iinv.vf[ug.fId] = (iinv.f1[ug.fId] * iinv.vl + iinv.f2[ug.fId] * iinv.vr) + iinv.Devn * iinv.Bpe;
	iinv.wf[ug.fId] = (iinv.f1[ug.fId] * iinv.wl + iinv.f2[ug.fId] * iinv.wr) + iinv.Dewn * iinv.Bpe;


	iinv.vnflow = (*ug.xfn)[ug.fId] * iinv.uf[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vf[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wf[ug.fId];

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * (*ug.farea)[ug.fId];  //��һʱ�̽���Ԥ��ͨ��


	//Real clr = MAX(0, iinv.fq[ug.fId]);  //�ӽ�����൥Ԫ�����Ҳ൥Ԫ�ĳ�ʼ��������

	//Real crl = clr - iinv.fq[ug.fId];   //�ӽ����Ҳ൥Ԫ������൥Ԫ�ĳ�ʼ��������

	//iinv.spc = 0;

	//iinv.ai[ug.fId][0] = crl;
	//iinv.ai[ug.fId][1] = clr;

	//iinv.spc[ug.lc] -= crl;
	//iinv.spc[ug.rc] -= clr;
}


void INsInvterm::CmpINsBcFaceflux()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	iinv.Vau = ((*ug.cvol1)[ug.lc] * gcom.xfn / (iinv.spc[ug.lc])); //Df*n������
	iinv.Vav = ((*ug.cvol1)[ug.lc] * gcom.yfn / (iinv.spc[ug.lc]));
	iinv.Vaw = ((*ug.cvol1)[ug.lc] * gcom.zfn / (iinv.spc[ug.lc]));

	iinv.dist = gcom.xfn * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + gcom.yfn * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) 
		+ gcom.zfn * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);

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
			iinv.rf[ug.fId] = iinv.rl;    //��ʼ�����ϵ�ֵ��u��v��w ��

			iinv.uf[ug.fId] = (*ug.vfx)[ug.fId];

			iinv.vf[ug.fId] = (*ug.vfy)[ug.fId];

			iinv.wf[ug.fId] = (*ug.vfz)[ug.fId];

			iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

			iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;
		}
		else
		{
			iinv.rf[ug.fId] = iinv.rl;    //��ʼ�����ϵ�ֵ��u��v��w ��

			iinv.uf[ug.fId] = (*inscom.bcflow)[IIDX::IIU];

			iinv.vf[ug.fId] = (*inscom.bcflow)[IIDX::IIV];

			iinv.wf[ug.fId] = (*inscom.bcflow)[IIDX::IIW];

			iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

			iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;
		}

	}

	else if (ug.bctype == BC::INFLOW)
	{
		iinv.rf[ug.fId] = inscom.inflow[IIDX::IIR];    //��ʼ�����ϵ�ֵ��u��v��w ��

		iinv.uf[ug.fId] = inscom.inflow[IIDX::IIU];

		iinv.vf[ug.fId] = inscom.inflow[IIDX::IIV];

		iinv.wf[ug.fId] = inscom.inflow[IIDX::IIW];

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea; //��ʼ�����ϵ�����ͨ��
	}

	else if (ug.bctype == BC::OUTFLOW)
	{

		iinv.rf[ug.fId] = iinv.rl;    //��ʼ�����ϵ�ֵ��u��v��w ��

		//iinv.uf[ug.fId] = iinv.ul + iinv.Deun * iinv.Bpe;

		//iinv.vf[ug.fId] = iinv.vl + iinv.Devn * iinv.Bpe;

		//iinv.wf[ug.fId] = iinv.wl + iinv.Dewn * iinv.Bpe;
		
		iinv.uf[ug.cId] = iinv.ul;

		iinv.vf[ug.cId] = iinv.vl;

		iinv.wf[ug.cId] = iinv.wl;

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea;
	}


	//Real clr = MAX(0, iinv.fq[ug.fId]);  //�ӽ�����൥Ԫ�����Ҳ൥Ԫ�ĳ�ʼ��������

	//Real crl = clr - iinv.fq[ug.fId];   //�ӽ����Ҳ൥Ԫ������൥Ԫ�ĳ�ʼ��������

	//iinv.ai[ug.fId][0] = crl;
	//iinv.ai[ug.fId][1] = crl;
	//iinv.spc[ug.lc] -= crl;

}

void INsInvterm::CmpINsFaceCorrectPresscoef()
{

	iinv.Vdvu[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));  // -Mf*n�����������ٶ�������
	iinv.Vdvv[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));
	iinv.Vdvw[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));

	iinv.dist = ((*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc])) * (*ug.farea)[ug.fId];

	//iinv.ajp[ug.fId] = iinv.rf[ug.fId] * (iinv.Vdvu[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + iinv.Vdvv[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Vdvw[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) * (*ug.farea)[ug.fId] / iinv.dist;

	iinv.spp[ug.lc] -= iinv.rf[ug.fId] * iinv.Vdvu[ug.fId] * (*ug.farea)[ug.fId] * (*ug.farea)[ug.fId] / iinv.dist;
	iinv.spp[ug.rc] -= iinv.rf[ug.fId] * iinv.Vdvu[ug.fId] * (*ug.farea)[ug.fId] * (*ug.farea)[ug.fId] / iinv.dist;

	iinv.ajp[ug.fId] += iinv.rf[ug.fId] * iinv.Vdvu[ug.fId] * (*ug.farea)[ug.fId] * (*ug.farea)[ug.fId] / iinv.dist;

	iinv.bp[ug.lc] += iinv.fq[ug.fId];
	iinv.bp[ug.rc] -= iinv.fq[ug.fId];

}

void INsInvterm::CmpINsBcFaceCorrectPresscoef()
{
	int bcType = ug.bcRecord->bcType[ug.fId];

	iinv.Vdvu[ug.fId] = (*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc]);
	iinv.Vdvv[ug.fId] = (*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc]);
	iinv.Vdvw[ug.fId] = (*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc]);

	iinv.spp[ug.lc] -= iinv.rf[ug.fId] * iinv.Vdvu[ug.fId] * (*ug.farea)[ug.fId] * (*ug.farea)[ug.fId] / iinv.dist;
	iinv.ajp[ug.fId] = 0;

	iinv.bp[ug.lc] += iinv.fq[ug.fId];
}

















EndNameSpace