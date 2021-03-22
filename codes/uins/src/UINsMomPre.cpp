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

#include "UINsMomPre.h"
#include "INsInvterm.h"
//#include "UINsVisterm.h"
//#include "UINsGrad.h"
#include "UGrad.h"
#include "BcData.h"
#include "Zone.h"
#include "Atmosphere.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "UCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Multigrid.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "FieldImp.h"
#include "Iteration.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "Ctrl.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace(ONEFLOW)

UINsMomPre::UINsMomPre()
{
	;
}

UINsMomPre::~UINsMomPre()
{
	;
}


void UINsMomPre::CmpConv()
{
	this->CmpConvTerm();
}


void UINsMomPre::CmpINsPreflux()
{
		iinv.Init();
		ug.Init();
		uinsf.Init();
        //this->Init();

		Real rl, ul, vl, wl, pl;
		Real rr, ur, vr, wr, pr;

		ug.nRegion = ug.bcRecord->bcInfo->bcType.size();
		for (int ir = 0; ir < ug.nRegion; ++ir)
		{
			ug.ir = ir;
			ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
			ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

			for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
			{
				BcInfo* bcInfo = ug.bcRecord->bcInfo;
				int fId = bcInfo->bcFace[ug.ir][ibc];
				ug.bcr = bcInfo->bcRegion[ug.ir][ibc];
				ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];
				int lc = (*ug.lcf)[fId];

				if (ug.bcr == -1) return; //interface
				int dd = bcdata.r2d[ug.bcr];
				if (dd != -1)
				{
					ug.bcdtkey = 1;
					inscom.bcflow = &bcdata.dataList[dd];
				}

				
				INsExtractl(*uinsf.q, rl, ul,vl,wl,pl,lc); 

				if (ug.bctype == BC::SOLID_SURFACE)
				{
					//if (ug.bcdtkey == 0)     
					//{
						/*iinv.rf = iinv.rl;   

						iinv.uf[fId] = (*ug.vfx)[fId];

						iinv.vf[fId] = (*ug.vfy)[fId];

						iinv.wf[fId] = (*ug.vfz)[fId];

						iinv.pf[fId] = iinv.pl;*/

						//iinv.flux[fId] = iinv.rf * ((*ug.a1)[fId] * iinv.uf[fId] + (*ug.a2)[fId] * iinv.vf[fId] + (*ug.a3)[fId] * iinv.wf[fId]);

						iinv.flux[fId] =0;
					//}
					//else
					//{
						/*iinv.rf = iinv.rl;   

						iinv.uf[fId] = (*inscom.bcflow)[1];

						iinv.vf[fId] = (*inscom.bcflow)[2];

						iinv.wf[fId] = (*inscom.bcflow)[3];

						iinv.pf[fId] = iinv.pl;*/

						//iinv.flux[fId] = iinv.rf * ((*ug.a1)[fId] * iinv.uf[fId] + (*ug.a2)[fId] * iinv.vf[fId] + (*ug.a3)[fId] * iinv.wf[fId]);

						//iinv.flux[fId] = rl * ((*ug.a1)[fId] * (*inscom.bcflow)[1] + (*ug.a2)[fId] * (*inscom.bcflow)[2] + (*ug.a3)[fId] * (*inscom.bcflow)[3]);

					//}

				}
				else if (ug.bctype == BC::OUTFLOW)
				{

					/*iinv.rf = iinv.rl;   

					iinv.uf[fId] = iinv.ul;

					iinv.vf[fId] = iinv.vl;

					iinv.wf[fId] = iinv.wl;

					iinv.pf[fId] = iinv.pl;*/

					//iinv.flux[fId] = iinv.rf * ((*ug.a1)[fId] * iinv.uf[fId] + (*ug.a2)[fId] * iinv.vf[fId] + (*ug.a3)[fId] * iinv.wf[fId]);

					iinv.flux[fId] = rl * ((*ug.a1)[fId] * ul + (*ug.a2)[fId] * vl + (*ug.a3)[fId] * wl);
				}

				else if (ug.bctype == BC::INFLOW)
				{
					/*iinv.rf = inscom.inflow[0];   

					iinv.uf[fId] = inscom.inflow[1];

					iinv.vf[fId] = inscom.inflow[2];

					iinv.wf[fId] = inscom.inflow[3];

					iinv.pf[fId] = inscom.inflow[4];*/

					//iinv.flux[fId] = iinv.rf * ((*ug.a1)[fId] * iinv.uf[fId] + (*ug.a2)[fId] * iinv.vf[fId] + (*ug.a3)[fId] * iinv.wf[fId]);

					iinv.flux[fId] = inscom.inflow[0] * ((*ug.a1)[fId] * inscom.inflow[1] + (*ug.a2)[fId] * inscom.inflow[2] + (*ug.a3)[fId] * inscom.inflow[3]);
				}

				else if (ug.bctype == BC::SYMMETRY)
				{
					/*iinv.rf = iinv.rl;

					iinv.uf[fId] = 0;

					iinv.vf[fId] = 0;

					iinv.wf[fId] = 0;

					iinv.pf[fId] = iinv.pl;*/

					//iinv.flux[fId] = iinv.rf * ((*ug.a1)[fId] * iinv.uf[fId] + (*ug.a2)[fId] * iinv.vf[fId] + (*ug.a3)[fId] * iinv.wf[fId]);

					iinv.flux[fId] =0;
				}

			}
		}

		Real uf, vf, wf;

		for (int fId = ug.nBFace; fId < ug.nFace; fId++)
		{
			int lc = (*ug.lcf)[fId];
			int rc = (*ug.rcf)[fId];

			INsExtractl(*uinsf.q, rl, ul, vl, wl, pl,lc);
			INsExtractr(*uinsf.q, rr, ur, vr, wr, pr,rc);

			//iinv.rf = iinv.rl * (*ug.fl)[fId] + iinv.rr * (*ug.fr)[fId];

			uf = ul * (*ug.fl)[fId] + ur * (*ug.fr)[fId];

			vf = vl * (*ug.fl)[fId] + vr * (*ug.fr)[fId];

			wf = wl * (*ug.fl)[fId] + wr * (*ug.fr)[fId];

		   //iinv.pf[fId] = iinv.pl * (*ug.fl)[fId] + iinv.pr * (*ug.fr)[fId];

			iinv.flux[fId] = rl * ((*ug.a1)[fId] * uf + (*ug.a2)[fId] * vf + (*ug.a3)[fId] * wf);
		}


		/*RealField massflux = 0;
		massflux.resize(ug.nCell);
		for (int cId = 0; cId < ug.nCell; cId++)
		{
			int fn = (*ug.c2f)[cId].size();
			for (int iFace = 0; iFace < fn; iFace++)
			{
				int fId = (*ug.c2f)[cId][iFace];
				massflux[cId] += iinv.fq[fId];
			}
			std::cout << "cId: " << cId << ", massflux[cId]: " << massflux[cId] << std::endl;
		}*/
}


void UINsMomPre::InitMomCoeff()
{
	iinv.spu = 0;
	iinv.bu = 0;
	iinv.bv = 0;
	iinv.bw = 0;

	iinv.ai[0] = 0;
	iinv.ai[1] = 0;
}

void UINsMomPre::CmpConvTerm()
{
	InitMomCoeff();

	RealField ub, vb, wb;
	ub.resize(ug.nBFace);
	vb.resize(ug.nBFace);
	wb.resize(ug.nBFace);

	BcVelocity(ub, vb, wb);

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		this->InConvCoff(fId);
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		Real ub1 = ub[fId];
		Real vb1 = vb[fId];
		Real wb1 = wb[fId];

		this->BcConvCoff(ub1, vb1, wb1, fId);
	}
}

void UINsMomPre::InConvCoff(int&fId)
{
	Real clr = MAX(0, iinv.flux[fId]);
	Real crl = clr - iinv.flux[fId];

	iinv.ai[0][fId] = crl;
	iinv.ai[1][fId] = clr;
}

void UINsMomPre::BcConvCoff(Real &ub1, Real &vb1, Real &wb1, int&fId)
{
	int lc = (*ug.lcf)[fId];
	int bctype = ug.bcRecord->bcType[fId];

	Real clr = MAX(0, iinv.flux[fId]);
	Real crl = clr - iinv.flux[fId];

	if (bctype == BC::SOLID_SURFACE || bctype == BC::SYMMETRY)
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

}

void UINsMomPre::BcVelocity(RealField& ub, RealField& vb, RealField& wb)
{
	ub = 0;
	vb = 0;
	wb = 0;

	ug.nRegion = ug.bcRecord->bcInfo->bcType.size();

	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
		{
			BcInfo* bcInfo = ug.bcRecord->bcInfo;
			int fId = bcInfo->bcFace[ug.ir][ibc];
			ug.bcr = bcInfo->bcRegion[ug.ir][ibc];
			ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];
			ug.lc = (*ug.lcf)[fId];
			ug.rc = (*ug.rcf)[fId];

			if (ug.bcr == -1) return; //interface
			int dd = bcdata.r2d[ug.bcr];
			if (dd != -1)
			{
				ug.bcdtkey = 1;
				inscom.bcflow = &bcdata.dataList[dd];
			}

			if (ug.bctype == BC::SOLID_SURFACE)
			{
				if (ug.bcdtkey == 0)
				{

					ub[fId] = (*ug.vfx)[fId];

					vb[fId] = (*ug.vfy)[fId];

					wb[fId] = (*ug.vfz)[fId];
				}
				else
				{
					ub[fId] = (*inscom.bcflow)[1];

					vb[fId] = (*inscom.bcflow)[2];

					wb[fId] = (*inscom.bcflow)[3];
				}

			}
			else if (ug.bctype == BC::OUTFLOW)
			{
				ub[fId] = (*uinsf.q)[IIDX::IIU][ug.lc];

				vb[fId] = (*uinsf.q)[IIDX::IIV][ug.lc];

				wb[fId] = (*uinsf.q)[IIDX::IIW][ug.lc];
			}

			else if (ug.bctype == BC::INFLOW)
			{
				ub[fId] = inscom.inflow[1];

				vb[fId] = inscom.inflow[2];

				wb[fId] = inscom.inflow[3];
			}

			else if (ug.bctype == BC::SYMMETRY)
			{
				ub[fId] = 0;

				vb[fId] = 0;

				wb[fId] = 0;
			}
		}
	}
}

void UINsMomPre::CmpDiffus()
{
	this->CmpDiffusTerm();
}

void UINsMomPre::CmpDiffusTerm()
{
	RealField dudx, dudy, dudz;
	RealField dvdx, dvdy, dvdz;
	RealField dwdx, dwdy, dwdz;

	RealField ub, vb, wb;
	RealField uf, vf, wf;

	dudx.resize(ug.nCell);
	dudy.resize(ug.nCell);
	dudz.resize(ug.nCell);
	dvdx.resize(ug.nCell);
	dvdy.resize(ug.nCell);
	dvdz.resize(ug.nCell);
	dwdx.resize(ug.nCell);
	dwdy.resize(ug.nCell);
	dwdz.resize(ug.nCell);

	ub.resize(ug.nBFace);
	vb.resize(ug.nBFace);
	wb.resize(ug.nBFace);

	uf.resize(ug.nFace);
	vf.resize(ug.nFace);
	wf.resize(ug.nFace);

    BcVelocity(ub, vb, wb);
	FaceVelocity(ub, vb, wb, uf, vf, wf);

	ONEFLOW::CmpINsGrad(uf, dudx, dudy, dudz);
	ONEFLOW::CmpINsGrad(vf, dvdx, dvdy, dvdz);
	ONEFLOW::CmpINsGrad(wf, dwdx, dwdy, dwdz);

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		this->InDiffusCoff(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, fId);
	}

	//Direchlet Boundary Condition
	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		Real ub1 = ub[fId];
		Real vb1 = vb[fId];
		Real wb1 = wb[fId];

		this->BcDiffusCoff(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, ub1, vb1, wb1, fId);
	}

}

void UINsMomPre::InDiffusCoff(RealField& dudx, RealField& dudy, RealField& dudz, RealField& dvdx, RealField& dvdy, RealField& dvdz, RealField& dwdx, RealField& dwdy, RealField& dwdz, int& fId)
{
	int lc = (*ug.lcf)[fId];
	int rc = (*ug.rcf)[fId];

	Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
	Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
	Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];

	Real vis = GetDataValue< Real >("k_vis");

	Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

	Real Fn = (*ug.a1)[fId] * (*ug.a1)[fId] + (*ug.a2)[fId] * (*ug.a2)[fId] + (*ug.a3)[fId] * (*ug.a3)[fId];

	Fn = Fn / dist;

	Real T1 = (*ug.a1)[fId] - l2rdx * Fn;
	Real T2 = (*ug.a2)[fId] - l2rdy * Fn;
	Real T3 = (*ug.a3)[fId] - l2rdz * Fn;

	Real fdudx = (*ug.fl)[fId] * dudx[lc] + (*ug.fr)[fId] * dudx[rc];
	Real fdudy = (*ug.fl)[fId] * dudy[lc] + (*ug.fr)[fId] * dudy[rc];
	Real fdudz = (*ug.fl)[fId] * dudz[lc] + (*ug.fr)[fId] * dudz[rc];
	Real fdvdx = (*ug.fl)[fId] * dvdx[lc] + (*ug.fr)[fId] * dvdx[rc];
	Real fdvdy = (*ug.fl)[fId] * dvdy[lc] + (*ug.fr)[fId] * dvdy[rc];
	Real fdvdz = (*ug.fl)[fId] * dvdz[lc] + (*ug.fr)[fId] * dvdz[rc];
	Real fdwdx = (*ug.fl)[fId] * dwdx[lc] + (*ug.fr)[fId] * dwdx[rc];
	Real fdwdy = (*ug.fl)[fId] * dwdy[lc] + (*ug.fr)[fId] * dwdy[rc];
	Real fdwdz = (*ug.fl)[fId] * dwdz[lc] + (*ug.fr)[fId] * dwdz[rc];

	iinv.ai[0][fId] += vis * Fn;
	iinv.ai[1][fId] += vis * Fn;

	iinv.bu[lc] += vis * (fdudx * T1 + fdudy * T2 + fdudz * T3);
	iinv.bu[rc] -= vis * (fdudx * T1 + fdudy * T2 + fdudz * T3);

	iinv.bv[lc] += vis * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);
	iinv.bv[rc] -= vis * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);

	iinv.bw[lc] += vis * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);
	iinv.bw[rc] -= vis * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);
}

void UINsMomPre::BcDiffusCoff(RealField& dudx, RealField& dudy, RealField& dudz, RealField& dvdx, RealField& dvdy, RealField& dvdz, RealField& dwdx, RealField& dwdy, RealField& dwdz, Real& ub1, Real& vb1, Real& wb1, int& fId)
{
	int lc = (*ug.lcf)[fId];

	Real l2rdx = (*ug.xfc)[fId] - (*ug.xcc)[lc];
	Real l2rdy = (*ug.yfc)[fId] - (*ug.ycc)[lc];
	Real l2rdz = (*ug.zfc)[fId] - (*ug.zcc)[lc];

	Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

	Real Fn = (*ug.a1)[fId] * (*ug.a1)[fId] + (*ug.a2)[fId] * (*ug.a2)[fId] + (*ug.a3)[fId] * (*ug.a3)[fId];

	Fn = Fn / dist;

	Real vis = GetDataValue< Real >("k_vis");

	Real T1 = (*ug.a1)[fId] - l2rdx * Fn;
	Real T2 = (*ug.a2)[fId] - l2rdy * Fn;
	Real T3 = (*ug.a3)[fId] - l2rdz * Fn;

	Real fdudx = dudx[lc];
	Real fdudy = dudy[lc];
	Real fdudz = dudz[lc];
	Real fdvdx = dvdx[lc];
	Real fdvdy = dvdy[lc];
	Real fdvdz = dvdz[lc];
	Real fdwdx = dwdx[lc];
	Real fdwdy = dwdy[lc];
	Real fdwdz = dwdz[lc];

	iinv.spu[lc] += vis * Fn;

	iinv.bu[lc] += vis * Fn * ub1 + vis * (fdudx * T1 + fdudy * T2 + fdudz * T3);
	iinv.bv[lc] += vis * Fn * vb1 + vis * (fdvdx * T1 + fdvdy * T2 + fdvdz * T3);
	iinv.bw[lc] += vis * Fn * wb1 + vis * (fdwdx * T1 + fdwdy * T2 + fdwdz * T3);

}

void UINsMomPre::BcPressure(RealField& pb)
{
	pb = 0;

	ug.nRegion = ug.bcRecord->bcInfo->bcType.size();

	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
		{
			BcInfo* bcInfo = ug.bcRecord->bcInfo;
			int fId = bcInfo->bcFace[ug.ir][ibc];
			ug.lc = (*ug.lcf)[fId];
			ug.rc = (*ug.rcf)[fId];

			if (ug.bctype == BC::SOLID_SURFACE|| ug.bctype == BC::OUTFLOW|| ug.bctype == BC::SYMMETRY)
			{
				pb[fId] = (*uinsf.q)[IIDX::IIP][ug.lc];
			}

			else if (ug.bctype == BC::INFLOW)
			{
				pb[fId] = inscom.inflow[4];
			}
		}
	}
}

void UINsMomPre::FaceVelocity(RealField& ub, RealField& vb, RealField& wb, RealField& uf, RealField& vf, RealField& wf)
{
	uf = 0;
	vf = 0;
	wf = 0;

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		uf[fId] = ub[fId];
		vf[fId] = vb[fId];
		wf[fId] = wb[fId];
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		uf[fId] = (*ug.fl)[fId] * (*uinsf.q)[IIDX::IIU][lc] + (*ug.fr)[fId] * (*uinsf.q)[IIDX::IIU][rc];
		vf[fId] = (*ug.fl)[fId] * (*uinsf.q)[IIDX::IIV][lc] + (*ug.fr)[fId] * (*uinsf.q)[IIDX::IIV][rc];
		wf[fId] = (*ug.fl)[fId] * (*uinsf.q)[IIDX::IIW][lc] + (*ug.fr)[fId] * (*uinsf.q)[IIDX::IIW][rc];
	}

}

void UINsMomPre::FacePressure(RealField& pb, RealField& pf)
{
	pf = 0;

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		pf[fId] = pb[fId];
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		pf[fId] = (*ug.fl)[fId] * (*uinsf.q)[IIDX::IIP][lc] + (*ug.fr)[fId] * (*uinsf.q)[IIDX::IIP][rc];
	}

}

void UINsMomPre::CmpTranst()
{
	Real timestep = GetDataValue< Real >("global_dt");

	for (int cId = 0; cId < ug.nCell; ++cId)
	{

		iinv.spu[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][cId] / timestep;  //矩阵对角线元素的非稳态项

		iinv.bu[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * (*uinsf.q)[IIDX::IIU][cId] / timestep; //源项的非稳态项
		iinv.bv[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * (*uinsf.q)[IIDX::IIV][cId] / timestep;
		iinv.bw[cId] += (*ug.cvol)[cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * (*uinsf.q)[IIDX::IIW][cId] / timestep;
	}

}

void UINsMomPre::CmpSrc()
{
	RealField dpdx, dpdy, dpdz;
	RealField pb, pf;

	dpdx.resize(ug.nCell);
	dpdy.resize(ug.nCell);
	dpdz.resize(ug.nCell);

	pb.resize(ug.nBFace);
	pf.resize(ug.nFace);

	BcPressure(pb);
	FacePressure(pb, pf);

	ONEFLOW::CmpINsGrad(pf, dpdx, dpdy, dpdz);

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		Real vol = (*ug.cvol)[cId];

		iinv.bu[cId] += -vol * dpdx[cId];
		iinv.bv[cId] += -vol * dpdy[cId];
		iinv.bw[cId] += -vol * dpdz[cId];
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		iinv.spu[lc] += iinv.ai[0][fId];
		iinv.spu[rc] += iinv.ai[1][fId];
	}
}

void UINsMomPre::MomEquCoeff()
{
	iinv.remax_up = 0;
	iinv.remax_vp = 0;
	iinv.remax_wp = 0;

	for (int fId = ug.nBFace; fId < ug.nFace; fId++)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		iinv.bu[lc] += iinv.ai[0][fId] * (*uinsf.q)[IIDX::IIU][rc];
		iinv.bv[lc] += iinv.ai[0][fId] * (*uinsf.q)[IIDX::IIV][rc];
		iinv.bw[lc] += iinv.ai[0][fId] * (*uinsf.q)[IIDX::IIW][rc];

		iinv.bu[rc] += iinv.ai[1][fId] * (*uinsf.q)[IIDX::IIU][lc];
		iinv.bv[rc] += iinv.ai[1][fId] * (*uinsf.q)[IIDX::IIV][lc];
		iinv.bw[rc] += iinv.ai[1][fId] * (*uinsf.q)[IIDX::IIW][lc];
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		iinv.bu[cId] -= iinv.spu[cId] * (*uinsf.q)[IIDX::IIU][cId];
		iinv.bv[cId] -= iinv.spu[cId] * (*uinsf.q)[IIDX::IIV][cId];
		iinv.bw[cId] -= iinv.spu[cId] * (*uinsf.q)[IIDX::IIW][cId];

		/*iinv.remax_up = MAX(abs(iinv.remax_up), abs(iinv.buc[cId]));
		iinv.remax_vp = MAX(abs(iinv.remax_vp), abs(iinv.bvc[cId]));
		iinv.remax_wp = MAX(abs(iinv.remax_wp), abs(iinv.bwc[cId]));*/

		iinv.remax_up += abs(iinv.bu[cId]);
		iinv.remax_vp += abs(iinv.bv[cId]);
		iinv.remax_wp += abs(iinv.bw[cId]);
	}
}

void UINsMomPre::RelaxMom(Real a)
{
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.spu[cId] = iinv.spu[cId] * (1 + a);
	}

}

void UINsMomPre::Solveuvw()
{
	RealField uCorrect, vCorrect, wCorrect;
	uCorrect.resize(ug.nCell);
	vCorrect.resize(ug.nCell);
	wCorrect.resize(ug.nCell);
	this->CmpINsMomRes();
	int MaxIter = GetDataValue< int >("EquaMomIter");
	Real Tol = GetDataValue< Real >("EquaMomTol");
	int mRestarts = GetDataValue< int >("EquaMomRestarts");
	string gt = GetDataValue< string >("EquaMomMethod");
	SolveEqua(iinv.spu, iinv.ai, iinv.bu, uCorrect, iinv.res_u, gt, MaxIter, Tol, mRestarts);
	SolveEqua(iinv.spu, iinv.ai, iinv.bv, vCorrect, iinv.res_v, gt, MaxIter, Tol, mRestarts);
	SolveEqua(iinv.spu, iinv.ai, iinv.bw, wCorrect, iinv.res_w, gt, MaxIter, Tol, mRestarts);

	for (int cId = 0; cId < ug.nCell; cId++)
	{
		(*uinsf.q)[IIDX::IIU][cId] += uCorrect[cId];
		(*uinsf.q)[IIDX::IIV][cId] += vCorrect[cId];
		(*uinsf.q)[IIDX::IIW][cId] += wCorrect[cId];
	}
}


void UINsMomPre::CmpFaceflux()
{
	RealField dpdx, dpdy, dpdz;
	RealField pb, pf;

	dpdx.resize(ug.nCell);
	dpdy.resize(ug.nCell);
	dpdz.resize(ug.nCell);

	pb.resize(ug.nBFace);
	pf.resize(ug.nFace);

	BcPressure(pb);
	FacePressure(pb, pf);

	ONEFLOW::CmpINsGrad(pf, dpdx, dpdy, dpdz);

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		Real dpdx1 = dpdx[lc];
		Real dpdx2 = dpdx[rc];
		Real dpdy1 = dpdy[lc];
		Real dpdy2 = dpdy[rc];
		Real dpdz1 = dpdz[lc];
		Real dpdz2 = dpdz[rc];

		CmpINsFaceflux(dpdx1, dpdx2, dpdy1, dpdy2, dpdz1, dpdz2, fId);
	}

	ug.nRegion = ug.bcRecord->bcInfo->bcType.size();
	BcInfo* bcInfo = ug.bcRecord->bcInfo;

	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
		{
			BcInfo* bcInfo = ug.bcRecord->bcInfo;
			int fId = bcInfo->bcFace[ug.ir][ibc];
			ug.bcr = bcInfo->bcRegion[ug.ir][ibc];
			ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

			if (ug.bcr == -1) return; //interface
			int dd = bcdata.r2d[ug.bcr];
			if (dd != -1)
			{
				ug.bcdtkey = 1;
				inscom.bcflow = &bcdata.dataList[dd];
			}

			int lc = (*ug.lcf)[fId];

			Real dpdx1 = dpdx[lc];
			Real dpdy1 = dpdy[lc];
			Real dpdz1 = dpdz[lc];

			Real pb1 = pb[fId];

			CmpINsBcFaceflux(dpdx1, dpdy1, dpdz1, pb1,fId);
		}
	}


	/*RealField massflux = 0;
	massflux.resize(ug.nCell);
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		int fn = (*ug.c2f)[cId].size();
		for (int iFace = 0; iFace < fn; iFace++)
		{
			int fId = (*ug.c2f)[cId][iFace];
			massflux[cId] += iinv.fq[fId];
		}
		std::cout << "cId: " << cId << ", massflux[cId]: " << massflux[cId] << std::endl;
	}*/


}

void UINsMomPre::CmpINsFaceflux(Real & dpdx1, Real & dpdx2, Real & dpdy1, Real & dpdy2, Real & dpdz1, Real & dpdz2, int& fId)
{
	int lc = (*ug.lcf)[fId];
	int rc = (*ug.rcf)[fId];

	Real rl, ul, vl, wl, pl;
	Real rr, ur, vr, wr, pr;

	Real VdU1, VdU2, Vdvu;

	Real uf, vf, wf;

	Real vnflow;

	INsExtractl(*uinsf.q, rl, ul, vl, wl, pl, lc);
	INsExtractr(*uinsf.q, rr, ur, vr, wr, pr, rc);

	Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
	Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
	Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];
	Real rurf = 0.8 / (1 + 0.8);

	VdU1 = (*ug.cvol)[lc] / iinv.spu[lc];
	VdU2 = (*ug.cvol)[rc] / iinv.spu[rc];

	Vdvu = (*ug.fl)[fId] * VdU1 + (*ug.fr)[fId] * VdU2;

	Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;

	Real Df1 = Vdvu * (*ug.a1)[fId] / dist;
	Real Df2 = Vdvu * (*ug.a2)[fId] / dist;
	Real Df3 = Vdvu * (*ug.a3)[fId] / dist;

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

	vnflow = (*ug.a1)[fId] * (uf + fdpdx * Df1) + (*ug.a2)[fId] * (vf + fdpdy * Df2) + (*ug.a3)[fId] * (wf + fdpdz * Df3) + rurf * iinv.dun[fId];
	iinv.flux[fId] = rl * vnflow;

}

void UINsMomPre::CmpINsBcFaceflux(Real& dpdx1, Real& dpdy1, Real& dpdz1, Real& pb1, int& fId)
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

		uf = ul + fdpdx * Df1;
		vf = vl + fdpdy * Df2;
		wf = wl + fdpdz * Df3;

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

}

void UINsMomPre::CmpINsMomRes()
{
	iinv.res_u = 0;
	iinv.res_v = 0;
	iinv.res_w = 0;
}


EndNameSpace