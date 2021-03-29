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

//#include "UINsInvterm.h"
#include "UINsPressCorrect.h"
#include "UINsMomPre.h"
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

UINsPressCorrect::UINsPressCorrect()
{
	;
}

UINsPressCorrect::~UINsPressCorrect()
{
	;
}

//UINsPressCorrect NonZero;

void UINsPressCorrect::InitPressCoeff()
{
	iinv.spp = 0;
	iinv.bp = 0;
}

void UINsPressCorrect::PresEquCoeff()
{
	//this->CmpNewMomCoe();
	InitPressCoeff();

	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.dup[cId] = iinv.spu[cId];
	}
	for (int fId = ug.nBFace; fId < ug.nFace; fId++)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];
		iinv.dup[lc] = iinv.dup[lc] - iinv.ai[0][fId];
		iinv.dup[rc] = iinv.dup[rc] - iinv.ai[1][fId];
	}

	iinv.ai[0] = 0;
	iinv.ai[1] = 0;

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		this->CmpInPressCoeff(fId);
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		this->CmpBcPressCoeff(fId);
	}

	iinv.remax_pp = 0;
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		//iinv.remax_pp = MAX(abs(iinv.remax_pp), abs(iinv.bp[cId]));
		iinv.remax_pp += abs(iinv.bp[cId]);
	}
}

void UINsPressCorrect::CmpInPressCoeff(int& fId)
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

	iinv.spp[lc] += (*uinsf.r)[0][lc] * Sfarea / dist;
	iinv.spp[rc] += (*uinsf.r)[0][lc] * Sfarea / dist;
	iinv.ai[0][fId] = (*uinsf.r)[0][lc] * Sfarea / dist;
	iinv.ai[1][fId] = (*uinsf.r)[0][lc] * Sfarea / dist;

	iinv.bp[lc] -= iinv.flux[fId];
	iinv.bp[rc] += iinv.flux[fId];
}

void UINsPressCorrect::CmpBcPressCoeff(int& fId)
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
		iinv.spp[lc] += (*uinsf.r)[0][lc] * Sfarea / dist;
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
}

void UINsPressCorrect::CmpPressCorrectEqu()
{
	int MaxIter = GetDataValue< int >("EquaPressIter");
	Real Tol = GetDataValue< Real >("EquaPressTol");
	int mRestarts = GetDataValue< int >("EquaPressRestarts");
	string gt = GetDataValue< string >("EquaPressMethod");
	SolveEqua(iinv.spp, iinv.ai, iinv.bp, iinv.pp, iinv.res_p, gt, MaxIter, Tol, mRestarts);
	//std::cout << iinv.res_p << std::endl;

	//边界单元
	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		int bcType = ug.bcRecord->bcType[fId];

		if (bcType == BC::OUTFLOW)
		{
			iinv.ppf[fId] = 0;
		}

		else if (bcType == BC::SOLID_SURFACE)
		{
			iinv.ppf[fId] = iinv.pp[lc];
		}

		else if (bcType == BC::INFLOW)
		{
			iinv.ppf[fId] = iinv.pp[lc];
		}

		else if (ug.bctype == BC::SYMMETRY)
		{
			iinv.ppf[fId] = iinv.pp[lc];
		}
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];

		iinv.ppf[fId] = (*ug.fl)[fId] * iinv.pp[lc] + (*ug.fr)[fId] * iinv.pp[rc];
	}

	Real press_relax = GetDataValue< Real >("press_relax");

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		(*uinsf.p)[0][cId] = (*uinsf.p)[0][cId] + press_relax * (iinv.pp[cId]);
	}


	for (int fId = 0; fId < ug.nBFace; fId++)
	{
		int bcType = ug.bcRecord->bcType[fId];
		int lc = (*ug.lcf)[fId];

		if (bcType == BC::SOLID_SURFACE)
		{
			iinv.pb[fId] = (*uinsf.p)[0][lc];

		}
		else if (bcType == BC::INFLOW)
		{
			iinv.pb[fId] = (*uinsf.p)[0][lc];
		}
		else if (bcType == BC::OUTFLOW)
		{
			iinv.pb[fId] += 0;
		}
		else if (ug.bctype == BC::SYMMETRY)
		{
			iinv.pb[fId] = (*uinsf.p)[0][lc];
		}
	}

	/*for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];
		iinv.pf[fId] = (*ug.fl)[fId] * (*uinsf.q)[IIDX::IIP][lc] + (*ug.fr)[fId] * (*uinsf.q)[IIDX::IIP][rc];
	}

	for (int fId = 0; fId < ug.nBFace; fId++)
	{
		int bcType = ug.bcRecord->bcType[fId];
		int lc = (*ug.lcf)[fId];

		if (bcType == BC::SOLID_SURFACE)
		{
			iinv.pf[fId] = (*uinsf.q)[IIDX::IIP][lc];
			
		}
		else if (bcType == BC::INFLOW)
		{
			iinv.pf[fId] = (*uinsf.q)[IIDX::IIP][lc];
		}
		else if (bcType == BC::OUTFLOW)
		{
			iinv.pf[fId] += 0;
		}
		else if (ug.bctype == BC::SYMMETRY)
		{
			iinv.pf[fId] = (*uinsf.q)[IIDX::IIP][lc];
		}
	}*/

	for (int fId = 0; fId < ug.nBFace; fId++)
	{
		int rc = (*ug.rcf)[fId];
		(*uinsf.p)[0][rc] = iinv.pb[fId];
	}
}

void UINsPressCorrect::UpdateFaceflux()
{

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		CmpUpdateINsFaceflux(fId);
	    
		CmpDun(fId);
	}

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

			CmpUpdateINsBcFaceflux(fId);
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

void UINsPressCorrect::CmpUpdateINsBcFaceflux(int& fId)
{
	int lc = (*ug.lcf)[fId];

	if (ug.bctype == BC::SOLID_SURFACE)
	{
		iinv.flux[fId] = 0;
	}

	else if (ug.bctype == BC::INFLOW)
	{
		iinv.flux[fId] += 0;
	}

	else if (ug.bctype == BC::OUTFLOW)
	{
		Real dupf, dvpf, dwpf;
		dupf = (*ug.cvol)[lc] / iinv.dup[lc];
		dvpf = (*ug.cvol)[lc] / iinv.dup[lc];
		dwpf = (*ug.cvol)[lc] / iinv.dup[lc];
		Real Df1 = dupf * (*ug.a1)[fId];
		Real Df2 = dvpf * (*ug.a2)[fId];
		Real Df3 = dwpf * (*ug.a3)[fId];

		Real l2rdx = (*ug.xfc)[fId] - (*ug.xcc)[lc];
		Real l2rdy = (*ug.yfc)[fId] - (*ug.ycc)[lc];
		Real l2rdz = (*ug.zfc)[fId] - (*ug.zcc)[lc];

		Real Df = Df1 * (*ug.a1)[fId] + Df2 * (*ug.a2)[fId] + Df3 * (*ug.a3)[fId];

		Real dist = l2rdx * (*ug.a1)[fId] + l2rdy * (*ug.a2)[fId] + l2rdz * (*ug.a3)[fId];

		Real fux = (*uinsf.r)[0][lc] * Df / dist * (iinv.pp[lc] - iinv.ppf[fId]);
		iinv.flux[fId] = iinv.flux[fId] + fux;
	}

	else if (ug.bctype == BC::SYMMETRY)
	{
		iinv.flux[fId] = 0;
	}

}

void UINsPressCorrect::CmpUpdateINsFaceflux(int& fId)
{
	int lc = (*ug.lcf)[fId];
	int rc = (*ug.rcf)[fId];

	Real dupf;
	dupf = (*ug.fl)[fId] * ((*ug.cvol)[lc] / iinv.dup[lc]) + (*ug.fr)[fId] * ((*ug.cvol)[rc] / iinv.dup[rc]);
	Real Df1 = dupf * (*ug.a1)[fId];
	Real Df2 = dupf * (*ug.a2)[fId];
	Real Df3 = dupf * (*ug.a3)[fId];

	Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
	Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
	Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];

	Real Df = Df1 * (*ug.a1)[fId] + Df2 * (*ug.a2)[fId] + Df3 * (*ug.a3)[fId];

	Real dist = l2rdx * (*ug.a1)[fId] + l2rdy * (*ug.a2)[fId] + l2rdz * (*ug.a3)[fId];

	Real fux = (*uinsf.r)[0][lc] * Df / dist * (iinv.pp[lc] - iinv.pp[rc]);
	iinv.flux[fId] = iinv.flux[fId] + fux;

}

void UINsPressCorrect::CmpDun(int& fId)
{
	int lc = (*ug.lcf)[fId];
	int rc = (*ug.rcf)[fId];

	//Real ub1, vb1, wb1;

	//if (fId > ug.nBFace - 1)
	//{
		Real uf1 = (*ug.fl)[fId] * (*uinsf.u)[0][lc] + (*ug.fr)[fId] * (*uinsf.u)[0][rc];
		Real vf1 = (*ug.fl)[fId] * (*uinsf.v)[0][lc] + (*ug.fr)[fId] * (*uinsf.v)[0][rc];
		Real wf1 = (*ug.fl)[fId] * (*uinsf.w)[0][lc] + (*ug.fr)[fId] * (*uinsf.w)[0][rc];

		/*iinv.uf[ug.fId] = (*ug.fl)[ug.fId] * (*uinsf.q)[IIDX::IIU][ug.lc] + (*ug.fr)[ug.fId] * (*uinsf.q)[IIDX::IIU][ug.rc];
		iinv.vf[ug.fId] = (*ug.fl)[ug.fId] * (*uinsf.q)[IIDX::IIV][ug.lc] + (*ug.fr)[ug.fId] * (*uinsf.q)[IIDX::IIV][ug.rc];
		iinv.wf[ug.fId] = (*ug.fl)[ug.fId] * (*uinsf.q)[IIDX::IIW][ug.lc] + (*ug.fr)[ug.fId] * (*uinsf.q)[IIDX::IIW][ug.rc];*/
		Real un = uf1 * (*ug.a1)[fId] + vf1 * (*ug.a2)[fId] + wf1 * (*ug.a3)[fId];

		iinv.dun[fId] = iinv.flux[fId] / ((*uinsf.r)[0][lc] + SMALL) - un;
	//}

	/*else 
	{
		if (ug.bctype == BC::INFLOW)
		{
			iinv.uf[ug.fId] = inscom.inflow[1];

			iinv.vf[ug.fId] = inscom.inflow[2];

			iinv.wf[ug.fId] = inscom.inflow[3];

			 ub1 = inscom.inflow[1];

			 vb1 = inscom.inflow[2];

			 wb1 = inscom.inflow[3];


		}

		else if (ug.bctype == BC::SOLID_SURFACE)
		{
			if (ug.bcdtkey == 0)     //静止流动状态时固壁边界面的速度应该为零
			{
				inv.uf[ug.fId] = (*ug.vfx)[ug.fId];

				iinv.vf[ug.fId] = (*ug.vfy)[ug.fId];

				iinv.wf[ug.fId] = (*ug.vfz)[ug.fId];

				 ub1 = (*ug.vfx)[fId];

				 vb1 = (*ug.vfy)[fId];

				 wb1 = (*ug.vfz)[fId];
			}
			else
			{
				iinv.uf[ug.fId] = (*inscom.bcflow)[1];

				iinv.vf[ug.fId] = (*inscom.bcflow)[2];

				iinv.wf[ug.fId] = (*inscom.bcflow)[3];

				 ub1 = (*inscom.bcflow)[1];

				 vb1 = (*inscom.bcflow)[2];

				 wb1 = (*inscom.bcflow)[3];
			}
		}

		else if (ug.bctype == BC::OUTFLOW)
		{
			iinv.uf[ug.fId] = (*uinsf.q)[IIDX::IIU][ug.lc] -(*ug.cvol)[ug.lc] / iinv.dup[ug.lc] * iinv.dppbdx[ug.fId];

			iinv.vf[ug.fId] = (*uinsf.q)[IIDX::IIV][ug.lc] -(*ug.cvol)[ug.lc] / iinv.dup[ug.lc] * iinv.dppbdy[ug.fId];

			iinv.wf[ug.fId] = (*uinsf.q)[IIDX::IIW][ug.lc] -(*ug.cvol)[ug.lc] / iinv.dup[ug.lc] * iinv.dppbdz[ug.fId];

			 ub1 = (*uinsf.q)[IIDX::IIU][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * iinv.dppbdx[fId];

			 vb1 = (*uinsf.q)[IIDX::IIV][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * iinv.dppbdy[fId];

			 wb1 = (*uinsf.q)[IIDX::IIW][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * iinv.dppbdz[fId];
		}

		else if (ug.bctype == BC::SYMMETRY)
		{
			 ub1 = 0;

			 vb1 = 0;

			 wb1 = 0;
		}

		(*uinsf.q)[IIDX::IIU][rc] = ub1;
		(*uinsf.q)[IIDX::IIV][rc] = vb1;
		(*uinsf.q)[IIDX::IIW][rc] = wb1;
	}*/
}

void UINsPressCorrect::UpdateSpeed()
{
	RealField dppdx, dppdy, dppdz;
	dppdx.resize(ug.nCell);
	dppdy.resize(ug.nCell);
	dppdz.resize(ug.nCell);
	ONEFLOW::CmpINsGrad(iinv.ppf, dppdx, dppdy, dppdz);

	/*for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		int lc = (*ug.lcf)[fId];

		iinv.dppbdx[fId] = dqqdx[lc];
		iinv.dppbdy[fId] = dqqdy[lc];
		iinv.dppbdz[fId] = dqqdz[lc];
	}*/

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		(*uinsf.u)[0][cId] -= (*ug.cvol)[cId] / iinv.dup[cId] * dppdx[cId];
		(*uinsf.v)[0][cId] -= (*ug.cvol)[cId] / iinv.dup[cId] * dppdy[cId];
		(*uinsf.w)[0][cId] -= (*ug.cvol)[cId] / iinv.dup[cId] * dppdz[cId];
	}

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
			int rc = (*ug.rcf)[fId];

			if (ug.bcr == -1) return; //interface
			int dd = bcdata.r2d[ug.bcr];
			if (dd != -1)
			{
				ug.bcdtkey = 1;
				inscom.bcflow = &bcdata.dataList[dd];
			}

			if (ug.bctype == BC::INFLOW)
			{
				/*ub1 = inscom.inflow[IIDX::IIU];

				vb1 = inscom.inflow[IIDX::IIV];

				wb1 = inscom.inflow[IIDX::IIW];*/

				;
			}

			else if (ug.bctype == BC::SOLID_SURFACE)
			{
				/*if (ug.bcdtkey == 0)     //静止流动状态时固壁边界面的速度应该为零
				{
					ub1 = (*ug.vfx)[fId];

					vb1 = (*ug.vfy)[fId];

					wb1 = (*ug.vfz)[fId];
				}
				else
				{
					ub1 = (*inscom.bcflow)[IIDX::IIU];

					vb1 = (*inscom.bcflow)[IIDX::IIV];

					wb1 = (*inscom.bcflow)[IIDX::IIW];
				}*/

				;
			}

			else if (ug.bctype == BC::OUTFLOW)
			{

				/*ub1 = (*uinsf.u)[0][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdx[lc];

				vb1 = (*uinsf.v)[0][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdy[lc];

				wb1 = (*uinsf.w)[0][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdz[lc];*/

				iinv.ub[fId] = (*uinsf.u)[0][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdx[lc];

				iinv.vb[fId] = (*uinsf.v)[0][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdy[lc];

				iinv.wb[fId] = (*uinsf.w)[0][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdz[lc];

			}

			else if (ug.bctype == BC::SYMMETRY)
			{
				/*iinv.ub[fId] = 0;

				iinv.vb[fId] = 0;

				iinv.wb[fId] = 0;*/
				;
			}

			(*uinsf.u)[0][rc] = iinv.ub[fId];
			(*uinsf.v)[0][rc] = iinv.vb[fId];
			(*uinsf.w)[0][rc] = iinv.wb[fId];

		}
	}


	/*for (int cId = 0; cId < ug.nCell; ++cId)
	{
		iinv.uu[cId] = iinv.VdU[cId] * dqqdx[cId]; 
		iinv.vv[cId] = iinv.VdV[cId] * dqqdy[cId];
		iinv.ww[cId] = iinv.VdW[cId] * dqqdz[cId];

		(*uinsf.q)[IIDX::IIU][cId] -= iinv.uu[cId];
		(*uinsf.q)[IIDX::IIV][cId] -= iinv.vv[cId];
		(*uinsf.q)[IIDX::IIW][cId] -= iinv.ww[cId];

	}*/
}


EndNameSpace