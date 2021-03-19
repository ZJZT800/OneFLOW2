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

#include "UINsInvterm.h"
#include "INsInvterm.h"
#include "UINsVisterm.h"
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

UINsInvterm::UINsInvterm()
{
	;
}

UINsInvterm::~UINsInvterm()
{
	;
}


void UINsInvterm::CmpConv()
{
	this->CmpConvTerm();
}


void UINsInvterm::CmpINsPreflux()
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

/*void UINsInvterm::Init()
{
	//iinv.f1.resize(ug.nFace);
	//iinv.f2.resize(ug.nFace);
	
	//iinv.uf.resize(ug.nFace);
	//iinv.vf.resize(ug.nFace);
	//iinv.wf.resize(ug.nFace);
	//iinv.Vdvu.resize(ug.nFace);
	//iinv.VdU.resize(ug.nCell);
	iinv.bu.resize(ug.nCell);
	iinv.bv.resize(ug.nCell);
	iinv.bw.resize(ug.nCell);
	iinv.bp.resize(ug.nCell);
	iinv.flux.resize(ug.nFace);
	iinv.spu.resize(ug.nCell);
	iinv.ai.resize(ug.nFace, 2);
	//iinv.ajp.resize(ug.nFace, 2);
	iinv.spp.resize(ug.nCell);
	iinv.pp.resize(ug.nCell);
	//iinv.pf.resize(ug.nFace);
	iinv.ppf.resize(ug.nFace);
	iinv.dup.resize(ug.nCell);
	iinv.dun.resize(ug.nFace);
	//iinv.dppbdx.resize(ug.nBFace);
	//iinv.dppbdy.resize(ug.nBFace);
	//iinv.dppbdz.resize(ug.nBFace);

	iinv.bu = 0;
	iinv.bv = 0;
	iinv.bw = 0;
	iinv.bp = 0;
	iinv.pp = 0;
}*/

void UINsInvterm::InitMomCoeff()
{
	iinv.spu = 0;
	iinv.bu = 0;
	iinv.bv = 0;
	iinv.bw = 0;

	iinv.ai[0] = 0;
	iinv.ai[1] = 0;
}

void UINsInvterm::CmpConvTerm()
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

void UINsInvterm::BcVelocity(RealField& ub, RealField& vb, RealField& wb)
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

void UINsInvterm::BcPressure(RealField& pb)
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

void UINsInvterm::FaceVelocity(RealField& ub, RealField& vb, RealField& wb, RealField& uf, RealField& vf, RealField& wf)
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

void UINsInvterm::FacePressure(RealField& pb, RealField& pf)
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

UINsInvterm NonZero;

void UINsInvterm::Solveuvw()
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


void UINsInvterm::CmpFaceflux()
{
	RealField dpdx, dpdy, dpdz;
	RealField pb, pf;

	dpdx.resize(ug.nCell);
	dpdy.resize(ug.nCell);
	dpdz.resize(ug.nCell);

	pb.resize(ug.nBFace);
	pf.resize(ug.nFace);

	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->BcPressure(pb);
	uINsInvterm->FacePressure(pb, pf);

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

void UINsInvterm::CmpINsMomRes()
{
	iinv.res_u = 0;
	iinv.res_v = 0;
	iinv.res_w = 0;
}

void UINsInvterm::AddFlux()
{
	UnsGrid* grid = Zone::GetUnsGrid();
	MRField* res = GetFieldPointer< MRField >(grid, "res");
	int nEqu = res->GetNEqu();
	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];
		//if ( ug.lc == 0 ) cout << fId << endl;

		for (int iEqu = 0; iEqu < nEqu; ++iEqu)
		{
			(*res)[iEqu][ug.lc] -= (*iinvflux)[iEqu][ug.fId];
		}
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//if ( ug.lc == 0 || ug.rc == 0 ) cout << fId << endl;

		for (int iEqu = 0; iEqu < nEqu; ++iEqu)
		{
			(*res)[iEqu][ug.lc] -= (*iinvflux)[iEqu][ug.fId];
			(*res)[iEqu][ug.rc] += (*iinvflux)[iEqu][ug.fId];
		}
	}

	//ONEFLOW::AddF2CField(res, iinvflux);

}

void UINsInvterm::InitPressCoeff()
{
	iinv.spp = 0;
	iinv.bp = 0;
}

void UINsInvterm::PresEquCoeff()
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

void UINsInvterm::maxmin(RealField& a, Real& max_a, Real& min_a)
{
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		max_a = MAX(max_a, a[cId]);
		min_a = MIN(min_a, a[cId]);
	}
}

void UINsInvterm::CmpPressCorrectEqu()
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

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		(*uinsf.q)[IIDX::IIP][cId] = (*uinsf.q)[IIDX::IIP][cId] + 0.4 * (iinv.pp[cId]);
	}

	RealField pb;

	pb.resize(ug.nBFace);

	for (int fId = 0; fId < ug.nBFace; fId++)
	{
		int bcType = ug.bcRecord->bcType[fId];
		int lc = (*ug.lcf)[fId];

		if (bcType == BC::SOLID_SURFACE)
		{
			pb[fId] = (*uinsf.q)[IIDX::IIP][lc];

		}
		else if (bcType == BC::INFLOW)
		{
			pb[fId] = (*uinsf.q)[IIDX::IIP][lc];
		}
		else if (bcType == BC::OUTFLOW)
		{
			pb[fId] = 0;
		}
		else if (ug.bctype == BC::SYMMETRY)
		{
			pb[fId] = (*uinsf.q)[IIDX::IIP][lc];
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
		(*uinsf.q)[IIDX::IIP][rc] = pb[fId];
	}
}


void UINsInvterm::CmpINsPreRes()
{
	//iinv.res_p = 0;


	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		if (ug.cId == 0)
		{
			iinv.res_p = abs(iinv.bp[ug.cId]);
		}
		else
		{
			iinv.res_p = MAX(abs(iinv.bp[ug.cId]), abs(iinv.bp[ug.cId - 1]));
		}

		//iinv.res_p += (iinv.bp[ug.cId]+iinv.mp[ug.cId] - iinv.pp1[ug.cId]* (0.01+iinv.spp[ug.cId]))*(iinv.bp[ug.cId]+iinv.mp[ug.cId] - iinv.pp1[ug.cId]* (0.01+iinv.spp[ug.cId]));
	}

	//iinv.res_p = sqrt(iinv.res_p);
	//iinv.res_p = 0;
}


void UINsInvterm::UpdateFaceflux()
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

void UINsInvterm::CmpUpdateINsBcFaceflux(int& fId)
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

		Real fux = (*uinsf.q)[IIDX::IIR][lc] * Df / dist * (iinv.pp[lc] - iinv.ppf[fId]);
		iinv.flux[fId] = iinv.flux[fId] + fux;
	}

	else if (ug.bctype == BC::SYMMETRY)
	{
		iinv.flux[fId] = 0;
	}

}


void UINsInvterm::CmpUpdateINsFaceflux(int& fId)
{
	int lc = (*ug.lcf)[fId];
	int rc = (*ug.rcf)[fId];

	Real dupf;
	dupf = (*ug.fl)[fId] * ((*ug.cvol)[lc] / iinv.dup[lc]) + (*ug.fr)[fId] * ((*ug.cvol)[lc] / iinv.dup[rc]);
	Real Df1 = dupf * (*ug.a1)[fId];
	Real Df2 = dupf * (*ug.a2)[fId];
	Real Df3 = dupf * (*ug.a3)[fId];

	Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
	Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
	Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];

	Real Df = Df1 * (*ug.a1)[fId] + Df2 * (*ug.a2)[fId] + Df3 * (*ug.a3)[fId];

	Real dist = l2rdx * (*ug.a1)[fId] + l2rdy * (*ug.a2)[fId] + l2rdz * (*ug.a3)[fId];

	Real fux = (*uinsf.q)[IIDX::IIR][lc] * Df / dist * (iinv.pp[lc] - iinv.pp[rc]);
	iinv.flux[fId] = iinv.flux[fId] + fux;

}

void UINsInvterm::CmpDun(int& fId)
{
	int lc = (*ug.lcf)[fId];
	int rc = (*ug.rcf)[fId];

	Real ub1, vb1, wb1;

	//if (fId > ug.nBFace - 1)
	//{
		Real uf1 = (*ug.fl)[fId] * (*uinsf.q)[IIDX::IIU][lc] + (*ug.fr)[fId] * (*uinsf.q)[IIDX::IIU][rc];
		Real vf1 = (*ug.fl)[fId] * (*uinsf.q)[IIDX::IIV][lc] + (*ug.fr)[fId] * (*uinsf.q)[IIDX::IIV][rc];
		Real wf1 = (*ug.fl)[fId] * (*uinsf.q)[IIDX::IIW][lc] + (*ug.fr)[fId] * (*uinsf.q)[IIDX::IIW][rc];

		/*iinv.uf[ug.fId] = (*ug.fl)[ug.fId] * (*uinsf.q)[IIDX::IIU][ug.lc] + (*ug.fr)[ug.fId] * (*uinsf.q)[IIDX::IIU][ug.rc];
		iinv.vf[ug.fId] = (*ug.fl)[ug.fId] * (*uinsf.q)[IIDX::IIV][ug.lc] + (*ug.fr)[ug.fId] * (*uinsf.q)[IIDX::IIV][ug.rc];
		iinv.wf[ug.fId] = (*ug.fl)[ug.fId] * (*uinsf.q)[IIDX::IIW][ug.lc] + (*ug.fr)[ug.fId] * (*uinsf.q)[IIDX::IIW][ug.rc];*/
		Real un = uf1 * (*ug.a1)[fId] + vf1 * (*ug.a2)[fId] + wf1 * (*ug.a3)[fId];

		iinv.dun[ug.fId] = iinv.flux[fId] / ((*uinsf.q)[IIDX::IIR][lc] + SMALL) - un;
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

void UINsInvterm::UpdateSpeed()
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
		(*uinsf.q)[IIDX::IIU][cId] -= (*ug.cvol)[cId] / iinv.dup[cId] * dppdx[cId];
		(*uinsf.q)[IIDX::IIV][cId] -= (*ug.cvol)[cId] / iinv.dup[cId] * dppdy[cId];
		(*uinsf.q)[IIDX::IIW][cId] -= (*ug.cvol)[cId] / iinv.dup[cId] * dppdz[cId];
	}

	Real ub1, vb1, wb1;



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
				ub1 = inscom.inflow[1];

				vb1 = inscom.inflow[2];

				wb1 = inscom.inflow[3];
			}

			else if (ug.bctype == BC::SOLID_SURFACE)
			{
				if (ug.bcdtkey == 0)     //静止流动状态时固壁边界面的速度应该为零
				{
					ub1 = (*ug.vfx)[fId];

					vb1 = (*ug.vfy)[fId];

					wb1 = (*ug.vfz)[fId];
				}
				else
				{
					ub1 = (*inscom.bcflow)[1];

					vb1 = (*inscom.bcflow)[2];

					wb1 = (*inscom.bcflow)[3];
				}
			}

			else if (ug.bctype == BC::OUTFLOW)
			{

				ub1 = (*uinsf.q)[IIDX::IIU][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdx[lc];

				vb1 = (*uinsf.q)[IIDX::IIV][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdy[lc];

				wb1 = (*uinsf.q)[IIDX::IIW][lc] - (*ug.cvol)[lc] / iinv.dup[lc] * dppdz[lc];
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

void UINsInvterm::UpdateINsRes()
{

	std::cout << "iinv.remax_up:" << iinv.remax_up << std::endl;
	std::cout << "iinv.remax_vp:" << iinv.remax_vp << std::endl;
	std::cout << "iinv.remax_wp:" << iinv.remax_wp << std::endl;
	std::cout << "iinv.remax_pp:" << iinv.remax_pp << std::endl;

	ofstream fileres_up("residual_up.txt", ios::app);
	//fileres_p << "residual_p:" <<residual_p << endl;
	fileres_up << iinv.remax_up << endl;
	fileres_up.close();

	ofstream fileres_vp("residual_vp.txt", ios::app);
	//fileres_p << "residual_p:" <<residual_p << endl;
	fileres_vp << iinv.remax_vp << endl;
	fileres_vp.close();

	ofstream fileres_wp("residual_wp.txt", ios::app);
	//fileres_p << "residual_p:" <<residual_p << endl;
	fileres_wp << iinv.remax_wp << endl;
	fileres_wp.close();

	ofstream fileres_pp("residual_pp.txt", ios::app);
	//fileres_p << "residual_p:" <<residual_p << endl;
	fileres_pp << iinv.remax_pp << endl;
	fileres_pp.close();
}


EndNameSpace