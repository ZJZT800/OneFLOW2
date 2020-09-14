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
#include "UINsGrad.h"
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
#include "UINsLimiter.h"
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
	limiter = new INsLimiter();
	limf = limiter->limf;
}

UINsInvterm::~UINsInvterm()
{
	delete limiter;
}

void UINsInvterm::CmpLimiter()
{
	limiter->CmpLimiter();
}

void UINsInvterm::CmpInvFace()  //单元数据重构
{
	this->CmpLimiter();   //不改

	this->GetQlQrField();  //不改

	this->BoundaryQlQrFixField();  //不改
}

void UINsInvterm::GetQlQrField()
{
	limf->GetQlQr();
}

void UINsInvterm::ReconstructFaceValueField()
{
	limf->CmpFaceValue();
	//limf->CmpFaceValueWeighted();
}

void UINsInvterm::BoundaryQlQrFixField()
{
	limf->BcQlQrFix();
}

void UINsInvterm::CmpInvcoff()
{
	if (inscom.icmpInv == 0) return;
	iinv.Init();
	ug.Init();
	uinsf.Init();
	
	this->CmpInvMassFlux();  //需要改动

}

void UINsInvterm::CmpINsTimestep()
{
	iinv.timestep = GetDataValue< Real >("global_dt");
}
void UINsInvterm::CmpINsPreflux()
{
	if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
	{
		if (inscom.icmpInv == 0) return;
		iinv.Init();
		ug.Init();
		uinsf.Init();

		this->CmpInvFace();
		this->INsPreflux();

	}

}

void UINsInvterm::INsPreflux()
{
	this->Initflux();


	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareFaceValue();

		this->CmpINsinvFlux();

	}

	ug.nRegion = ug.bcRecord->bcInfo->bcType.size();
	BcInfo * bcInfo = ug.bcRecord->bcInfo;

	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
		{
			ug.bcfId = ibc;

			BcInfo * bcInfo = ug.bcRecord->bcInfo;

			ug.fId = bcInfo->bcFace[ug.ir][ibc];
			ug.bcr = bcInfo->bcRegion[ug.ir][ibc];

			ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			inscom.bcdtkey = 0;
			if (ug.bcr == -1) return; //interface
			int dd = bcdata.r2d[ug.bcr];
			if (dd != -1)
			{
				inscom.bcdtkey = 1;
				inscom.bcflow = &bcdata.dataList[dd];
			}

			this->PrepareFaceValue();

			this->CmpINsBcinvFlux();

		}
	}

}
void UINsInvterm::Initflux()
{
	iinv.f1.resize(ug.nFace);
	iinv.f2.resize(ug.nFace);
	iinv.rf.resize(ug.nFace);
	iinv.uf.resize(ug.nFace);
	iinv.vf.resize(ug.nFace);
	iinv.wf.resize(ug.nFace);
	iinv.Vdvu.resize(ug.nFace);
	iinv.Vdvv.resize(ug.nFace);
	iinv.Vdvw.resize(ug.nFace);
	iinv.aju.resize(ug.nFace);
	iinv.ajv.resize(ug.nFace);
	iinv.ajw.resize(ug.nFace);
	iinv.VdU.resize(ug.nTCell);
	iinv.VdV.resize(ug.nTCell);
	iinv.VdW.resize(ug.nTCell);
	iinv.buc.resize(ug.nTCell);
	iinv.bvc.resize(ug.nTCell);
	iinv.bwc.resize(ug.nTCell);
	iinv.bp.resize(ug.nTCell);
	iinv.ajp.resize(ug.nFace);
	iinv.sju.resize(ug.nTCell);
	iinv.sjv.resize(ug.nTCell);
	iinv.sjw.resize(ug.nTCell);
	iinv.fq.resize(ug.nFace);
	iinv.spc.resize(ug.nTCell);
	iinv.ai.resize(ug.nFace,2);
	//iinv.biu.resize(ug.nFace,2);
	//iinv.biv.resize(ug.nFace,2);
	//iinv.biw.resize(ug.nFace,2);
	//iinv.sj.resize(ug.nTCell, 4);
	//iinv.sd.resize(ug.nTCell, 4);
	//iinv.sjp.resize(ug.nCell, ug.nCell);
	//iinv.sjd.resize(ug.nCell, ug.nCell);
	iinv.spp.resize(ug.nTCell);
	iinv.pp.resize(ug.nTCell);
	iinv.uu.resize(ug.nTCell);
	iinv.vv.resize(ug.nTCell);
	iinv.ww.resize(ug.nTCell);
	iinv.uuj.resize(ug.nFace);
	iinv.vvj.resize(ug.nFace);
	iinv.wwj.resize(ug.nFace);
	iinv.muc.resize(ug.nTCell);
	iinv.mvc.resize(ug.nTCell);
	iinv.mwc.resize(ug.nTCell);
	iinv.mp.resize(ug.nTCell);
	iinv.uc.resize(ug.nTCell);
	iinv.vc.resize(ug.nTCell);
	iinv.wc.resize(ug.nTCell);
	iinv.up.resize(ug.nTCell);
	iinv.vp.resize(ug.nTCell);
	iinv.wp.resize(ug.nTCell);
	iinv.spt.resize(ug.nTCell);
	iinv.but.resize(ug.nTCell);
	iinv.bvt.resize(ug.nTCell);
	iinv.bwt.resize(ug.nTCell);
	iinv.dqqdx.resize(ug.nTCell);
	iinv.dqqdy.resize(ug.nTCell);
	iinv.dqqdz.resize(ug.nTCell);
	iinv.Fn.resize(ug.nFace);
	iinv.Fnu.resize(ug.nFace);
	iinv.Fnv.resize(ug.nFace);
	iinv.Fnw.resize(ug.nFace);
	iinv.Fpu.resize(ug.nFace);
	iinv.Fpv.resize(ug.nFace);
	iinv.Fpw.resize(ug.nFace);
	iinv.dsrl.resize(ug.nFace);
	iinv.elrn.resize(ug.nFace);
	//iinv.value.resize(ug.nFace);
	iinv.mu.resize(ug.nCell);
	iinv.mv.resize(ug.nCell);
	iinv.mw.resize(ug.nCell);
	iinv.mua.resize(ug.nCell);
	iinv.mva.resize(ug.nCell);
	iinv.mwa.resize(ug.nCell);
	iinv.res_pp.resize(ug.nCell);
	iinv.res_up.resize(ug.nCell);
	iinv.res_vp.resize(ug.nCell);
	iinv.res_wp.resize(ug.nCell);
	iinv.op.resize(ug.nBFace);
	iinv.dj.resize(ug.nCell);
	iinv.pf.resize(ug.nFace);
	iinv.ppf.resize(ug.nFace);
	iinv.uuf.resize(ug.nBFace);
	iinv.vvf.resize(ug.nBFace);
	iinv.wwf.resize(ug.nBFace);

	iinv.ai1 = 0;
	iinv.ai2 = 0;
	iinv.spu1 = 1;
	iinv.spv1 = 1;
	iinv.spw1 = 1;
	iinv.spu2 = 1;
	iinv.spv2 = 1;
	iinv.spw2 = 1;

	iinv.buc = 0;
	iinv.bvc = 0;
	iinv.bwc = 0;
	iinv.sp1 = 0;
	iinv.sp2 = 0;
	iinv.spj = 0;
	iinv.spp = 0;
	iinv.sppu = 0;
	iinv.sppv = 0;
	iinv.sppw = 0;

	iinv.bpu = 0;
	iinv.bpv = 0;
	iinv.bpw = 0;
	iinv.bp = 0;
	iinv.pp = 0;

	iinv.muc = 0;
	iinv.mvc = 0;
	iinv.mwc = 0;

	iinv.uc = 0;
	iinv.vc = 0;
	iinv.wc = 0;

	iinv.uu = 0;
	iinv.vv = 0;
	iinv.ww = 0;
}

void UINsInvterm::CmpInvMassFlux()
{
	iinv.buc = 0;
	iinv.bvc = 0;
	iinv.bwc = 0;
	iinv.spc = 0;

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->CmpINsinvTerm();
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->CmpINsBcinvTerm();
	}
}

void UINsInvterm::PrepareFaceValue()
{
	gcom.xfn = (*ug.xfn)[ug.fId];
	gcom.yfn = (*ug.yfn)[ug.fId];
	gcom.zfn = (*ug.zfn)[ug.fId];
	gcom.vfn = (*ug.vfn)[ug.fId];
	gcom.farea = (*ug.farea)[ug.fId];

	inscom.gama1 = (*uinsf.gama)[0][ug.lc];
	inscom.gama2 = (*uinsf.gama)[0][ug.rc];

	iinv.gama1 = inscom.gama1;
	iinv.gama2 = inscom.gama2;

	for (int iEqu = 0; iEqu < limf->nEqu; ++iEqu)
	{
		iinv.prim1[iEqu] = (*limf->q)[iEqu][ug.lc];
		iinv.prim2[iEqu] = (*limf->q)[iEqu][ug.rc];
	}
}

void UINsInvterm::PrepareProFaceValue()
{
	gcom.xfn = (*ug.xfn)[ug.fId];
	gcom.yfn = (*ug.yfn)[ug.fId];
	gcom.zfn = (*ug.zfn)[ug.fId];
	gcom.vfn = (*ug.vfn)[ug.fId];
	gcom.farea = (*ug.farea)[ug.fId];

	iinv.prim1[IIDX::IIR] = (*uinsf.q)[IIDX::IIR][ug.lc];
	iinv.prim1[IIDX::IIU] = iinv.uc[ug.lc];
	iinv.prim1[IIDX::IIV] = iinv.vc[ug.lc];
	iinv.prim1[IIDX::IIW] = iinv.wc[ug.lc];
	iinv.prim1[IIDX::IIP] = (*uinsf.q)[IIDX::IIP][ug.lc];

	iinv.prim2[IIDX::IIR] = (*uinsf.q)[IIDX::IIR][ug.rc];
	iinv.prim2[IIDX::IIU] = iinv.uc[ug.rc];
	iinv.prim2[IIDX::IIV] = iinv.vc[ug.rc];
	iinv.prim2[IIDX::IIW] = iinv.wc[ug.rc];
	iinv.prim2[IIDX::IIP] = (*uinsf.q)[IIDX::IIP][ug.rc];

}

UINsInvterm NonZero;
void UINsInvterm::Init()
{
	int Number = 0;
}

void UINsInvterm::MomPre()
{
	this->CmpINsMomRes();

	//BGMRES求解
	NonZero.Number = 0;
	for (int cId = 0; cId < ug.nCell; ++cId)
	{                                          
		int fn = (*ug.c2f)[cId].size();                             //相邻单元的个数                                    
		NonZero.Number += fn;                                          //非对角线上非零元的个数
	}
	NonZero.Number = NonZero.Number + ug.nCell;                     //非零元的总个数         
	Rank.RANKNUMBER = ug.nCell;                                     // 矩阵的行列大小
	Rank.NUMBER = NonZero.Number;                                    // 矩阵非零元素个数传到计算程序中
	Rank.COLNUMBER = 1;                                              //右端项个数
	Rank.Init();                                                     //传入GMRES计算程序的中间变量
	double residual_u, residual_v, residual_w;
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		Rank.TempIA[0] = 0;
		int n = Rank.TempIA[cId];
		int fn = (*ug.c2f)[cId].size();
		Rank.TempIA[cId + 1] = Rank.TempIA[cId] + iinv.dj[cId] + 1;                                                  // 前n+1行非零元素的个数
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];                                                            // 相邻面的编号
			ug.lc = (*ug.lcf)[fId];                                                                     // 面左侧单元
			ug.rc = (*ug.rcf)[fId];                                                                     // 面右侧单元

			if (fId > ug.nBFace - 1)
			{
				if (cId == ug.lc)
				{
					Rank.TempA[n + iFace] = -iinv.ai[fId][0];
					Rank.TempJA[n + iFace] = ug.rc;
				}
				else if (cId == ug.rc)
				{
					Rank.TempA[n + iFace] = -iinv.ai[fId][1];
					Rank.TempJA[n + iFace] = ug.lc;
				}
			}
			else
			{
				continue;
			}
		}

		int fj = iinv.dj[cId];
		Rank.TempA[n + fj] = iinv.spc[cId];                          //主对角线元素值
		Rank.TempJA[n + fj] = cId;                                      //主对角线纵坐标

	}
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		Rank.TempB[cId][0] = iinv.buc[cId];
	}
	bgx.BGMRES();
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.uc[cId] = Rank.TempX[cId][0];                       // 解的输出
	}
	residual_u = Rank.residual;
	iinv.res_u = residual_u;

	Rank.Deallocate();
	//cout << "residual_u:" << residual_u << endl;


	NonZero.Number = 0;
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		int fn = (*ug.c2f)[cId].size();                             //相邻单元的个数                                    
		NonZero.Number += fn;                                          //非对角线上非零元的个数
	}
	NonZero.Number = NonZero.Number + ug.nCell;                     //非零元的总个数         
	Rank.RANKNUMBER = ug.nCell;                                     // 矩阵的行列大小
	Rank.NUMBER = NonZero.Number;                                    // 矩阵非零元素个数传到计算程序中
	Rank.COLNUMBER = 1;                                              //右端项个数
	Rank.Init();                                                     //传入GMRES计算程序的中间变量
	//double residual_u, residual_v, residual_w;
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		Rank.TempIA[0] = 0;
		int n = Rank.TempIA[cId];
		int fn = (*ug.c2f)[cId].size();
		Rank.TempIA[cId + 1] = Rank.TempIA[cId] + iinv.dj[cId] + 1;                                                  // 前n+1行非零元素的个数
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];                                                            // 相邻面的编号
			ug.lc = (*ug.lcf)[fId];                                                                     // 面左侧单元
			ug.rc = (*ug.rcf)[fId];                                                                     // 面右侧单元

			if (fId > ug.nBFace - 1)
			{
				if (cId == ug.lc)
				{
					Rank.TempA[n + iFace] = -iinv.ai[fId][0];
					Rank.TempJA[n + iFace] = ug.rc;
				}
				else if (cId == ug.rc)
				{
					Rank.TempA[n + iFace] = -iinv.ai[fId][1];
					Rank.TempJA[n + iFace] = ug.lc;
				}
			}
			else
			{
				continue;
			}
		}
		int fj = iinv.dj[cId];
		Rank.TempA[n + fj] = iinv.spc[cId];                          //主对角线元素值
		Rank.TempJA[n + fj] = cId;                                      //主对角线纵坐标

	}


	for (int cId = 0; cId < ug.nCell; cId++)
	{
		Rank.TempB[cId][0] = iinv.bvc[cId];
	}	
	bgx.BGMRES();
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.vc[cId] = Rank.TempX[cId][0];
	}
	residual_v = Rank.residual;
	iinv.res_v = residual_v;

	Rank.Deallocate();

	//cout << "residual_v:" << residual_v << endl;


	NonZero.Number = 0;
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		int fn = (*ug.c2f)[cId].size();                             //相邻单元的个数                                    
		NonZero.Number += fn;                                          //非对角线上非零元的个数
	}
	NonZero.Number = NonZero.Number + ug.nCell;                     //非零元的总个数         
	Rank.RANKNUMBER = ug.nCell;                                     // 矩阵的行列大小
	Rank.NUMBER = NonZero.Number;                                    // 矩阵非零元素个数传到计算程序中
	Rank.COLNUMBER = 1;                                              //右端项个数
	Rank.Init();                                                     //传入GMRES计算程序的中间变量
	//double residual_u, residual_v, residual_w;
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		Rank.TempIA[0] = 0;
		int n = Rank.TempIA[cId];
		int fn = (*ug.c2f)[cId].size();
		Rank.TempIA[cId + 1] = Rank.TempIA[cId] + iinv.dj[cId] + 1;                                                  // 前n+1行非零元素的个数
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];                                                            // 相邻面的编号
			ug.lc = (*ug.lcf)[fId];                                                                     // 面左侧单元
			ug.rc = (*ug.rcf)[fId];                                                                     // 面右侧单元

			if (fId > ug.nBFace - 1)
			{
				if (cId == ug.lc)
				{
					Rank.TempA[n + iFace] = -iinv.ai[fId][0];
					Rank.TempJA[n + iFace] = ug.rc;
				}
				else if (cId == ug.rc)
				{
					Rank.TempA[n + iFace] = -iinv.ai[fId][1];
					Rank.TempJA[n + iFace] = ug.lc;
				}
			}
			else
			{
				continue;
			}

		}
		int fj = iinv.dj[cId];
		Rank.TempA[n + fj] = iinv.spc[cId];                          //主对角线元素值
		Rank.TempJA[n + fj] = cId;                                      //主对角线纵坐标

	}

	for (int cId = 0; cId < ug.nCell; cId++)
	{
		Rank.TempB[cId][0] = iinv.bwc[cId];
	}
	bgx.BGMRES();
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		iinv.wc[cId] = Rank.TempX[cId][0];
	}
	residual_w = Rank.residual;
	iinv.res_w = residual_w;

	Rank.Deallocate();

	//cout << "residual_w:" << residual_w << endl;
}

void UINsInvterm::CmpFaceflux()
{

	iinv.Init();
	ug.Init();
	uinsf.Init();

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareProFaceValue();

		this->CmpINsFaceflux();
	}

	ug.nRegion = ug.bcRecord->bcInfo->bcType.size();
	BcInfo * bcInfo = ug.bcRecord->bcInfo;

	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
		{
			ug.bcfId = ibc;

			BcInfo * bcInfo = ug.bcRecord->bcInfo;

			ug.fId = bcInfo->bcFace[ug.ir][ibc];
			ug.bcr = bcInfo->bcRegion[ug.ir][ibc];

			ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			inscom.bcdtkey = 0;
			if (ug.bcr == -1) return; //interface
			int dd = bcdata.r2d[ug.bcr];
			if (dd != -1)
			{
				inscom.bcdtkey = 1;
				inscom.bcflow = &bcdata.dataList[dd];
			}

			this->PrepareProFaceValue();

			this->CmpINsBcFaceflux();
		}
	}

}

void UINsInvterm::CmpINsMomRes()
{
	iinv.res_u = 0;
	iinv.res_v = 0;
	iinv.res_w = 0;

	//判别迭代收敛的条件
	//double phiscale, temp;
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	phiscale = iinv.uc[0];
	//	if (phiscale < iinv.uc[cId])
	//	{
	//		phiscale = iinv.uc[cId];
	//	}
	//}
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	if (iinv.spc[cId] * phiscale - 0.0 > 1e-6)
	//	{
	//		temp = iinv.buc[cId]/(iinv.spc[cId]*phiscale);
	//		iinv.res_u += temp * temp;
	//	}

	//}
	//iinv.res_u = sqrt(iinv.res_u);

	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	phiscale = iinv.vc[0];
	//	if (phiscale < iinv.vc[cId])
	//	{
	//		phiscale = iinv.vc[cId];
	//	}
	//}
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	if (iinv.spc[cId] * phiscale - 0.0 > 1e-6)
	//	{
	//		temp = iinv.bvc[cId] / (iinv.spc[cId] * phiscale);
	//		iinv.res_v += temp * temp;
	//	}

	//}
	//iinv.res_v = sqrt(iinv.res_v);

	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	phiscale = iinv.wc[0];
	//	if (phiscale < iinv.wc[cId])
	//	{
	//		phiscale = iinv.wc[cId];
	//	}
	//}
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	if (iinv.spc[cId] * phiscale - 0.0 > 1e-6)
	//	{
	//		temp = iinv.bwc[cId] / (iinv.spc[cId] * phiscale);
	//		iinv.res_w += temp * temp;
	//	}

	//}
	//iinv.res_w = sqrt(iinv.res_w);


	/*for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.res_u += (iinv.buc[ug.cId]+iinv.muc[ug.cId] - iinv.ump[ug.cId]* (iinv.spu[ug.cId]))*(iinv.buc[ug.cId]+iinv.muc[ug.cId]  - iinv.ump[ug.cId] * (iinv.spu[ug.cId]));
		iinv.res_v += (iinv.bvc[ug.cId]+iinv.mvc[ug.cId] - iinv.vmp[ug.cId] * (iinv.spv[ug.cId]))*(iinv.bvc[ug.cId]+iinv.mvc[ug.cId] - iinv.vmp[ug.cId] * (iinv.spv[ug.cId]));
		iinv.res_w += (iinv.bwc[ug.cId]+iinv.mwc[ug.cId] - iinv.wmp[ug.cId] * (iinv.spw[ug.cId]))*(iinv.bwc[ug.cId]+iinv.mwc[ug.cId] - iinv.wmp[ug.cId] * (iinv.spw[ug.cId]));

	}

	iinv.res_u = sqrt(iinv.res_u);
	iinv.res_v = sqrt(iinv.res_v);
	iinv.res_w = sqrt(iinv.res_w);*/

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

void UINsInvterm::CmpCorrectPresscoef()
{
	//this->CmpNewMomCoe();
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->CmpINsFaceCorrectPresscoef();
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->CmpINsBcFaceCorrectPresscoef();
	}

	iinv.spp = 0;
	iinv.bp = 0;

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.spp[ug.lc] += iinv.ajp[ug.fId];
		iinv.spp[ug.rc] += iinv.ajp[ug.fId];

		iinv.bp[ug.lc] += -iinv.fq[ug.fId];
		iinv.bp[ug.rc] += iinv.fq[ug.fId];
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		int bcType = ug.bcRecord->bcType[ug.fId];

		if (bcType == BC::SOLID_SURFACE)
		{
			;
		}

		else if (bcType == BC::INFLOW)
		{
			iinv.bp[ug.lc] += -iinv.fq[ug.fId];
		}

		else if (bcType == BC::OUTFLOW)
		{
			//iinv.spp[ug.lc] += iinv.ajp[ug.fId];
			iinv.bp[ug.lc] += -iinv.fq[ug.fId];
		}
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.VdU[ug.cId] = -(*ug.cvol)[ug.cId] / ((1+1)*iinv.spc[ug.cId]); //用于求单元修正速度量;
		iinv.VdV[ug.cId] = -(*ug.cvol)[ug.cId] / ((1+1)*iinv.spc[ug.cId]);
		iinv.VdW[ug.cId] = -(*ug.cvol)[ug.cId] / ((1+1)*iinv.spc[ug.cId]);

		int fn = (*ug.c2f)[ug.cId].size();
		iinv.dj[ug.cId] = fn;

		if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
		{
			iinv.sjp.resize(ug.nCell, fn);
			iinv.sjd.resize(ug.nCell, fn);
		}
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			if (fId>ug.nBFace-1)
			{
				if (ug.cId == ug.lc)
				{
					iinv.sjp[ug.cId][iFace] = -iinv.ajp[ug.fId]; //求解压力修正方程的非零系数
					iinv.sjd[ug.cId][iFace] = ug.rc;
				}
				else if (ug.cId == ug.rc)
				{
					iinv.sjp[ug.cId][iFace] = -iinv.ajp[ug.fId];
					iinv.sjd[ug.cId][iFace] = ug.lc;
				}
	
			}
			else
			{
				iinv.dj[ug.cId] -= 1;
			}
		}
	}
}

void UINsInvterm::CmpNewMomCoe()
{
	iinv.spc = 0;

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.ai[ug.fId][1] += iinv.Fn[ug.fId];
		iinv.ai[ug.fId][0] += iinv.Fn[ug.fId];

		iinv.spc[ug.lc] += iinv.ai[ug.fId][1];
		iinv.spc[ug.rc] += iinv.ai[ug.fId][0];
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		int bcType = ug.bcRecord->bcType[ug.fId];

		iinv.ai[ug.fId][1] += iinv.Fn[ug.fId];
		iinv.ai[ug.fId][0] += iinv.Fn[ug.fId];

		if (bcType == BC::SOLID_SURFACE)
		{
			//iinv.spc[ug.lc] += iinv.ai[ug.fId][1];
			;
		}

		/*if(bcType == BC::INFLOW|| bcType == BC::OUTFLOW)
		{
			if (iinv.fq[ug.fId] < 0)
			{
				iinv.spc[ug.lc] += iinv.ai[ug.fId][0];
			}
		}*/

		else if (bcType == BC::INFLOW)
		{
			if (iinv.fq[ug.fId] < 0)
			{
				iinv.spc[ug.lc] += iinv.ai[ug.fId][0];
			}
		}

		else if (bcType == BC::OUTFLOW)
		{
			iinv.spc[ug.lc] += iinv.ai[ug.fId][1];
		}

	}

	
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.spc[ug.cId] += iinv.spt[ug.cId];
	}
}

void UINsInvterm::CmpPressCorrectEqu()
{

	//BGMRES求解
	NonZero.Number = 0;

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		//ug.cId = cId;                                                                  // 主单元编号
		int fn = (*ug.c2f)[cId].size();                                                                 // 单元相邻面的个数
		//NonZero.Number += iinv.dj[cId];
		NonZero.Number += fn;
	}

	NonZero.Number = NonZero.Number + ug.nCell;                                                        // 非零元素的计数
	Rank.RANKNUMBER = ug.nCell;                                                                        // 矩阵的行列
	Rank.COLNUMBER = 1;
	Rank.NUMBER = NonZero.Number;                                                                      // 矩阵非零元素个数
	Rank.Init();
	double residual_p;

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		iinv.ppd = iinv.pp[cId];
		Rank.TempIA[0] = 0;
		int n = Rank.TempIA[cId];
		int fn = (*ug.c2f)[cId].size();
		Rank.TempIA[cId + 1] = Rank.TempIA[cId] + iinv.dj[cId] + 1;                  // 前n+1行非零元素的个数
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];                           // 相邻面的编号
			ug.fId = fId;
			ug.lc = (*ug.lcf)[fId];                                    // 面左侧单元
			ug.rc = (*ug.rcf)[fId];                                    // 面右侧单元

			if (fId > ug.nBFace - 1)
			{
				if (cId == ug.lc)
				{
					Rank.TempA[n + iFace] = iinv.sjp[cId][iFace];          //非对角线元素值
					Rank.TempJA[n + iFace] = ug.rc;                           //非对角线元素纵坐标
				}
				else if (cId == ug.rc)
				{
					Rank.TempA[n + iFace] = iinv.sjp[cId][iFace];          //非对角线元素值
					Rank.TempJA[n + iFace] = ug.lc;                           //非对角线元素纵坐标
				}

			}
			else
			{
				continue;
			}
		}

		int fj = iinv.dj[cId];
		Rank.TempA[n + fj] = iinv.spp[cId];                            //主对角线元素
		Rank.TempJA[n + fj] = cId;                                        //主对角线纵坐标

		Rank.TempB[cId][0] = iinv.bp[cId];                             //右端项
	}
	bgx.BGMRES();
	residual_p = Rank.residual;
	//cout << "residual_p:" << residual_p << endl;
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		//ug.cId = cId;
		iinv.pp[cId] = Rank.TempX[cId][0]; //当前时刻的压力修正值
	}

	Rank.Deallocate();

	//iinv.res_p = 0;
	//iinv.res_p = MAX(iinv.res_p, abs(iinv.ppd - iinv.pp[ug.cId]));

	//边界单元
	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		int bcType = ug.bcRecord->bcType[ug.fId];

		if (bcType == BC::OUTFLOW)
		{
			iinv.ppf[ug.fId] = 0;//Dirichlet
		}

		else if (ug.bctype == BC::SOLID_SURFACE)
		{
			iinv.ppf[ug.fId] = iinv.pp[ug.lc];
		}

		else if (ug.bctype == BC::INFLOW)
		{
			iinv.ppf[ug.fId] = iinv.pp[ug.lc];//Neumann
		}
	}

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real dxl = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dyl = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dzl = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real dxr = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
		Real dyr = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
		Real dzr = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

		Real delt1 = DIST(dxl, dyl, dzl);
		Real delt2 = DIST(dxr, dyr, dzr);
		Real delta = 1.0 / (delt1 + delt2);

		Real cl = delt2 * delta;
		Real cr = delt1 * delta;

		iinv.ppf[ug.fId] = cl * iinv.pp[ug.lc] + cr * iinv.pp[ug.rc];
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		(*uinsf.q)[IIDX::IIP][ug.cId] = (*uinsf.q)[IIDX::IIP][ug.cId] + 0.8*iinv.pp[ug.cId];
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		int bcType = ug.bcRecord->bcType[ug.fId];

		if (bcType == BC::SOLID_SURFACE|| ug.bctype == BC::INFLOW)
		{
			iinv.pf[ug.fId] = (*uinsf.q)[IIDX::IIP][ug.lc];
		}

		else if (bcType == BC::OUTFLOW)
		{
			iinv.pf[ug.fId] += 0;
		}

		//iinv.pf[ug.fId] = iinv.pf[ug.fId] + 0.8*iinv.ppf[ug.fId];
	}


	ug.nRegion = ug.bcRecord->bcInfo->bcType.size();
	BcInfo * bcInfo = ug.bcRecord->bcInfo;

	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		if (ug.bctype < 0)
		{
			false;
		}

		else if (ug.bctype == BC::SOLID_SURFACE)
		{
			for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
			{
				ug.bcfId = ibc;

				BcInfo * bcInfo = ug.bcRecord->bcInfo;

				ug.fId = bcInfo->bcFace[ug.ir][ibc];
				ug.bcr = bcInfo->bcRegion[ug.ir][ibc];

				ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

				ug.lc = (*ug.lcf)[ug.fId];
				ug.rc = (*ug.rcf)[ug.fId];
				inscom.bcdtkey = 0;
				if (ug.bcr == -1) return; //interface
				int dd = bcdata.r2d[ug.bcr];
				if (dd != -1)
				{
					inscom.bcdtkey = 1;
					inscom.bcflow = &bcdata.dataList[dd];
				}

				//(*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.rc]+0.8*iinv.pp[ug.lc];
				(*uinsf.q)[IIDX::IIP][ug.rc] = -(iinv.f1[ug.fId]/ iinv.f2[ug.fId])*(*uinsf.q)[IIDX::IIP][ug.lc] + (1/ iinv.f2[ug.fId]) * iinv.pf[ug.fId];
			}
		}

		else if (ug.bctype == BC::INFLOW)
		{
			for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
			{
				ug.bcfId = ibc;

				BcInfo * bcInfo = ug.bcRecord->bcInfo;

				ug.fId = bcInfo->bcFace[ug.ir][ibc];
				ug.bcr = bcInfo->bcRegion[ug.ir][ibc];

				ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

				ug.lc = (*ug.lcf)[ug.fId];
				ug.rc = (*ug.rcf)[ug.fId];

				inscom.bcdtkey = 0;
				if (ug.bcr == -1) return; //interface
				int dd = bcdata.r2d[ug.bcr];
				if (dd != -1)
				{
					inscom.bcdtkey = 1;
					inscom.bcflow = &bcdata.dataList[dd];
				}

				//(*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.rc] + 0.8*iinv.pp[ug.lc];
				(*uinsf.q)[IIDX::IIP][ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*(*uinsf.q)[IIDX::IIP][ug.lc] + (1 / iinv.f2[ug.fId])  * iinv.pf[ug.fId];
			}
		}
		
		else if (ug.bctype == BC::OUTFLOW)
		{
			for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
			{
				ug.bcfId = ibc;

				BcInfo * bcInfo = ug.bcRecord->bcInfo;

				ug.fId = bcInfo->bcFace[ug.ir][ibc];
				ug.bcr = bcInfo->bcRegion[ug.ir][ibc];

				ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

				ug.lc = (*ug.lcf)[ug.fId];
				ug.rc = (*ug.rcf)[ug.fId];

				inscom.bcdtkey = 0;
				if (ug.bcr == -1) return; //interface
				int dd = bcdata.r2d[ug.bcr];
				if (dd != -1)
				{
					inscom.bcdtkey = 1;
					inscom.bcflow = &bcdata.dataList[dd];
				}

				//(*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.rc]- 0.8*iinv.pp[ug.lc];
				(*uinsf.q)[IIDX::IIP][ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*(*uinsf.q)[IIDX::IIP][ug.lc] + (1 / iinv.f2[ug.fId])  * iinv.pf[ug.fId];
			}
		}

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
	iinv.Init();
	ug.Init();
	uinsf.Init();
	//Alloc();
	//this->CmpInvFace();  //边界处理
	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//this->PrepareFaceValue();

		this->CmpUpdateINsFaceflux();

	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//this->PrepareFaceValue();

		this->CmpUpdateINsBcFaceflux();
	}

}

void UINsInvterm::CmpUpdateINsBcFaceflux()
{
	int bcType = ug.bcRecord->bcType[ug.fId];

	if (bcType == BC::OUTFLOW)
	{
		iinv.fux = iinv.rf[ug.fId] * ((*ug.xfn)[ug.fId] * iinv.uuf[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vvf[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wwf[ug.fId]) * (*ug.farea)[ug.fId];
		iinv.fq[ug.fId] = iinv.fq[ug.fId]+iinv.fux;
	}

	else if (bcType == BC::INFLOW)
	{
		iinv.fux = iinv.rf[ug.fId] * ((*ug.xfn)[ug.fId] * iinv.uuf[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vvf[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wwf[ug.fId]) * (*ug.farea)[ug.fId];
		iinv.fq[ug.fId] = iinv.fq[ug.fId] +iinv.fux;
	}

	else if (bcType == BC::SOLID_SURFACE)
	{
		iinv.fux = iinv.rf[ug.fId] * ((*ug.xfn)[ug.fId] * iinv.uuf[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vvf[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wwf[ug.fId]) * (*ug.farea)[ug.fId];
		iinv.fq[ug.fId] = iinv.fq[ug.fId] +iinv.fux;
	}

	/*Real clr = MAX(0, iinv.fq[ug.fId]);  //从界面左侧单元流入右侧单元的初始质量流量

	Real crl = clr - iinv.fq[ug.fId];   //从界面右侧单元流入左侧单元的初始质量流量

	iinv.ai[ug.fId][0] = crl;
	iinv.ai[ug.fId][1] = clr;*/

}


void UINsInvterm::CmpUpdateINsFaceflux()
{

	iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	iinv.uuj[ug.fId] = iinv.Vdvu[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) * (*ug.xfn)[ug.fId] / iinv.dist; //面速度修正量
	iinv.vvj[ug.fId] = iinv.Vdvv[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) * (*ug.yfn)[ug.fId] / iinv.dist;
	iinv.wwj[ug.fId] = iinv.Vdvw[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) * (*ug.zfn)[ug.fId] / iinv.dist;

	iinv.uf[ug.fId] = iinv.f1[ug.fId]* (*uinsf.q)[IIDX::IIU][ug.lc]+ iinv.f2[ug.fId] * (*uinsf.q)[IIDX::IIU][ug.rc];
	iinv.vf[ug.fId] = iinv.f1[ug.fId]* (*uinsf.q)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.q)[IIDX::IIV][ug.rc];
	iinv.wf[ug.fId] = iinv.f1[ug.fId]* (*uinsf.q)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.q)[IIDX::IIW][ug.rc];

	iinv.fux = iinv.rf[ug.fId] * ((*ug.xfn)[ug.fId] * iinv.uuj[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vvj[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wwj[ug.fId]) * (*ug.farea)[ug.fId];
	iinv.fq[ug.fId] = iinv.fq[ug.fId]+iinv.fux;

	/*Real clr = MAX(0, iinv.fq[ug.fId]);  //从界面左侧单元流入右侧单元的初始质量流量

	Real crl = clr - iinv.fq[ug.fId];   //从界面右侧单元流入左侧单元的初始质量流量

	iinv.ai[ug.fId][0] = crl;
	iinv.ai[ug.fId][1] = clr;*/
}

void UINsInvterm::UpdateSpeed()
{
	ONEFLOW::CmpINsGrad(iinv.ppf, iinv.dqqdx, iinv.dqqdy, iinv.dqqdz);

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.uu[ug.cId] = iinv.VdU[ug.cId] * iinv.dqqdx[ug.cId]*0.8; //速度修正量
		iinv.vv[ug.cId] = iinv.VdV[ug.cId] * iinv.dqqdy[ug.cId]*0.8;
		iinv.ww[ug.cId] = iinv.VdW[ug.cId] * iinv.dqqdz[ug.cId]*0.8;

		iinv.up[ug.cId] = iinv.uc[cId] + iinv.uu[ug.cId];  //下一时刻的速度值
		iinv.vp[ug.cId] = iinv.vc[cId] + iinv.vv[ug.cId];
		iinv.wp[ug.cId] = iinv.wc[cId] + iinv.ww[ug.cId];

		(*uinsf.q)[IIDX::IIU][ug.cId] = iinv.up[ug.cId];
		(*uinsf.q)[IIDX::IIV][ug.cId] = iinv.vp[ug.cId];
		(*uinsf.q)[IIDX::IIW][ug.cId] = iinv.wp[ug.cId];

	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		int bcType = ug.bcRecord->bcType[ug.fId];

		if (bcType == BC::SOLID_SURFACE)
		{
			iinv.uuf[ug.fId] = 0;
			iinv.vvf[ug.fId] = 0;
			iinv.wwf[ug.fId] = 0;

			iinv.uf[ug.fId] = iinv.uf[ug.fId] + iinv.uuf[ug.fId];
			iinv.vf[ug.fId] = iinv.vf[ug.fId] + iinv.vvf[ug.fId];
			iinv.wf[ug.fId] = iinv.wf[ug.fId] + iinv.wwf[ug.fId];

		}

		else if (bcType == BC::INFLOW)
		{
			iinv.uuf[ug.fId] = 0;
			iinv.vvf[ug.fId] = 0;
			iinv.vvf[ug.fId] = 0;

			iinv.uf[ug.fId] = iinv.uf[ug.fId] + iinv.uuf[ug.fId];
			iinv.vf[ug.fId] = iinv.vf[ug.fId] + iinv.vvf[ug.fId];
			iinv.wf[ug.fId] = iinv.wf[ug.fId] + iinv.wwf[ug.fId];
		}

		else if (bcType == BC::OUTFLOW)
		{
			iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);
			
			iinv.uuf[ug.fId] = 0.8*iinv.Vdvu[ug.fId] * (iinv.pp[ug.lc] - iinv.ppf[ug.fId]) * (*ug.xfn)[ug.fId] / iinv.dist;
			iinv.vvf[ug.fId] = 0.8*iinv.Vdvv[ug.fId] * (iinv.pp[ug.lc] - iinv.ppf[ug.fId]) * (*ug.yfn)[ug.fId] / iinv.dist;
			iinv.wwf[ug.fId] = 0.8*iinv.Vdvw[ug.fId] * (iinv.pp[ug.lc] - iinv.ppf[ug.fId]) * (*ug.zfn)[ug.fId] / iinv.dist;

			iinv.uf[ug.fId] = iinv.uf[ug.fId] + iinv.uuf[ug.fId];//iinv.VdU[ug.lc] * iinv.dqqdx[ug.lc] * 0.8;
			iinv.vf[ug.fId] = iinv.vf[ug.fId] + iinv.vvf[ug.fId];//iinv.VdV[ug.lc] * iinv.dqqdy[ug.lc] * 0.8;
			iinv.wf[ug.fId] = iinv.wf[ug.fId] + iinv.wwf[ug.fId];//iinv.VdW[ug.lc] * iinv.dqqdz[ug.lc] * 0.8;


		}
	}


	ug.nRegion = ug.bcRecord->bcInfo->bcType.size();
	BcInfo * bcInfo = ug.bcRecord->bcInfo;

	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		if (ug.bctype < 0)
		{
			false;
		}

		else if (ug.bctype == BC::SOLID_SURFACE)
		{
			for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
			{
				ug.bcfId = ibc;

				BcInfo * bcInfo = ug.bcRecord->bcInfo;

				ug.fId = bcInfo->bcFace[ug.ir][ibc];
				ug.bcr = bcInfo->bcRegion[ug.ir][ibc];

				ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

				ug.lc = (*ug.lcf)[ug.fId];
				ug.rc = (*ug.rcf)[ug.fId];

				inscom.bcdtkey = 0;
				if (ug.bcr == -1) return; //interface
				int dd = bcdata.r2d[ug.bcr];
				if (dd != -1)
				{
					inscom.bcdtkey = 1;
					inscom.bcflow = &bcdata.dataList[dd];
				}
				(*uinsf.q)[IIDX::IIP][ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*(*uinsf.q)[IIDX::IIP][ug.lc] + (1 / iinv.f2[ug.fId]) * iinv.pf[ug.fId];
				if (inscom.bcdtkey == 0)
				{
					iinv.up[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.up[ug.lc] + (1 / iinv.f2[ug.fId])  * gcom.vfx;
					iinv.vp[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.vp[ug.lc] + (1 / iinv.f2[ug.fId])  * gcom.vfy;
					iinv.wp[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.wp[ug.lc] + (1 / iinv.f2[ug.fId])  * gcom.vfz;

					(*uinsf.q)[IIDX::IIU][ug.rc] = iinv.up[ug.rc];
					(*uinsf.q)[IIDX::IIV][ug.rc] = iinv.vp[ug.rc];
					(*uinsf.q)[IIDX::IIW][ug.rc] = iinv.wp[ug.rc];
				}
				else
				{
					iinv.up[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.up[ug.lc] + (1 / iinv.f2[ug.fId]) * (*inscom.bcflow)[IIDX::IIU];
					iinv.vp[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.vp[ug.lc] + (1 / iinv.f2[ug.fId]) * (*inscom.bcflow)[IIDX::IIV];
					iinv.wp[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.wp[ug.lc] + (1 / iinv.f2[ug.fId]) * (*inscom.bcflow)[IIDX::IIW];

					(*uinsf.q)[IIDX::IIU][ug.rc] = iinv.up[ug.rc];
					(*uinsf.q)[IIDX::IIV][ug.rc] = iinv.vp[ug.rc];
					(*uinsf.q)[IIDX::IIW][ug.rc] = iinv.wp[ug.rc];
				}
			}
		}

		else if (ug.bctype == BC::INFLOW)
		{
			for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
			{
				ug.bcfId = ibc;

				BcInfo * bcInfo = ug.bcRecord->bcInfo;

				ug.fId = bcInfo->bcFace[ug.ir][ibc];
				ug.bcr = bcInfo->bcRegion[ug.ir][ibc];

				ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

				ug.lc = (*ug.lcf)[ug.fId];
				ug.rc = (*ug.rcf)[ug.fId];

				inscom.bcdtkey = 0;
				if (ug.bcr == -1) return; //interface
				int dd = bcdata.r2d[ug.bcr];
				if (dd != -1)
				{
					inscom.bcdtkey = 1;
					inscom.bcflow = &bcdata.dataList[dd];
				}
				iinv.up[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.up[ug.lc]+ (1 / iinv.f2[ug.fId]) *inscom.inflow[IIDX::IIU];
				iinv.vp[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.vp[ug.lc]+ (1 / iinv.f2[ug.fId]) *inscom.inflow[IIDX::IIV];
				iinv.wp[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.wp[ug.lc]+ (1 / iinv.f2[ug.fId]) *inscom.inflow[IIDX::IIW];

				(*uinsf.q)[IIDX::IIU][ug.rc] = iinv.up[ug.rc];
				(*uinsf.q)[IIDX::IIV][ug.rc] = iinv.vp[ug.rc];
				(*uinsf.q)[IIDX::IIW][ug.rc] = iinv.wp[ug.rc];
			}
		}

		else if (ug.bctype == BC::OUTFLOW)
		{
			for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
			{
				ug.bcfId = ibc;

				BcInfo * bcInfo = ug.bcRecord->bcInfo;

				ug.fId = bcInfo->bcFace[ug.ir][ibc];
				ug.bcr = bcInfo->bcRegion[ug.ir][ibc];

				ug.bcdtkey = bcInfo->bcdtkey[ug.ir][ibc];

				ug.lc = (*ug.lcf)[ug.fId];
				ug.rc = (*ug.rcf)[ug.fId];

				inscom.bcdtkey = 0;
				if (ug.bcr == -1) return; //interface
				int dd = bcdata.r2d[ug.bcr];
				if (dd != -1)
				{
					inscom.bcdtkey = 1;
					inscom.bcflow = &bcdata.dataList[dd];
				}

				iinv.up[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.up[ug.lc]+ (1 / iinv.f2[ug.fId])*iinv.uf[ug.fId];
				iinv.vp[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.vp[ug.lc]+ (1 / iinv.f2[ug.fId])*iinv.vf[ug.fId];
				iinv.wp[ug.rc] = -(iinv.f1[ug.fId] / iinv.f2[ug.fId])*iinv.wp[ug.lc]+ (1 / iinv.f2[ug.fId])*iinv.wf[ug.fId];

				(*uinsf.q)[IIDX::IIU][ug.rc] = iinv.up[ug.rc];
		        (*uinsf.q)[IIDX::IIV][ug.rc] = iinv.vp[ug.rc];
		        (*uinsf.q)[IIDX::IIW][ug.rc] = iinv.wp[ug.rc];
			}
		}

	}
}

void UINsInvterm::UpdateINsRes()
{
	/*iinv.remax_V = 0;
	iinv.remax_pp = 0;

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.bp[ug.lc] += -iinv.fq[ug.fId];
		iinv.bp[ug.rc] += iinv.fq[ug.fId];
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		iinv.res_V[ug.cId] = 10*iinv.bp[ug.cId];

		iinv.remax_V = MAX(iinv.remax_V, abs(iinv.res_V[ug.cId]));
		iinv.remax_pp = MAX(iinv.remax_pp, abs(iinv.pp[ug.cId]));

	}
	cout << "iinv.remax_V:" << iinv.remax_V << endl;
	cout << "iinv.remax_pp:" << iinv.remax_pp << endl;
	cout <<"innerSteps:"<< Iteration::innerSteps<< endl;
	//cout << "outerSteps:" << Iteration::outerSteps << endl;

	ofstream fileres_vv("residual_vv.txt", ios::app);
	//fileres_p << "residual_p:" <<residual_p << endl;
	fileres_vv << iinv.remax_V << endl;
	fileres_vv.close();
	

	ofstream fileres_pp("residual_pp.txt", ios::app);
	//fileres_p << "residual_p:" <<residual_p << endl;
	fileres_pp << iinv.remax_pp << endl;
	fileres_pp.close();*/


	/*iinv.remax_up = 0;
	iinv.remax_vp = 0;
	iinv.remax_wp = 0;
	iinv.remax_pp = 0;

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.bp[ug.lc] += -iinv.fq[ug.fId];
		iinv.bp[ug.rc] += iinv.fq[ug.fId];
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.res_up[ug.cId] = (iinv.buc[ug.cId]- iinv.spc[ug.cId]*iinv.up[ug.cId])* (iinv.buc[ug.cId] - iinv.spc[ug.cId] * iinv.up[ug.cId]);
		iinv.res_vp[ug.cId] = (iinv.bvc[ug.cId]- iinv.spc[ug.cId] * iinv.vp[ug.cId]) * (iinv.bvc[ug.cId] - iinv.spc[ug.cId] * iinv.vp[ug.cId]);
		iinv.res_wp[ug.cId] = (iinv.bwc[ug.cId] - iinv.spc[ug.cId] * iinv.wp[ug.cId]) * (iinv.bwc[ug.cId] - iinv.spc[ug.cId] * iinv.wp[ug.cId]);

		iinv.res_pp[ug.cId] = iinv.bp[ug.cId] * iinv.bp[ug.cId];

		iinv.remax_up += iinv.res_up[ug.cId];
		iinv.remax_vp += iinv.res_vp[ug.cId];
		iinv.remax_wp += iinv.res_wp[ug.cId];
		iinv.remax_pp += iinv.res_pp[ug.cId];
	}

	iinv.remax_up = sqrt(iinv.remax_up);
	iinv.remax_vp = sqrt(iinv.remax_vp);
	iinv.remax_wp = sqrt(iinv.remax_wp);
	iinv.remax_pp = sqrt(iinv.remax_pp);*/



	iinv.remax_up = 0;
	iinv.remax_vp = 0;
	iinv.remax_wp = 0;
	iinv.remax_pp = 0;

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.bp[ug.lc] += -iinv.fq[ug.fId];
		iinv.bp[ug.rc] += iinv.fq[ug.fId];
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		int fn = (*ug.c2f)[ug.cId].size();

		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			if (fId > ug.nBFace - 1)
			{
				if (ug.cId == ug.lc)
				{
					iinv.mu[ug.cId] += -iinv.ai[ug.fId][0] * (iinv.up[ug.rc] - iinv.uc[ug.rc]);  //矩阵非零系数，动量方程中与主单元相邻的单元面通量
					iinv.mv[ug.cId] += -iinv.ai[ug.fId][0] * (iinv.vp[ug.rc] - iinv.vc[ug.rc]);
					iinv.mw[ug.cId] += -iinv.ai[ug.fId][0] * (iinv.wp[ug.rc] - iinv.wc[ug.rc]);
				}
				else if (ug.cId == ug.rc)
				{
					iinv.mu[ug.cId] += -iinv.ai[ug.fId][1] * (iinv.up[ug.lc] - iinv.uc[ug.lc]);  //矩阵非零系数，动量方程中与主单元相邻的单元面通量
					iinv.mv[ug.cId] += -iinv.ai[ug.fId][1] * (iinv.vp[ug.lc] - iinv.vc[ug.lc]);
					iinv.mw[ug.cId] += -iinv.ai[ug.fId][1] * (iinv.wp[ug.lc] - iinv.wc[ug.lc]);
				}
			}
			else
			{
				iinv.mu[ug.cId] += -iinv.ai[ug.fId][0]*iinv.uuf[ug.fId];
				iinv.mv[ug.cId] += -iinv.ai[ug.fId][0]*iinv.vvf[ug.fId];
				iinv.mu[ug.cId] += -iinv.ai[ug.fId][0]*iinv.wwf[ug.fId];
			}
		}

		iinv.mua[ug.cId] = iinv.spc[ug.cId] * (iinv.up[ug.cId]- iinv.uc[ug.cId]) + iinv.mu[ug.cId];
		iinv.mva[ug.cId] = iinv.spc[ug.cId] * (iinv.vp[ug.cId]- iinv.vc[ug.cId]) + iinv.mv[ug.cId];
		iinv.mwa[ug.cId] = iinv.spc[ug.cId] * (iinv.wp[ug.cId]- iinv.wc[ug.cId]) + iinv.mw[ug.cId];

		iinv.res_up[ug.cId] = (iinv.mua[ug.cId]) * (iinv.mua[ug.cId]);
		iinv.res_vp[ug.cId] = (iinv.mva[ug.cId]) * (iinv.mva[ug.cId]);
		iinv.res_wp[ug.cId] = (iinv.mwa[ug.cId]) * (iinv.mwa[ug.cId]);
		iinv.res_pp[ug.cId] = iinv.bp[ug.cId] * iinv.bp[ug.cId];

		iinv.remax_up += iinv.res_up[ug.cId];
		iinv.remax_vp += iinv.res_vp[ug.cId];
		iinv.remax_wp += iinv.res_wp[ug.cId];
		iinv.remax_pp += iinv.res_pp[ug.cId];
	}

	iinv.remax_up = sqrt(iinv.remax_up);
	iinv.remax_vp = sqrt(iinv.remax_vp);
	iinv.remax_wp = sqrt(iinv.remax_wp);
	iinv.remax_pp = sqrt(iinv.remax_pp);


	cout << "iinv.remax_up:" << iinv.remax_up << endl;
	cout << "iinv.remax_vp:" << iinv.remax_vp << endl;
	cout << "iinv.remax_wp:" << iinv.remax_wp << endl;
	cout << "iinv.remax_pp:" << iinv.remax_pp << endl;


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


void UINsInvterm::CmpPreGrad()
{
	iinv.dqqdx = 0;
	iinv.dqqdy = 0;
	iinv.dqqdz = 0;

	for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real dxl = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dyl = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dzl = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real dxr = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
		Real dyr = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
		Real dzr = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

		Real delt1 = DIST(dxl, dyl, dzl);
		Real delt2 = DIST(dxr, dyr, dzr);
		Real delta = 1.0 / (delt1 + delt2 + SMALL);

		Real cl = delt2 * delta;
		Real cr = delt1 * delta;

		iinv.value = cl * iinv.pp[ug.lc] + cr * iinv.pp[ug.rc];

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		iinv.dqqdx[ug.lc] += fnxa * iinv.value;
		iinv.dqqdy[ug.lc] += fnya * iinv.value;
		iinv.dqqdz[ug.lc] += fnza * iinv.value;

		//if (ug.fId < ug.nBFace) continue;

		iinv.dqqdx[ug.rc] += -fnxa * iinv.value;
		iinv.dqqdy[ug.rc] += -fnya * iinv.value;
		iinv.dqqdz[ug.rc] += -fnza * iinv.value;
	}

	for (int fId =0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.value = iinv.ppf[ug.fId];

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		iinv.dqqdx[ug.lc] += fnxa * iinv.value;
		iinv.dqqdy[ug.lc] += fnya * iinv.value;
		iinv.dqqdz[ug.lc] += fnza * iinv.value;
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		Real ovol = one / (*ug.cvol)[ug.cId];
		iinv.dqqdx[ug.cId] *= ovol;
		iinv.dqqdy[ug.cId] *= ovol;
		iinv.dqqdz[ug.cId] *= ovol;
	}

	/*for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//if (ug.rc > ug.nCell)
		//{
		iinv.dqqdx[ug.rc] = iinv.dqqdx[ug.lc];
		iinv.dqqdy[ug.rc] = iinv.dqqdy[ug.lc];
		iinv.dqqdz[ug.rc] = iinv.dqqdz[ug.lc];
		//}

	}*/

}


void UINsInvterm::Alloc()
{
	uinsf.qf = new MRField(inscom.nEqu, ug.nFace);
}

void UINsInvterm::DeAlloc()
{
	delete uinsf.qf;
}


void UINsInvterm::ReadTmp()
{
	static int iii = 0;
	if (iii) return;
	iii = 1;
	fstream file;
	file.open("nsflow.dat", ios_base::in | ios_base::binary);
	if (!file) exit(0);

	uinsf.Init();

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		for (int iEqu = 0; iEqu < 5; ++iEqu)
		{
			file.read(reinterpret_cast<char*>(&(*uinsf.q)[iEqu][cId]), sizeof(double));
		}
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		file.read(reinterpret_cast<char*>(&(*uinsf.visl)[0][cId]), sizeof(double));
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		file.read(reinterpret_cast<char*>(&(*uinsf.vist)[0][cId]), sizeof(double));
	}

	vector< Real > tmp1(ug.nTCell), tmp2(ug.nTCell);

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		tmp1[cId] = (*uinsf.timestep)[0][cId];
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		file.read(reinterpret_cast<char*>(&(*uinsf.timestep)[0][cId]), sizeof(double));
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		tmp2[cId] = (*uinsf.timestep)[0][cId];
	}

	turbcom.Init();
	uturbf.Init();
	for (int iCell = 0; iCell < ug.nTCell; ++iCell)
	{
		for (int iEqu = 0; iEqu < turbcom.nEqu; ++iEqu)
		{
			file.read(reinterpret_cast<char*>(&(*uturbf.q)[iEqu][iCell]), sizeof(double));
		}
	}
	file.close();
	file.clear();
}



EndNameSpace