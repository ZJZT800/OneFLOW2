
#include "HXMath.h"
#include "DataBase.h"
#include "UCom.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "UINsCom.h"
#include "EquaDiscrete.h"
#include "SIMPLEC.h"
#include "FaceValue.h"
#include "CmpGrad.h"
#include "Direchlet.h"

using namespace std;

BeginNameSpace(ONEFLOW)

EquaDisc::EquaDisc(string &EquaVary, RealField &vary, RealField &vary_b, RealField &vary_old, RealField &flux, RealField &p, RealField &pb, RealField &fdiffus_cof, string &conv_ischeme, string &diffus_ischeme, Real &relax, int &transt,RealField &spu, RealField2D &ai, RealField &bu, Real &resmax)
{
	ConvDiscrete(vary, flux, conv_ischeme, ai);

	DiffusDiscrete(EquaVary, vary, vary_b, fdiffus_cof, diffus_ischeme, ai, bu);

	BcDiscrete(EquaVary, vary, vary_b, flux, fdiffus_cof, conv_ischeme, diffus_ischeme, spu, ai, bu);

	SpecDiscrete(EquaVary, vary, vary_b, p, pb, bu);

	if (transt != 0) TranstDiscrete(vary, vary_old, spu, bu, ai);

	SrcDiscrete(vary, spu, ai, bu, resmax);

	Relax(relax, spu);
}

EquaDisc::~EquaDisc()
{
	;
}

void EquaDisc::ConvDiscrete(RealField &vary, RealField &flux, string &conv_ischeme, RealField2D &ai)
{
	if (conv_ischeme == "FOU")
	{
		for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
		{
			Real clr = MAX(0, flux[fId]);
			Real crl = clr - flux[fId];

			iinv.ai[0][fId] = crl;
			iinv.ai[1][fId] = clr;
		}
	}

	else if (conv_ischeme == "CENTRAL")
	{
		;
	}
	else if (conv_ischeme == "SOU")
	{
		;
	}
}

void EquaDisc::DiffusDiscrete(string &EquaVary,RealField &vary, RealField &vary_b, RealField &fdiffus_cof, string &diffus_ischeme, RealField2D &ai, RealField &bu)
{
	RealField dqdx, dqdy, dqdz;
	RealField qf;

	dqdx.resize(ug.nCell);
	dqdy.resize(ug.nCell);
	dqdz.resize(ug.nCell);

	qf.resize(ug.nFace);

	FaceValue(vary, vary_b, qf);

	CmpUnsGrad(qf, dqdx, dqdy, dqdz);

	if (EquaVary == "vary_h")
	{
		;
	}
	else
	{
		for (int fId = ug.nBFace; fId < ug.nFace; ++fId)
		{
			int lc = (*ug.lcf)[fId];
			int rc = (*ug.rcf)[fId];
			Real l2rdx = (*ug.xcc)[rc] - (*ug.xcc)[lc];
			Real l2rdy = (*ug.ycc)[rc] - (*ug.ycc)[lc];
			Real l2rdz = (*ug.zcc)[rc] - (*ug.zcc)[lc];
			Real dist = (*ug.a1)[fId] * l2rdx + (*ug.a2)[fId] * l2rdy + (*ug.a3)[fId] * l2rdz;
			Real Fn = (*ug.a1)[fId] * (*ug.a1)[fId] + (*ug.a2)[fId] * (*ug.a2)[fId] + (*ug.a3)[fId] * (*ug.a3)[fId];
			Fn = Fn / dist;
			Real T1 = (*ug.a1)[fId] - l2rdx * Fn;
			Real T2 = (*ug.a2)[fId] - l2rdy * Fn;
			Real T3 = (*ug.a3)[fId] - l2rdz * Fn;

			Real fdqdx = (*ug.fl)[fId] * dqdx[lc] + (*ug.fr)[fId] * dqdx[rc];
			Real fdqdy = (*ug.fl)[fId] * dqdy[lc] + (*ug.fr)[fId] * dqdy[rc];
			Real fdqdz = (*ug.fl)[fId] * dqdz[lc] + (*ug.fr)[fId] * dqdz[rc];

			ai[0][fId] += fdiffus_cof[fId] * Fn;
			ai[1][fId] += fdiffus_cof[fId] * Fn;
			bu[lc] += fdiffus_cof[fId] * (fdqdx * T1 + fdqdy * T2 + fdqdz * T3);
			bu[rc] -= fdiffus_cof[fId] * (fdqdx * T1 + fdqdy * T2 + fdqdz * T3);
		}
	}
}

void EquaDisc::BcDiscrete(string&EquaVary,RealField &vary, RealField &vary_b, RealField &flux, RealField &fdiffus_cof, string &conv_ischeme, string &diffus_ischeme, RealField &spu, RealField2D &ai, RealField &bu)
{
	RealField dqdx, dqdy, dqdz;
	RealField qf;

	dqdx.resize(ug.nCell);
	dqdy.resize(ug.nCell);
	dqdz.resize(ug.nCell);

	qf.resize(ug.nFace);

	FaceValue(vary, vary_b, qf);

	CmpUnsGrad(qf, dqdx, dqdy, dqdz);

	for (int ir = 0; ir < ug.nRegion; ++ir)
	{
		ug.ir = ir;
		ug.bctype = ug.bcRecord->bcInfo->bcType[ir];
		ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ir].size();

		for (int ibc = 0; ibc < ug.nRBFace; ++ibc)
		{
			BcInfo* bcInfo = ug.bcRecord->bcInfo;
			int fId = bcInfo->bcFace[ug.ir][ibc];
			int lc = (*ug.lcf)[fId];

			Real clr = MAX(0, flux[fId]);
			Real crl = clr - flux[fId];

			if (ug.bctype == BC::SOLID_SURFACE)
			{
				Real qb = vary_b[fId];

				DirechletBc(EquaVary, dqdx, dqdy, dqdz, qb, fdiffus_cof, diffus_ischeme, fId, spu, bu);
			}
			else if (ug.bctype == BC::INFLOW)
			{
				if (conv_ischeme == "FOU")
				{
					if (iinv.flux[fId] < 0)
					{
						spu[lc] += crl;
						bu[lc] += crl * vary_b[fId];
					}
				}
				else if (conv_ischeme == "CENTRAL")
				{
					;
				}
				else if (conv_ischeme == "SOU")
				{
					;
				}

				Real qb = vary_b[fId];

				DirechletBc(EquaVary, dqdx, dqdy, dqdz, qb, fdiffus_cof, diffus_ischeme,fId, spu, bu);
			}

			else if (ug.bctype == BC::OUTFLOW)
			{
				if (conv_ischeme == "FOU")
				{
					if (iinv.flux[fId] < 0)
					{
						spu[lc] += crl;
						bu[lc] += crl * vary_b[fId];
					}
				}
				else if (conv_ischeme == "CENTRAL")
				{
					;
				}
				else if (conv_ischeme == "SOU")
				{
					;
				}

				Real qb = vary_b[fId];

				DirechletBc(EquaVary, dqdx, dqdy, dqdz, qb, fdiffus_cof, diffus_ischeme, fId, spu, bu);
			}
		}
	}
}

void EquaDisc::SpecDiscrete(string&EquaVary, RealField &vary, RealField &vary_b, RealField &p, RealField &pb, RealField &bu)
{
	RealField dqdx, dqdy, dqdz;
	RealField qf;

	dqdx.resize(ug.nCell);
	dqdy.resize(ug.nCell);
	dqdz.resize(ug.nCell);

	qf.resize(ug.nFace);

	FaceValue(p, pb, qf);

	CmpUnsGrad(qf, dqdx, dqdy, dqdz);

	if (EquaVary == "vary_u")
	{
		for (int cId = 0; cId < ug.nCell; ++cId)
		{
			Real vol = (*ug.cvol)[cId];

			bu[cId] += -vol * dqdx[cId];
		}
	}

	else if(EquaVary == "vary_v")
	{
		for (int cId = 0; cId < ug.nCell; ++cId)
		{
			Real vol = (*ug.cvol)[cId];

			bu[cId] += -vol * dqdy[cId];
		}
	}
	
	else if (EquaVary == "vary_w")
	{
		for (int cId = 0; cId < ug.nCell; ++cId)
		{
			Real vol = (*ug.cvol)[cId];

			bu[cId] += -vol * dqdz[cId];
		}
	}

	else if (EquaVary == "vary_h")
	{
		;
	}

}

void EquaDisc::TranstDiscrete(RealField &vary, RealField &vary_old, RealField &spu, RealField &bu, RealField2D &ai)
{
		Real timestep = GetDataValue< Real >("global_dt");

		for (int cId = 0; cId < ug.nCell; ++cId)
		{
			spu[cId] += (*ug.cvol)[cId] * (*uinsf.rho)[0][cId] / timestep;  
			bu[cId] += (*ug.cvol)[cId] * (*uinsf.rho)[0][cId] * vary_old[cId] / timestep; 
		}
}

void EquaDisc::SrcDiscrete(RealField &vary, RealField &spu, RealField2D &ai, RealField &bu, Real &resmax)
{
	    resmax = 0;

	for (int fId = ug.nBFace; fId < ug.nFace; fId++)
	{
		int lc = (*ug.lcf)[fId];
		int rc = (*ug.rcf)[fId];
		spu[lc] += ai[0][fId];
		spu[rc] += ai[1][fId];
		bu[lc] += ai[0][fId] * vary[rc];
		bu[rc] += ai[1][fId] * vary[lc];
	}
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		bu[cId] -= spu[cId] * vary[cId];
		resmax += pow(bu[cId], 2);
	}

	    resmax = sqrt(resmax);
}

void EquaDisc::Relax(Real&relax, RealField &spu)
{
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		spu[cId] = spu[cId] * (1 + relax);
	}
}


EndNameSpace