#pragma once
#include "HXDefine.h"

BeginNameSpace(ONEFLOW)

class EquaDisc
{
public:
	EquaDisc(string &EquaVary, RealField &vary, RealField &vary_b, RealField &vary_old, RealField &flux, RealField &p, RealField &pb, RealField &fdiffus_cof, string &conv_ischeme, string &diffus_ischeme, Real &relax, int &transt,RealField &spu, RealField2D &ai, RealField &bu,Real &resmax);
	~EquaDisc();
public:
	void ConvDiscrete(RealField &vary, RealField &flux, string &conv_ischeme, RealField2D &ai);
	void DiffusDiscrete(string &EquaVary,RealField &vary, RealField &vary_b, RealField &fdiffus_cof, string &diffus_ischeme,RealField2D &ai, RealField &bu);
	void BcDiscrete(string&EquaVary,RealField &vary, RealField &vary_b, RealField &flux, RealField &fdiffus_cof, string &conv_ischeme, string &diffus_ischeme, RealField &spu, RealField2D &ai, RealField &bu);
	void SpecDiscrete(string&EquaVary, RealField &vary, RealField &vary_b, RealField &p, RealField &pb, RealField &bu);
	void TranstDiscrete(RealField &vary, RealField &vary_old,RealField &spu, RealField &bu, RealField2D &ai);
	void SrcDiscrete(RealField &vary, RealField &spu, RealField2D &ai, RealField &bu, Real &resmax);
	void Relax(Real&relax, RealField &spu);
};

EndNameSpace