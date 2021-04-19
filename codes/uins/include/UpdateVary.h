#pragma once
#include "HXDefine.h"

BeginNameSpace(ONEFLOW)

class UpateVar
{
public:
	UpateVar(RealField &vary, RealField&diff_value);
	~UpateVar();
};

EndNameSpace