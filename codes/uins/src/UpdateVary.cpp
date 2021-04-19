#include "UpdateVary.h"
#include "UCom.h"

using namespace std;

BeginNameSpace(ONEFLOW)

UpateVar::UpateVar(RealField &vary, RealField&diff_value)
{
	for (int cId = 0; cId < ug.nCell; cId++)
	{
		vary[cId] += diff_value[cId];
	}
}

UpateVar::~UpateVar()
{

}

EndNameSpace