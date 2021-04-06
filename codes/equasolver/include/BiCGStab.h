
#include "BasicComputation.h"

BeginNameSpace(ONEFLOW)

class BiCGStab
{
public:
    BiCGStab(RealField& sp, RealField2D& aii, RealField& rhs, RealField& x, int unknowns, int Iter, Real Tol, bool ifPrecond);

    ~BiCGStab();

private:
    BasicCompute basic;
    RealField r;
    RealField p;
    RealField v;
    RealField t;
    RealField s;
    RealField rHat;
    RealField r_temp;
    RealField apTmp;

private:
    void InitPre(int unknowns);
    void Init(int unknowns);
    void Solve(RealField& sp, RealField2D& aii, RealField& rhs, RealField& x, int unknowns, int Iter, Real Tol);
    void SolvePrecond(RealField& sp, RealField2D& aii, RealField& rhs, RealField& x, int unknowns, int Iter, Real Tol);
    void PrecondILU(RealField& sp, RealField2D& aii, RealField& apTmp, int unknowns);
    void PrecondRU(RealField2D& aii, RealField& apTmp, RealField& b, int unknowns);
};

EndNameSpace