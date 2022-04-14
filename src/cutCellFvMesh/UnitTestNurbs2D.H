#ifndef UNITTESTNURBS2D_H
#define UNITTESTNURBS2D_H

#include "Nurbs2D.H"
#include "UnitTestNurbs1D.H"
#include <math.h>
#include <cassert>
#include <random>

namespace Foam
{
class TESTNURBS2D : public TESTNURBS1D
{
    public:
        TESTNURBS2D()
        {
            testBSplineBasis();
            testControlPointDerivative();
            testCurveDerivative();
            testA();
            Info<<"Nurbs2D Test done"<<endl;
        }
        
    protected:
        void testControlPointDerivative();
        void testA();
        void testCurveDerivative();
};
}
#endif
