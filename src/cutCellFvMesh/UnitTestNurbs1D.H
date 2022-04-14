#ifndef UNITTESTNURBS1D_H
#define UNITTESTNURBS1D_H

#include "Nurbs1D.H"
#include <math.h>
#include <cassert>
#include <random>

namespace Foam
{
class TESTNURBS1D
{
    public:
        TESTNURBS1D()
        {
            testBSplineBasis();
            testControlPointDerivative();
            testCurveDerivative();
            testA();
            Info<<"Nurbs1D Test done"<<endl;
        }
        
    protected:
        bool approxEqual(vector actual, vector target, scalar percErrorMargin = 0.0001);
        bool approxEqual(scalar actual, scalar target, scalar percErrorMargin = 0.0001);
        void testBSplineBasis();
        void testControlPointDerivative();
        void testA();
        void testCurveDerivative();
};
}
#endif
