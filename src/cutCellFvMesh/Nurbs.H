#ifndef NURBS_H
#define NURBS_H

#include "fvCFD.H"
#include <memory>

namespace Foam
{   
struct BoundingBox
{
    vector Min={    std::numeric_limits<scalar>::max(),
                    std::numeric_limits<scalar>::max(),
                    std::numeric_limits<scalar>::max()
                };
    vector Max={
                    std::numeric_limits<scalar>::lowest(),
                    std::numeric_limits<scalar>::lowest(),
                    std::numeric_limits<scalar>::lowest()
                };
    bool isInside(vector point) const;
};
    
class Nurbs
{
    private:
        const scalarList knots;
        const List<vector> controlPoints;
        const scalarList weights;
        const int m;
        const int n;
        const int p;
        List<vector> weightedControlPoints;
        const scalar _min_U;
        const scalar _max_U;
        const scalar diameter;
        const scalar deltaX;
        
    public:
        explicit Nurbs
        (
            scalarList knots,
            List<vector> controlPoints,
            scalarList weights,
            int degree,
            scalar diameter=0,
            scalar deltaX=0
        );
        
        scalar getBoundingBoxOverhang() const {return diameter+deltaX;}
        
        vector Curve_Derivative(int k,scalar u) const;
        
        BoundingBox computeBoundingBox(scalar start_u, scalar end_u) const;
        
        BoundingBox computeBoundingBox() const;
        
        scalar min_U() const {return _min_U;}
        
        scalar max_U() const {return std::nexttoward(_max_U,_min_U);}
        
        int nbrKnots() const {return knots.size();}
        
        int degree() const {return p;}
        
        scalar newtonIterateNearestNeighbour(scalar u_0,vector point) const;
        scalar newtonIterateNearestNeighbour_alt(scalar u_0,vector point) const;
        
        scalar distanceToNurbsSurface(scalar para, vector point) const;
                
        void print() const {Info<<controlPoints[0]<<"->"<<controlPoints[1]<<endl;}
        void printBox() const {Info<<computeBoundingBox().Min<<"->"<<computeBoundingBox().Max<<endl;}
                
    private:
        scalar B_Spline_Basis(int i, int p, scalar u) const;
        template <typename T> T Control_Point_Derivative(int k, int i, const List<T>&  controlPoints) const;
        scalar Weights_B_Spline_Derivative(int k, scalar u) const;
        vector A(int k, scalar u) const;
        int binomial(int n, int k) const;
        scalar supremum_Derivative2(scalar start_u, scalar end_u) const;
        vector vectorCurveToPoint_Derivative(int k,scalar para,vector point) const;
        scalar distCurveToPoint_Deriv0(scalar para,vector point) const;
        scalar distCurveToPoint_Deriv1(scalar para,vector point) const;
        scalar distCurveToPoint_Deriv2(scalar para,vector point) const;
        
        friend void UnitTest_BsTree();
};

void UnitTest_Nurbs();
inline scalar euklidianNorm(vector vec) {return sqrt(vec && vec);}
    
}
#endif
