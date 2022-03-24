#ifndef NURBS2D_H
#define NURBS2D_H

#include "fvCFD.H"
#include <memory>
#include "Nurbs.H"

namespace Foam
{    
class Nurbs2D : public Nurbs1D
{
    private:      

    public:
        Nurbs2D
        (
            List<scalarList> knots,
            List<List<vector>> controlPoints,
            scalarListList weights,
            int degree_p,
            int degree_q
            scalar diameter=0.2,
            scalar deltaX=0.2
        );
    
        
        //scalar getBoundingBoxOverhang() const {return diameter+deltaX;}
        
        vector Curve_Derivative(int k,scalar u,nurbsStatus state=curr) const = delete;        
        vector Surface_Derivative(int k,int l, scalar u, scalar v, nurbsStatus state=curr) const;
        
        BoundingBox computeBoundingBox(scalar start_u, scalar end_u) const = delete;
        BoundingBox computeBoundingBox(scalar start_u, scalar start_v, scalar end_u, scalar end_v) const;
                
        scalar min_V() const {return minPara[v];}    
        scalar max_V() const const {return std::nexttoward(maxPara[v],minPara[v]);}
        
        int nbrKnotsV() const {return knots[v].size();}
        
        int degreeV() const {return p_q[v];}
        
        scalar newtonIterateNearestNeighbour(scalar u_0,vector point) const;
        scalar newtonIterateNearestNeighbour_alt(scalar u_0,vector point) const;
        
        scalar distanceToNurbsSurface(scalar para, vector point) const;
        
        void moveNurbs(List<vector> controlPoints);
        vector movementVector(scalar u);
        
        scalar radius() const {return diameter/2;};
                
        void print() const {Info<<controlPoints[0]<<"->"<<controlPoints[1]<<endl;}
        void printBox() const {Info<<computeBoundingBox().Min<<"->"<<computeBoundingBox().Max<<endl;}
                
    private:
        scalar Weights_B_Spline_Derivative(int k, int l, scalar u, scalar v, nurbsStatus state=curr) const;
        vector A(int k,int l, scalar u, scalar v, nurbsStatus state=curr) const;
        scalar supremum_Derivative2(scalar start_u, scalar end_u) const = delete;
        scalar supremum_Derivative2(scalar start_u, scalar start_v, scalar end_u, scalar end_v) const;
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
