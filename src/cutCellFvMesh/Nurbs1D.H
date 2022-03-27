#ifndef NURBS_H
#define NURBS_H

#include "fvCFD.H"
#include <memory>
//#include "Nurbs.H"

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
    
class Nurbs1D
{
    protected:      
        // Outer List[initial,previous,current] - Inner List[u,v]
        List<List<scalarList>> knots;
        
        // Outer List[initial,previous,current] - Inner Matrix[dim(u),dim(v)]
        List<List<List<vector>>> controlPoints;
        List<List<List<scalar>>> weights;
        List<List<List<vector>>> weightedControlPoints;
        
        const labelList n_m;
        const labelList cPdim;
        const labelList p_q;
        
        const scalarList minPara;
        const scalarList maxPara;
        const scalar diameter;
        const scalar deltaX;
        
        bool nurbsMoved = false;
        
        enum timeState {initial=0,previous=1,current=2};
        enum dimState {u=0,v=1};
        
    public:
        Nurbs1D
        (
            scalarList knots,
            List<vector> controlPoints,
            scalarList weights,
            int degree,
            scalar diameter=0.2,
            scalar deltaX=0.2
        );
        
        enum nurbsStatus
        {
            prev = 1,
            curr = 2
        };
        
        scalar getBoundingBoxOverhang() const {return diameter+deltaX;}
        
        vector Curve_Derivative(int k,scalar u,nurbsStatus state=curr) const;
        
        BoundingBox computeBoundingBox(scalar start_u, scalar end_u) const;
        
        BoundingBox computeBoundingBox() const;
        
        scalar min_U() const {return minPara[u];}
        
        scalar max_U() const {return std::nexttoward(maxPara[u],minPara[u]);}
        
        int nbrKnotsU() const {return knots[u].size();}
        
        int degreeU() const {return p_q[u];}
        
        scalar newtonIterateNearestNeighbour(scalar u_0,vector point) const;
        scalar newtonIterateNearestNeighbour_alt(scalar u_0,vector point) const;
        
        scalar distanceToNurbsSurface(scalar para, vector point) const;
        
        void moveNurbs(List<List<vector>> controlPoints);
        vector movementVector(scalar u);
        
        scalar radius() const {return diameter/2;};
                
        void print() const {Info<<controlPoints[0]<<"->"<<controlPoints[1]<<endl;}
        void printBox() const {Info<<computeBoundingBox().Min<<"->"<<computeBoundingBox().Max<<endl;}
                
    protected:
        Nurbs1D
        (
            List<scalarList> knots,
            List<List<vector>> controlPoints,
            scalarListList weights,
            int degree_p,
            int degree_q,
            scalar diameter,
            scalar deltaX
        ); 
        
        scalar B_Spline_Basis(int i, int p, scalar u, dimState dState, nurbsStatus state=curr) const;
        template <typename T> T Control_Point_Derivative(List<int> k_l, List<int> i_j, const List<List<T>>& controlPoints, nurbsStatus state=curr) const;
        scalar Weights_B_Spline_Derivative(int k, scalar u, nurbsStatus state=curr) const;
        vector A(int k, scalar u, nurbsStatus state=curr) const;
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
