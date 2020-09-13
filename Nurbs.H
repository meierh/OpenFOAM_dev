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
        const std::unique_ptr<scalarList> knots;
        const std::unique_ptr<List<vector>> controlPoints;
        const std::unique_ptr<scalarList> weights;
        const int m;
        const int n;
        const int p;
        std::unique_ptr<List<vector>> weightedControlPoints;
        const scalar _min_U;
        const scalar _max_U;
        const scalar diameter;
        const scalar deltaX;
        
    public:
        explicit Nurbs
        (
            std::unique_ptr<scalarList> knots,
            std::unique_ptr<List<vector>> controlPoints,
            std::unique_ptr<scalarList> weights,
            int degree,
            scalar diameter=0,
            scalar deltaX=0
        );
        
        vector Curve_Derivative(int k,scalar u) const;
        
        BoundingBox computeBoundingBox(scalar start_u, scalar end_u) const;
        
        BoundingBox computeBoundingBox() const;
        
        scalar min_U() const {return _min_U;}
        
        scalar max_U() const {return _max_U;}
                
    private:
        scalar B_Spline_Basis(int i, int p, scalar u) const;
        template <typename T> T Control_Point_Derivative(int k, int i, const List<T>*  controlPoints) const;
        scalar Weights_B_Spline_Derivative(int k, scalar u) const;
        vector A(int k, scalar u) const;
        int binomial(int n, int k) const;
        inline scalar euklidianNorm(vector vec) const;
        scalar supremum_Derivative(scalar start_u, scalar end_u) const;
};

void UnitTest_Nurbs();
}
#endif
