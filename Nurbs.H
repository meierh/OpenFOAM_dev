#include "fvCFD.H"
#include <memory>

namespace Foam
{   
struct BoundingBox
{
    vector Min;
    vector Max;
};
    
class Nurbs
{
    private:
        std::unique_ptr<scalarList> knots;
        int m;
        int n;
        int p;
        std::unique_ptr<List<vector>> controlPoints;
        std::unique_ptr<scalarList> weights;
        std::unique_ptr<List<vector>> weightedControlPoints;
        scalar _min_U;
        scalar _max_U;
        
    public:
        explicit Nurbs
        (
            std::unique_ptr<scalarList> knots,
            std::unique_ptr<List<vector>> controlPoints,
            std::unique_ptr<scalarList> weights,
            int degree
        );
        vector Curve_Derivative(int k, scalar u);
        BoundingBox computeBoundingBox(scalar start_u, scalar end_u);
        BoundingBox computeBoundingBox();
        scalar min_U() const {return _min_U;}
        scalar max_U() const {return _max_U;}
        
    private:
        scalar B_Spline_Basis(int i, int p, scalar u);
        template <typename T> T Control_Point_Derivative(int k, int i, const List<T>*  controlPoints);
        scalar Weights_B_Spline_Derivative(int k, scalar u);
        vector A(int k, scalar u);
        int binomial(int n, int k);
        inline scalar euklidianNorm(vector vec);
        scalar supremum_Derivative(scalar start_u, scalar end_u);
};
}
