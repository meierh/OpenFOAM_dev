#include "NURBS.H"
#include <math.h> 

Foam::Nurbs::Nurbs
(
    std::unique_ptr<scalarList> knots,
    std::unique_ptr<List<vector>> controlPoints,
    std::unique_ptr<scalarList> weights,
    int degree
)
{
    if(weights->size() != controlPoints->size())
    {
        FatalErrorInFunction
        << " A Nurbs Curve must have an equal amount of control Points and weights!"<<endl
        << " Currently there are "<<weights->size()<<" weights and "<<controlPoints->size()<<" control Points!"
        << abort(FatalError);
    }
    m = knots->size();
    n = controlPoints->size();
    p = degree;
    if(m < n+p+1)
    {
        FatalErrorInFunction
        << " A Nurbs Curve must have a number of knots greater or equal than the number of control Points plus the degree plus one"<<endl
        << " Currently there are "<<m<<" knots and "<<n<<" control Points and degree "<<p
        << abort(FatalError);
    }
    this->knots = std::unique_ptr<scalarList>(knots);
    this->controlPoints = std::unique_ptr<List<vector>>(controlPoints);
    this->weights = std::unique_ptr<scalarList>(weights);    
}

scalar Foam::Nurbs::B_Spline_Basis  // The Nurbs Book Equation 2.5 S.50
(
    int i,
    int p,
    scalar u
)
{
    if(p==0)
    {
        if(knots[i] <= u && u < knots[i+1])
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        scalar factor1_Z = u - knots[i];
        scalar factor1_N = knots[i+p] - knots[i];
        scalar factor1 = 0;
        if(factor1_N != 0)
        {
            factor1 = factor1_Z / factor1_N;
        }
        scalar factor2_Z = knots[i+p+1] - u;
        scalar factor2_N = knots[i+p+1] - knots[i+1]
        scalar factor2 = 0;
        if(factor2_N != 0)
        {
            factor2 = factor2_Z / factor2_N;
        }
        return factor1*B_Spline_Basis(i,p-1,u) + factor2*B_Spline_Basis(i+1,p-1,u);
    }
}

template <typename T>
T Foam::Nurbs::Control_Point_Derivative //The Nurbs Book Equation 3.8 S.97
(
    int k,
    int i,
    const unique_ptr<List<T>>&  controlPoints
)
{
    if(m < controlPoints->size()+p+1)
    {
        FatalErrorInFunction
        << " A Nurbs Curve must have a number of knots greater or equal than the number of control Points plus the degree plus one"<<endl
        << " Currently there are "<<m<<" knots and "<<controlPoints->size()<<" control Points and degree "<<p
        << abort(FatalError);
    }
    if(m <= i+k)
    {
        FatalErrorInFunction
        << " Procedure must not be called with an i+k greater or equal than the number of knot points"<<endl
        << " Currently there are "<<m<<" knots and i:"<<i<<" k:"<<k
        << abort(FatalError);
    }
    if(k==0)
    {
        return controlPoints[i]
    }
    else
    {
        scalar koeff = 0;
        scalar nenner = knots[i+p+1] - knots[i+k]
        if(nenner == 0)
        {
            koeff = 0;
        }
        else
        {
            koeff = (p-k+1)/nenner;
        }
        T P_ip1_km1 = Control_Point_Derivative(k-1,i+1,controlPoints);
        T P_i_km1 = Control_Point_Derivative(k-1,i,controlPoints);
        
        return koeff*(P_ip1_km1 - P_i_km1);
    }
}

scalar Foam::Nurbs::Weights_B_Spline_Derivative //The Nurbs Book Equations 3.8,3.4 S.97
(
    int k,
    scalar u
)
{
    scalar res = 0;
    for(int i = 0;i<n-k;i++)
    {
        res += B_Spline_Basis(i+k,p-k,u) * Control_Point_Derivative<scalar>(k,i,this->weights);
    }
    return res;
}

vector Foam::Nurbs::A   //The Nurbs Book Equation 4.8 S.125 and Equation 3.8,3.4 S.97
(
    int k,
    scalar u
)
{
    unique_ptr<List<vector>> weightedControlPoints = unique_ptr<List<vector>>();
    for(int i=0;i<controlPoints->size();i++)
    {
        weightedControlPoints.append(weights[i] * controlPoints[i]);
    }
    
    vector res;
    for(int i=0;i<controlPoints->size();i++)
    {
        res += B_Spline_Basis(i+k,p-k,u) * Control_Point_Derivative<vector>(k,i,weightedControlPoints);
    }
    return res;
}

int Foam::Nurbs::binomial
(
    int n,
    int k
)
{
    if(k==0 || n==k)
    {
        return 1;
    }
    else
    {
        return binomial(n-1,k-1)+binomial(n-1,k);
    }
}

vector Foam::Nurbs::Curve_Derivative
(
    int k,
    int u
)
{
    vector A = A(k,u)
    w = Weights_B_Spline_Derivative(0,u);
    vector B;
    for(int i=1;i<k+1;i++)
    {
        scalar koeff = binomial(n,k);
        scalar w_i = Weights_B_Spline_Derivative(i,u);
        vector C_kmi = Curve_Derivative(k-i,u);
        B = koeff * w_i * C_kmi;
    }    
    return (A-B)/w;
}

inline scalar Foam::Nurbs::euklidianNorm
(
    vector vec
)
{
        return sqrt(vec.x()*vec.x()+vec.y()*vec.y()+vec.z()*vec.z());
}


}
