#include "Nurbs.H"
#include <math.h> 

Foam::Nurbs::Nurbs
(
)
{
}

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
    this->knots = std::unique_ptr<scalarList>(std::move(knots));
    this->controlPoints = std::unique_ptr<List<vector>>(std::move(controlPoints));
    this->weights = std::unique_ptr<scalarList>(std::move(weights));
    
    _min_U = this->knots->first();
    Info<<_min_U<<endl;
    _max_U = this->knots->last();
    Info<<_max_U<<endl;
    Info<<"Constructed"<<endl;
    
    this->weightedControlPoints = std::unique_ptr<List<vector>>(new List<vector>(n));
    for(int i=0;i<n;i++)
    {
        Info<<i<<endl;
        (*(this->weightedControlPoints))[i] = (*(this->weights))[i] * (*(this->controlPoints))[i];
    }
    Info<<"Constructed"<<endl;
}

scalar Foam::Nurbs::B_Spline_Basis  // The Nurbs Book Equation 2.5 S.50
(
    int i,
    int p,
    scalar u
)
{
    if(i+p+1 >= knots->size())
    {
        FatalErrorInFunction
        << " B_Spline_Basis is called by i:"<<i<<" p:"<<p<<" although the Nurbs Curve does only have "<<knots->size()<<"knots"<<endl
        << abort(FatalError);
    }
    if(p==0)
    {
        if((*knots)[i] <= u && u < (*knots)[i+1])
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
        scalar factor1_Z = u - (*knots)[i];
        scalar factor1_N = (*knots)[i+p] - (*knots)[i];
        scalar factor1 = 0;
        if(factor1_N != 0)
        {
            factor1 = factor1_Z / factor1_N;
        }
        scalar factor2_Z = (*knots)[i+p+1] - u;
        scalar factor2_N = (*knots)[i+p+1] - (*knots)[i+1];
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
    const List<T>*  controlPoints
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
        if(i>=controlPoints->size())
        {
            FatalErrorInFunction
            << " Something is wrong here. i:"<<i<<" while controlPoints.size():"<<controlPoints->size()<<endl
            << abort(FatalError);
        }
        return (*controlPoints)[i];
    }
    else
    {
        scalar koeff = 0;
        scalar nenner = (*knots)[i+p+1] - (*knots)[i+k];
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
        res += B_Spline_Basis(i+k,p-k,u) * Control_Point_Derivative<scalar>(k,i,weights.get());
    }
    return res;
}

vector Foam::Nurbs::A   //The Nurbs Book Equation 4.8 S.125 and Equation 3.8,3.4 S.97
(
    int k,
    scalar u
)
{    
    vector res(0,0,0);
    for(int i=0;i<n-k;i++)
    {
        Info<<"B_Spline_Basis("<<i+k<<","<<p-k<<","<<u<<")"<<endl;
        res += B_Spline_Basis(i+k,p-k,u) * Control_Point_Derivative<vector>(k,i,weightedControlPoints.get());
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
    scalar u
)
{
    if(k > p)
    {
        FatalErrorInFunction
        << " Called the "<<k<<"-th Derivative of a "<<p<<"-th order Nurbs. This will not work!"<<endl
        << abort(FatalError);
    }
    Info<<"Start"<<endl;
    vector A_res = A(k,u);
    Info<<"Computed A"<<endl;
    scalar w = Weights_B_Spline_Derivative(0,u);
    Info<<"Computed w"<<endl;
    vector B(0,0,0);
    for(int i=1;i<k+1;i++)
    {
        scalar koeff = binomial(n,k);
        scalar w_i = Weights_B_Spline_Derivative(i,u);
        vector C_kmi = Curve_Derivative(k-i,u);
        B += koeff * w_i * C_kmi;
    }
    Info<<"Computed B"<<endl;
    return (A_res-B)/w;
}

inline scalar Foam::Nurbs::euklidianNorm
(
    vector vec
)
{
        return sqrt(vec.x()*vec.x()+vec.y()*vec.y()+vec.z()*vec.z());
}

Foam::BoundingBox Foam::Nurbs::computeBoundingBox()
{
    if(controlPoints->size() == 0)
    {
        FatalErrorInFunction
        << " Nurbs Curve has no controlPoints. Therefore no bounds can be computed!"<<endl
        << abort(FatalError);
    }
    
    BoundingBox MinMaxBox;
    MinMaxBox.Min = controlPoints->first();
    MinMaxBox.Max = controlPoints->first();
    
    for(int i=0;i<n;i++)
    {
        for(int d=0;d<3;d++)
        {
            if((*controlPoints)[i][d] > MinMaxBox.Max[d])
                MinMaxBox.Max[d] = (*controlPoints)[i][d];
                
            if((*controlPoints)[i][d] < MinMaxBox.Min[d])
                MinMaxBox.Min[d] = (*controlPoints)[i][d];
        }
    }
    
    return MinMaxBox;
}

Foam::BoundingBox Foam::Nurbs::computeBoundingBox(scalar start, scalar end)
{
    if(start < _min_U || start >= _max_U)
    {
        FatalErrorInFunction
        << " Start value of Bounding Box has to be in ["<<_min_U<<","<<_max_U<<")"<<endl
        << abort(FatalError);
    }
    if(end < _min_U || end >= _max_U)
    {
        FatalErrorInFunction
        << " End value of Bounding Box has to be in ["<<_min_U<<","<<_max_U<<")"<<endl
        << abort(FatalError);
    }
    BoundingBox empty;
    scalar sup_D2 = 0;
    return empty;    
}
