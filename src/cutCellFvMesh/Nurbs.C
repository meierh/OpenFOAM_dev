#include "Nurbs.H"
#include <math.h> 

bool Foam::BoundingBox::isInside(vector point) const
{
    bool inside = true;
    for(int d=0;d<3;d++)
    {
        if(Min[d]>point[d])
            inside = false;
        if(Max[d]<point[d])
            inside = false;
    }
    return inside;
}

Foam::Nurbs1D::Nurbs1D
(
    scalarList knots,
    List<vector> controlPoints,
    scalarList weights,
    int degree,
    scalar diameter,
    scalar deltaX
):
controlPoints(List<List<List<vector>>>(3,List<List<vector>>(1))),
knots(List<List<scalarList>>(3,List<scalarList>(1))),
weights(List<List<scalarList>>(3,List<scalarList>(1))),
weightedControlPoints(List<List<List<vector>>>(3,List<List<vector>>(1))),
m(knots.size()),
n(controlPoints.size()),
p(degree),
minPara(scalarList(1,knots.first())),
maxPara(scalarList(1,knots.last())),
diameter(diameter),
deltaX(deltaX)
{
    controlPoints[initial][u] = controlPoints;
    controlPoints[previous][u] = controlPoints;
    controlPoints[current][u] = controlPoints;
    
    knots[initial][u] = knots;
    knots[previous][u] = knots;
    knots[current][u] = knots;

    weights[initial][u] = weights;
    weights[previous][u] = weights;
    weights[current][u] = weights;
    
    //Info<<"Construct Nurbs"<<endl;
    if(weights.size() != controlPoints.size())
    {
        FatalErrorInFunction
        << " A Nurbs Curve must have an equal amount of control Points and weights!"<<endl
        << " Currently there are "<<weights.size()<<" weights and "<<controlPoints.size()<<" control Points!"
        << exit(FatalError);
    }
    //m = knots->size();
    //n = controlPoints->size();
    //p = degree;
    //this->diameter = diameter;
    //this->deltaX = deltaX;
    if(m < n+p+1)
    {
        FatalErrorInFunction
        << " A Nurbs Curve must have a number of knots equal to the number of control Points plus the degree plus one"<<endl
        << " Currently there are "<<m<<" knots and "<<n<<" control Points and degree "<<p
        << exit(FatalError);
    }
    //this->knots = std::unique_ptr<scalarList>(std::move(knots));
    //this->controlPoints = std::unique_ptr<List<vector>>(std::move(controlPoints));
    //this->weights = std::unique_ptr<scalarList>(std::move(weights));
    
    //_min_U = this->knots->first();
    //Info<<_min_U<<endl;
    //_max_U = this->knots->last();
    //Info<<_max_U<<endl;
    //Info<<"Constructed"<<endl;
    
    weightedControlPoints[initial][u] = List<vector>(n);
    weightedControlPoints[previous][u] = List<vector>(n);
    weightedControlPoints[current][u] = List<vector>(n);
    for(int i=0;i<n;i++)
    {
        weightedControlPoints[initial][u][i] = weights[initial][u][i] * controlPoints[initial][u][i];
        weightedControlPoints[previous][u][i] = weights[previous][u][i] * controlPoints[previous][u][i];
        weightedControlPoints[current][u][i] = weights[current][u][i] * controlPoints[current][u][i];
    }
}

scalar Foam::Nurbs::B_Spline_Basis // The Nurbs Book Equation 2.5 S.50
(
    int i,
    int p,
    scalar u,
    nurbsStatus state,
    dimState dState
) const
{
    scalarList& knots = this->knots[state][dState];
    
    if(i+p+1 >= knots.size())
    {
        FatalErrorInFunction
        << " B_Spline_Basis is called by i:"<<i<<" p:"<<p<<" although the Nurbs Curve does only have "<<knots.size()<<"knots"<<endl
        << exit(FatalError);
    }
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
        scalar factor2_N = knots[i+p+1] - knots[i+1];
        scalar factor2 = 0;
        if(factor2_N != 0)
        {
            factor2 = factor2_Z / factor2_N;
        }
        return factor1*B_Spline_Basis(i,p-1,u,state,dState) + factor2*B_Spline_Basis(i+1,p-1,u,state,dState);
    }
}

template <typename T>
T Foam::Nurbs::Control_Point_Derivative //The Nurbs Book Equation 3.8 S.97
(
    int k,
    int i,
    const List<T>&  controlPoints,
    nurbsStatus state
) const
{
    scalarList knots = (state==curr)*this->knots+(state==prev)*this->knots_prev;

    if(m < controlPoints.size()+p+1)
    {
        FatalErrorInFunction
        << " A Nurbs Curve must have a number of knots greater or equal than the number of control Points plus the degree plus one"<<endl
        << " Currently there are "<<m<<" knots and "<<controlPoints.size()<<" control Points and degree "<<p
        << exit(FatalError);
    }
    if(m <= i+k)
    {
        FatalErrorInFunction
        << " Procedure must not be called with an i+k greater or equal than the number of knot points"<<endl
        << " Currently there are "<<m<<" knots and i:"<<i<<" k:"<<k
        << exit(FatalError);
    }
    if(k==0)
    {
        if(i>=controlPoints.size())
        {
            FatalErrorInFunction
            << " Something is wrong here. i:"<<i<<" while controlPoints.size():"<<controlPoints.size()<<endl
            << exit(FatalError);
        }
        return controlPoints[i];
    }
    else
    {
        scalar koeff = 0;
        scalar nenner = knots[i+p+1] - knots[i+k];
        if(nenner == 0)
        {
            koeff = 0;
        }
        else
        {
            koeff = (p-k+1)/nenner;
        }
        T P_ip1_km1 = Control_Point_Derivative(k-1,i+1,controlPoints,state);
        T P_i_km1 = Control_Point_Derivative(k-1,i,controlPoints,state);
        
        return koeff*(P_ip1_km1 - P_i_km1);
    }
}

scalar Foam::Nurbs::Weights_B_Spline_Derivative //The Nurbs Book Equations 3.8,3.4 S.97
(
    int k,
    scalar u,
    nurbsStatus state
) const
{
    scalarList weights = (state==curr)*this->weights+(state==prev)*this->weights_prev;
    scalar res = 0;
    for(int i = 0;i<n-k;i++)
    {
        res += B_Spline_Basis(i+k,p-k,u,state) * Control_Point_Derivative<scalar>(k,i,weights,state);
    }
    return res;
}

vector Foam::Nurbs::A //The Nurbs Book Equation 4.8 S.125 and Equation 3.8,3.4,3.3 S.97
(
    int k,
    scalar u,
    nurbsStatus state
) const
{
    List<vector> weightedControlPoints = (state==curr)*this->weightedControlPoints+(state==prev)*this->weightedControlPoints_prev;

    //Info<<"n-k:"<<n-k<<endl;
    vector res(0,0,0);
    for(int i=0;i<n-k;i++)
    {
        //Info<<"B_Spline_Basis("<<i+k<<","<<p-k<<","<<u<<")"<<endl;
        res += B_Spline_Basis(i+k,p-k,u,state) * Control_Point_Derivative<vector>(k,i,weightedControlPoints,state);
        //Info<<"Done"<<endl;
    }
    //Info<<"Return "<<res<<endl;
    return res;
}

int Foam::Nurbs::binomial
(
    int n,
    int k
) const
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
    scalar u,
    nurbsStatus state
) const
{
    if(state==prev && !nurbsMoved)
        FatalErrorInFunction<< "Prev nurbs status can only be used if Nurbs has been moved"<< exit(FatalError);
    
    if(u>=_max_U || u<_min_U)
    {
        FatalErrorInFunction
        << " Parameter u has to be within ["<<_min_U<<","<<_max_U<<") but is "<<u<<"!"<<endl
        << exit(FatalError);
    }
    if(k > p)
    {
        FatalErrorInFunction
        << " Called the "<<k<<"-th Derivative of a "<<p<<"-th order Nurbs. This will not work!"<<endl
        << exit(FatalError);
    }
    //Info<<"-----------------------"<<endl;
    vector A_res = A(k,u,state);
    scalar w = Weights_B_Spline_Derivative(0,u,state);
    //Info<<"A:"<<A_res<<" w:"<<w<<endl;
    vector B(0,0,0);
    for(int i=1;i<k+1;i++)
    {
        scalar koeff = binomial(k,i);
        scalar w_i = Weights_B_Spline_Derivative(i,u,state);
        vector C_kmi = Curve_Derivative(k-i,u,state);
        B += koeff * w_i * C_kmi;
        //Info<<koeff<<"  "<<w_i<<"   "<<C_kmi<<endl;
    }
    
    //Info<<"A_"<<k<<"_:"<<A_res<<endl;
    //Info<<"w_"<<k<<"_:"<<w<<endl;
    //Info<<"B_"<<k<<"_:"<<B<<endl;
    //Info<<"-----------------------"<<endl;
    
    if(w == 0)
        w = 1;
    
    return (A_res-B)/w;
}

Foam::BoundingBox Foam::Nurbs::computeBoundingBox() const
{
    if(controlPoints.size() == 0)
    {
        FatalErrorInFunction
        << " Nurbs Curve has no controlPoints. Therefore no bounds can be computed!"<<endl
        << exit(FatalError);
    }
    
    BoundingBox MinMaxBox;
    MinMaxBox.Min = controlPoints.first();
    MinMaxBox.Max = controlPoints.first();
    
    for(int i=0;i<n;i++)
    {
        for(int d=0;d<3;d++)
        {
            if(controlPoints[i][d] > MinMaxBox.Max[d])
                MinMaxBox.Max[d] = controlPoints[i][d];
                
            if(controlPoints[i][d] < MinMaxBox.Min[d])
                MinMaxBox.Min[d] = controlPoints[i][d];
        }
    }
    
    for(int d=0;d<3;d++)
    {
        MinMaxBox.Min[d] = MinMaxBox.Min[d]-(2*deltaX+diameter);
        MinMaxBox.Max[d] = MinMaxBox.Max[d]+(2*deltaX+diameter);
    }
    
    return MinMaxBox;
}

Foam::BoundingBox Foam::Nurbs::computeBoundingBox
(
    scalar start,
    scalar end
) const
{
    //Info<<"Compute Bounding Box"<<endl;
    if(start < _min_U || start >= _max_U)
    {
        FatalErrorInFunction
        << " Start value of Bounding Box has to be in ["<<_min_U<<","<<_max_U<<")"<<endl
        << exit(FatalError);
    }
    if(end < _min_U || end >= _max_U)
    {
        FatalErrorInFunction
        << " End value of Bounding Box has to be in ["<<_min_U<<","<<_max_U<<") but is "<<end<<endl
        << exit(FatalError);
    }
    
    vector startP = Curve_Derivative(0,start);
    //Info<<"StartP: "<<startP<<endl;
    vector endP = Curve_Derivative(0,end);
    //Info<<"EndP: "<<endP<<endl;
    BoundingBox Box;
    for(int d=0;d<3;d++)
    {
        if(Box.Min[d] > startP[d])
            Box.Min[d] = startP[d];
        if(Box.Max[d] < startP[d])
            Box.Max[d] = startP[d];
        if(Box.Min[d] > endP[d])
            Box.Min[d] = endP[d];
        if(Box.Max[d] < endP[d])
            Box.Max[d] = endP[d];
    }
    //Info<<"Initial Bounding Box"<<endl;
    scalar koeff = (1.0/8.0)*(end-start)*(end-start)*supremum_Derivative2(start,end);
    for(int d=0;d<3;d++)
    {
        Box.Min[d] -= (koeff+diameter+deltaX);
        Box.Max[d] += (koeff+diameter+deltaX);
    }
    return Box;    
}

scalar Foam::Nurbs::supremum_Derivative2
(
    scalar start_u,
    scalar end_u
) const
{
    vector maxD2(   std::numeric_limits<scalar>::lowest(),
                    std::numeric_limits<scalar>::lowest(),
                    std::numeric_limits<scalar>::lowest()
                );
    
    vector res;
    scalar dist = end_u-start_u;
    label NUM_POINTS = 100;
    scalar p;
    for(int i=0;i<NUM_POINTS;i++)
    {
        p = (static_cast<double>(i)/static_cast<double>(NUM_POINTS))*dist+start_u;
        res = Curve_Derivative(2,p);
        for(int d=0;d<3;d++)
        {
            res[d] = std::abs(res[d]);
            maxD2[d] = std::max(res[d],maxD2[d]);
        }        
    }
    scalar sup = 0;
    for(int d=0;d<3;d++)
    {
        sup += maxD2[d] * maxD2[d];
    }
    return sup;
}

vector Foam::Nurbs::vectorCurveToPoint_Derivative
(
    int k,
    scalar para,
    vector point
) const
{
    vector C = Curve_Derivative(k,para);
    if(k==0)
        return C-point;
    else
        return C;
}

scalar Foam::Nurbs::distCurveToPoint_Deriv0
(
    scalar para,
    vector point
) const
{
    return euklidianNorm(vectorCurveToPoint_Derivative(0,para,point));
}

scalar Foam::Nurbs::distCurveToPoint_Deriv1
(
    scalar para,
    vector point
) const
{
    /*
    vector PointToCurve = vectorCurveToPoint_Derivative(0,para,point);
    vector PointToCurve_D1 = vectorCurveToPoint_Derivative(1,para,point);
    scalar N = euklidianNorm(PointToCurve);
    scalar Z = 0;
    for(int d=0;d<3;d++)
        Z += PointToCurve[d]*PointToCurve_D1[d];
    return N/Z;
    */
    
    vector C = Curve_Derivative(0,para);
    vector C_D1 = Curve_Derivative(1,para);
    vector C_P = C-point;
    
    return (C_P && C_D1)/std::sqrt(C_P && C_P);
}

scalar Foam::Nurbs::distCurveToPoint_Deriv2
(
    scalar para,
    vector point
) const
{
    /*
    vector PointToCurve = vectorCurveToPoint_Derivative(0,para,point);
    vector PointToCurve_D1 = vectorCurveToPoint_Derivative(1,para,point);
    vector PointToCurve_D2 = vectorCurveToPoint_Derivative(2,para,point);
    
    //Info<<endl<<"Computed all vectors"<<PointToCurve<<PointToCurve_D1<<PointToCurve_D2<<endl;
    
    scalar v_i_2 = 0;
    for(int d=0;d<3;d++)
        v_i_2 += PointToCurve[d]*PointToCurve[d];
    
    scalar v_i_D0_D1 = 0;
    for(int d=0;d<3;d++)
        v_i_D0_D1 += PointToCurve[d]*PointToCurve_D1[d];
    
    scalar v_i_D0_D2 = 0;
    for(int d=0;d<3;d++)
        v_i_D0_D2 += PointToCurve[d]*PointToCurve_D2[d];
    
    scalar sqrt_v_i_2 = std::sqrt(v_i_2);
    scalar sum1 = -(1/(sqrt_v_i_2*sqrt_v_i_2*sqrt_v_i_2))*v_i_D0_D1*v_i_D0_D1;
    scalar sum2 = (1/(sqrt_v_i_2))*(v_i_2+v_i_D0_D2);
    
    return sum1+sum2;
    */
    
    vector C = Curve_Derivative(0,para);
    vector C_D1 = Curve_Derivative(1,para);
    vector C_D2 = Curve_Derivative(2,para);
    vector C_P = C-point;
    
    return ((C_D2 && C_P)+(C_D1 && C_D1))/(C_P && C_P);    
}

scalar Foam::Nurbs::newtonIterateNearestNeighbour_alt
(
    scalar u_0,
    vector point
) const
{
    scalar change = distCurveToPoint_Deriv1(u_0,point)/distCurveToPoint_Deriv2(u_0,point);
    //bool posSign = (change > 0);
    int i=0;
    scalar f = distCurveToPoint_Deriv1(u_0,point);
    while(abs(f) > 1e-6)
    {   
        /*
        Info<<"u:"<<u_0<<"    "<<"change:"<<change<<endl;
        Info<<vectorCurveToPoint_Derivative(0,u_0,point)<<
        " d:"<<distCurveToPoint_Deriv0(u_0,point)<<
        " d_1:"<<distCurveToPoint_Deriv1(u_0,point)<<
        " d_2:"<<distCurveToPoint_Deriv2(u_0,point)<<endl;
        */
        u_0 = u_0 - change;
        f = distCurveToPoint_Deriv1(u_0,point);
        change = f/distCurveToPoint_Deriv2(u_0,point);
        if(i++ > 1000000)
            break;
    }
    /*
    Info<<"Count: "<<i<<endl;
    Info<<"u:"<<u_0<<"    "<<"change:"<<change<<endl;
    Info<<vectorCurveToPoint_Derivative(0,u_0,point)<<
    " d:"<<distCurveToPoint_Deriv0(u_0,point)<<
    " d_1:"<<distCurveToPoint_Deriv1(u_0,point)<<
    " d_2:"<<distCurveToPoint_Deriv2(u_0,point)<<endl;
    */
    return u_0;
}

scalar Foam::Nurbs::newtonIterateNearestNeighbour
(
    scalar u_0,
    vector point
) const
{
    vector C = Curve_Derivative(0,u_0);
    vector C_D1 = Curve_Derivative(1,u_0);
    vector C_D2 = Curve_Derivative(2,u_0);
    scalar f = (point-C) && C_D1;
    scalar f_D1 = -(C_D1 && C_D1)+((point-C)&&C_D2);
    scalar epsilon = 1e-10;
    int iterations = 0;
    int maxIterations = 100;
    //Info<<"Point: "<<point<<endl;
    //Info<<"It:0 f("<<u_0<<")="<<f<<endl;
    scalar min_U = this->min_U();
    scalar max_U = this->max_U();
    int hitMaxOrMinCounter = 0;
    while(abs(f) > epsilon)
    {
        iterations++;
        C = Curve_Derivative(0,u_0);
        //Info<<"\tC("<<u_0<<"): "<<C<<endl;
        C_D1 = Curve_Derivative(1,u_0);
        //Info<<"\tC_D1("<<u_0<<"): "<<C_D1<<endl;
        C_D2 = Curve_Derivative(2,u_0);
        //Info<<"\tC_D2("<<u_0<<"): "<<C_D2<<endl;
        f = (point-C) && C_D1;
        //Info<<"\tf("<<u_0<<"): "<<f<<endl;
        f_D1 = -(C_D1 && C_D1)+((point-C)&&C_D2);
        //Info<<"\tf_D1("<<u_0<<"): "<<f_D1<<endl;
        if(f_D1 == 0 || iterations > maxIterations)
            break;
        u_0 = u_0 - (f/f_D1);
        //Info<<"It:"<<iterations<<" f("<<u_0<<")="<<f<<endl;
        if(u_0 > max_U)
        {
            u_0 = max_U;
            //Info<<"Break: f("<<u_0<<")"<<endl;
            if(hitMaxOrMinCounter >= 2)
                break;
            hitMaxOrMinCounter++;

        }
        if(u_0 < min_U)
        {
            u_0 = min_U;
            //Info<<"Break: f("<<u_0<<")"<<endl;
            if(hitMaxOrMinCounter >= 2)
                break;
            hitMaxOrMinCounter++;
        }
    }
    //Info<<"Res U: "<<u_0<<endl;
    return u_0;
}

scalar Foam::Nurbs::distanceToNurbsSurface
(
    scalar para,
    vector point
) const
{
    vector C = Curve_Derivative(0,para);
    scalar distToCentre = euklidianNorm(C-point);
    return distToCentre - diameter;
}

void Foam::Nurbs::moveNurbs
(
    List<vector> controlPoints
)
{
    if(controlPoints.size() != this->controlPoints.size())
        FatalErrorInFunction<<"Number of nurbs curves must stay the same at all times"<<exit(FatalError);
    
    knots_prev = this->knots;
    controlPoints_prev = this->controlPoints;
    weights_prev = this->weights;
    weightedControlPoints_prev = this->weightedControlPoints;
    
    this->controlPoints = controlPoints;    
    this->weightedControlPoints = List<vector>(n);
    for(int i=0;i<n;i++)
        this->weightedControlPoints[i] = this->weights[i] * this->controlPoints[i];
    
    nurbsMoved = true;
}

vector Foam::Nurbs::movementVector
(
    scalar u
)
{
    vector C_curr = Curve_Derivative(0,u,curr);
    vector C_prev = Curve_Derivative(0,u,prev);
    //Info<<"C_curr:"<<C_curr<<endl;
    //Info<<"C_prev:"<<C_prev<<endl;
    //FatalErrorInFunction<<"Temp Stop"<< exit(FatalError);            
    return C_curr-C_prev;
}


