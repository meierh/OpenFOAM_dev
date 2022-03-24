#include "Nurbs.H"
#include <math.h> 

Foam::Nurbs2D::Nurbs2D
(
    List<scalarList> knots,
    List<List<vector>> controlPoints,
    scalarListList weights,
    int degree_p,
    int degree_q
    scalar diameter,
    scalar deltaX
):
testContructorParameters(knots,controlPoints),
controlPoints(List<List<List<vector>>>(3,List<List<vector>>(1))),
knots(List<List<scalarList>>(3,List<scalarList>(1))),
weights(List<List<scalarList>>(3,List<scalarList>(1))),
weightedControlPoints(List<List<List<vector>>>(3,List<List<vector>>(1))),
n_m({knots[0].size(),knots[1].size()}), // former m
cPdim({controlPoints.size(),controlPoints[0].size()}),
p_q({degree_p,degree_q}),
minPara(scalarList({knots[u].first(),knots[v].first()})),
maxPara(scalarList({knots[u].last(),knots[v].last()})),
diameter(diameter),
deltaX(deltaX)
{
    this->controlPoints[initial] = controlPoints;
    this->controlPoints[previous] = controlPoints;
    this->controlPoints[current] = controlPoints;
    
    this->knots[initial] = knots;
    this->knots[previous] = knots;
    this->knots[current] = knots;

    this->weights[initial] = weights;
    this->weights[previous] = weights;
    this->weights[current] = weights;
    
    //Info<<"Construct Nurbs"<<endl;
    if(weights.size() != controlPoints.size())
    {
        FatalErrorInFunction
        << " A Nurbs Curve must have an equal amount of control Points and weights!"<<endl
        << " Currently there are "<<weights.size()<<" weights and "<<controlPoints.size()<<" control Points!"
        << exit(FatalError);
    }
    for(int i=0;i<controlPoints.size();i++)
    {
        if(weights[i].size() != controlPoints[i].size())
        {
            FatalErrorInFunction
            << " A Nurbs Curve must have an equal amount of control Points and weights!"<<endl
            << " Currently there are "<<weights[i].size()<<" weights and "<<controlPoints[i].size()<<" control Points!"
            << exit(FatalError);
        }
    }
    //m = knots->size();
    //n = controlPoints->size();
    //p = degree;
    //this->diameter = diameter;
    //this->deltaX = deltaX;
    if(n_m[u] < cPdim[u]+p_q[u]+1)
    {
        FatalErrorInFunction
        << " A Nurbs Curve must have a number of knots equal to the number of control Points plus the degree plus one"<<endl
        << " Currently there are "<<n_m[u]<<" knots and "<<cPdim[u]<<" control Points and degree "<<p_q[u]<<" in u direction"
        << exit(FatalError);
    }
    if(n_m[v] < cPdim[v]+p_q[v]+1)
    {
        FatalErrorInFunction
        << " A Nurbs Curve must have a number of knots equal to the number of control Points plus the degree plus one"<<endl
        << " Currently there are "<<n_m[v]<<" knots and "<<cPdim[v]<<" control Points and degree "<<p_q[v]<<" in v direction"
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
    
    this->weightedControlPoints[initial] = List<List<vector>>(cPdim[u]);
    this->weightedControlPoints[previous] = List<List<vector>>(cPdim[u]);
    this->weightedControlPoints[current] = List<List<vector>>(cPdim[u]);
    for(int i=0;i<cPdim[u];i++)
    {
        this->weightedControlPoints[initial][i] = List<vector>>(cPdim[v]);
        for(int j=0;j<cPdim[v];j++)
        {
            this->weightedControlPoints[initial][i][j] =
            this->weightedControlPoints[previous][i][j] =
            this->weightedControlPoints[current][i][j] = this->weights[current][i][j] * this->controlPoints[current][i][j];
        }
    }
}

void Foam::Nurbs2D::testContructorParameters(List<scalarList> knots, List<List<vector>> controlPoints)
{
    if(knots.size() != 2)
        FatalErrorInFunction<<"Knots vector is not two-dimensional"<<exit(FatalError);
    if(controlPoints.size()<1)
        FatalErrorInFunction<<"Control Point vector is not two-dimensional"<<exit(FatalError);
}

scalar Foam::Nurbs2D::Weights_B_Spline_Derivative //The Nurbs Book Equations 3.8,3.4 S.97
(
    int k,
    int l,
    scalar u,
    scalar v,
    nurbsStatus state
) const
{
    scalarListList& weights = this->weights[state];
    scalar res = 0;    
    for(int i=0; i<n_m[u]-k; i++)
    {
        for(int j=0; j<n_m[v]-l; j++)
        {
            res += B_Spline_Basis(i,p_q[u]-k,u,state) *
                   B_Spline_Basis(j,p_q[v]-l,v,state) * 
                   Control_Point_Derivative<scalar>({k,l},{i,j},weights,state);
        }
    }
    return res;
}

vector Foam::Nurbs2D::A //The Nurbs Book Equation 4.8 S.125 and Equation 3.8,3.4,3.3 S.97
(
    int k,
    int l,
    scalar u,
    scalar v,
    nurbsStatus state
) const
{
    List<List<vector>>& weightedControlPoints = this->weightedControlPoints[state];

    //Info<<"n-k:"<<n-k<<endl;
    vector res(0,0,0);
    for(int i=0; i<n_m[u]-k; i++)
    {
        for(int j=0; j<n_m[v]-l; j++)
        {
            res += B_Spline_Basis(i,p_q[u]-k,u,state) *
                   B_Spline_Basis(j,p_q[v]-l,v,state) * 
                   Control_Point_Derivative<vector>({k,l},{i,j},weightedControlPoints,state);
        }
    }
    //Info<<"Return "<<res<<endl;
    return res;
}

vector Foam::Nurbs2D::Surface_Derivative
(
    int k,
    int l,
    scalar u,
    scalar v,
    nurbsStatus state=curr
) const
{
    if(state==prev && !nurbsMoved)
        FatalErrorInFunction<< "Prev nurbs status can only be used if Nurbs has been moved"<< exit(FatalError);
    
    if(u>=maxPara[dimState.u] || u<minPara[dimState.u] || v>=maxPara[dimState.v] || v<minPara[dimState.v])
    {
        FatalErrorInFunction
        << " Parameter u:"<<u<<" and v:"<<v<<" has to be within ["<<minPara[dimState.u]<<","<<maxPara[dimState.u]<<"),["<<minPara[dimState.v]<<","<<maxPara[dimState.v]<<") !"<<endl
        << exit(FatalError);
    }
    if(k > p_q[dimState.u] || l > p_q[dimState.v])
    )
    {
        FatalErrorInFunction
        << " Called the ("<<k<<","<<l<<")-th Derivative of a ("<<p_q[dimState.u]<<","<<p_q[dimState.v]<<")-th order Nurbs. This will not work!"<<endl
        << exit(FatalError);
    }
    
    vector A_res = A(k,l,u,v,state);
    scalar w = Weights_B_Spline_Derivative(0,0,u,v,state);

    vector B_i_0(0,0,0);
    for(int i=1;i<k+1;i++)
    {
        scalar koeff = binomial(k,i);
        scalar w_i_0 = Weights_B_Spline_Derivative(i,0,u,v,state);
        vector S_kmi_l = Surface_Derivative(k-i,l,u,v,state);
        B_i_0 += koeff * w_i_0 * S_kmi_l;
    }
    vector B_0_j(0,0,0);
    for(int j=1;j<l+1;j++)
    {
        scalar koeff = binomial(l,j);
        scalar w_0_j = Weights_B_Spline_Derivative(0,j,u,v,state);
        vector S_k_lmj = Surface_Derivative(k,l-j,u,v,state);
        B_0_j += koeff * w_0_j * S_k_lmj;
    }
    vector B_i_j(0,0,0);
    for(int i=1;i<k+1;i++)
    {
        vector inter(0,0,0);
        for(int j=1;j<l+1;j++)
        {
            scalar koeff = binomial(l,j);
            scalar w_i_j = Weights_B_Spline_Derivative(i,j,u,v,state);
            vector S_kmi_lmj = Surface_Derivative(k-i,l-j,u,v,state);
            inter += koeff * w_i_j * S_kmi_lmj;
        }
        B_i_j += binomial(k,i) * inter;
    }
    
    if(w == 0)
        w = 1;
    
    return (A_res-B_i_0-B_0_j-B_i_j)/w;    
}

Foam::BoundingBox Foam::Nurbs::computeBoundingBox
(
    scalar start_u,
    scalar start_v,
    scalar end_u,
    scalar end_v
)
) const
{
    //Info<<"Compute Bounding Box"<<endl;
    if(start_u < minPara[u] || start_u >= maxPara[u])
    {
        FatalErrorInFunction
        << " Start value of Bounding Box has to be in ["<<minPara[u]<<","<<maxPara[u]<<")"<<endl
        << exit(FatalError);
    }
    if(end_u < minPara[u] || end_u >= maxPara[u])
    {
        FatalErrorInFunction
        << " End value of Bounding Box has to be in ["<<minPara[u]<<","<<maxPara[u]<<") but is "<<end<<endl
        << exit(FatalError);
    }
    if(start_v < minPara[v] || start_v >= maxPara[v])
    {
        FatalErrorInFunction
        << " Start value of Bounding Box has to be in ["<<minPara[v]<<","<<maxPara[v]<<")"<<endl
        << exit(FatalError);
    }
    if(end_v < minPara[v] || end_v >= maxPara[v])
    {
        FatalErrorInFunction
        << " End value of Bounding Box has to be in ["<<minPara[v]<<","<<maxPara[v]<<") but is "<<end<<endl
        << exit(FatalError);
    }
    
    FatalErrorInFunction<< " Implementation not complete"<< exit(FatalError);
    
    vector startP = Surface_Derivative(0,0,start_u,start_v);
    //Info<<"StartP: "<<startP<<endl;
    vector endP = Curve_Derivative(0,end_u);
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
    scalar koeff = (1.0/8.0)*(end_u-start_u)*(end_u-start_u)*supremum_Derivative2(start_u,end_u);
    for(int d=0;d<3;d++)
    {
        Box.Min[d] -= (koeff+diameter+deltaX);
        Box.Max[d] += (koeff+diameter+deltaX);
    }
    return Box;    
}

scalar Foam::Nurbs2D::supremum_Derivative2
(
    scalar start_u,
    scalar start_v,
    scalar end_u,
    scalar end_v
) const
{
    vector maxD2(   std::numeric_limits<scalar>::lowest(),
                    std::numeric_limits<scalar>::lowest(),
                    std::numeric_limits<scalar>::lowest()
                );
    
    vector res;
    scalar distU = end_u-start_u;
    scalar distV = end_v-start_v;
    label NUM_POINTS = 100;
    scalar Ux,Vx;
    for(int i=0;i<NUM_POINTS;i++)
    {
        Ux = (static_cast<double>(i)/static_cast<double>(NUM_POINTS))*distU+start_u;
        for(int j=0;j<NUM_POINTS;j++)
        {
            Vx = (static_cast<double>(j)/static_cast<double>(NUM_POINTS))*distV+start_v;
            vector du2 = Surface_Derivative(2,0,Ux,Vx);
            vector dv2 = Surface_Derivative(0,2,Ux,Vx);
            vector dudv = Surface_Derivative(1,1,Ux,Vx);

        }
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


