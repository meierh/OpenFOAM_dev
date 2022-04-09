#include "Nurbs1D.H"
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
knots(List<List<scalarList>>(3,List<scalarList>(1))),
controlPoints(List<List<List<vector>>>(3,List<List<vector>>(1))),
weights(List<List<scalarList>>(3,List<scalarList>(1))),
weightedControlPoints(List<List<List<vector>>>(3,List<List<vector>>(1))),
n_m({knots.size()}), // former m
cPdim({controlPoints.size()}),
p_q({degree}),
minPara(scalarList(1,knots.first())),
maxPara(scalarList(1,knots.last())),
diameter(diameter),
deltaX(deltaX)
{
    this->controlPoints[initial][u] = controlPoints;
    this->controlPoints[previous][u] = controlPoints;
    this->controlPoints[current][u] = controlPoints;
    
    this->knots[initial][u] = knots;
    this->knots[previous][u] = knots;
    this->knots[current][u] = knots;

    this->weights[initial][u] = weights;
    this->weights[previous][u] = weights;
    this->weights[current][u] = weights;
    
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
    if(n_m[u] < cPdim[u]+p_q[u]+1)
    {
        FatalErrorInFunction
        << " A Nurbs Curve must have a number of knots equal to the number of control Points plus the degree plus one"<<endl
        << " Currently there are "<<n_m[u]<<" knots and "<<cPdim[u]<<" control Points and degree "<<p_q[u]<<" in u direction"
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
    
    this->weightedControlPoints[initial][u] = List<vector>(cPdim[u]);
    this->weightedControlPoints[previous][u] = List<vector>(cPdim[u]);
    this->weightedControlPoints[current][u] = List<vector>(cPdim[u]);
    for(int i=0;i<cPdim[u];i++)
    {
        this->weightedControlPoints[initial][u][i] =
        this->weightedControlPoints[previous][u][i] =
        this->weightedControlPoints[current][u][i] = this->weights[current][u][i] * this->controlPoints[current][u][i];
    }
    //Info<<"Construct Nurbs done"<<endl;
}

Foam::Nurbs1D::Nurbs1D
(
    List<scalarList> knots,
    List<List<vector>> controlPoints,
    scalarListList weights,
    int degree_p,
    int degree_q,
    scalar diameter,
    scalar deltaX
):
knots(List<List<scalarList>>(3,List<scalarList>(1))),
controlPoints(List<List<List<vector>>>(3,List<List<vector>>(1))),
weights(List<List<scalarList>>(3,List<scalarList>(1))),
weightedControlPoints(List<List<List<vector>>>(3,List<List<vector>>(1))),
n_m(knots.size()==2 ? labelList({knots[0].size(),knots[1].size()}) : labelList({})), // former m
cPdim(controlPoints.size()==2 ? labelList({controlPoints.size(),controlPoints[0].size()}) : labelList({})),
p_q({degree_p,degree_q}),
minPara(knots.size()==2 ? scalarList({knots[u].first(),knots[v].first()}) : scalarList({})),
maxPara(knots.size()==2 ? scalarList({knots[u].last(),knots[v].last()}) : scalarList({})),
diameter(diameter),
deltaX(deltaX)
{
    if(knots.size()!=2 || controlPoints.size()!=2)
        FatalErrorInFunction<<"Illformed 2D Nurbs"<< exit(FatalError);
}

scalar Foam::Nurbs1D::B_Spline_Basis // The Nurbs Book Equation 2.5 S.50
(
    int i,
    int p,
    scalar u,
    dimState dState,
    nurbsStatus state
) const
{
    const scalarList& knots = this->knots[state][dState];
    
    //Info<<"i:"<<i<<"  p:"<<p<<endl;
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
        return factor1*B_Spline_Basis(i,p-1,u,dState,state) + factor2*B_Spline_Basis(i+1,p-1,u,dState,state);
    }
}

template <typename T>
T Foam::Nurbs1D::Control_Point_Derivative //The Nurbs Book Equation 3.8 S.97
(
    List<int> k_l,
    List<int> i_j,
    const List<List<T>>&  controlPoints,
    nurbsStatus state
) const
{
    //Info<<"k_l:"<<k_l<<"  i_j:"<<i_j<<endl;
    dimState treatedDim=dimState::u;
    if(k_l.size()==1)
        treatedDim=dimState::u;
    else if(k_l.size()==2)
    {
        if(k_l[0] >= k_l[1])
            treatedDim=dimState::u;
        else
            treatedDim=dimState::v;
    }
    else
        FatalErrorInFunction<<"Parameter k_l must be of size 1 or 2!"<<exit(FatalError);
        
    const scalarListList& knots = this->knots[state];
    
    //Info<<"treatedDim:"<<treatedDim<<endl;
    
    for(int uv=0;uv<k_l.size();uv++)
    {
        if(n_m[uv] < controlPoints.size()+p_q[uv]+1)
        {
            FatalErrorInFunction
            << " A Nurbs Curve must have a number of knots greater or equal than the number of control Points plus the degree plus one"<<endl
            << " Currently there are "<<n_m[uv]<<" knots and "<<controlPoints.size()<<" control Points and degree "<<p_q[uv]
            << exit(FatalError);
        }
        if(n_m[uv] <= i_j[uv]+k_l[uv])
        {
            FatalErrorInFunction
            << " Procedure must not be called with an i+k greater or equal than the number of knot points"<<endl
            << " Currently there are "<<n_m[uv]<<" knots and i_j[uv]:"<<i_j[uv]<<" k_l[uv]:"<<k_l[uv]
            << exit(FatalError);
        }
    }
    
    label derivTotal=0;
    for(int x=0;x<k_l.size();x++)
        derivTotal+=k_l[x];
    
    //Info<<"derivTotal:"<<derivTotal<<endl;
    
    if(derivTotal==0)
    {
        if(i_j.size()==1)
        {
            //Info<<"controlPoints:"<<controlPoints<<endl;
            if(controlPoints.size()!=1 || i_j[0]>=controlPoints[0].size())
            {
                FatalErrorInFunction
                << " Something is wrong here. controlPoints.size()="<<controlPoints.size()<<"  i:"<<i_j[0]<<" while controlPoints.size():"<<controlPoints[0].size()<<endl
                << exit(FatalError);
            }
            return controlPoints[0][i_j[0]];
        }
        else
        {
            if(i_j[0]>=controlPoints.size() || i_j[1]>=controlPoints[0].size())
            {
                FatalErrorInFunction
                << " Something is wrong here. i:"<<i_j[0]<<" while controlPoints.size():"<<controlPoints[0].size()<<" or "
                <<"j:"<<i_j[1]<<" while controlPoints.size():"<<controlPoints.size() <<endl
                << exit(FatalError);
            }
            return controlPoints[i_j[0]][i_j[1]];
        }
    }
    else
    {
        scalar koeff = 0;
        scalar nenner = knots[treatedDim][i_j[treatedDim]+p_q[treatedDim]+1]
                        - knots[treatedDim][i_j[treatedDim]+k_l[treatedDim]];
        if(nenner == 0)
        {
            koeff = 0;
        }
        else
        {
            koeff = (p_q[treatedDim]-k_l[treatedDim]+1)/nenner;
        }
        List<int> k_l_next(k_l);
        k_l_next[treatedDim]--;
        
        List<int> i_j_next(i_j);
        i_j_next[treatedDim]--;
        
        
        T P_ip1_km1 = Control_Point_Derivative(k_l_next,i_j_next,controlPoints,state);
        T P_i_km1 = Control_Point_Derivative(k_l_next,i_j,controlPoints,state);
        
        return koeff*(P_ip1_km1 - P_i_km1);
    }
}

scalar Foam::Nurbs1D::Weights_B_Spline_Derivative //The Nurbs Book Equations 3.8,3.4 S.97
(
    int k,
    scalar u,
    nurbsStatus state
) const
{
    const scalarListList& weights = this->weights[state];
    scalar res = 0;    
    for(int i=0; i<cPdim[dimState::u]-k; i++)
    {
        res += B_Spline_Basis(i,p_q[dimState::u]-k,u,dimState::u,state) * Control_Point_Derivative<scalar>({k},{i},weights,state);
    }    
    return res;
}

vector Foam::Nurbs1D::A //The Nurbs Book Equation 4.8 S.125 and Equation 3.8,3.4,3.3 S.97
(
    int k,
    scalar u,
    nurbsStatus state
) const
{
    const List<List<vector>>& weightedControlPoints = this->weightedControlPoints[state];

    //Info<<"n_m:"<<n_m<<endl;
    //Info<<"cPdim:"<<cPdim<<endl;
    //Info<<"n_m[u]-k:"<<n_m[dimState::u]-k<<endl;
    vector res(0,0,0);
    for(int i=0;i<cPdim[dimState::u]-k;i++)
    {
        //Info<<"B_Spline_Basis("<<i+k<<","<<p-k<<","<<u<<")"<<endl;
        //Info<<"B_Spline_Basis("<<i<<","<<p_q[dimState::u]-k<<","<<u<<")"<<endl;
        Info<<"B_Spline_Basis("<<i+k<<","<<p_q[dimState::u]-k<<","<<u<<","<<dimState::u<<","<<state<<"):"<<B_Spline_Basis(i+k,p_q[dimState::u]-k,u,dimState::u,state)<<endl;
        Info<<"Control_Point_Derivative<vector>("<<k<<","<<i<<","<<weightedControlPoints<<","<<state<<"):"<<Control_Point_Derivative<vector>({k},{i},weightedControlPoints,state)<<endl;
        res += B_Spline_Basis(i+k,p_q[dimState::u]-k,u,dimState::u,state) * Control_Point_Derivative<vector>({k},{i},weightedControlPoints,state);
        //Info<<"Done"<<endl;
    }
    //Info<<"Return "<<res<<endl;
    return res;
}

int Foam::Nurbs1D::binomial
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

vector Foam::Nurbs1D::Curve_Derivative
(
    int k,
    scalar u,
    nurbsStatus state
) const
{
    if(state==prev && !nurbsMoved)
        FatalErrorInFunction<< "Prev nurbs status can only be used if Nurbs has been moved"<< exit(FatalError);
    
    if(u>=maxPara[u] || u<minPara[u])
    {
        FatalErrorInFunction
        << " Parameter u has to be within ["<<minPara[u]<<","<<maxPara[u]<<") but is "<<u<<"!"<<endl
        << exit(FatalError);
    }
    if(k > p_q[u])
    {
        FatalErrorInFunction
        << " Called the "<<k<<"-th Derivative of a "<<p_q[u]<<"-th order Nurbs. This will not work!"<<endl
        << exit(FatalError);
    }
    Info<<"-----------------------"<<endl;
    vector A_res = A(k,u,state);
    Info<<"\tA:"<<A_res;
    scalar w = Weights_B_Spline_Derivative(0,u,state);
    Info<<"  w:"<<w<<endl;
    vector B(0,0,0);
    for(int i=1;i<k+1;i++)
    {
        scalar koeff = binomial(k,i);
        scalar w_i = Weights_B_Spline_Derivative(i,u,state);
        vector C_kmi = Curve_Derivative(k-i,u,state);
        B += koeff * w_i * C_kmi;
        Info<<"\tkoeff:"<<koeff<<"  w_i:"<<w_i<<"  C_kmi:"<<C_kmi<<endl;
    }
    
    //Info<<"A_"<<k<<"_:"<<A_res<<endl;
    //Info<<"w_"<<k<<"_:"<<w<<endl;
    //Info<<"B_"<<k<<"_:"<<B<<endl;
    //Info<<"-----------------------"<<endl;
    
    if(w == 0)
        w = 1;
    
    return (A_res-B)/w;
}

Foam::BoundingBox Foam::Nurbs1D::computeBoundingBox() const
{
    const List<List<vector>>& controlPoints = this->controlPoints[current];
    if(controlPoints.size() == 0 || controlPoints[0].size() == 0)
    {
        FatalErrorInFunction
        << " Nurbs Curve has no controlPoints. Therefore no bounds can be computed!"<<endl
        << exit(FatalError);
    }
    
    BoundingBox MinMaxBox;
    MinMaxBox.Min = MinMaxBox.Max = controlPoints.first().first();
    
    for(int i=0;i<cPdim[0];i++)
    {
        for(int j=0;j<cPdim[1];j++)
        {
            for(int d=0;d<3;d++)
            {
                if(controlPoints[i][j][d] > MinMaxBox.Max[d])
                    MinMaxBox.Max[d] = controlPoints[i][j][d];
                    
                if(controlPoints[i][j][d] < MinMaxBox.Min[d])
                    MinMaxBox.Min[d] = controlPoints[i][j][d];
            }
        }
    }    
    for(int d=0;d<3;d++)
    {
        MinMaxBox.Min[d] = MinMaxBox.Min[d]-(2*deltaX+diameter);
        MinMaxBox.Max[d] = MinMaxBox.Max[d]+(2*deltaX+diameter);
    }    
    return MinMaxBox;
}

Foam::BoundingBox Foam::Nurbs1D::computeBoundingBox
(
    scalar start_u,
    scalar end_u
) const
{
    //Info<<"Compute Bounding Box"<<endl;
    if(start_u < minPara[u] || start_u >= maxPara[u])
    {
        FatalErrorInFunction
        << " Start value of Bounding Box has to be in ["<<minPara[u]<<","<<maxPara[u]<<") but is"<<start_u<<endl
        << exit(FatalError);
    }
    if(end_u < minPara[u] || end_u >= maxPara[u])
    {
        FatalErrorInFunction
        << " End value of Bounding Box has to be in ["<<minPara[u]<<","<<maxPara[u]<<") but is "<<end_u<<endl
        << exit(FatalError);
    }
    
    vector startP = Curve_Derivative(0,start_u);
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

scalar Foam::Nurbs1D::supremum_Derivative2
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
    return std::sqrt(sup);
}

vector Foam::Nurbs1D::vectorCurveToPoint_Derivative
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

scalar Foam::Nurbs1D::distCurveToPoint_Deriv0
(
    scalar para,
    vector point
) const
{
    return euklidianNorm(vectorCurveToPoint_Derivative(0,para,point));
}

scalar Foam::Nurbs1D::distCurveToPoint_Deriv1
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

scalar Foam::Nurbs1D::distCurveToPoint_Deriv2
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

scalar Foam::Nurbs1D::newtonIterateNearestNeighbour_alt
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

scalar Foam::Nurbs1D::newtonIterateNearestNeighbour
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

scalar Foam::Nurbs1D::distanceToNurbsSurface
(
    scalar para,
    vector point
) const
{
    vector C = Curve_Derivative(0,para);
    scalar distToCentre = euklidianNorm(C-point);
    return distToCentre - diameter;
}

void Foam::Nurbs1D::moveNurbs
(
    List<List<vector>> controlPoints
)
{
    if(controlPoints.size() != this->controlPoints.size())
        FatalErrorInFunction<<"Number of nurbs curves must stay the same at all times"<<exit(FatalError);
    for(int i=0;i<controlPoints.size();i++)
        if(controlPoints[i].size() != this->controlPoints[i].size())
            FatalErrorInFunction<<"Number of nurbs curves must stay the same at all times"<<exit(FatalError);
    
    knots[previous] = knots[current];
    this->controlPoints[previous] = this->controlPoints[current];
    weights[previous] = weights[current];
    weightedControlPoints[previous] = weightedControlPoints[current];
    
    this->controlPoints[current] = controlPoints;
    
    for(int i=0;i<weightedControlPoints[current].size();i++)
    {
        for(int j=0;j<weightedControlPoints[current][i].size();j++)
        {
            weightedControlPoints[current][i][j] = weights[current][i][j] * this->controlPoints[current][i][j];
        }
    }   
    nurbsMoved = true;
}

vector Foam::Nurbs1D::movementVector
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


