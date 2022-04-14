#include "Nurbs2D.H"
#include <math.h> 

Foam::Nurbs2D::Nurbs2D
(
    List<scalarList> knots,
    List<List<vector>> controlPoints,
    scalarListList weights,
    int degree_p,
    int degree_q,
    scalar diameter,
    scalar deltaX
):
Nurbs1D(knots,controlPoints,weights,degree_p,degree_q,diameter,deltaX)
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
    for(int i=0;i<p_q[u]+1;i++)
    {
        if(knots[u][i]!=minPara[u])
        {
            FatalErrorInFunction
            << "First "<<p_q[u]+1<<" knots must be equal to minPara[u]:"<<minPara[u]<<" but is "<<knots[u][i]<<endl
            << exit(FatalError);
        }
    }
    for(int i=0;i<p_q[u]+1;i++)
    {
        if(knots[u][knots[u].size()-1-i]!=maxPara[u])
        {
            FatalErrorInFunction
            << "Last "<<p_q[u]+1<<" knots must be equal to maxPara[u]:"<<maxPara[u]<<" but is "<<knots[u][knots[u].size()-1-i]<<endl
            << exit(FatalError);
        }
    }
    for(int i=0;i<p_q[v]+1;i++)
    {
        if(knots[v][i]!=minPara[v])
        {
            FatalErrorInFunction
            << "First "<<p_q[v]+1<<" knots must be equal to minPara[v]:"<<minPara[v]<<" but is "<<knots[v][i]<<endl
            << exit(FatalError);
        }
    }
    for(int i=0;i<p_q[v]+1;i++)
    {
        if(knots[v][knots[v].size()-1-i]!=maxPara[v])
        {
            FatalErrorInFunction
            << "Last "<<p_q[v]+1<<" knots must be equal to maxPara[v]:"<<maxPara[v]<<" but is "<<knots[v][knots[v].size()-1-i]<<endl
            << exit(FatalError);
        }
    }

    if(n_m[u] != cPdim[u]+p_q[u]+1)
    {
        FatalErrorInFunction
        << " A Nurbs Surface must have a number of knots equal to the number of control Points plus the degree plus one"<<endl
        << " Currently there are "<<n_m[u]<<" knots and "<<cPdim[u]<<" control Points and degree "<<p_q[u]<<" in u direction"
        << exit(FatalError);
    }
    if(n_m[v] != cPdim[v]+p_q[v]+1)
    {
        FatalErrorInFunction
        << " A Nurbs Surface must have a number of knots equal to the number of control Points plus the degree plus one"<<endl
        << " Currently there are "<<n_m[v]<<" knots and "<<cPdim[v]<<" control Points and degree "<<p_q[v]<<" in v direction"
        << exit(FatalError);
    }

    //Info<<"Constructed"<<endl;
    
    this->weightedControlPoints[initial] = List<List<vector>>(cPdim[u]);
    this->weightedControlPoints[previous] = List<List<vector>>(cPdim[u]);
    this->weightedControlPoints[current] = List<List<vector>>(cPdim[u]);
    //Info<<"cPdim:"<<cPdim<<endl;
    for(int i=0;i<cPdim[u];i++)
    {
        this->weightedControlPoints[initial][i] = List<vector>(cPdim[v]);
        this->weightedControlPoints[previous][i] = List<vector>(cPdim[v]);
        this->weightedControlPoints[current][i] = List<vector>(cPdim[v]);
        for(int j=0;j<cPdim[v];j++)
        {
            //Info<<"i,j:"<<i<<","<<j<<endl;
            this->weightedControlPoints[initial][i][j] =
            this->weightedControlPoints[previous][i][j] =
            this->weightedControlPoints[current][i][j] = this->weights[current][i][j] * this->controlPoints[current][i][j];
        }
    }
    //Info<<"Done"<<endl;
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
    const scalarListList& weights = this->weights[state];
    scalar res = 0;    
    for(int i=0; i<cPdim[dimState::u]-k; i++)
    {
        for(int j=0; j<cPdim[dimState::v]-l; j++)
        {
            //Change04/10 +k and +l
            res += B_Spline_Basis(i+k,p_q[dimState::u]-k,u,dimState::u,state) *
                   B_Spline_Basis(j+l,p_q[dimState::v]-l,v,dimState::v,state) * 
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
    const List<List<vector>>& weightedControlPoints = this->weightedControlPoints[state];

    //Info<<"A("<<k<<","<<l<<","<<u<<","<<v<<")"<<endl;
    //Info<<"cPdim:"<<cPdim<<endl;
    vector res(0,0,0);
    for(int i=0; i<cPdim[dimState::u]-k; i++)
    {
        for(int j=0; j<cPdim[dimState::v]-l; j++)
        {
            //Info<<"----i,j:"<<i<<","<<j<<"/"<<cPdim[dimState::v]-l<<endl;
            //Change04/10 +k and +l
            res += B_Spline_Basis(i+k,p_q[dimState::u]-k,u,dimState::u,state) *
                   B_Spline_Basis(j+l,p_q[dimState::v]-l,v,dimState::v,state) * 
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
    nurbsStatus state
) const
{
    if(state==prev && !nurbsMoved)
        FatalErrorInFunction<< "Prev nurbs status can only be used if Nurbs has been moved"<< exit(FatalError);
    
    if(u>=maxPara[dimState::u] || u<minPara[dimState::u] || v>=maxPara[dimState::v] || v<minPara[dimState::v])
    {
        FatalErrorInFunction
        << " Parameter u:"<<u<<" and v:"<<v<<" has to be within ["<<minPara[dimState::u]<<","<<maxPara[dimState::u]<<"),["<<minPara[dimState::v]<<","<<maxPara[dimState::v]<<") !"<<endl
        << exit(FatalError);
    }
    if(k > p_q[dimState::u] || l > p_q[dimState::v])
    {
        FatalErrorInFunction
        << " Called the ("<<k<<","<<l<<")-th Derivative of a ("<<p_q[dimState::u]<<","<<p_q[dimState::v]<<")-th order Nurbs. This will not work!"<<endl
        << exit(FatalError);
    }
    
    //Info<<"Start:"<<endl;
    vector A_res = A(k,l,u,v,state);
    //Info<<"A:"<<A_res<<endl;
    scalar w = Weights_B_Spline_Derivative(0,0,u,v,state);
    //Info<<"A:"<<A_res<<" w:"<<w<<endl;

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

Foam::BoundingBox Foam::Nurbs2D::computeBoundingBox
(
    scalar start_u,
    scalar start_v,
    scalar end_u,
    scalar end_v
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
    if(start_v < minPara[v] || start_v >= maxPara[v])
    {
        FatalErrorInFunction
        << " Start value of Bounding Box has to be in ["<<minPara[v]<<","<<maxPara[v]<<") but is "<<start_v<<endl
        << exit(FatalError);
    }
    if(end_v < minPara[v] || end_v >= maxPara[v])
    {
        FatalErrorInFunction
        << " End value of Bounding Box has to be in ["<<minPara[v]<<","<<maxPara[v]<<") but is "<<end_v<<endl
        << exit(FatalError);
    }
    
    FatalErrorInFunction<< " Implementation not complete"<< exit(FatalError);
    
    vector startP = Surface_Derivative(0,0,start_u,start_v);
    //Info<<"StartP: "<<startP<<endl;
    vector endP = Surface_Derivative(0,0,end_u,end_v);
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
    scalar koeff = (1.0/8.0)*(end_u-start_u)*(end_u-start_u)*supremum_Derivative2(start_u,start_v,end_u,end_v);
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

            vector maxVec(0,0,0);
            for(int d=0;d<3;d++)
            {
                maxVec[d] = std::max({std::abs(du2[d]),std::abs(dv2[d]),std::abs(dudv[d])});
                maxD2[d] = std::max(maxVec[d],maxD2[d]);
            }    
        }  
    }
    scalar sup = 0;
    for(int d=0;d<3;d++)
    {
        sup += maxD2[d];
    }
    return sup;
}

/*
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
*/

/*
scalar Foam::Nurbs::distCurveToPoint_Deriv0
(
    scalar para,
    vector point
) const
{
    return euklidianNorm(vectorCurveToPoint_Derivative(0,para,point));
}
*/

/*
scalar Foam::Nurbs::distCurveToPoint_Deriv1
(
    scalar para,
    vector point
) const
{
    vector C = Curve_Derivative(0,para);
    vector C_D1 = Curve_Derivative(1,para);
    vector C_P = C-point;
    
    return (C_P && C_D1)/std::sqrt(C_P && C_P);
}
*/

/*
scalar Foam::Nurbs::distCurveToPoint_Deriv2
(
    scalar para,
    vector point
) const
{    
    vector C = Curve_Derivative(0,para);
    vector C_D1 = Curve_Derivative(1,para);
    vector C_D2 = Curve_Derivative(2,para);
    vector C_P = C-point;
    
    return ((C_D2 && C_P)+(C_D1 && C_D1))/(C_P && C_P);    
}
*/

FixedList<scalar,2> Foam::Nurbs2D::newtonIterateNearestNeighbour
(
    scalar u_0,
    scalar v_0,
    vector point
) const
{
    vector S = Surface_Derivative(0,0,u_0,v_0);
    
    vector S_D1u = Surface_Derivative(1,0,u_0,v_0);
    vector S_D1v = Surface_Derivative(0,1,u_0,v_0);
    
    vector S_D2u2 = Surface_Derivative(2,0,u_0,v_0);
    vector S_D2v2 = Surface_Derivative(0,2,u_0,v_0);
    vector S_D2uv = Surface_Derivative(1,1,u_0,v_0);
    
    scalar f1 = S_D1u && (S-point);
    scalar f2 = S_D1v && (S-point);
    
    scalar f1_u = (S_D2u2 && (S-point)) + (S_D1u && S_D1u);
    scalar f1_v,f2_u;
    f1_v = f2_u = (S_D2uv && (S-point)) + (S_D1u && S_D1v);
    scalar f2_v = (S_D2v2 && (S-point)) + (S_D1v && S_D1v);
    
    scalar epsilon = 1e-10;
    int iterations = 0;
    int maxIterations = 100;
    //Info<<"Point: "<<point<<endl;
    //Info<<"It:0 f("<<u_0<<")="<<f<<endl;
    scalar min_U = this->min_U();
    scalar max_U = this->max_U();
    scalar min_V = this->min_V();
    scalar max_V = this->max_V();
    int hitMaxOrMinCounter = 0;
    while((abs(f1)+abs(f2)) > epsilon)
    {
        iterations++;
        
        S = Surface_Derivative(0,0,u_0,v_0);
    
        S_D1u = Surface_Derivative(1,0,u_0,v_0);
        S_D1v = Surface_Derivative(0,1,u_0,v_0);
    
        S_D2u2 = Surface_Derivative(2,0,u_0,v_0);
        S_D2v2 = Surface_Derivative(0,2,u_0,v_0);
        S_D2uv = Surface_Derivative(1,1,u_0,v_0);
    
        f1 = S_D1u && (S-point);
        f2 = S_D1v && (S-point);
    
        f1_u = (S_D2u2 && (S-point)) + (S_D1u && S_D1u);
        f1_v = f2_u = (S_D2uv && (S-point)) + (S_D1u && S_D1v);
        f2_v = (S_D2v2 && (S-point)) + (S_D1v && S_D1v);

        //Test for singularity and maxIteration
        scalar detJ = f1_u*f2_v-f1_v*f2_u;
        //Prepare for cramer rule
        scalar detJ_1 = (-f1)*f2_v-f1_v*(-f2);
        scalar detJ_2 = f1_u*(-f2)-(-f1)*f2_u;        
        
        if((f1_u*f2_v-f1_v*f2_u) == 0 || iterations > maxIterations)
            break;
        
        u_0 = u_0 + detJ_1/detJ;
        v_0 = v_0 + detJ_2/detJ;
        
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
        if(v_0 > max_V)
        {
            v_0 = max_V;
            //Info<<"Break: f("<<u_0<<")"<<endl;
            if(hitMaxOrMinCounter >= 2)
                break;
            hitMaxOrMinCounter++;
        }
        if(v_0 < min_V)
        {
            v_0 = min_V;
            //Info<<"Break: f("<<u_0<<")"<<endl;
            if(hitMaxOrMinCounter >= 2)
                break;
            hitMaxOrMinCounter++;
        }
    }
    //Info<<"Res U: "<<u_0<<endl;
    return FixedList<scalar,2>{u_0,v_0};
}

scalar Foam::Nurbs2D::distanceToNurbsSurface
(
    scalar paraU,
    scalar paraV,
    vector point
) const
{
    vector S = Surface_Derivative(0,0,paraU,paraV);
    scalar distToCentre = euklidianNorm(S-point);
    return distToCentre - diameter;
}

vector Foam::Nurbs2D::movementVector
(
    scalar u,
    scalar v
)
{
    vector S_curr = Surface_Derivative(0,0,u,v,curr);
    vector S_prev = Surface_Derivative(0,0,u,v,prev);           
    return S_curr-S_prev;
}


