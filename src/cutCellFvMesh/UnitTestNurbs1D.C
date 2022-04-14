#include "UnitTestNurbs1D.H"

bool Foam::TESTNURBS1D::approxEqual(vector actual, vector target, scalar percErrorMargin)
{
    vector diff = actual-target;
    scalar diffNorm = std::sqrt(diff.x()*diff.x() + diff.y()*diff.y() + diff.z()*diff.z());
    scalar targetNorm = std::sqrt(target.x()*target.x() + target.y()*target.y() + target.z()*target.z());
    if(targetNorm==0.0)
    {
        if(diffNorm < percErrorMargin)
            return true;
        else
        {
            Info<<"diff:"<<diff<<endl;
            Info<<"diffNorm:"<<diffNorm<<endl;
            Info<<"Fail:  target:"<<target<<"\t"<<"actual:"<<actual<<endl;
            return false;
        }
    }
    scalar perc = diffNorm/targetNorm;
    bool eq = perc<percErrorMargin;
    if(eq) return true;
    else
    {
        Info<<"diffNorm:"<<diffNorm<<endl;
        Info<<"targetNorm:"<<targetNorm<<endl;
        Info<<"perc:"<<perc<<endl;
        Info<<"diff:"<<diff<<endl;
        Info<<"diffNorm:"<<diffNorm<<endl;
        Info<<"Fail:  target:"<<target<<"\t"<<"actual:"<<actual<<endl;
        return false;
    }
}

bool Foam::TESTNURBS1D::approxEqual(scalar actual, scalar target, scalar percErrorMargin)
{
    scalar diff = actual-target;
    if(target==0.0)
    {
        if(std::abs(diff) < percErrorMargin)
            return true;
        else
        {
            Info<<"Fail:  target:"<<target<<"\t"<<"actual:"<<actual<<endl;
            return false;
        }
    }
    scalar perc = std::abs(diff)/std::abs(target);
    bool eq = perc<percErrorMargin;
    if(eq) return true;
    else
    {
        Info<<"Fail:  target:"<<target<<"\t"<<"actual:"<<actual<<endl;
        return false;
    }
}

void Foam::TESTNURBS1D::testBSplineBasis() // The Nurbs Book by Prof Les Piegl S.52 Ex2.1 & Ex2.2
{
    Info<<"testBSplineBasis";
    Nurbs1D testNurbs1({0,0,0,1,1,1},{vector(1,0,0),vector(1,1,0),vector(0,1,0)},{1,1,1},2);
    Nurbs1D testNurbs2({0,0,0,1,2,3,4,4,5,5,5},{vector(1,0,0),vector(1,1,0),vector(0,1,0),vector(0,1,1),vector(0,1,1),vector(0,1,1),vector(0,1,1),vector(0,1,1)},{1,1,1,1,1,1,1,1},2);

    
    std::random_device devInf;
    std::mt19937 rngInf(devInf());
    std::uniform_real_distribution<> distUinf(std::numeric_limits<scalar>::min(),std::numeric_limits<scalar>::max());
    
    std::random_device dev01;
    std::mt19937 rng01(dev01());
    std::uniform_real_distribution<> distU01(0,1);
    
    std::random_device dev12;
    std::mt19937 rng12(dev12());
    std::uniform_real_distribution<> distU12(1,2);
    
    std::random_device dev23;
    std::mt19937 rng23(dev23());
    std::uniform_real_distribution<> distU23(2,3);
    
    std::random_device dev34;
    std::mt19937 rng34(dev34());
    std::uniform_real_distribution<> distU34(3,4);
    
    std::random_device dev45;
    std::mt19937 rng45(dev45());
    std::uniform_real_distribution<> distU45(4,5);
    
    label nbrTestRuns = 100;
    scalar u=0;
    for(int i=0;i<nbrTestRuns;i++) //Ex2.1
    {  
        //N0,0 N1,0
        assert(approxEqual(testNurbs1.B_Spline_Basis(0,0,distUinf(rngInf),Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs1.B_Spline_Basis(1,0,distUinf(rngInf),Nurbs1D::dimState::u),0));
        
        //N2,0
        assert(approxEqual(testNurbs1.B_Spline_Basis(2,0,distU01(rng01),Nurbs1D::dimState::u),1));
        assert(approxEqual(testNurbs1.B_Spline_Basis(2,0,-1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs1.B_Spline_Basis(2,0,2,Nurbs1D::dimState::u),0));

        //N3,0 N4,0
        assert(approxEqual(testNurbs1.B_Spline_Basis(3,0,distUinf(rngInf),Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs1.B_Spline_Basis(4,0,distUinf(rngInf),Nurbs1D::dimState::u),0));
        
        //N0,1
        assert(approxEqual(testNurbs1.B_Spline_Basis(0,1,distUinf(rngInf),Nurbs1D::dimState::u),0));
        
        //N1,1
        u = distU01(rng01);
        assert(approxEqual(testNurbs1.B_Spline_Basis(1,1,u,Nurbs1D::dimState::u),1-u));
        assert(approxEqual(testNurbs1.B_Spline_Basis(1,1,-1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs1.B_Spline_Basis(1,1,2,Nurbs1D::dimState::u),0));
        
        //N2,1
        u = distU01(rng01);
        assert(approxEqual(testNurbs1.B_Spline_Basis(2,1,u,Nurbs1D::dimState::u),u));
        assert(approxEqual(testNurbs1.B_Spline_Basis(2,1,-1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs1.B_Spline_Basis(2,1,2,Nurbs1D::dimState::u),0));

        //N3,1
        assert(approxEqual(testNurbs1.B_Spline_Basis(3,1,distUinf(rngInf),Nurbs1D::dimState::u),0));
        
        //N0,2
        u = distU01(rng01);
        assert(approxEqual(testNurbs1.B_Spline_Basis(0,2,u,Nurbs1D::dimState::u),(1-u)*(1-u)));
        assert(approxEqual(testNurbs1.B_Spline_Basis(0,2,-1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs1.B_Spline_Basis(0,2,2,Nurbs1D::dimState::u),0));
        
        //N1,2
        u = distU01(rng01);
        assert(approxEqual(testNurbs1.B_Spline_Basis(1,2,u,Nurbs1D::dimState::u),2*u*(1-u)));
        assert(approxEqual(testNurbs1.B_Spline_Basis(1,2,-1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs1.B_Spline_Basis(1,2,2,Nurbs1D::dimState::u),0));
        
        //N2,2
        u = distU01(rng01);
        assert(approxEqual(testNurbs1.B_Spline_Basis(2,2,u,Nurbs1D::dimState::u),u*u));
        assert(approxEqual(testNurbs1.B_Spline_Basis(2,2,-1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs1.B_Spline_Basis(2,2,2,Nurbs1D::dimState::u),0));       
    }
    for(int i=0;i<nbrTestRuns;i++) //Ex2.2
    {  
        //N0,0 N1,0
        assert(approxEqual(testNurbs2.B_Spline_Basis(0,0,distUinf(rngInf),Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(1,0,distUinf(rngInf),Nurbs1D::dimState::u),0));
        
        //N2,0
        assert(approxEqual(testNurbs2.B_Spline_Basis(2,0,distU01(rng01),Nurbs1D::dimState::u),1));
        assert(approxEqual(testNurbs2.B_Spline_Basis(2,0,-1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(2,0,2,Nurbs1D::dimState::u),0));

        //N3,0
        assert(approxEqual(testNurbs2.B_Spline_Basis(3,0,distU12(rng12),Nurbs1D::dimState::u),1));
        assert(approxEqual(testNurbs2.B_Spline_Basis(3,0,0,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(3,0,3,Nurbs1D::dimState::u),0));
        
        //N4,0
        assert(approxEqual(testNurbs2.B_Spline_Basis(4,0,distU23(rng23),Nurbs1D::dimState::u),1));
        assert(approxEqual(testNurbs2.B_Spline_Basis(4,0,1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(4,0,4,Nurbs1D::dimState::u),0));
        
        //N5,0
        assert(approxEqual(testNurbs2.B_Spline_Basis(5,0,distU34(rng34),Nurbs1D::dimState::u),1));
        assert(approxEqual(testNurbs2.B_Spline_Basis(5,0,2,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(5,0,5,Nurbs1D::dimState::u),0));
        
        //N6,0
        assert(approxEqual(testNurbs2.B_Spline_Basis(6,0,distUinf(rngInf),Nurbs1D::dimState::u),0));
        
        //N7,0
        assert(approxEqual(testNurbs2.B_Spline_Basis(7,0,distU45(rng45),Nurbs1D::dimState::u),1));
        assert(approxEqual(testNurbs2.B_Spline_Basis(7,0,3,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(7,0,6,Nurbs1D::dimState::u),0));

        //N8,0 N9,0
        assert(approxEqual(testNurbs2.B_Spline_Basis(8,0,distUinf(rngInf),Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(9,0,distUinf(rngInf),Nurbs1D::dimState::u),0));
        
        
        //N0,1
        assert(approxEqual(testNurbs2.B_Spline_Basis(0,1,distUinf(rngInf),Nurbs1D::dimState::u),0));
        
        //N1,1
        u = distU01(rng01);
        assert(approxEqual(testNurbs2.B_Spline_Basis(1,1,u,Nurbs1D::dimState::u),1-u));
        assert(approxEqual(testNurbs2.B_Spline_Basis(1,1,-1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(1,1,2,Nurbs1D::dimState::u),0));

        //N2,1
        u = distU01(rng01);
        assert(approxEqual(testNurbs2.B_Spline_Basis(2,1,u,Nurbs1D::dimState::u),u));
        u = distU12(rng12);
        assert(approxEqual(testNurbs2.B_Spline_Basis(2,1,u,Nurbs1D::dimState::u),2-u));        assert(approxEqual(testNurbs2.B_Spline_Basis(2,1,-1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(2,1,3,Nurbs1D::dimState::u),0));
        
        //N3,1
        u = distU12(rng12);
        assert(approxEqual(testNurbs2.B_Spline_Basis(3,1,u,Nurbs1D::dimState::u),u-1));
        u = distU23(rng23);
        assert(approxEqual(testNurbs2.B_Spline_Basis(3,1,u,Nurbs1D::dimState::u),3-u));        assert(approxEqual(testNurbs2.B_Spline_Basis(3,1,0,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(3,1,4,Nurbs1D::dimState::u),0));
        
        //N4,1
        u = distU23(rng23);
        assert(approxEqual(testNurbs2.B_Spline_Basis(4,1,u,Nurbs1D::dimState::u),u-2));
        u = distU34(rng34);
        assert(approxEqual(testNurbs2.B_Spline_Basis(4,1,u,Nurbs1D::dimState::u),4-u));        assert(approxEqual(testNurbs2.B_Spline_Basis(4,1,1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(4,1,5,Nurbs1D::dimState::u),0));
        
        //N5,1
        u = distU34(rng34);
        assert(approxEqual(testNurbs2.B_Spline_Basis(5,1,u,Nurbs1D::dimState::u),u-3));      assert(approxEqual(testNurbs2.B_Spline_Basis(5,1,1,Nurbs1D::dimState::u),0));
        assert(approxEqual(testNurbs2.B_Spline_Basis(5,1,5,Nurbs1D::dimState::u),0));
        
        //N2,2
        u = distU01(rng01);
        assert(approxEqual(testNurbs2.B_Spline_Basis(2,2,u,Nurbs1D::dimState::u),0.5*u*u));
        u = distU12(rng12);
        assert(approxEqual(testNurbs2.B_Spline_Basis(2,2,u,Nurbs1D::dimState::u),-3.0/2.0+3*u-u*u));
        u = distU23(rng23);
        assert(approxEqual(testNurbs2.B_Spline_Basis(2,2,u,Nurbs1D::dimState::u),0.5*(3-u)*(3-u)));
        
        //N3,2
        u = distU12(rng12);
        assert(approxEqual(testNurbs2.B_Spline_Basis(3,2,u,Nurbs1D::dimState::u),0.5*(u-1)*(u-1)));
        u = distU23(rng23);
        assert(approxEqual(testNurbs2.B_Spline_Basis(3,2,u,Nurbs1D::dimState::u),-11.0/2.0+5*u-u*u));
        u = distU34(rng34);
        assert(approxEqual(testNurbs2.B_Spline_Basis(3,2,u,Nurbs1D::dimState::u),0.5*(4-u)*(4-u)));
        
        //N4,2
        u = distU23(rng23);
        assert(approxEqual(testNurbs2.B_Spline_Basis(4,2,u,Nurbs1D::dimState::u),0.5*(u-2)*(u-2)));
        u = distU34(rng34);
        assert(approxEqual(testNurbs2.B_Spline_Basis(4,2,u,Nurbs1D::dimState::u),-16+10*u-u*u*3.0/2.0));
        
        //N5,2
        u = distU34(rng34);
        assert(approxEqual(testNurbs2.B_Spline_Basis(5,2,u,Nurbs1D::dimState::u),(u-3)*(u-3)));
        u = distU45(rng45);
        assert(approxEqual(testNurbs2.B_Spline_Basis(5,2,u,Nurbs1D::dimState::u),(5-u)*(5-u)));
        
        //N6,2
        u = distU45(rng45);
        assert(approxEqual(testNurbs2.B_Spline_Basis(6,2,u,Nurbs1D::dimState::u),2*(u-4)*(5-u)));
        
        //N7,2
        u = distU45(rng45);
        assert(approxEqual(testNurbs2.B_Spline_Basis(7,2,u,Nurbs1D::dimState::u),(u-4)*(u-4)));
    }
    Info<<" done"<<endl;
}

void Foam::TESTNURBS1D::testControlPointDerivative() // The Nurbs Book by Prof Les Piegl S.94ff
{
    Info<<"testControlPointDerivative ";
    Nurbs1D testNurbs1
    (
        {0,0,0,2.0/5.0,3.0/5.0,1,1,1},
        {vector(1,0,0),vector(0,1,0),vector(0,0,1),vector(-10,100,0),vector(0,10,-100)},
        {1,1,1,1,1},
        2
    );
    List<List<vector>>& cP = testNurbs1.controlPoints[Nurbs1D::nurbsStatus::curr];
    List<vector>& P = cP[0];    
    assert(approxEqual(testNurbs1.Control_Point_Derivative<vector>({1},{0},cP),5*(P[1]-P[0])));
    assert(approxEqual(testNurbs1.Control_Point_Derivative<vector>({1},{1},cP),(10.0/3.0)*(P[2]-P[1])));    
    assert(approxEqual(testNurbs1.Control_Point_Derivative<vector>({1},{2},cP),(10.0/3.0)*(P[3]-P[2])));    
    assert(approxEqual(testNurbs1.Control_Point_Derivative<vector>({1},{3},cP),5*(P[4]-P[3])));
    
    Nurbs1D testNurbs2
    (
        //Error in book corrected knots here
        {0,0,0,0,1.0/3.0,2.0/3.0,2.0/3.0,1,1,1,1},
        {vector(1,0,0),vector(0,1,0),vector(0,0,1),vector(-10,100,0),vector(0,10,-100),vector(-100,0,10),vector(10,-100,1000)},
        {1,1,1,1,1,1,1},
        3
    );
    cP = testNurbs2.controlPoints[Nurbs1D::nurbsStatus::curr];
    assert(approxEqual(testNurbs2.Control_Point_Derivative<vector>({1},{0},cP),9*(P[1]-P[0])));
    assert(approxEqual(testNurbs2.Control_Point_Derivative<vector>({1},{1},cP),(9.0/2.0)*(P[2]-P[1])));    
    assert(approxEqual(testNurbs2.Control_Point_Derivative<vector>({1},{2},cP),(9.0/2.0)*(P[3]-P[2])));    
    assert(approxEqual(testNurbs2.Control_Point_Derivative<vector>({1},{3},cP),(9.0/2.0)*(P[4]-P[3])));
    assert(approxEqual(testNurbs2.Control_Point_Derivative<vector>({1},{4},cP),9*(P[5]-P[4])));
    assert(approxEqual(testNurbs2.Control_Point_Derivative<vector>({1},{5},cP),9*(P[6]-P[5])));
    
    Info<<"done"<<endl;
}

void Foam::TESTNURBS1D::testA()
{
    Info<<"testA";
    Nurbs1D testNurbs1
    (
        {0,0,0,0,2.0/5.0,3.0/5.0,3.0/5.0,1,1,1,1},
        {vector(1,0,0),vector(0,1,0),vector(0,0,1),vector(-10,100,0),vector(0,10,-100),vector(-100,0,10),vector(10,-100,1000)},
        {1,1,1,1,1,1,1},
        3
    );
    List<List<vector>>& cP = testNurbs1.controlPoints[Nurbs1D::nurbsStatus::curr];
    List<vector>& P = cP[0];    
    assert(approxEqual(testNurbs1.A(1,0),(3.0/(2.0/5.0))*(P[1]-P[0])));
    assert(approxEqual(testNurbs1.A(1,std::nexttoward(1,0)),(3.0/(1-(3.0/5.0)))*(P[6]-P[5])));
    
    vector C_2_0 = (6/(2.0/5.0))*(P[0]/(2.0/5.0) - P[1]/(6.0/25.0) + P[2]/(3.0/5.0));
    assert(approxEqual(testNurbs1.A(2,0),C_2_0));
    
    vector C_2_1 = (6/(1-3.0/5.0))*(P[6]/(1-3.0/5.0) - (4.0/5.0)*P[5]/((1-3.0/5.0)*(1-3.0/5.0)) + P[4]/(1-3.0/5.0));
    assert(approxEqual(testNurbs1.A(2,std::nexttoward(1,0)),C_2_1));
    
    Info<<" done"<<endl;
}

void Foam::TESTNURBS1D::testCurveDerivative() // The Nurbs Book by Prof Les Piegl S.125ff
{
    Info<<"testCurveDerivative";
    Nurbs1D testNurbs1
    (
        {0,0,0,1,2,3,3,3},
        {vector(0,0,0),vector(1,1,0),vector(3,2,0),vector(4,1,0),vector(5,-1,0)},
        {1,4,1,1,1},
        2
    );
    assert(approxEqual(testNurbs1.Curve_Derivative(0,1),vector(7.0/5.0,6.0/5.0,0)));

    Nurbs1D testNurbs2
    (
        {0,0,0,1,1,1},
        {vector(1,0,0),vector(1,1,0),vector(0,1,0)},
        {1,1,2},
        2
    );
    assert(approxEqual(testNurbs2.Curve_Derivative(1,0),vector(0,2,0)));
    assert(approxEqual(testNurbs2.Curve_Derivative(1,std::nexttoward(1,0)),vector(-1,0,0)));
    assert(approxEqual(testNurbs2.Curve_Derivative(2,0),vector(-4,0,0)));
    
    Nurbs1D testNurbs3
    (
        {0,0,0,0.4,0.6,1,1,1},
        {vector(0,0,0),vector(1,1,0),vector(3,2,0),vector(4,1,0),vector(5,-1,0)},
        {1,4,1,1,1},
        2
    );
    List<List<vector>>& cP = testNurbs3.controlPoints[Nurbs1D::nurbsStatus::curr];
    List<vector>& P = cP[0];    
    assert(approxEqual(testNurbs3.Curve_Derivative(1,0),(2.0/0.4)*4*(P[1]-P[0])));
    assert(approxEqual(testNurbs3.Curve_Derivative(1,std::nexttoward(1,0)),(2.0/0.4)*(P[4]-P[3])));
    
    Info<<" done"<<endl;
}
