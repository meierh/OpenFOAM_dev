#include "UnitTestNurbs2D.H"

void Foam::TESTNURBS2D::testControlPointDerivative() // The Nurbs Book by Prof Les Piegl S.94ff
{    
    label p=2,q=2;
    scalarList u = {0,0,0,1,2,3,4,4,5,5,5};
    scalarList v = {0,0,0,1.0/5.0,1.0/2.0,4.0/5.0,1,1,1};
    Nurbs2D testNurbs1
    (
        {u,v},
        {
         {vector(0,0,0),vector(0,1,-2),vector(0,2,-4),vector(0,3,-4),vector(0,2,-2),vector(0,3,0)},
         {vector(1,0,0),vector(1,1,-1),vector(1,2,-2),vector(1,3,-2), vector(1,2,2), vector(1,3,0)},
         {vector(2,0,0),vector(2,1,0), vector(2,2,-1),vector(2,3,0), vector(2,2,2), vector(2,3,0)}, 
         {vector(3,0,0),vector(3,1,1), vector(3,2,1), vector(3,3,2), vector(3,2,2), vector(3,3,0)},
         {vector(4,0,0),vector(4,1,1), vector(4,2,1), vector(4,3,2), vector(4,2,2), vector(4,3,0)},
         {vector(5,0,0),vector(5,1,0), vector(5,2,-1),vector(5,3,0), vector(5,2,2), vector(5,3,0)},
         {vector(6,0,0),vector(6,1,-1),vector(6,2,-2),vector(6,3,-2), vector(6,2,2), vector(6,3,0)},
         {vector(7,0,0),vector(7,1,-2),vector(7,2,-4),vector(7,3,-4),vector(7,2,-2),vector(7,3,0)}
        },
        {
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1}
        },
        p,q
    );
    List<List<vector>>& cP = testNurbs1.controlPoints[Nurbs2D::nurbsStatus::curr];
    List<List<vector>>& P = cP;
    assert(approxEqual(testNurbs1.Control_Point_Derivative<vector>({1,0},{0,0},cP),(p/u[p+1])*(P[1][0]-P[0][0])));
    assert(approxEqual(testNurbs1.Control_Point_Derivative<vector>({0,1},{0,0},cP),(q/v[q+1])*(P[0][1]-P[0][0])));    
    assert(approxEqual(testNurbs1.Control_Point_Derivative<vector>({1,1},{0,0},cP),((p*q)/(u[p+1]*v[q+1]))*(P[1][1]-P[0][1]-P[1][0]+P[0][0])));
    Info<<"testControlPointDerivative done"<<endl;
}

void Foam::TESTNURBS2D::testA()
{
    Info<<"testA ";
    label p=2,q=2;
    scalarList u = {0,0,0,1,2,3,4,4,5,5,5};
    scalarList v = {0,0,0,1.0/5.0,1.0/2.0,4.0/5.0,1,1,1};
    Nurbs2D testNurbs1
    (
        {u,v},
        {
         {vector(0,0,0),vector(0,1,-2),vector(0,2,-4),vector(0,3,-4),vector(0,2,-2),vector(0,3,0)},
         {vector(1,0,0),vector(1,1,-1),vector(1,2,-2),vector(1,3,-2), vector(1,2,2), vector(1,3,0)},
         {vector(2,0,0),vector(2,1,0), vector(2,2,-1),vector(2,3,0), vector(2,2,2), vector(2,3,0)}, 
         {vector(3,0,0),vector(3,1,1), vector(3,2,1), vector(3,3,2), vector(3,2,2), vector(3,3,0)},
         {vector(4,0,0),vector(4,1,1), vector(4,2,1), vector(4,3,2), vector(4,2,2), vector(4,3,0)},
         {vector(5,0,0),vector(5,1,0), vector(5,2,-1),vector(5,3,0), vector(5,2,2), vector(5,3,0)},
         {vector(6,0,0),vector(6,1,-1),vector(6,2,-2),vector(6,3,-2), vector(6,2,2), vector(6,3,0)},
         {vector(7,0,0),vector(7,1,-2),vector(7,2,-4),vector(7,3,-4),vector(7,2,-2),vector(7,3,0)}
        },
        {
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1}
        },
        p,q
    );
    List<List<vector>>& cP = testNurbs1.controlPoints[Nurbs2D::nurbsStatus::curr];
    List<List<vector>>& P = cP;
    assert(approxEqual(testNurbs1.A(1,0,0,0),(p/u[p+1])*(P[1][0]-P[0][0])));
    assert(approxEqual(testNurbs1.A(0,1,0,0),(q/v[q+1])*(P[0][1]-P[0][0])));    
    assert(approxEqual(testNurbs1.A(1,1,0,0),((p*q)/(u[p+1]*v[q+1]))*(P[1][1]-P[0][1]-P[1][0]+P[0][0])));
    Info<<"done"<<endl;
}

void Foam::TESTNURBS2D::testCurveDerivative() // The Nurbs Book by Prof Les Piegl S.125ff
{
    Info<<"testCurveDerivative ";
    label p=2,q=2;
    scalarList u = {0,0,0,1,2,3,4,4,5,5,5};
    scalarList v = {0,0,0,1.0/5.0,1.0/2.0,4.0/5.0,1,1,1};
    Nurbs2D testNurbs1
    (
        {u,v},
        {
         {vector(0,0,0),vector(0,1,-2),vector(0,2,-4),vector(0,3,-4),vector(0,2,-2),vector(0,3,0)},
         {vector(1,0,0),vector(1,1,-1),vector(1,2,-2),vector(1,3,-2), vector(1,2,2), vector(1,3,0)},
         {vector(2,0,0),vector(2,1,0), vector(2,2,-1),vector(2,3,0), vector(2,2,2), vector(2,3,0)}, 
         {vector(3,0,0),vector(3,1,1), vector(3,2,1), vector(3,3,2), vector(3,2,2), vector(3,3,0)},
         {vector(4,0,0),vector(4,1,1), vector(4,2,1), vector(4,3,2), vector(4,2,2), vector(4,3,0)},
         {vector(5,0,0),vector(5,1,0), vector(5,2,-1),vector(5,3,0), vector(5,2,2), vector(5,3,0)},
         {vector(6,0,0),vector(6,1,-1),vector(6,2,-2),vector(6,3,-2), vector(6,2,2), vector(6,3,0)},
         {vector(7,0,0),vector(7,1,-2),vector(7,2,-4),vector(7,3,-4),vector(7,2,-2),vector(7,3,0)}
        },
        {
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1},
         {1,1,1,1,1,1}
        },
        p,q
    );
    List<List<vector>>& cP = testNurbs1.controlPoints[Nurbs2D::nurbsStatus::curr];
    List<List<vector>>& P = cP;
    assert(approxEqual(testNurbs1.Surface_Derivative(1,0,0,0),(p/u[p+1])*(P[1][0]-P[0][0])));
    assert(approxEqual(testNurbs1.Surface_Derivative(0,1,0,0),(q/v[q+1])*(P[0][1]-P[0][0])));    
    assert(approxEqual(testNurbs1.Surface_Derivative(1,1,0,0),((p*q)/(u[p+1]*v[q+1]))*(P[1][1]-P[0][1]-P[1][0]+P[0][0])));
    Info<<"done"<<endl;
}
