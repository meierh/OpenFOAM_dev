#include "Nurbs.H"
#include <math.h> 

bool vectorIsEqual(vector vec1, vector vec2, scalar errorMargin)
{
    vector diff = vec1-vec2;
    scalar dist = diff.x()*diff.x() + diff.y()*diff.y() + diff.z()*diff.z();
    dist = std::sqrt(dist);
    if(dist<errorMargin)
        return true;
    else
        return false;
}

void Foam::UnitTest_Nurbs()
{
    Foam::Info<<"NURBS"<<Foam::endl;    
    std::unique_ptr<scalarList> knots(new scalarList(6));
    (*knots)[0] = 0;    (*knots)[1] = 0;    (*knots)[2] = 0;    (*knots)[3] = 1;
    (*knots)[4] = 1;    (*knots)[5] = 1;    
    //Info<<"Knoten"<<endl;    
    int testdegree = 2;    
    std::unique_ptr<scalarList> weights(new scalarList(3));
    (*weights)[0] = 1;    (*weights)[1] = 1;    (*weights)[2] = 2;
    //Info<<"Gewichte"<<endl;    
    std::unique_ptr<List<vector>> controlPoints(new List<vector>(3));
    (*controlPoints)[0] = vector(1,0,0);    
    (*controlPoints)[1] = vector(1,1,0);    
    (*controlPoints)[2] = vector(0,1,0);
    Info<<"Kontrollpunkte"<<endl;    
    Nurbs QuarterCircle(std::move(knots),std::move(controlPoints),std::move(weights),testdegree);
    
    int correctRes = 0;
    
    vector X0_D0 = QuarterCircle.Curve_Derivative(0,0);
    vector X0_D0_soll(1,0,0);
    if(vectorIsEqual(X0_D0,X0_D0_soll,1e-2))
        correctRes++;
    else
        Info<<"C(0) Wrong: "<<X0_D0<<"!="<<X0_D0_soll<<endl;
    
    vector X1_D0 = QuarterCircle.Curve_Derivative(0,0.99999999999);
    vector X1_D0_soll(0,1,0);
    if(vectorIsEqual(X1_D0,X1_D0_soll,1e-2))
        correctRes++;
    else
        Info<<"C(1) Wrong: "<<X1_D0<<"!="<<X1_D0_soll<<endl;
    
    vector X0_D1 = QuarterCircle.Curve_Derivative(1,0);
    vector X0_D1_soll(0,2,0);
    if(vectorIsEqual(X0_D1,X0_D1_soll,1e-2))
        correctRes++;
    else
        Info<<"C'(0) Wrong: "<<X0_D1<<"!="<<X0_D1_soll<<endl;
    
    vector X1_D1 = QuarterCircle.Curve_Derivative(1,0.999);
    vector X1_D1_soll(-1,0,0);
    if(vectorIsEqual(X1_D1,X1_D1_soll,1e-2))
        correctRes++;
    else
        Info<<"C'(1) Wrong: "<<X1_D1<<"!="<<X1_D1_soll<<endl;
    
    vector X0_D2 = QuarterCircle.Curve_Derivative(2,0);
    vector X0_D2_soll(-4,0,0);
    if(vectorIsEqual(X0_D2,X0_D2_soll,1e-2))
        correctRes++;
    else
        Info<<"C''(0) Wrong: "<<X0_D2<<"!="<<X0_D2_soll<<endl;
    
    vector X1_D2 = QuarterCircle.Curve_Derivative(2,0.999);
    vector X1_D2_soll(1,-1,0);
    if(vectorIsEqual(X1_D2,X1_D2_soll,1e-2))
        correctRes++;
    else
        Info<<"C''(1) Wrong: "<<X1_D2<<"!="<<X1_D2_soll<<endl;
    
    Info<<QuarterCircle.max_U()<<" "<<(QuarterCircle.max_U()==1.0)<<endl;
    Info<<QuarterCircle.Curve_Derivative(0,QuarterCircle.max_U())<<endl;
   

    Info<<"UnitTest Nurbs Quarter Circle Done:"<<correctRes<<"/6 correct"<<endl;
    
    /*
    Foam::Info<<"Test Quarter Circle Nurbs Curve"<<Foam::endl;    
    std::unique_ptr<scalarList> knots(new scalarList(12));
    (*knots)[0] = 0;    (*knots)[1] = 0;    (*knots)[2] = 0;    (*knots)[3] = 1;
    (*knots)[4] = 1;    (*knots)[5] = 2;    (*knots)[6] = 2;    (*knots)[7] = 3;
    (*knots)[8] = 3;    (*knots)[9] = 4;    (*knots)[10] = 4;   (*knots)[11] = 4;
    Info<<"Knoten"<<endl;    
    int testdegree = 2;
    std::unique_ptr<scalarList> weights(new scalarList(9));
    (*weights)[0] = 1;    (*weights)[1] = sqrt(2)/2;    (*weights)[2] = 1;
    (*weights)[3] = sqrt(2)/2;    (*weights)[4] = 1;    (*weights)[5] = sqrt(2)/2;
    (*weights)[6] = 1;    (*weights)[7] = sqrt(2)/2;    (*weights)[8] = 1;
    Info<<"Gewichte"<<endl;
    std::unique_ptr<List<vector>> controlPoints(new List<vector>(9));
    (*controlPoints)[0] = vector(1,0,0);    (*controlPoints)[1] = vector(1,1,0);    (*controlPoints)[2] = vector(0,1,0);
    (*controlPoints)[3] = vector(-1,1,0);   (*controlPoints)[4] = vector(-1,0,0);   (*controlPoints)[5] = vector(-1,-1,0);
    (*controlPoints)[6] = vector(0,-1,0);   (*controlPoints)[7] = vector(1,-1,0);   (*controlPoints)[8] = vector(1,0,0);
    Info<<"Kontrollpunkte"<<endl;
    Nurbs Circle(std::move(knots),std::move(controlPoints),std::move(weights),testdegree);
    */
}
