#include "Nurbs.H"
#include "BsTree.H"

void Foam::UnitTest_BsTree()
{
    Foam::Info<<"BSTREE"<<Foam::endl;    
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
    std::unique_ptr<Nurbs> QuarterCircle(new Nurbs(
        std::move(knots),std::move(controlPoints),std::move(weights),testdegree,0.05,0.05));
    
    BsTree Tree(std::move(QuarterCircle),1);
    
    vector testPointOutside(-2,-2,0);
    
    vector testPointInside(0.75,0.8,0);
    
    scalarList resOut = Tree.nearestPoint(testPointOutside);
    Info<<resOut.size()<<endl;
    
    scalarList resIn = Tree.nearestPoint(testPointInside);
    Info<<resIn.size()<<endl;
    for(int i=0;i<resIn.size();i++)
        Info<<resIn[i]<<" "<<Tree.Curve->vectorCurveToPoint_Derivative(0,resIn[i],testPointInside)<<
        " d:"<<Tree.Curve->distCurveToPoint_Deriv0(resIn[i],testPointInside)<<
        " d_1:"<<Tree.Curve->distCurveToPoint_Deriv1(resIn[i],testPointInside)<<
        " d_2:"<<Tree.Curve->distCurveToPoint_Deriv2(resIn[i],testPointInside)<<endl;
    
    Info<<"-----------------Nearest-------------------"<<endl;
    scalar nearest = Tree.Curve->newtonIterateNearestNeighbour(resIn.last(),testPointInside);
    Info<<"Result:"<<nearest<<endl;
    
    Info<<"-----------------Nearest-------------------"<<endl;
    nearest = Tree.Curve->newtonIterateNearestNeighbour_alt(resIn.last(),testPointInside);
    Info<<"Result:"<<nearest<<endl;

    Info<<"-----------------Nearest-------------------"<<endl;
    nearest = Tree.distanceToNurbs(testPointInside);
    Info<<"Result:"<<nearest<<endl;
}
