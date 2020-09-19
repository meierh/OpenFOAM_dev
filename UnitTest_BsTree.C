#include "Nurbs.H"
#include "BsTree.H"

void Foam::UnitTest_BsTree()
{
    Foam::Info<<"BSTREE"<<Foam::endl;    
    scalarList knots(6);
    knots[0] = 0;    knots[1] = 0;    knots[2] = 0;    knots[3] = 1;
    knots[4] = 1;    knots[5] = 1;    
    //Info<<"Knoten"<<endl;    
    int testdegree = 2;    
    scalarList weights(3);
    weights[0] = 1;    weights[1] = 1;    weights[2] = 2;
    //Info<<"Gewichte"<<endl;    
    List<vector> controlPoints(3);
    controlPoints[0] = vector(1,0,0);    
    controlPoints[1] = vector(1,1,0);    
    controlPoints[2] = vector(0,1,0);
    //Info<<"Kontrollpunkte"<<endl;
    std::shared_ptr<Nurbs> QuarterCircle(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.05,0.05)));
    
    BsTree Tree(QuarterCircle);
    
    vector testPointOutside(-2,-2,0);
    
    vector testPointInside(0.2,0.9,0);
    
    scalarList resOut = Tree.nearestPoints(testPointOutside);
    //Info<<resOut.size()<<endl;
    
    
    scalarList resIn = Tree.nearestPoints(testPointInside);
    /*
    Info<<resIn.size()<<endl;
    for(int i=0;i<resIn.size();i++)
        Info<<resIn[i]<<" "<<Tree.Curve->vectorCurveToPoint_Derivative(0,resIn[i],testPointInside)<<
        " d:"<<Tree.Curve->distCurveToPoint_Deriv0(resIn[i],testPointInside)<<
        " d_1:"<<Tree.Curve->distCurveToPoint_Deriv1(resIn[i],testPointInside)<<
        " d_2:"<<Tree.Curve->distCurveToPoint_Deriv2(resIn[i],testPointInside)<<endl;
    */
    scalar nearest = Tree.Curve->newtonIterateNearestNeighbour(resIn.last(),testPointInside);
    Info<<"UnitTest BsTree nearest Point on Nurbs found 1/3"<<endl;
    
    nearest = Tree.Curve->newtonIterateNearestNeighbour_alt(resIn.last(),testPointInside);
    Info<<"UnitTest BsTree nearest Point on Nurbs found 2/3"<<endl;

    nearest = Tree.closestParaOnNurbsToPoint(testPointInside);
    Info<<"UnitTest BsTree nearest Point on Nurbs found 3/3"<<endl;
    
    
    knots = scalarList(6);
    knots[0]=0; knots[1]=0; knots[2]=0; knots[3]=1; knots[4]=1; knots[5]=1;    
    weights = scalarList(2);
    weights[0] = 1;    weights[1] = 1;    
    controlPoints = List<vector>(2);
    controlPoints[0]=vector(0,0,0); controlPoints[1]=vector(0,0,1);    
    std::shared_ptr<Nurbs> Line(std::shared_ptr<Nurbs>(new Nurbs(knots,controlPoints,weights,testdegree,0.3,1)));
    
    BsTree Tree2(Line);
    vector pointOnNurbs(0,0,0.6);
    Info<<"Compute closest Parameter: ";
    scalar closestPara = Tree2.closestParaOnNurbsToPoint(pointOnNurbs);
    Info<<closestPara<<" "<<Tree2.Curve->Curve_Derivative(0,closestPara)<<endl;
    Info<<"Compute closest Distance ";
    scalar closestDist = Tree2.Curve->distanceToNurbsSurface(closestPara,pointOnNurbs);
    Info<<closestDist<<endl;
}
