#include "UnitTestNurbs2D.H"

DynamicList<scalar> Foam::TESTBSTREE::collectLeafCenter(DynamicList<scalar>& leafCenters, const Node* const node)
{
    if(node->left!=_nil)
        collectLeafCenter(leafCenters,node->left);
    else
        leafCenters.append(node->divideBound);
    if(node->right!=_nil)
        collectLeafCenter(leafCenters,node->right);
    else
        leafCenters.append(node->divideBound);
}

DynamicList<scalar> Foam::TESTBSTREE::checkNearestPoints(const BsTree& tree, scalarList nearestPoints, vector point)
{
    DynamicList<scalar>& leafCenters;
    collectLeafCenter(leafCenters, tree.node)
    
    std::vector<scalar> nearestPointsDist;
    for(int i=0;i<nearestPoints.size();i++)
    {
        vector pntToNurbs = point - tree.Curve.Curve_Derivative(0,nearestPoints[i]);
        nearestPointsDist.push_back(std::sqrt(pntToNurbs&&pntToNurbs));
    }
    std::sort(nearestPointsDist.begin(),nearestPointsDist.end());
    
    std::vector<scalar> leafCentersDist;
    for(int i=0;i<leafCenters.size();i++)
    {
        vector pntToNurbs = point - tree.Curve.Curve_Derivative(0,leafCenters[i]);
        leafCentersDist.push_back(std::sqrt(pntToNurbs&&pntToNurbs));
    }
    std::sort(leafCentersDist.begin(),leafCentersDist.end());
    
    if(nearestPointsDist.empty())
    {
        Info<<"No nearest point"<<endl;
        return true;
    }
    // Cont here
    
}

void Foam::TESTBSTREE::nearestLeaf()
{    
    Foam::Info<<"BSTREE"<<Foam::endl;
    Nurbs1D QuarterCircle
    (
        {0,0,0,1,1,1},
        {vector(1,0,0),vector(1,1,0),vector(0,1,0)},
        {1,1,2},
        2
    );    
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
    scalar nearest = Tree.Curve.newtonIterateNearestNeighbour(resIn.last(),testPointInside);
    Info<<"UnitTest BsTree nearest Point on Nurbs found 1/5"<<endl;
    
    nearest = Tree.Curve.newtonIterateNearestNeighbour_alt(resIn.last(),testPointInside);
    Info<<"UnitTest BsTree nearest Point on Nurbs found 2/5"<<endl;

    nearest = Tree.closestParaOnNurbsToPoint(testPointInside);
    Info<<"UnitTest BsTree nearest Point on Nurbs found 3/5"<<endl;
    
    
    knots = scalarList(6);
    knots[0]=0; knots[1]=0; knots[2]=0; knots[3]=1; knots[4]=1; knots[5]=1;    
    weights = scalarList(2);
    weights[0] = 1;    weights[1] = 1;    
    controlPoints = List<vector>(2);
    controlPoints[0]=vector(0,0,0); controlPoints[1]=vector(0,0,1);    
    Nurbs1D Line(knots,controlPoints,weights,testdegree,0.3,1);
    
    BsTree Tree2(Line);
    vector pointOnNurbs(0,0,0.6);
    //Info<<"Compute closest Parameter: ";
    scalar closestPara = Tree2.closestParaOnNurbsToPoint(pointOnNurbs);
    //Info<<closestPara<<" "<<Tree2.Curve->Curve_Derivative(0,closestPara)<<endl;
    Info<<"UnitTest BsTree nearest Point on Nurbs found 4/5"<<endl;
    //Info<<"Compute closest Distance ";
    scalar closestDist = Tree2.Curve.distanceToNurbsSurface(closestPara,pointOnNurbs);
    Info<<"UnitTest BsTree nearest Point on Nurbs found 5/5"<<endl;
    //Info<<closestDist<<endl;
}
