#include "UnitTestBsTree.H"

void Foam::TESTBSTREE::collectLeafCenter(DynamicList<scalar>& leafCenters, const BsTree::Node* const node, const BsTree& tree)
{
    if(node->left!=tree._nil)
        collectLeafCenter(leafCenters,node->left,tree);
    if(node->right!=tree._nil)
        collectLeafCenter(leafCenters,node->right,tree);
    
    if(node->left==tree._nil && node->right==tree._nil)
        leafCenters.append(node->divideBound);
    else if(node->left!=tree._nil && node->right!=tree._nil)
        return;
    else
        FatalErrorInFunction<<"Invalid BsTree"<< exit(FatalError);
}

bool Foam::TESTBSTREE::checkNearestPoints(const BsTree& tree, scalarList nearestPoints, vector point)
{
    DynamicList<scalar> leafCenters;
    collectLeafCenter(leafCenters, tree.root, tree);
    
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
    
    label nbrNearestPD = nearestPointsDist.size();
    scalar lastElementNPD = nearestPointsDist.back();
    scalar nthElementleafCenterDist = leafCentersDist[nbrNearestPD];
    if(lastElementNPD<=nthElementleafCenterDist)
        return true;
    else
        return false;   
}

bool Foam::TESTBSTREE::checkNearestNeighbor(const BsTree& tree, scalar checkNearestNeighbor, vector point)
{
    DynamicList<scalar> leafCenters;
    collectLeafCenter(leafCenters, tree.root, tree);
    std::vector<scalar> leafCentersDist;
    for(int i=0;i<leafCenters.size();i++)
    {
        vector pntToNurbs = point - tree.Curve.Curve_Derivative(0,leafCenters[i]);
        leafCentersDist.push_back(std::sqrt(pntToNurbs&&pntToNurbs));
    }
    std::sort(leafCentersDist.begin(),leafCentersDist.end());
    
    vector pntTocheckNearestNeighbor = point - tree.Curve.Curve_Derivative(0,checkNearestNeighbor);
    scalar pntTocheckNearestNeighborDist = std::sqrt(pntTocheckNearestNeighbor&&pntTocheckNearestNeighbor);
    
    if(pntTocheckNearestNeighborDist<=leafCentersDist.front())
        return true;
    else
        return false; 
}

void Foam::TESTBSTREE::nearestLeaf()
{    
    Nurbs1D QuarterCircle
    (
        {0,0,0,1,1,1},
        {vector(1,0,0),vector(1,1,0),vector(0,1,0)},
        {1,1,2},
        2
    );    
    BsTree Tree(QuarterCircle);
    
    vector testPointOutside(-2,-2,0);    
    scalarList resOut = Tree.nearestPoints(testPointOutside);
    assert(checkNearestPoints(Tree,resOut,testPointOutside));
    
    vector testPointInside(0.2,0.9,0);
    scalarList resIn = Tree.nearestPoints(testPointInside);
    assert(checkNearestPoints(Tree,resIn,testPointInside));

    scalar nearest = Tree.Curve.newtonIterateNearestNeighbour(resIn.last(),testPointInside);
    assert(checkNearestNeighbor(Tree,nearest,testPointInside));
    
    nearest = Tree.Curve.newtonIterateNearestNeighbour_alt(resIn.last(),testPointInside);
    assert(checkNearestNeighbor(Tree,nearest,testPointInside));

    nearest = Tree.closestParaOnNurbsToPoint(testPointInside);
    assert(checkNearestNeighbor(Tree,nearest,testPointInside));
    
    Nurbs1D Line
    (
        {0,0,0,1,1,1},
        {vector(1,0,0),vector(0.5,0,0.5),vector(0,0,1)},
        {1,1,1},
        2,
        0.3,
        1
    );    
    BsTree Tree2(Line);
    vector pointOnNurbs(0,0,0.6);
    //Info<<"Compute closest Parameter: ";
    scalar closestPara = Tree2.closestParaOnNurbsToPoint(pointOnNurbs);
    assert(checkNearestNeighbor(Tree,closestPara,pointOnNurbs));
}
