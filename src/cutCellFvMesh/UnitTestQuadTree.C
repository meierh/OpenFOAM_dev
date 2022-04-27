#include "UnitTestQuadTree.H"

void Foam::TESTQUADTREE::collectLeafCenter(DynamicList<FixedList<scalar,2>>& leafCenters, const QuadTree::Node* const node, const QuadTree& tree)
{
    if(node->leftU_leftV!=tree._nil)
        collectLeafCenter(leafCenters,node->leftU_leftV,tree);
    if(node->rightU_leftV!=tree._nil)
        collectLeafCenter(leafCenters,node->rightU_leftV,tree);
    if(node->leftU_rightV!=tree._nil)
        collectLeafCenter(leafCenters,node->leftU_rightV,tree);
    if(node->rightU_rightV!=tree._nil)
        collectLeafCenter(leafCenters,node->rightU_rightV,tree);
    
    if(node->leftU_leftV==tree._nil && node->rightU_leftV==tree._nil && node->leftU_rightV==tree._nil && node->rightU_rightV==tree._nil)
        leafCenters.append(FixedList<scalar,2>({node->divideBoundU,node->divideBoundV}));
    else if(node->leftU_leftV!=tree._nil && node->rightU_leftV!=tree._nil && node->leftU_rightV!=tree._nil && node->rightU_rightV!=tree._nil)
        return;
    else
        FatalErrorInFunction<<"Invalid QuadTree"<< exit(FatalError);
}

void Foam::TESTQUADTREE::collectLeafCenterWithPointInBox(DynamicList<FixedList<scalar,2>>& leafCenters, const QuadTree::Node* const node, const QuadTree& tree, const vector point)
{
    if(node->leftU_leftV!=tree._nil && node->leftU_leftV->MinMaxBox.isInside(point))
        collectLeafCenter(leafCenters,node->leftU_leftV,tree);
    if(node->rightU_leftV!=tree._nil && node->rightU_leftV->MinMaxBox.isInside(point))
        collectLeafCenter(leafCenters,node->rightU_leftV,tree);
    if(node->leftU_rightV!=tree._nil && node->leftU_rightV->MinMaxBox.isInside(point))
        collectLeafCenter(leafCenters,node->leftU_rightV,tree);
    if(node->rightU_rightV!=tree._nil && node->rightU_rightV->MinMaxBox.isInside(point))
        collectLeafCenter(leafCenters,node->rightU_rightV,tree);
    
    if(node->leftU_leftV==tree._nil && node->rightU_leftV==tree._nil && node->leftU_rightV==tree._nil && node->rightU_rightV==tree._nil && node->MinMaxBox.isInside(point))
        leafCenters.append(FixedList<scalar,2>({node->divideBoundU,node->divideBoundV}));
    else if(node->leftU_leftV!=tree._nil && node->rightU_leftV!=tree._nil && node->leftU_rightV!=tree._nil && node->rightU_rightV!=tree._nil)
        return;
    else
        FatalErrorInFunction<<"Invalid QuadTree"<< exit(FatalError);
}

bool Foam::TESTQUADTREE::checkNearestPoints(const QuadTree& tree, DynamicList<FixedList<scalar,2>> nearestPoints, vector point)
{
    DynamicList<FixedList<scalar,2>> leafCenters;
    collectLeafCenter(leafCenters, tree.root, tree);
    
    DynamicList<FixedList<scalar,2>> leafCentersInBox;
    collectLeafCenterWithPointInBox(leafCentersInBox, tree.root, tree, point);
    
    std::vector<scalar> nearestPointsDist;
    for(int i=0;i<nearestPoints.size();i++)
    {
        vector pntToNurbs = point - tree.Surface.Surface_Derivative(0,0,nearestPoints[i][0],nearestPoints[i][1]);
        nearestPointsDist.push_back(std::sqrt(pntToNurbs&&pntToNurbs));
    }
    std::sort(nearestPointsDist.begin(),nearestPointsDist.end());
    
    std::vector<scalar> leafCentersDist;
    for(int i=0;i<leafCenters.size();i++)
    {
        vector pntToNurbs = point - tree.Surface.Surface_Derivative(0,0,leafCenters[i][0],leafCenters[i][1]);
        leafCentersDist.push_back(std::sqrt(pntToNurbs&&pntToNurbs));
    }
    std::sort(leafCentersDist.begin(),leafCentersDist.end());
    
    if(leafCentersDist.empty())
        Info<<"No nearest point"<<endl;
    else
    {
        if(!(*(nearestPointsDist.begin())<=*(leafCentersDist.begin())))
            return false;
    }    
    
    std::vector<scalar> leafCentersInBoxDist;
    for(int i=0;i<leafCentersInBox.size();i++)
    {
        vector pntToNurbs = point - tree.Surface.Surface_Derivative(0,0,leafCentersInBox[i][0],leafCentersInBox[i][1]);
        leafCentersInBoxDist.push_back(std::sqrt(pntToNurbs&&pntToNurbs));
    }
    std::sort(leafCentersInBoxDist.begin(),leafCentersInBoxDist.end());
    
    if(leafCentersInBoxDist.empty())
    {
        Info<<"No nearest point"<<endl;
        return true;
    }
    
    if(leafCentersInBoxDist.empty())
        Info<<"No nearest point"<<endl;
    else
    {
        if(!(*(nearestPointsDist.begin())<=*(leafCentersInBoxDist.begin())))
            return false;
    }
    return true;
}

bool Foam::TESTQUADTREE::checkNearestNeighbor(const QuadTree& tree, FixedList<scalar,2> nearestPoint, vector point)
{
    DynamicList<FixedList<scalar,2>> leafCenters;
    collectLeafCenter(leafCenters, tree.root, tree);
    std::vector<scalar> leafCentersDist;
    for(int i=0;i<leafCenters.size();i++)
    {
        vector pntToNurbs = point - tree.Surface.Surface_Derivative(0,0,leafCenters[i][0],leafCenters[i][1]);
        leafCentersDist.push_back(std::sqrt(pntToNurbs&&pntToNurbs));
    }
    std::sort(leafCentersDist.begin(),leafCentersDist.end());
    
    vector pntTocheckNearestNeighbor = point - tree.Surface.Surface_Derivative(0,0,nearestPoint[0],nearestPoint[1]);
    scalar pntTocheckNearestNeighborDist = std::sqrt(pntTocheckNearestNeighbor&&pntTocheckNearestNeighbor);
    
    if(pntTocheckNearestNeighborDist<=leafCentersDist.front())
        return true;
    else
    {
        
        return false;
    }
}

void Foam::TESTQUADTREE::nearestLeaf()
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
    //Info<<"Built Nurbs2D"<<endl;
    QuadTree Tree(testNurbs1);
    //Info<<"Built Quadtree"<<endl;
    
    vector testPointOutside(-2,-2,0);    
    DynamicList<FixedList<scalar,2>> resOut = Tree.nearestPoints(testPointOutside);
    assert(checkNearestPoints(Tree,resOut,testPointOutside));
    //Info<<"1"<<endl;
    
    vector testPointInside(0.2,0.9,0);
    DynamicList<FixedList<scalar,2>> resIn = Tree.nearestPoints(testPointInside);
    assert(checkNearestPoints(Tree,resIn,testPointInside));
    //Info<<"2"<<endl;

    //Info<<"resIn:"<<resIn<<endl;
    //Info<<"resIn.last():"<<resIn.last()<<endl;
    
    FixedList<scalar,2> nearestTotal = {0,0};
    scalar  nearestDistTotal = std::numeric_limits<scalar>::max();
    for(FixedList<scalar,2>& nearesPoint : resIn)
    {
        FixedList<scalar,2> nearest = Tree.Surface.newtonIterateNearestNeighbour(nearesPoint[0],nearesPoint[1],testPointInside);
        scalar nearestDist = Tree.Surface.distanceToNurbsSurface(nearest[0],nearest[1],testPointInside);
        if(nearestDistTotal>=nearestDist)
        {
            nearestDistTotal = nearestDist;
            nearestTotal = nearest;
        }
    }
    assert(checkNearestNeighbor(Tree,nearestTotal,testPointInside));
    //Info<<"3"<<endl;
}
