#include "QuadTree.H"

Foam::QuadTree::QuadTree
(
    Nurbs2D Surface,
    label nbrSplitsBetweenCPs
):
BsTree(std::max(nbrSplitsBetweenCPs*Surface.nbrKnotsU()*Surface.degreeU(),nbrSplitsBetweenCPs*Surface.nbrKnotsV()*Surface.degreeV())),
Surface(Surface)
{
    Info<<"Construction Tree start"<<endl;
    _nil =  new Node();
    root = newNode(_nil,this->Surface.min_U(),this->Surface.max_U(),this->Surface.min_V(),this->Surface.max_V());
    Info<<"Inter"<<endl;
    constructTree(root);
    Info<<"Construction Tree done"<<endl;
}

/*
Foam::QuadTree::~QuadTree()
{
    recursiveNodeDeleter(root);
    delete _nil;
}
*/

void Foam::QuadTree::recursiveNodeDeleter
(
    Node*
    thisNode
)
{
    if(thisNode->leftU_leftV != _nil)
        recursiveNodeDeleter(thisNode->leftU_leftV);
    if(thisNode->rightU_leftV != _nil)
        recursiveNodeDeleter(thisNode->rightU_leftV);
    if(thisNode->leftU_rightV != _nil)
        recursiveNodeDeleter(thisNode->leftU_rightV);
    if(thisNode->rightU_rightV != _nil)
        recursiveNodeDeleter(thisNode->rightU_rightV);
    delete thisNode;     
}

Foam::QuadTree::Node *Foam::QuadTree::newNode
(
    Node* parent,
    scalar minU,
    scalar maxU,
    scalar minV,
    scalar maxV
)
{
    Node* newNodeItem = new Node();

    newNodeItem->leftU_leftV = _nil;
    newNodeItem->rightU_leftV = _nil;
    newNodeItem->leftU_rightV = _nil;
    newNodeItem->rightU_rightV = _nil;
    newNodeItem->parent = parent;
    newNodeItem->minU = minU;
    newNodeItem->maxU = maxU;
    newNodeItem->minV = minV;
    newNodeItem->maxV = maxV;
    newNodeItem->divideBoundU = 0.5*maxU+0.5*minU;
    newNodeItem->divideBoundV = 0.5*maxV+0.5*minV;

    return newNodeItem;
}

void Foam::QuadTree::constructTree
(
    Node* thisNode,
    int height
)
{
    Info<<"Construct Tree at "<<height<<"/"<<maxHeight<<endl;
    thisNode->MinMaxBox = Surface.computeBoundingBox(thisNode->minU,thisNode->minV,thisNode->maxU,thisNode->maxV);

    Info<<"Computed Box"<<endl;
    scalar maxWidth = 0;
    for(int d=0;d<3;d++)
    {
        maxWidth = std::max(thisNode->MinMaxBox.Max[d]-thisNode->MinMaxBox.Min[d],maxWidth);
    }
    Info<<"Computed maxWidth "<<maxWidth<<endl;
    maxWidth = maxWidth-2*Surface.getBoundingBoxOverhang();
    if(maxWidth > Surface.getBoundingBoxOverhang() && height<maxHeight)
    {
        thisNode->leftU_leftV =  newNode(thisNode,thisNode->minU,thisNode->divideBoundU,thisNode->minV,thisNode->divideBoundV);
        thisNode->rightU_leftV = newNode(thisNode,thisNode->divideBoundU,thisNode->maxU,thisNode->minV,thisNode->divideBoundV);
        thisNode->leftU_rightV = newNode(thisNode,thisNode->minU,thisNode->divideBoundU,thisNode->divideBoundV,thisNode->maxV);
        thisNode->rightU_rightV = newNode(thisNode,thisNode->divideBoundU,thisNode->maxU,thisNode->divideBoundV,thisNode->maxV);
        
        constructTree(thisNode->leftU_leftV,height+1);
        constructTree(thisNode->rightU_leftV,height+1);
        constructTree(thisNode->leftU_rightV,height+1);
        constructTree(thisNode->rightU_rightV,height+1);
    }
}

DynamicList<FixedList<scalar,2>> Foam::QuadTree::nearestPoints
(
    vector point
) const
{
    DynamicList<FixedList<scalar,2>> coordNurbs;
    traverseBsTree(root,point,coordNurbs);
    return coordNurbs;
}

void Foam::QuadTree::traverseBsTree
(
    Node* currentNode,
    vector point,
    DynamicList<FixedList<scalar,2>>& coordNurbs
) const
{
    if(!currentNode->MinMaxBox.isInside(point))
    {
        //Info<<"Point "<<point<<" not in Box "<<currentNode->MinMaxBox.Min<<" "<<currentNode->MinMaxBox.Max<<endl;
        return;
    }
    if(currentNode->leftU_leftV == _nil && currentNode->rightU_leftV == _nil && currentNode->leftU_rightV == _nil && currentNode->rightU_rightV == _nil)
    {
        coordNurbs.append({currentNode->divideBoundU,currentNode->divideBoundV});
    }
    else if(currentNode->leftU_leftV != _nil && currentNode->rightU_leftV != _nil && currentNode->leftU_rightV != _nil && currentNode->rightU_rightV != _nil)
    {
        traverseBsTree(currentNode->leftU_leftV,point,coordNurbs);
        traverseBsTree(currentNode->rightU_leftV,point,coordNurbs);
        traverseBsTree(currentNode->leftU_rightV,point,coordNurbs);
        traverseBsTree(currentNode->rightU_rightV,point,coordNurbs);
    }
    else
        FatalErrorInFunction<<"Unbalanced Quadtree"<<exit(FatalError);
}

FixedList<scalar,2> Foam::QuadTree::closestParaOnNurbsToPoint
(
    vector point
) const
{
    //Info<<"\tPoint: "<<point<<endl;
    DynamicList<FixedList<scalar,2>> testUV = nearestPoints(point);
    if(testUV.size() == 0)
        return {Surface.min_U()-1,Surface.min_V()-1};
    
    /*
    Info<<"BsTree nearest Points: ";
    for(int i=0;i<testU.size();i++)
        Info<<testU[i]<<endl;
    */
    
    DynamicList<FixedList<scalar,2>> uv_min_List;
    uv_min_List.setCapacity(testUV.size());
    FixedList<scalar,2> uv_min;
    for(int i=0;i<testUV.size();i++)
    {
        uv_min = Surface.newtonIterateNearestNeighbour(testUV[i][0],testUV[i][1],point);
        uv_min_List.append(uv_min);
        
        /*
         * deprecated because of instability
        while(i<testU.size() && u_min > testU[i])
        {
            i++;
        }
        */
    }
    //Info<<"U_min: "<<u_min<<endl;
    DynamicList<scalar> uv_min_Dist;
    uv_min_Dist.setCapacity(uv_min_List.size());
    for(int i=0;i<uv_min_List.size();i++)
    {
        uv_min_Dist.append(euklidianNorm(Surface.Surface_Derivative(0,0,uv_min_List[i][0],uv_min_List[i][1])-point));
    }
    scalar uv_min_Dist_min = uv_min_Dist[0];
    FixedList<scalar,2> uv_min_min = uv_min_List[0];
    for(int i=0;i<uv_min_Dist.size();i++)
    {
        if(uv_min_Dist[i] < uv_min_Dist_min)
        {
            uv_min_Dist_min = uv_min_Dist[i];
            uv_min_min = uv_min_List[i];
        }
    }
    //Info<<"\tRes: "<<u_min_min<<endl;
    return uv_min_min;    
}
