#include "QuadTree.H"

Foam::QuadTree::QuadTree
(
    Nurbs2D Surface,
    label nbrSplitsBetweenCPs
):
Surface(Surface),
BsTree(std::max(nbrSplitsBetweenCPs*this->Surface.nbrKnotsU()*this->Surface.degreeU(),nbrSplitsBetweenCPs*this->Surface.nbrKnotsV()*this->Surface.degreeV()))
{
    _nil =  new Node();
    root = newNode(_nil,this->Surface.min_U(),this->Surface.max_U(),this->Surface.min_V(),this->Surface.max_V());
    //Info<<"Create first box from "<<this->Curve->min_U()<<" to "<<this->Curve->max_U()<<endl;
    //Info<<"Create first box from "<<root->min<<" to "<<root->max<<"//"<<(root->max == 1)<<endl;
    constructTree(root);
    //Info<<"Construction Tree done"<<endl;
}

Foam::QuadTree::~QuadTree()
{
    recursiveNodeDeleter(root);
    delete _nil;
}

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
    //Info<<"Construct Tree at "<<height<<endl;
    thisNode->MinMaxBox = Surface.computeBoundingBox(thisNode->minU,thisNode->minV,thisNode->maxU,thisNode->maxV);

    //Info<<"Computed Box"<<endl;
    scalar maxWidth = 0;
    for(int d=0;d<3;d++)
    {
        maxWidth = std::max(thisNode->MinMaxBox.Max[d]-thisNode->MinMaxBox.Min[d],maxWidth);
    }
    //Info<<"Computed maxWidth "<<maxWidth<<endl;
    maxWidth = maxWidth-2*Surface.getBoundingBoxOverhang();
    if(maxWidth > Surface.getBoundingBoxOverhang() && height<maxHeight)
    {
        thisNode->leftU_leftV =  newNode(thisNode,thisNode->minU,thisNode->divideBoundU,thisNode->minV,thisNode->divideBoundV);
        thisNode->rightU_leftV = newNode(thisNode,thisNode->divideBoundU,thisNode->maxU,thisNode->minV,thisNode->divideBoundV);
        thisNode->leftU_rightV = newNode(thisNode,thisNode->minU,thisNode->divideBoundU,thisNode->divideBoundV,thisNode->maxV);
        thisNode->rightU_rightV =newNode(thisNode,thisNode->divideBoundU,thisNode->maxU,thisNode->divideBoundV,thisNode->maxV);
        
        constructTree(thisNode->leftU_leftV,height+1);
        constructTree(thisNode->rightU_leftV,height+1);
        constructTree(thisNode->leftU_rightV,height+1);
        constructTree(thisNode->rightU_rightV,height+1);
    }
}

DynamicList<scalar> Foam::BsTree::nearestPoints
(
    vector point
) const
{
    DynamicList<scalar> coordNurbs;
    traverseBsTree(root,point,coordNurbs);
    return coordNurbs;
}

void Foam::BsTree::traverseBsTree
(
    Node* currentNode,
    vector point,
    DynamicList<scalar>& coordNurbs
) const
{
    if(!currentNode->MinMaxBox.isInside(point))
    {
        //Info<<"Point "<<point<<" not in Box "<<currentNode->MinMaxBox.Min<<" "<<currentNode->MinMaxBox.Max<<endl;
        return;
    }
    if(currentNode->left == _nil && currentNode->right == _nil)
    {
        coordNurbs.append(currentNode->divideBound);
    }
    else if(currentNode->left != _nil && currentNode->right != _nil)
    {
        traverseBsTree(currentNode->left,point,coordNurbs);
        traverseBsTree(currentNode->right,point,coordNurbs);
    }
}

scalar Foam::BsTree::closestParaOnNurbsToPoint
(
    vector point
) const
{
    //Info<<"\tPoint: "<<point<<endl;
    DynamicList<scalar> testU = nearestPoints(point);
    if(testU.size() == 0)
        return Curve.min_U()-1;
    
    /*
    Info<<"BsTree nearest Points: ";
    for(int i=0;i<testU.size();i++)
        Info<<testU[i]<<endl;
    */
    
    DynamicList<scalar> u_min_List;
    u_min_List.setCapacity(testU.size());
    scalar u_min;
    for(int i=0;i<testU.size();i++)
    {
        u_min = Curve.newtonIterateNearestNeighbour(testU[i],point);
        u_min_List.append(u_min);
        
        /*
         * deprecated because of instability
        while(i<testU.size() && u_min > testU[i])
        {
            i++;
        }
        */
    }
    //Info<<"U_min: "<<u_min<<endl;
    DynamicList<scalar> u_min_Dist;
    u_min_Dist.setCapacity(u_min_List.size());
    for(int i=0;i<u_min_List.size();i++)
    {
        u_min_Dist.append(euklidianNorm(Curve.Curve_Derivative(0,u_min_List[i])-point));
    }
    scalar u_min_Dist_min = u_min_Dist[0];
    scalar u_min_min = u_min_List[0];
    for(int i=0;i<u_min_Dist.size();i++)
    {
        if(u_min_Dist[i] < u_min_Dist_min)
        {
            u_min_Dist_min = u_min_Dist[i];
            u_min_min = u_min_List[i];
        }
    }
    //Info<<"\tRes: "<<u_min_min<<endl;
    return u_min_min;    
}
