#include "BsTree.H"

Foam::BsTree::BsTree
(
    Nurbs Curve,
    label nbrSplitsBetweenCPs
):
Curve(Curve),
maxHeight(nbrSplitsBetweenCPs*this->Curve.nbrKnots()*this->Curve.degree())
{
    //Info<<"Constructe the tree"<<endl;
    _nil =  new Node();
    root = newNode(_nil,this->Curve.min_U(),this->Curve.max_U());
    //Info<<"Create first box from "<<this->Curve->min_U()<<" to "<<this->Curve->max_U()<<endl;
    //Info<<"Create first box from "<<root->min<<" to "<<root->max<<"//"<<(root->max == 1)<<endl;
    constructTree(root);
    //Info<<"Construction Tree done"<<endl;
}

Foam::BsTree::~BsTree()
{
    recursiveNodeDeleter(root);
    delete _nil;
}

void Foam::BsTree::recursiveNodeDeleter
(
    Node* 
    thisNode
)
{
    if(thisNode->left != _nil)
        recursiveNodeDeleter(thisNode->left);
    if(thisNode->right != _nil)
        recursiveNodeDeleter(thisNode->right);
    delete thisNode;     
}

Foam::BsTree::Node *Foam::BsTree::newNode
(
    Node* parent,
    scalar min,
    scalar max
)
{
    Node* newNodeItem = new Node();

    newNodeItem->left = _nil;
    newNodeItem->right = _nil;
    newNodeItem->parent = parent;
    newNodeItem->min = min;
    newNodeItem->max = max;
    newNodeItem->divideBound = 0.5*max+0.5*min;

    return newNodeItem;
}

void Foam::BsTree::constructTree
(
    Node* thisNode,
    int height
)
{
    //Info<<"Construct Tree at "<<height<<endl;
    thisNode->MinMaxBox = Curve.computeBoundingBox(thisNode->min,thisNode->max);
    //Info<<"Computed Box"<<endl;
    scalar maxWidth = 0;
    for(int d=0;d<3;d++)
    {
        maxWidth = std::max(thisNode->MinMaxBox.Max[d]-thisNode->MinMaxBox.Min[d],maxWidth);
    }
    //Info<<"Computed maxWidth "<<maxWidth<<endl;
    maxWidth = maxWidth-2*Curve.getBoundingBoxOverhang();
    if(maxWidth > Curve.getBoundingBoxOverhang() && height<maxHeight)
    {
        thisNode->left = newNode(thisNode,thisNode->min,thisNode->divideBound);
        thisNode->right = newNode(thisNode,thisNode->divideBound,thisNode->max);
        constructTree(thisNode->left,height+1);
        constructTree(thisNode->right,height+1);
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
