#include "BsTree.H"

Foam::BsTree::BsTree
(
    std::shared_ptr<Nurbs> Curve,
    label nbrSplitsBetweenCPs
):
Curve(Curve),
maxHeight(nbrSplitsBetweenCPs*this->Curve->nbrKnots()*this->Curve->degree())
{
    _nil =  new Node();
    root = newNode(_nil,this->Curve->min_U(),this->Curve->max_U());
    Info<<"Create first box from "<<this->Curve->min_U()<<" to "<<this->Curve->max_U()<<endl;
    Info<<"Create first box from "<<root->min<<" to "<<root->max<<"//"<<(root->max == 1)<<endl;
    constructTree(root);
}

Foam::BsTree::~BsTree()
{
    recursiveNodeDeleter(root);
    delete _nil;
}

void Foam::BsTree::recursiveNodeDeleter(Node* thisNode)
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

void Foam::BsTree::constructTree(Node* thisNode,int height)
{
    thisNode->MinMaxBox = Curve->computeBoundingBox(thisNode->min,thisNode->max);
    scalar maxWidth = 0;
    for(int d=0;d<3;d++)
    {
        maxWidth = std::max(thisNode->MinMaxBox.Max[d]-thisNode->MinMaxBox.Min[d],maxWidth);
    }
    maxWidth = maxWidth-2*Curve->getBoundingBoxOverhang();
    if(maxWidth > Curve->getBoundingBoxOverhang() && height<maxHeight)
    {
        thisNode->left = newNode(thisNode,thisNode->min,thisNode->divideBound);
        thisNode->right = newNode(thisNode,thisNode->divideBound,thisNode->max);
        constructTree(thisNode->left,height+1);
        constructTree(thisNode->right,height+1);
    }
}

scalarList Foam::BsTree::nearestPoint(vector point) const
{
    scalarList coordNurbs(0);
    traverseBsTree(root,point,coordNurbs);
    return coordNurbs;
}

void Foam::BsTree::traverseBsTree
(
    Node* currentNode,
    vector point,
    scalarList& coordNurbs
) const
{
    if(!currentNode->MinMaxBox.isInside(point))
    {
        Info<<"Point "<<point<<" not in Box "<<currentNode->MinMaxBox.Min<<" "<<currentNode->MinMaxBox.Max<<endl;
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

scalar Foam::BsTree::closestParaOnNurbsToPoint(vector point) const
{
    scalarList testPoints = nearestPoint(point);
    for(int i=0;i<testPoints.size();i++)
        Info<<testPoints[i]<<endl;
    scalarList u_min_List(0);
    scalar u_min;
    for(int i=0;i<testPoints.size();i++)
    {
        u_min = Curve->newtonIterateNearestNeighbour(testPoints[i],point);
        u_min_List.append(u_min);
        while(u_min > testPoints[i])
            i++;
    }
    scalarList u_min_Dist(0);
    for(int i=0;i<u_min_List.size();i++)
    {
        u_min_Dist.append(euklidianNorm(Curve->Curve_Derivative(0,u_min_List[i])-point));
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
    return u_min_min;    
}

scalar minDistanceToPoint(vector point) const
{
    return euklidianNorm(Curve->Curve_Derivative(0,closestParaOnNurbsToPoint(point))-point); 
}
