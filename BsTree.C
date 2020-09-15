#include "BsTree.H"

Foam::BsTree::BsTree
(
    std::unique_ptr<Nurbs> Curve,
    label nbrSplitsBetweenCPs
):
Curve(std::move(Curve)),
maxHeight(nbrSplitsBetweenCPs*this->Curve->nbrKnots()*this->Curve->degree())
{
    _nil =  new Node();
    root = newNode(_nil,this->Curve->min_U(),this->Curve->max_U());
    Info<<"Create first box from "<<this->Curve->min_U()<<" to "<<this->Curve->max_U()<<endl;
    Info<<"Create first box from "<<root->min<<" to "<<root->max<<endl;
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
    newNodeItem->divideBound = 0.5*(max-min);

    return newNodeItem;
}

void Foam::BsTree::constructTree(Node* thisNode,int height)
{
    thisNode->MinMaxBox = Curve->computeBoundingBox(thisNode->min,thisNode->max);
    if(height<maxHeight)
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
