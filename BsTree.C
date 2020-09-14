#include "BsTree.H"

Foam::BsTree::BsTree
(
    std::unique_ptr<Nurbs> Curve,
    label nbrSplitsBetweenCPs
):
Curve(std::move(Curve)),
nbrSplitsBetweenCPs(nbrSplitsBetweenCPs)
{
    _nil =  new Node();
    root = newNode(_nil);
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
    Node* parent
)
{
    Node* newNodeItem = new Node();
    
    newNodeItem->left = _nil;
    newNodeItem->right = _nil;
    newNodeItem->parent = parent;

    return newNodeItem;
}

void Foam::BsTree::constructTree(Node* thisNode)
{
}

