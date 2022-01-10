#ifndef BSTREE_H
#define BSTREE_H

#include "fvCFD.H"
#include <memory>
#include "Nurbs.H"

namespace Foam
{
class BsTree
{
    private:
        struct Node
        {
            Node* parent;
            Node* left;
            Node* right;
            scalar min;
            scalar max;
            scalar divideBound;
            BoundingBox MinMaxBox;
        };
                
        Node* _nil;
        Node* root;
                
        const Nurbs Curve;
        const label maxHeight;
        
    public:
        BsTree
        (
            Nurbs Curve,
            label nbrSplitsBetweenCPs = 1
        );
        ~BsTree();
        
        DynamicList<scalar> nearestPoints(vector point) const;
        scalar closestParaOnNurbsToPoint(vector point) const;
                
    private:
        Node *newNode(Node* parent = nullptr,scalar min=0, scalar max=0);
        void constructTree(Node* thisNode, int height=0);
        
        void recursiveNodeDeleter(Node* thisNode);
        
        void traverseBsTree
        (
            Node* currentNode,
            vector point,
            DynamicList<scalar>& coordNurbs
        ) const;
        
        friend void UnitTest_BsTree();
};

void UnitTest_BsTree();
}

#endif
