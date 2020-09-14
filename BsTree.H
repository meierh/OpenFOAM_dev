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
        };
                
        Node* _nil;
        Node* root;
                
        const std::unique_ptr<Nurbs> Curve;
        const label nbrSplitsBetweenCPs;
        
    public:
        explicit BsTree
        (
            std::unique_ptr<Nurbs> Curve,
            label nbrSplitsBetweenCPs
        );
        ~BsTree();
        
        scalar nearestPoint(vector point);
                
    private:
        Node *newNode(Node* parent);
        void constructTree
        (
            Node* thisNode
        );
        
        void recursiveNodeDeleter(Node* thisNode);
        
        friend void UnitTest_BsTree();
};
}

#endif
