#ifndef BSTREE_H
#define BSTREE_H

#include "fvCFD.H"
#include <memory>
#include "Nurbs1D.H"

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
        const Nurbs1D Curve;

        Node* _nil;
        Node* root;
        
    protected:              
        const label maxHeight;
        
        BsTree(label maxHeight):
            Curve(Nurbs1D(scalarList(),List<vector>(),scalarList(),0)),
            maxHeight(maxHeight){};
        
    public:
        BsTree
        (
            Nurbs1D Curve,
            label nbrSplitsBetweenCPs = 1
        );
        ~BsTree();
        
        DynamicList<scalar> nearestPoints(vector point) const;
        scalar closestParaOnNurbsToPoint(vector point) const;
        
    protected:
        void constructTree(Node* thisNode, int height=0);
        void recursiveNodeDeleter(Node* thisNode);

    private:
        Node *newNode(Node* parent = nullptr,scalar min=0, scalar max=0);        
        
        void traverseBsTree
        (
            Node* currentNode,
            vector point,
            DynamicList<scalar>& coordNurbs
        ) const;
        
        friend void UnitTest_BsTree();
        friend class TESTBSTREE;

};

void UnitTest_BsTree();
}

#endif
