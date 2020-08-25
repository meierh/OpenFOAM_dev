#ifndef KDTREE_H
#define KDTREE_H

#include "fvCFD.H"
#include <memory>
#include "Nurbs.H"

namespace Foam
{   

class KdTree
{
    private:
        struct Node
        {
            Node* parent;
            Node* left;
            Node* right;
            BoundingBox MinMaxBox;
            labelList nurbsCurves;
        };
        
        Node* _nil;
        Node* root;
        
    private:
        std::unique_ptr<List<Nurbs*>> Items;
        List<Foam::BoundingBox> listMinMaxBoxes;
        
    public:
        explicit KdTree(std::unique_ptr<List<Nurbs*>> Items);
        
    private:
        Node *newNode(labelList nurbsCurves,Node* parent);
        void constructTree(Node* thisNode, labelList nurbsCurves, label treeHeight);
};
}

#endif
